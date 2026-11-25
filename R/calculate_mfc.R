#' Calculate metabolite-feature candidates from spectral match evidence
#'
#' Run the locked metabolite-feature candidate scorer and retain both all scored
#' candidates and the registered embedding-library filtered view. This entry
#' point accepts the result of match_spectral_library or a compatible step6
#' object. The object's pseudo-MS2 result is converted to the
#' scorer query format, and the input pixels are checked against the step6 object.
#' For features with measured centroids, saved step6 queries must agree with
#' those masses. Rebuild spectral matches if the result is unresolved.
#'
#' @param object A mass_dataset with ROI and generate_pseudo_ms2 history.
#' @param step6 A compatible legacy msi_dataset after step6; when NULL, use
#'   the most recent match_spectral_library result stored in object.
#' @param output_dir Existing writable directory for the scorer audit tables.
#' @param resources Optional named list or YAML file supplying `matching_config`, `ms1_library`, and `vector_ids`.
#'   Each entry is a path or a list with `path` and optional `sha256`.
#'   Relative paths in YAML resolve against the YAML file's directory.
#'   With NULL, use the package reference lock and verify its SHA-256 hashes.
#'   Set `MSIDATASET_REFERENCE_ROOT` to the root of the reference archive.
#'   An optional `ms1_project_root` defaults to the MS1 library directory.
#' @details Only negative ion mode is currently supported (`[M−H]⁻`, `[M+Cl]⁻`).
#'   If supplied, the resource declaration `ion_mode` must be `"negative"`.
#' @return The mass_dataset with candidates, excluded rows, and the full MFC
#'   result recorded in process_info.
#' @export
calculate_mfc <- function(object, step6 = NULL, output_dir, resources = NULL) {
  if (is.null(step6)) {
    keys <- grep('^match_spectral_library(_[0-9]+)?$', names(object@process_info), value=TRUE)
    if (!length(keys)) stop('No spectral library match result is available', call.=FALSE)
    result <- object@process_info[[tail(keys, 1L)]]$result
    if (!identical(result$status, 'OK')) return(.msidataset_append_step(object,
      'calculate_mfc', list(status='UNRESOLVED_STEP6', upstream=result), list(status='UNRESOLVED_STEP6')))
    step6 <- result$step6
  }
  .msidataset_validate(object, require_roi = TRUE)
  if (!methods::is(step6, 'msi_dataset')) stop('A compatible step6 object is required', call.=FALSE)
  if (length(output_dir) != 1L || !dir.exists(output_dir) || file.access(output_dir, 2L) != 0L)
    stop('output_dir must be an existing writable directory', call.=FALSE)
  pseudo_steps <- grep('^generate_pseudo_ms2(_[0-9]+)?$', names(object@process_info), value=TRUE)
  if (!length(pseudo_steps) || is.null(object@process_info[[tail(pseudo_steps, 1L)]]$result))
    stop('generate_pseudo_ms2 must run on the object first', call.=FALSE)
  pairs <- object@process_info[[tail(pseudo_steps, 1L)]]$result
  if (!is.data.frame(pairs)) return(.msidataset_append_step(object, 'calculate_mfc',
    list(status='UNRESOLVED_PSEUDO', upstream=pairs), list(status='UNRESOLVED_PSEUDO')))
  vars <- object@variable_info
  if (!'mz' %in% names(vars)) stop('Feature m/z values are required', call.=FALSE)
  masses <- .msidataset_feature_mz(vars)
  if (any(!is.finite(masses$mz))) return(.msidataset_append_step(object, 'calculate_mfc',
    list(status='UNRESOLVED_MZ', invalid_features=which(!is.finite(masses$mz))), list(status='UNRESOLVED_MZ')))
  old_vars <- step6@variable_info
  if (!identical(as.character(vars$variable_id), as.character(old_vars$variable_id)) ||
      !isTRUE(all.equal(as.numeric(vars$mz), as.numeric(old_vars$mz), tolerance=1e-12)))
    stop('step6 feature identifiers or m/z do not match the object', call.=FALSE)
  if (anyNA(object@sample_info$region))
    return(.msidataset_append_step(object, 'calculate_mfc',
      list(status='UNRESOLVED_ROI', missing_roi=which(is.na(object@sample_info$region))),
      list(status='UNRESOLVED_ROI')))
  use <- object@sample_info$region != 'VESSEL'
  s <- object@sample_info[use,,drop=FALSE]
  if (nrow(s) != nrow(step6@sample_info) ||
      !all(s$x == step6@sample_info$x & s$y == step6@sample_info$y))
    stop('step6 coordinates do not match non-vessel pixels', call.=FALSE)
  blocks <- split(seq_len(nrow(vars)), ceiling(seq_len(nrow(vars)) / 16L))
  for (j in blocks) {
    a <- as.matrix(object@expression_data[j, use, drop=FALSE])
    b <- as.matrix(step6@expression_data[j, , drop=FALSE])
    if (!isTRUE(all.equal(unname(a), unname(b), tolerance=0, check.attributes=FALSE)))
      stop('step6 intensities do not match the object', call.=FALSE)
  }
  query <- .mfc_prepare_pseudo(object, pairs)
  if (any(query$preparation_status != 'ok')) {
    return(.msidataset_append_step(object, 'calculate_mfc',
      list(status='UNRESOLVED_QUERY', query=query), list(status='UNRESOLVED_QUERY')))
  }
  original_query <- step6@pseudo_ms2$pseudo_ms2_table
  pair_key <- function(z) paste(z$precursor_feature_id, z$fragment_feature_id, sep='|')
  if (!is.data.frame(original_query) ||
      !all(c('precursor_feature_id','fragment_feature_id','SAS','intensity') %in% names(original_query)) ||
      anyDuplicated(pair_key(original_query)) || anyDuplicated(pair_key(query)))
    stop('step6 pseudo-MS2 evidence cannot be verified', call.=FALSE)
  ix <- match(pair_key(query), pair_key(original_query))
  if (nrow(query) != nrow(original_query) || anyNA(ix) ||
      !isTRUE(all.equal(query$SAS, original_query$SAS[ix], tolerance=1e-6)) ||
      !isTRUE(all.equal(query$intensity, original_query$intensity[ix], tolerance=1e-6)))
    stop('step6 pseudo-MS2 evidence does not match the object', call.=FALSE)
  full <- stats::setNames(masses$observed_mz_full, vars$variable_id)
  precursor_full <- is.finite(full[query$precursor_feature_id])
  fragment_full <- is.finite(full[query$fragment_feature_id])
  same_mass <- function(field, measured) {
    !any(measured) || isTRUE(all.equal(query[[field]][measured],
      original_query[[field]][ix[measured]], tolerance=1e-12))
  }
  if (!same_mass('anchor_mz', precursor_full) ||
      !same_mass('peak_mz', fragment_full) ||
      !same_mass('fragment_mz', fragment_full)) {
    status <- 'UNRESOLVED_STEP6_MZ'
    return(.msidataset_append_step(object, 'calculate_mfc',
      list(status=status, query=query, upstream=original_query), list(status=status)))
  }
  # The input object supplies measured masses even when a compatible saved
  # step6 omitted them for features absent from the pseudo-MS2 query.
  measured <- is.finite(masses$observed_mz_full)
  if (any(measured)) {
    if (!'observed_mz_full' %in% names(old_vars)) old_vars$observed_mz_full <- NA_real_
    old_vars$observed_mz_full[measured] <- masses$observed_mz_full[measured]
    step6@variable_info <- old_vars
  }
  step6@pseudo_ms2$pseudo_ms2_table <- query
  lock <- .msi_resource_entries(resources)
  required <- c('matching_config', 'ms1_library', 'ms1_project_root', 'vector_ids')
  provenance <- .msi_verify_resources(lock, required, locked = is.null(resources))
  paths <- lapply(lock[required], `[[`, 'path')
  cfg <- yaml::read_yaml(paths$matching_config)
  .msi_check_ion_mode(cfg$mfc$ion_mode)
  .msi_check_ion_mode(cfg$run$ion_mode)
  if (is.null(resources)) .msi_check_locked_ms1_path(cfg$mfc$ms1_library_path, lock$ms1_library)
  adducts <- c('[M-H]-','[M+Cl]-')
  fmc <- utils::modifyList(cfg$mfc,list(enabled=TRUE,ion_mode='negative',
    ms1_flat_project_root=paths$ms1_project_root,ms1_library_path=paths$ms1_library,
    sample_mode=cfg$run$sample_mode,analysis_domain=cfg$run$analysis_domain,
    ms1_mz_source=cfg$mfc$ms1_mz_source,allowed_adducts=adducts))
  ms1 <- list(enabled=TRUE,library_path=paths$ms1_library,ppm=cfg$mfc$ms1_ppm,
    ppm_sigma=cfg$mfc$ppm_sigma,max_candidates_per_feature=cfg$mfc$max_candidates_per_feature,
    allowed_adducts=adducts)
  scored <- withCallingHandlers(`7_annotate_compounds`(step6,reverse_masked_cosine_weight=.7,S_SAS_weight=.3,
    q_value_cutoff=cfg$annotation$q_value_cutoff,keep_all=TRUE,decoy_strategy='peak_shift',
    fdr_method='target_decoy',candidate_source_for_step7=cfg$annotation$candidate_source_for_step7,
    reverse_masked_cosine_threshold=0,S_SAS_threshold=cfg$annotation$candidate_sas_cutoff,
    primary_spectral_score='reverse_masked_cosine',ppm=15,ion_mode='negative',
    fmc_annotation=fmc,ms1_matching=ms1,fmc_ion_mode='negative',output_dir=output_dir),
    warning=function(w) {
      if (conditionMessage(w) %in% c(
        '7_annotate_compounds is retained for compatibility; use calculate_score for canonical scoring.',
        'FMC is now Module 8. Please call calculate_mfc() after calculate_score().'))
        invokeRestart('muffleWarning')
    })
  scored <- .legacy_calculate_mfc(scored, fmc_config=fmc, ion_mode='negative', output_dir=output_dir)
  # The reference workflow applied its filter after the scorer's CSV serialization. Preserve that
  # exact round trip: it can affect the equality of tied scores.
  all <- utils::read.csv(file.path(output_dir,'tables','mfc_final_candidates_all_ranked.csv'),check.names=FALSE)
  ids <- data.table::fread(paths$vector_ids,select='HMDB_ID',showProgress=FALSE)$HMDB_ID
  filtered <- .mfc_filter_candidates(all,ids)
  info <- list(status='OK',ion_mode='negative',adducts=adducts,source='Spectral match evidence',
               resource_source=if (is.null(resources)) 'reference_lock' else 'user',
               resources=provenance, time=as.character(Sys.time()))
  .msidataset_append_step(object,'calculate_mfc',list(status='OK',query=query,
      candidates_all=all,candidates=filtered$kept,removed=filtered$removed,
      mfc_details=scored@annotation_results$fmc),info)
}
