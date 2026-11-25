#' Match pseudo-MS2 spectra against the locked spectral library
#'
#' Run entropy retrieval and reverse spectral matching. The full input
#' object is retained; results and provenance are appended to process_info.
#' Spectra without a finite precursor cannot enter the registered precursor
#' search and are counted in the result. This function reads a large library.
#' Measured centroids are used for precursor and fragment matching when available.
#'
#' @param object A mass_dataset with ROI labels and generate_pseudo_ms2 history.
#' @param config A list containing an existing writable `output_dir`. Optionally
#'   `workers` sets the number of reverse matching workers (default 4). Resources
#'   and scoring parameters are read from the package's locked configuration.
#'   Set `MSIDATASET_D15_STUDY_ROOT` and `MSIDATASET_D15_REFERENCE_ROOT`
#'   to the directories containing the pinned spectral library and reference files.
#' @return The input mass_dataset with a match_spectral_library process record;
#'   its result contains `step6`, an msi_dataset usable by calculate_mfc, or an
#'   unresolved status with retained input evidence.
#' @export
match_spectral_library <- function(object, config) {
  .msidataset_validate(object, require_roi=TRUE)
  if (!is.list(config) || length(config$output_dir) != 1L ||
      !dir.exists(config$output_dir) || file.access(config$output_dir, 2L) != 0L)
    stop('config$output_dir must be an existing writable directory', call.=FALSE)
  unexpected <- setdiff(names(config), c('output_dir','workers'))
  if (length(unexpected)) stop('Unknown config fields: ', paste(unexpected, collapse=', '), call.=FALSE)
  workers <- if (is.null(config$workers)) 4L else config$workers
  if (length(workers)!=1L || is.na(workers) || !is.numeric(workers) || !is.finite(workers) || workers < 1 || workers != floor(workers))
    stop('config$workers must be a positive integer', call.=FALSE)
  keys <- grep('^generate_pseudo_ms2(_[0-9]+)?$', names(object@process_info), value=TRUE)
  if (!length(keys)) stop('generate_pseudo_ms2 must run first', call.=FALSE)
  pairs <- object@process_info[[tail(keys,1L)]]$result
  unresolved <- function(status, details) .msidataset_append_step(object, 'match_spectral_library',
    list(status=status, details=details), list(status=status))
  if (!is.data.frame(pairs)) return(unresolved('UNRESOLVED_PSEUDO', pairs))
  if (anyNA(object@sample_info$region)) return(unresolved('UNRESOLVED_ROI', which(is.na(object@sample_info$region))))
  v <- object@variable_info
  mz <- .msidataset_feature_mz(v)$mz
  if (!'mz' %in% names(v) || any(!is.finite(mz)))
    return(unresolved('UNRESOLVED_MZ', v))
  if (!nrow(pairs)) return(unresolved('UNRESOLVED_EMPTY_QUERY', pairs))
  bad_xy <- which(!is.finite(object@sample_info$x) | !is.finite(object@sample_info$y))
  if (length(bad_xy)) return(unresolved('UNRESOLVED_COORDINATES', bad_xy))
  expression <- object@expression_data
  bad_intensity <- which(vapply(seq_len(ncol(expression)), function(i)
    any(!is.finite(if (is.data.frame(expression)) expression[[i]] else expression[, i])), logical(1)))
  if (length(bad_intensity)) return(unresolved('UNRESOLVED_INTENSITY', bad_intensity))
  if (any(!is.finite(object@sample_info$tic) | object@sample_info$tic <= 0))
    return(unresolved('UNRESOLVED_TIC', which(!is.finite(object@sample_info$tic) | object@sample_info$tic <= 0)))
  if (!any(object@sample_info$region != 'VESSEL')) return(unresolved('UNRESOLVED_EMPTY_ROI', object@sample_info$region))
  query <- .mfc_prepare_pseudo(object, pairs)
  if (any(query$preparation_status != 'ok')) return(unresolved('UNRESOLVED_QUERY', query))
  if (anyDuplicated(query[c('precursor_feature_id','fragment_feature_id')]))
    return(unresolved('UNRESOLVED_DUPLICATE_QUERY', query))
  if (any(!is.finite(query$anchor_mz) | !is.finite(query$peak_mz)))
    return(unresolved('UNRESOLVED_QUERY_MZ', query))
  lockfile <- system.file('config','locked_d15.yaml',package='msidataset')
  if (length(lockfile)!=1L || !file.exists(lockfile)) stop('Locked configuration is missing', call.=FALSE)
  lock <- .msi_locked_resources(lockfile)
  required <- c('ms2_library','adp_config','ms1_library')
  for (key in required) {
    entry <- lock[[key]]
    if (is.null(entry$path) || is.null(entry$sha256) || !file.exists(entry$path))
      stop('Locked resource path or SHA-256 is missing: ', key,
           '. Check ', entry$root, '.', call.=FALSE)
    if (!identical(digest::digest(file=entry$path,algo='sha256'),entry$sha256))
      stop('Locked resource hash mismatch: ',key,call.=FALSE)
  }
  cfg <- yaml::read_yaml(lock$adp_config$path)
  .msi_check_locked_ms1_path(cfg$mfc$ms1_library_path, lock$ms1_library)
  use <- object@sample_info$region != 'VESSEL'
  s <- object@sample_info[use,,drop=FALSE]
  v$feature_id <- v$variable_id
  v$observed_mz_full <- mz
  legacy <- methods::new('msi_dataset', expression_data=as.matrix(object@expression_data[,use,drop=FALSE]),
    sample_info=as.data.frame(s), variable_info=as.data.frame(v), pseudo_ms2=list(pseudo_ms2_table=query))
  lib <- readRDS(lock$ms2_library$path)
  spectra <- lib$library_spectrum_table
  eligible <- rep(FALSE,nrow(spectra))
  for (mz in unique(query$anchor_mz))
    eligible <- eligible | (is.finite(spectra$precursor_mz) & abs(spectra$precursor_mz-mz)<=mz*10e-6)
  ids <- spectra$spectrum_id[eligible]
  unsearchable <- data.frame(spectrum_id=spectra$spectrum_id[!is.finite(spectra$precursor_mz)],
    status=rep('UNSEARCHABLE_MISSING_PRECURSOR',sum(!is.finite(spectra$precursor_mz))), stringsAsFactors=FALSE)
  counts <- list(total_spectra=nrow(spectra), no_precursor=sum(!is.finite(spectra$precursor_mz)), eligible_spectra=length(ids))
  view <- list(library_spectrum_table=spectra[eligible,,drop=FALSE],
    library_peak_table=lib$library_peak_table[lib$library_peak_table$spectrum_id %in% ids,,drop=FALSE],
    flash_entropy_index=lib$flash_entropy_index[lib$flash_entropy_index$spectrum_id %in% ids,,drop=FALSE])
  attr(view,'msidataset_standardized_library') <- TRUE
  attr(view,'msidataset_hmdb_processed_long') <- TRUE
  rm(lib,spectra); gc()
  contract <- list(fragment_retrieval=list(ppm=15,low_mz_da_floor=0),
    fragment_evidence=list(ppm=15,low_mz_da_floor=0),
    precursor_candidate=list(ppm=10,da_floor_mode='explicit',da_floor=0))
  legacy <- `5_calculate_flash_entropy`(legacy,spectral_library=view,ppm=15,mz_tol_da=NULL,
    top_n_entropy_candidates=cfg$flash_entropy$top_n_candidates,top_n_spectra=cfg$flash_entropy$top_n_spectra,
    top_n_compounds=cfg$flash_entropy$top_n_compounds,max_spectra_per_compound=cfg$flash_entropy$max_spectra_per_compound,
    min_library_peaks=cfg$library$min_peaks,precursor_filter=TRUE,precursor_tol=.05,polarity='negative',
    min_match_fe=cfg$flash_entropy$min_matches,min_s_fe=cfg$flash_entropy$min_score,
    output_dir=config$output_dir,tolerance_contract=contract)
  legacy <- `6_reverse_match`(legacy,use_flash_hits=TRUE,use_archive_hits=FALSE,fragment_tolerance_ppm=15,
    precursor_tol=.05,fragment_tolerance_mode='ppm',reverse_filter_mode='reverse_masked_cosine_evidence',
    min_matched_peaks=2,min_matched_peaks_evidence=2,
    min_reverse_masked_cosine_evidence=cfg$reverse_matching$min_reverse_masked_cosine,
    matched_library_intensity_fraction_min=0,moderate_min_matched_peaks=2,moderate_max_precursor_ppm=50,
    workers=as.integer(workers),output_dir=config$output_dir,tolerance_contract=contract)
  .msidataset_append_step(object,'match_spectral_library',list(status='OK',step6=legacy,query=query,library_counts=counts,
      unsearchable_spectra=unsearchable,
      queries_without_candidates=setdiff(unique(query$anchor_feature_id),
        unique(legacy@annotation_results$entropy_candidates$anchor_feature_id))),
    list(status='OK',resources=lock[required],workers=as.integer(workers),tolerance_contract=contract))
}
