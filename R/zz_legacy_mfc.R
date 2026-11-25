#' Internal metabolite-feature scoring implementation
#' @keywords internal
#' @noRd
NULL

# Internal helpers for candidate scoring and evidence aggregation.
# Runtime sampling and source-name attribution are in inst/design/runtime_trace_mfc.tsv.
# Other branches in this closure are retained conservatively. No functions here are exported.
# Scoring helpers retained for compatibility with saved evidence.
`.ann_add_annotation_confidence` <- 
function(x, confidence_thresholds = list()) {
  out <- as.data.frame(x, stringsAsFactors = FALSE, check.names = FALSE)
  out <- .ann_add_missing_columns(out, c(
    "reverse_masked_cosine", "S_SAS", "S_total", "rank_within_feature",
    "S_lib_coverage", "matched_peak_fraction", "precursor_error_ppm",
    "reverse_match_status", "match_mode", "n_matched_peaks"
  ))
  get_threshold <- function(name, default) {
    value <- confidence_thresholds[[name]]
    value <- suppressWarnings(as.numeric(value))
    if (length(value) == 0L || !is.finite(value[[1]])) {
      return(default)
    }
    value[[1]]
  }
  high_reverse <- get_threshold("high_reverse_masked_cosine", 0.05)
  high_ppm <- get_threshold("high_precursor_ppm", 50)
  medium_reverse <- get_threshold("medium_reverse_masked_cosine", 0)
  medium_ppm <- get_threshold("medium_precursor_ppm", 100)

  reverse_masked_cosine <- suppressWarnings(as.numeric(out$reverse_masked_cosine))
  sas <- suppressWarnings(as.numeric(out$S_SAS))
  s_total <- suppressWarnings(as.numeric(out$S_total))
  rank_within_feature <- suppressWarnings(as.integer(out$rank_within_feature))
  n_matched <- suppressWarnings(as.integer(out$n_matched_peaks))
  precursor_ppm <- abs(suppressWarnings(as.numeric(out$precursor_error_ppm)))
  status <- trimws(as.character(out$reverse_match_status))
  status[is.na(status)] <- ""
  matched <- status == "matched" | n_matched > 0
  unknown <- status == ""
  precursor_ok_high <- is.finite(precursor_ppm) & precursor_ppm <= high_ppm
  precursor_ok_medium <- is.finite(precursor_ppm) & precursor_ppm <= medium_ppm

  high <- matched &
    n_matched >= 3L &
    is.finite(reverse_masked_cosine) & reverse_masked_cosine >= high_reverse &
    is.finite(sas) & sas >= 0.60 &
    is.finite(s_total) & s_total >= 0.20 &
    (is.na(rank_within_feature) | rank_within_feature <= 3L) &
    (precursor_ok_high | !is.finite(precursor_ppm))

  medium <- !high & matched &
    n_matched >= 3L &
    is.finite(reverse_masked_cosine) & reverse_masked_cosine > medium_reverse &
    (precursor_ok_medium | !is.finite(precursor_ppm))

  review <- !high & !medium & matched & (
    (n_matched == 2L & is.finite(reverse_masked_cosine) & reverse_masked_cosine >= 0.05) |
      (n_matched == 1L & is.finite(reverse_masked_cosine) & reverse_masked_cosine >= 0.20 &
         is.finite(sas) & sas >= 0.80)
  )

  out$confidence_tier <- ifelse(
    high, "Tier1_pseudoMS2_spatial",
    ifelse(medium, "Tier2_pseudoMS2_candidate",
      ifelse(review, "Tier3_review_only", ifelse(unknown, "unknown", "low")))
  )
  out$annotation_confidence <- ifelse(high, "high", ifelse(medium, "medium", ifelse(review, "review_only", ifelse(unknown, "unknown", "low"))))
  out
}

`.ann_add_missing_columns` <- 
function(df, cols) {
  df <- as.data.frame(df, stringsAsFactors = FALSE, check.names = FALSE)
  for (col in cols) {
    if (!col %in% names(df)) {
      df[[col]] <- rep(NA, nrow(df))
    }
  }
  df
}

`.ann_add_rank_total` <- 
function(out) {
  out <- as.data.frame(out, stringsAsFactors = FALSE, check.names = FALSE)
  if (nrow(out) == 0L) {
    out$rank_total <- integer()
    return(out)
  }
  out <- .ann_add_missing_columns(out, c("query_id", "cluster_id", "anchor_feature_id", "S_total"))
  if (all(is.na(out$query_id))) {
    out$query_id <- as.character(out$cluster_id)
  }
  split_query <- .ann_split_key(out, c("query_id", "anchor_feature_id"))
  out$rank_total <- unsplit(lapply(split(seq_len(nrow(out)), split_query), function(idx) {
    score <- suppressWarnings(as.numeric(out$S_total[idx]))
    score[!is.finite(score)] <- -Inf
    rank(-score, ties.method = "first")
  }), split_query)
  out$rank_total <- suppressWarnings(as.integer(out$rank_total))
  out
}

`.ann_build_mfc_input_from_pseudo_evidence` <- 
function(object, lower_confidence_weight = 0.5) {
  final <- if (is.data.frame(object@annotation_results$final_annotation)) {
    object@annotation_results$final_annotation
  } else {
    data.frame()
  }
  candidates_source <- "final_annotation_pass"
  if (nrow(final) && "pass_annotation_filter" %in% names(final)) {
    source_rows <- final[final$pass_annotation_filter %in% TRUE, , drop = FALSE]
  } else {
    source_rows <- final
  }
  if (!nrow(source_rows)) {
    evidence <- if (is.data.frame(object@annotation_results$pseudo_ms2_candidate_evidence)) {
      object@annotation_results$pseudo_ms2_candidate_evidence
    } else {
      data.frame()
    }
    evidence_all <- if (is.data.frame(object@annotation_results$pseudo_ms2_candidate_evidence_all)) {
      object@annotation_results$pseudo_ms2_candidate_evidence_all
    } else {
      data.frame()
    }
    thresholds <- .ann_reverse_evidence_thresholds(object)
    source_rows <- extract_pseudo_ms2_evidence(
      evidence, evidence_all, evidence_weight = 50,
      min_matched_peaks_used = thresholds$min_matched_peaks,
      min_reverse_masked_cosine_used = thresholds$min_reverse_masked_cosine
    )
    candidates_source <- "pseudo_ms2_candidate_evidence_legacy"
  }
  if (nrow(source_rows) && !"source_row_id" %in% names(source_rows)) {
    source_rows$source_row_id <- paste0("step7_pass_row_", seq_len(nrow(source_rows)))
  }
  if (nrow(source_rows) && !"match_mode" %in% names(source_rows)) {
    source_rows$match_mode <- "pseudo_ms2"
  } else if (nrow(source_rows)) {
    missing_mode <- is.na(source_rows$match_mode) | !nzchar(as.character(source_rows$match_mode))
    source_rows$match_mode[missing_mode] <- "pseudo_ms2"
  }
  if (nrow(source_rows) && !"evidence_source" %in% names(source_rows)) {
    source_rows$evidence_source <- "pseudo_ms2_candidate_evidence"
  }
  candidates <- .ann_standardize_mfc_input_candidates(source_rows)
  candidates <- .mfc2_add_candidate_keys(candidates, feature_axis = .mfc2_feature_axis(object = object))
  n_pseudo <- sum(candidates$match_mode == "pseudo_ms2", na.rm = TRUE)
  n_ms1 <- sum(candidates$match_mode %in% c("ms1", "ms1_only"), na.rm = TRUE)
  candidate_hash <- .mfc2_candidate_hash(candidates)
  step7_pass_candidate_hash <- .mfc2_candidate_hash(candidates[candidates$match_mode == "pseudo_ms2", , drop = FALSE])
  summary <- data.frame(
    n_mfc_input_candidates = as.integer(nrow(candidates)),
    n_mfc_input_pseudo_ms2 = as.integer(n_pseudo),
    n_mfc_input_ms1_only = as.integer(n_ms1),
    ms1_only_status = if (n_ms1 > 0L) "generated" else "not_generated_or_no_eligible_features",
    lower_confidence_weight = as.numeric(lower_confidence_weight[[1]]),
    source = candidates_source,
    candidate_hash = candidate_hash,
    step7_pass_candidate_hash = step7_pass_candidate_hash,
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
  list(candidates = candidates, summary = summary)
}

`.ann_clip01_zero_na` <- 
function(x) {
  x <- suppressWarnings(as.numeric(x))
  x[is.na(x)] <- 0
  pmin(pmax(x, 0), 1)
}

`.ann_confidence_rank` <- 
function(confidence) {
  confidence <- as.character(confidence)
  out <- rep(1L, length(confidence))
  out[confidence == "medium"] <- 2L
  out[confidence == "high"] <- 3L
  out
}

`.ann_empty_decoy_library` <- 
function() {
  list(
    decoy_spectrum_table = data.frame(
      candidate_key = character(),
      decoy_candidate_key = character(),
      candidate_id = character(),
      spectrum_id = character(),
      candidate_name = character(),
      library_source = character(),
      database_source = character(),
      precursor_mz = numeric(),
      formula = character(),
      inchikey = character(),
      smiles = character(),
      adduct = character(),
      polarity = character(),
      is_decoy = logical(),
      stringsAsFactors = FALSE,
      check.names = FALSE
    ),
    decoy_peak_table = data.frame(
      decoy_candidate_key = character(),
      candidate_key = character(),
      spectrum_id = character(),
      peak_mz = numeric(),
      intensity_sum1 = numeric(),
      is_decoy = logical(),
      stringsAsFactors = FALSE,
      check.names = FALSE
    )
  )
}

`.ann_empty_decoy_scores` <- 
function() {
  data.frame(
    cluster_id = integer(),
    anchor_feature_id = character(),
    anchor_mz = numeric(),
    candidate_key = character(),
    decoy_candidate_key = character(),
    S_total = numeric(),
    S_FE = numeric(),
    reverse_masked_cosine = numeric(),
    primary_spectral_score = character(),
    S_SAS = numeric(),
    S_SAS_score_for_total = numeric(),
    S_SAS_missing = logical(),
    S_reverse_cos = numeric(),
    S_lib_coverage = numeric(),
    library_l2_coverage = numeric(),
    n_matched_peaks = integer(),
    query_match_ratio = numeric(),
    library_match_ratio = numeric(),
    S_spatial_reverse = numeric(),
    S_precursor = numeric(),
    S_chem = numeric(),
    s_fe = numeric(),
    s_arch = numeric(),
    s_rev = numeric(),
    u_spec = numeric(),
    s_chem = numeric(),
    s_spat = numeric(),
    S_FE_cal = numeric(),
    S_reverse_cos_cal = numeric(),
    S_lib_coverage_cal = numeric(),
    S_spatial_reverse_cal = numeric(),
    S_precursor_cal = numeric(),
    S_chem_cal = numeric(),
    s_final = numeric(),
    is_decoy = logical(),
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
}

`.ann_empty_final_annotation` <- 
function() {
  data.frame(
    cluster_id = integer(),
    query_id = character(),
    anchor_feature_id = character(),
    anchor_mz = numeric(),
    candidate_key = character(),
    candidate_id = character(),
    spectrum_id = character(),
    archive_id = character(),
    metabolite_name = character(),
    candidate_name = character(),
    library_source = character(),
    database_source = character(),
    precursor_mz = numeric(),
    formula = character(),
    inchikey = character(),
    smiles = character(),
    adduct = character(),
    polarity = character(),
    candidate_source = character(),
    entropy_S = numeric(),
    S_FE = numeric(),
    reverse_masked_cosine = numeric(),
    query_unmatched_peak_policy = character(),
    reference_unmatched_peak_policy = character(),
    primary_spectral_score = character(),
    S_SAS_raw = numeric(),
    S_SAS = numeric(),
    S_SAS_score_for_total = numeric(),
    S_SAS_missing = logical(),
    S_SAS_status = character(),
    S_total = numeric(),
    rank_total = integer(),
    rank_within_feature = integer(),
    rank_margin_to_next = numeric(),
    S_reverse_cos = numeric(),
    S_lib_coverage = numeric(),
    library_l2_coverage = numeric(),
    S_spatial_reverse = numeric(),
    S_precursor = numeric(),
    S_chem = numeric(),
    retrieval_rank = integer(),
    precursor_error_ppm = numeric(),
    matched_peak_fraction = numeric(),
    n_match_pairs = integer(),
    n_matched_peaks = integer(),
    n_reverse_match = integer(),
    unique_query_matched_peaks = integer(),
    unique_library_matched_peaks = integer(),
    n_query_peaks_total = integer(),
    n_library_peaks_total = integer(),
    query_match_ratio = numeric(),
    library_match_ratio = numeric(),
    query_peak_reuse_n = integer(),
    query_peak_reuse_fraction = numeric(),
    has_query_peak_reuse = logical(),
    n_SAS_matched_peaks = integer(),
    sum_matched_library_intensity = numeric(),
    matched_query_intensity = numeric(),
    matched_query_intensity_fraction = numeric(),
    matched_library_intensity = numeric(),
    matched_library_intensity_fraction = numeric(),
    reverse_match_status = character(),
    match_mode = character(),
    score_ms1 = numeric(),
    score_ms2 = numeric(),
    score_spatial = numeric(),
    score_total_pre_mfc = numeric(),
    pass_n_matched_peaks_ge3 = logical(),
    pass_reverse_masked_cosine = logical(),
    evidence_pass = logical(),
    evidence_failure_reason = character(),
    confidence_tier = character(),
    s_fe = numeric(),
    s_arch = numeric(),
    s_rev = numeric(),
    n_match = integer(),
    u_spec = numeric(),
    s_chem = numeric(),
    s_spat = numeric(),
    s_final = numeric(),
    annotation_q_value = numeric(),
    annotation_fdr_method = character(),
    rank_within_anchor = integer(),
    pass_annotation_filter = logical(),
    qvalue_enabled = logical(),
    passes_score_threshold = logical(),
    passes_q_value = logical(),
    passes_S_SAS_threshold = logical(),
    S_SAS_threshold = numeric(),
    pass_primary = logical(),
    pass_review = logical(),
    annotation_pass_class = character(),
    rule_fail_reason = character(),
    annotation_confidence = character(),
    S_FE_cal = numeric(),
    S_reverse_cos_cal = numeric(),
    S_lib_coverage_cal = numeric(),
    S_spatial_reverse_cal = numeric(),
    S_precursor_cal = numeric(),
    S_chem_cal = numeric(),
    s_fe_cal = numeric(),
    s_arch_cal = numeric(),
    s_rev_cal = numeric(),
    u_spec_cal = numeric(),
    s_chem_cal = numeric(),
    s_spat_cal = numeric(),
    chem_status = character(),
    chem_role = character(),
    chem_explanation = character(),
    role_delta_mz = numeric(),
    role_ppm_error = numeric(),
    role_evidence_n = integer(),
    spectral_shape_support = character(),
    passes_spectral_shape_filter = logical(),
    reverse_filter_mode = character(),
    threshold_validation_status = character(),
    selection_reason = character(),
    candidate_source_for_step7 = character(),
    diagnostic_only = logical(),
    annotation_status = character(),
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
}

`.ann_empty_mfc_input_candidates` <- 
function() {
  out <- .mfc2_empty_candidate()
  required <- .ann_mfc_input_required_columns()
  for (nm in setdiff(required, names(out))) {
    # Keep zero-row data frames zero-row.  Assigning scalar NA to a new
    # column is an error in base R when the data frame has no rows.
    out[[nm]] <- rep(NA, nrow(out))
  }
  out[, unique(c(required, setdiff(names(out), required))), drop = FALSE]
}

`.ann_empty_score_table` <- 
function() {
  data.frame(
    cluster_id = integer(),
    query_id = character(),
    anchor_feature_id = character(),
    anchor_mz = numeric(),
    candidate_key = character(),
    candidate_id = character(),
    spectrum_id = character(),
    archive_id = character(),
    metabolite_name = character(),
    candidate_name = character(),
    library_source = character(),
    database_source = character(),
    precursor_mz = numeric(),
    formula = character(),
    inchikey = character(),
    smiles = character(),
    adduct = character(),
    polarity = character(),
    candidate_source = character(),
    entropy_S = numeric(),
    S_FE = numeric(),
    reverse_masked_cosine = numeric(),
    query_unmatched_peak_policy = character(),
    reference_unmatched_peak_policy = character(),
    primary_spectral_score = character(),
    S_SAS_raw = numeric(),
    S_SAS = numeric(),
    S_SAS_score_for_total = numeric(),
    S_SAS_missing = logical(),
    S_SAS_status = character(),
    S_total = numeric(),
    rank_total = integer(),
    rank_within_feature = integer(),
    rank_margin_to_next = numeric(),
    S_reverse_cos = numeric(),
    S_lib_coverage = numeric(),
    library_l2_coverage = numeric(),
    S_spatial_reverse = numeric(),
    S_precursor = numeric(),
    S_chem = numeric(),
    retrieval_rank = integer(),
    precursor_error_ppm = numeric(),
    matched_peak_fraction = numeric(),
    n_match_pairs = integer(),
    n_matched_peaks = integer(),
    n_reverse_match = integer(),
    unique_query_matched_peaks = integer(),
    unique_library_matched_peaks = integer(),
    n_query_peaks_total = integer(),
    n_library_peaks_total = integer(),
    query_match_ratio = numeric(),
    library_match_ratio = numeric(),
    query_peak_reuse_n = integer(),
    query_peak_reuse_fraction = numeric(),
    has_query_peak_reuse = logical(),
    n_SAS_matched_peaks = integer(),
    sum_matched_library_intensity = numeric(),
    matched_query_intensity = numeric(),
    matched_query_intensity_fraction = numeric(),
    matched_library_intensity = numeric(),
    matched_library_intensity_fraction = numeric(),
    reverse_match_status = character(),
    match_mode = character(),
    score_ms1 = numeric(),
    score_ms2 = numeric(),
    score_spatial = numeric(),
    score_total_pre_mfc = numeric(),
    pass_n_matched_peaks_ge3 = logical(),
    pass_reverse_masked_cosine = logical(),
    evidence_pass = logical(),
    evidence_failure_reason = character(),
    confidence_tier = character(),
    s_fe = numeric(),
    s_arch = numeric(),
    s_rev = numeric(),
    n_match = integer(),
    u_spec = numeric(),
    s_chem = numeric(),
    s_spat = numeric(),
    S_FE_cal = numeric(),
    S_reverse_cos_cal = numeric(),
    S_lib_coverage_cal = numeric(),
    S_spatial_reverse_cal = numeric(),
    S_precursor_cal = numeric(),
    S_chem_cal = numeric(),
    s_fe_cal = numeric(),
    s_arch_cal = numeric(),
    s_rev_cal = numeric(),
    u_spec_cal = numeric(),
    s_chem_cal = numeric(),
    s_spat_cal = numeric(),
    annotation_confidence = character(),
    chem_status = character(),
    chem_role = character(),
    chem_explanation = character(),
    role_delta_mz = numeric(),
    role_ppm_error = numeric(),
    role_evidence_n = integer(),
    spectral_shape_support = character(),
    passes_spectral_shape_filter = logical(),
    reverse_filter_mode = character(),
    threshold_validation_status = character(),
    selection_reason = character(),
    candidate_source_for_step7 = character(),
    diagnostic_only = logical(),
    annotation_status = character(),
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
}

`.ann_empty_spatial_scores` <- 
function() {
  data.frame(
    cluster_id = integer(),
    anchor_feature_id = character(),
    anchor_mz = numeric(),
    candidate_key = character(),
    s_spat = numeric(),
    explained_spatial_peak_n = integer(),
    total_spatial_peak_n = integer(),
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
}

`.ann_ensure_evidence_pass` <- 
function(candidate_table, object) {
  candidate_table <- as.data.frame(candidate_table, stringsAsFactors = FALSE, check.names = FALSE)
  thresholds <- .ann_reverse_evidence_thresholds(object)
  if (!nrow(candidate_table)) {
    candidate_table$evidence_pass <- logical()
    candidate_table$min_matched_peaks_used <- integer()
    candidate_table$min_reverse_masked_cosine_used <- numeric()
    candidate_table$evidence_filter_reason <- character()
    return(candidate_table)
  }
  if ("evidence_pass" %in% names(candidate_table)) {
    evidence_pass <- candidate_table$evidence_pass
    if (!is.logical(evidence_pass)) {
      evidence_pass <- tolower(trimws(as.character(evidence_pass))) %in%
        c("true", "t", "1", "yes", "pass")
    }
    candidate_table$evidence_pass <- evidence_pass %in% TRUE
  } else {
    warning(
      "candidate table lacks evidence_pass; recomputing from module6_reverse_match thresholds for compatibility.",
      call. = FALSE
    )
    first_numeric <- function(cols, default = NA_real_) {
      out <- rep(default, nrow(candidate_table))
      for (col in cols) {
        if (!col %in% names(candidate_table)) next
        value <- suppressWarnings(as.numeric(candidate_table[[col]]))
        if (any(is.finite(value))) {
          return(value)
        }
      }
      out
    }
    n_matched <- suppressWarnings(as.integer(first_numeric(
      c("n_matched_peaks", "n_reverse_match", "n_match", "n_match_pairs"),
      default = 0
    )))
    n_matched[!is.finite(n_matched)] <- 0L
    rmc <- suppressWarnings(as.numeric(first_numeric(
      c("reverse_masked_cosine", "score_ms2", "S_reverse_cos"),
      default = 0
    )))
    rmc[!is.finite(rmc)] <- 0
    candidate_table$evidence_pass <- n_matched >= thresholds$min_matched_peaks &
      rmc > thresholds$min_reverse_masked_cosine
  }
  if (!"min_matched_peaks_used" %in% names(candidate_table)) {
    candidate_table$min_matched_peaks_used <- rep(thresholds$min_matched_peaks, nrow(candidate_table))
  }
  if (!"min_reverse_masked_cosine_used" %in% names(candidate_table)) {
    candidate_table$min_reverse_masked_cosine_used <- rep(thresholds$min_reverse_masked_cosine, nrow(candidate_table))
  }
  if (!"evidence_filter_reason" %in% names(candidate_table)) {
    candidate_table$evidence_filter_reason <- ifelse(candidate_table$evidence_pass, "pass", "below_step6_evidence_rule")
  }
  candidate_table
}

`.ann_first_numeric` <- 
function(df, cols) {
  out <- rep(NA_real_, nrow(df))
  for (col in cols) {
    if (!col %in% names(df)) next
    value <- suppressWarnings(as.numeric(df[[col]]))
    use <- is.na(out) & !is.na(value)
    out[use] <- value[use]
  }
  out
}

`.ann_fmc_requested` <- 
function(args) {
  if (!is.list(args) || !length(args)) return(FALSE)
  if (isTRUE(args$fmc_enabled)) return(TRUE)
  cfg <- args$fmc_annotation
  if (is.list(cfg) && isTRUE(cfg$enabled)) return(TRUE)
  non_enable <- args[setdiff(names(args), c("fmc_enabled", "fmc_annotation"))]
  any(vapply(non_enable, function(x) !is.null(x), logical(1)))
}

`.ann_match_query_reference` <- 
function(query, reference, ppm = 10,
                                       tolerance_da = NULL) {
  empty <- data.frame(
    query_peak_feature_id = character(),
    reference_peak_id = integer(),
    query_intensity = numeric(),
    reference_intensity = numeric(),
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
  query <- as.data.frame(query, stringsAsFactors = FALSE, check.names = FALSE)
  reference <- as.data.frame(reference, stringsAsFactors = FALSE, check.names = FALSE)
  if (nrow(query) == 0L || nrow(reference) == 0L) return(empty)
  query$peak_mz <- suppressWarnings(as.numeric(query$peak_mz))
  query$intensity <- suppressWarnings(as.numeric(query$intensity))
  reference$peak_mz <- suppressWarnings(as.numeric(reference$peak_mz))
  reference$intensity_sum1 <- suppressWarnings(as.numeric(reference$intensity_sum1))
  rows <- list()
  for (j in seq_len(nrow(reference))) {
    ref_mz <- reference$peak_mz[[j]]
    if (!is.finite(ref_mz) || ref_mz <= 0) next
    ppm_error <- (query$peak_mz - ref_mz) / ref_mz * 1e6
    ppm_tolerance_da <- ref_mz * suppressWarnings(as.numeric(ppm[[1]])) / 1e6
    absolute_tolerance_da <- if (is.null(tolerance_da)) {
      NA_real_
    } else {
      suppressWarnings(as.numeric(tolerance_da[[1]]))
    }
    effective_tolerance_da <- if (is.finite(ppm_tolerance_da) &&
        is.finite(absolute_tolerance_da)) {
      max(ppm_tolerance_da, absolute_tolerance_da)
    } else if (is.finite(absolute_tolerance_da)) {
      absolute_tolerance_da
    } else {
      ppm_tolerance_da
    }
    mz_error <- query$peak_mz - ref_mz
    hit_idx <- which(is.finite(mz_error) & is.finite(effective_tolerance_da) &
      abs(mz_error) <= effective_tolerance_da)
    if (length(hit_idx) == 0L) next
    candidates <- data.frame(
      idx = hit_idx,
      abs_error = abs(ppm_error[hit_idx]),
      query_intensity = query$intensity[hit_idx],
      stringsAsFactors = FALSE,
      check.names = FALSE
    )
    candidates$query_intensity[!is.finite(candidates$query_intensity)] <- 0
    candidates <- candidates[order(candidates$abs_error, -candidates$query_intensity), , drop = FALSE]
    i <- candidates$idx[[1]]
    rows[[length(rows) + 1L]] <- data.frame(
      query_peak_feature_id = as.character(query$peak_feature_id[[i]]),
      reference_peak_id = as.integer(j),
      query_intensity = as.numeric(query$intensity[[i]]),
      reference_intensity = as.numeric(reference$intensity_sum1[[j]]),
      stringsAsFactors = FALSE,
      check.names = FALSE
    )
  }
  if (length(rows) == 0L) return(empty)
  do.call(rbind, rows)
}

`.ann_mfc_input_required_columns` <- 
function() {
  c(
    "component", "feature_id", "anchor_feature_id", "mz", "metabolite_id",
    "hmdb_id", "metabolite_name", "adduct", "theoretical_mz", "mz_error_ppm",
    "match_mode", "evidence_source", "source", "lab_id", "n_matched_peaks",
    "matched_query_intensity_fraction", "matched_library_intensity_fraction",
    "library_l2_coverage", "reverse_masked_cosine", "S_SAS",
    "S_SAS_score_for_total", "S_SAS_missing", "S_total", "score_ms1",
    "score_ms2", "score_spatial", "score_total_pre_mfc",
    "confidence_tier", "mfc_entry_reason", "source_row_id",
    "annotation_pass_class", "pass_primary", "pass_review",
    "rule_fail_reason", "candidate_key", "dedup_key"
  )
}

`.ann_order_final_annotation` <- 
function(x) {
  out <- as.data.frame(x, stringsAsFactors = FALSE, check.names = FALSE)
  out <- .ann_add_missing_columns(out, c(
    "cluster_id", "anchor_feature_id", "annotation_confidence",
    "reverse_masked_cosine", "S_lib_coverage", "matched_peak_fraction", "n_reverse_match",
    "S_SAS", "precursor_error_ppm", "S_total", "s_final"
  ))
  confidence_rank <- .ann_confidence_rank(out$annotation_confidence)
  reverse_masked_cosine <- suppressWarnings(as.numeric(out$reverse_masked_cosine))
  lib_coverage <- suppressWarnings(as.numeric(out$S_lib_coverage))
  matched_fraction <- suppressWarnings(as.numeric(out$matched_peak_fraction))
  matched_n <- suppressWarnings(as.numeric(out$n_reverse_match))
  sas <- suppressWarnings(as.numeric(out$S_SAS))
  precursor_ppm <- abs(suppressWarnings(as.numeric(out$precursor_error_ppm)))
  final_score <- suppressWarnings(as.numeric(out[[.ann_primary_score_column(out)]]))
  reverse_masked_cosine[!is.finite(reverse_masked_cosine)] <- -Inf
  lib_coverage[!is.finite(lib_coverage)] <- -Inf
  matched_fraction[!is.finite(matched_fraction)] <- -Inf
  matched_n[!is.finite(matched_n)] <- -Inf
  sas[!is.finite(sas)] <- -Inf
  precursor_ppm[!is.finite(precursor_ppm)] <- Inf
  final_score[!is.finite(final_score)] <- -Inf
  out[order(
    out$cluster_id,
    out$anchor_feature_id,
    -final_score,
    -confidence_rank,
    -reverse_masked_cosine,
    -lib_coverage,
    -matched_n,
    -sas,
    -matched_fraction,
    precursor_ppm
  ), , drop = FALSE]
}

`.ann_primary_score_column` <- 
function(df) {
  if ("S_total" %in% names(df)) {
    value <- suppressWarnings(as.numeric(df$S_total))
    if (any(is.finite(value))) {
      return("S_total")
    }
  }
  "s_final"
}

`.ann_reference_peaks` <- 
function(candidate, object) {
  candidate <- as.data.frame(candidate, stringsAsFactors = FALSE, check.names = FALSE)
  empty <- data.frame(
    peak_mz = numeric(),
    intensity_sum1 = numeric(),
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
  source <- if ("candidate_source" %in% names(candidate)) as.character(candidate$candidate_source[[1]]) else "flash"
  if (source %in% c("flash", "both")) {
    peaks <- object@spectral_library$library_peak_table
    if (!is.null(peaks)) {
      peaks <- as.data.frame(peaks, stringsAsFactors = FALSE, check.names = FALSE)
      ref <- peaks[as.character(peaks$spectrum_id) == as.character(candidate$spectrum_id[[1]]), , drop = FALSE]
      if (nrow(ref) == 0L) {
        ref <- peaks[as.character(peaks$candidate_id) == as.character(candidate$candidate_id[[1]]), , drop = FALSE]
      }
      if (nrow(ref) > 0L) {
        return(data.frame(
          peak_mz = suppressWarnings(as.numeric(ref$peak_mz)),
          intensity_sum1 = suppressWarnings(as.numeric(ref$intensity_sum1)),
          stringsAsFactors = FALSE,
          check.names = FALSE
        ))
      }
    }
  }
  if (source %in% c("archive", "both") && length(object@archive_index) > 0L) {
    peaks <- object@archive_index$archive_peak_table
    if (!is.null(peaks)) {
      peaks <- as.data.frame(peaks, stringsAsFactors = FALSE, check.names = FALSE)
      ref <- peaks[as.character(peaks$archive_id) == as.character(candidate$archive_id[[1]]), , drop = FALSE]
      if (nrow(ref) == 0L) {
        ref <- peaks[as.character(peaks$spectrum_id) == as.character(candidate$spectrum_id[[1]]), , drop = FALSE]
      }
      if (nrow(ref) > 0L) {
        return(data.frame(
          peak_mz = suppressWarnings(as.numeric(ref$peak_mz)),
          intensity_sum1 = suppressWarnings(as.numeric(ref$intensity_sum1)),
          stringsAsFactors = FALSE,
          check.names = FALSE
        ))
      }
    }
  }
  empty
}

`.ann_reverse_cosine` <- 
function(query, reference, matches) {
  if (nrow(matches) == 0L || nrow(reference) == 0L) return(0)
  ref_int <- suppressWarnings(as.numeric(reference$intensity_sum1))
  q_int <- rep(0, length(ref_int))
  idx <- suppressWarnings(as.integer(matches$reference_peak_id))
  keep <- is.finite(idx) & idx >= 1L & idx <= length(q_int)
  q_int[idx[keep]] <- suppressWarnings(as.numeric(matches$query_intensity[keep]))
  numerator <- sum(ref_int * q_int, na.rm = TRUE)
  denom <- sqrt(sum(ref_int^2, na.rm = TRUE) * sum(q_int^2, na.rm = TRUE))
  if (!is.finite(denom) || denom <= 0) return(0)
  min(max(numerator / denom, 0), 1)
}

`.ann_reverse_evidence_thresholds` <- 
function(object) {
  info <- if (methods::is(object, "msi_dataset") &&
      is.list(object@process_info$module6_reverse_match)) {
    object@process_info$module6_reverse_match
  } else {
    list()
  }
  min_peaks <- suppressWarnings(as.integer(info$min_matched_peaks_evidence %||%
    info$min_matched_peaks_used %||% info$min_matched_peaks %||% NA_integer_))
  min_rmc <- suppressWarnings(as.numeric(info$min_reverse_masked_cosine_evidence %||%
    info$min_reverse_masked_cosine_used %||% info$min_reverse_masked_cosine %||% NA_real_))
  if (!is.finite(min_peaks) || min_peaks < 1L) {
    min_peaks <- 3L
  }
  if (!is.finite(min_rmc) || min_rmc < 0) {
    min_rmc <- 0
  }
  list(min_matched_peaks = min_peaks, min_reverse_masked_cosine = min_rmc)
}

`.ann_S_SAS_from_object` <- 
function(object, pairs = NULL) {
  if (methods::is(object, "msi_dataset")) {
    if (is.data.frame(pairs)) {
      pairs <- pairs
    } else if (is.data.frame(object@ms2_data$matched_peak_pairs)) {
      pairs <- object@ms2_data$matched_peak_pairs
    } else if (is.data.frame(object@annotation_results$matched_peak_pairs)) {
      pairs <- object@annotation_results$matched_peak_pairs
    } else if (is.data.frame(object@annotation_results$reverse_matched_peaks)) {
      pairs <- object@annotation_results$reverse_matched_peaks
    }
  }
  if (is.null(pairs)) {
    return(calculate_S_SAS(data.frame()))
  }
  calculate_S_SAS(pairs)
}

`.ann_split_key` <- 
function(df, cols) {
  df <- as.data.frame(df, stringsAsFactors = FALSE, check.names = FALSE)
  n <- nrow(df)
  if (length(cols) == 0L) {
    return(rep("", n))
  }
  parts <- lapply(cols, function(col) {
    value <- if (col %in% names(df)) as.character(df[[col]]) else rep(NA_character_, n)
    value[is.na(value)] <- "<NA>"
    paste0(nchar(value, type = "bytes"), ":", value)
  })
  do.call(paste, c(parts, sep = "|"))
}

`.ann_standardize_mfc_input_candidates` <- 
function(x) {
  out <- as.data.frame(x %||% data.frame(), stringsAsFactors = FALSE, check.names = FALSE)
  required <- .ann_mfc_input_required_columns()
  if (!nrow(out)) {
    return(.ann_empty_mfc_input_candidates())
  }
  if (!"source_row_id" %in% names(out)) {
    out$source_row_id <- paste0("step7_row_", seq_len(nrow(out)))
  } else {
    source_row_id <- as.character(out$source_row_id)
    missing_source <- is.na(source_row_id) | !nzchar(source_row_id)
    source_row_id[missing_source] <- paste0("step7_row_", which(missing_source))
    out$source_row_id <- source_row_id
  }
  for (nm in setdiff(required, names(out))) out[[nm]] <- rep(NA, nrow(out))
  forbidden <- c(
    paste0("reverse", "_dp"),
    paste0("rank_reverse", "_dp"),
    paste0("passes_reverse", "_dp_filter"),
    paste0("reverse", "_dp_weight"),
    paste0("old_pass_reverse", "_dp")
  )
  out <- out[, setdiff(names(out), forbidden), drop = FALSE]
  pseudo <- out$match_mode %in% "pseudo_ms2"
  if (!"evidence_type" %in% names(out)) out$evidence_type <- NA_character_
  missing_evidence_type <- is.na(out$evidence_type) | !nzchar(as.character(out$evidence_type))
  out$evidence_type[pseudo & missing_evidence_type] <- "pseudo_ms2_strict"
  out$evidence_source[pseudo] <- "pseudo_ms2_candidate_evidence"
  out$score_ms2[pseudo] <- suppressWarnings(as.numeric(out$reverse_masked_cosine[pseudo]))
  out$score_spatial[pseudo] <- suppressWarnings(as.numeric(out$S_SAS_score_for_total[pseudo]))
  out$score_total_pre_mfc[pseudo] <- suppressWarnings(as.numeric(out$S_total[pseudo]))
  out$mfc_entry_reason[pseudo] <- "pseudo_ms2_evidence_pass"
  if ("annotation_q_value" %in% names(out) && !"q_value" %in% names(out)) {
    out$q_value <- suppressWarnings(as.numeric(out$annotation_q_value))
  }
  out <- .mfc2_add_candidate_keys(out)
  out <- out[, unique(c(required, setdiff(names(out), required))), drop = FALSE]
  rownames(out) <- NULL
  out
}

`.ann_total_spatial_peak_counts` <- 
function(object) {
  empty <- data.frame(
    cluster_id = integer(),
    query_id = character(),
    anchor_feature_id = character(),
    total_spatial_peak_n = integer(),
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
  pseudo <- object@pseudo_ms2$pseudo_ms2_table
  if (is.null(pseudo)) {
    return(empty)
  }
  pseudo <- as.data.frame(pseudo, stringsAsFactors = FALSE, check.names = FALSE)
  if (nrow(pseudo) == 0L) {
    return(empty)
  }
  pseudo <- .ann_add_missing_columns(pseudo, c(
    "cluster_id", "anchor_feature_id", "fragment_feature_id", "peak_feature_id"
  ))
  feature_id <- as.character(pseudo$fragment_feature_id)
  fallback <- as.character(pseudo$peak_feature_id)
  missing_feature <- is.na(feature_id) | !nzchar(feature_id)
  feature_id[missing_feature] <- fallback[missing_feature]
  feature_id[is.na(feature_id)] <- ""

  key <- .ann_split_key(pseudo, c("cluster_id", "anchor_feature_id"))
  key_levels <- unique(key)
  first_idx <- match(key_levels, key)
  total <- vapply(key_levels, function(k) {
    values <- feature_id[key == k]
    values <- values[nzchar(values)]
    length(unique(values))
  }, integer(1))

  data.frame(
    cluster_id = pseudo$cluster_id[first_idx],
    anchor_feature_id = as.character(pseudo$anchor_feature_id[first_idx]),
    total_spatial_peak_n = total,
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
}

`.ann_within_anchor_percentile` <- 
function(df, col) {
  value <- suppressWarnings(as.numeric(df[[col]]))
  value[is.na(value)] <- 0
  out <- numeric(length(value))
  split_anchor <- .ann_split_key(df, c("cluster_id", "anchor_feature_id"))
  for (idx in split(seq_along(value), split_anchor)) {
    v <- value[idx]
    if (length(v) == 0L) next
    if (all(v == 0)) {
      out[idx] <- 0
    } else if (length(unique(v)) == 1L) {
      out[idx] <- 1
    } else {
      out[idx] <- (rank(v, ties.method = "average") - 1) / (length(v) - 1)
    }
  }
  pmin(pmax(out, 0), 1)
}

`.fmc_annotation_input_feature_table` <- 
function(object,
                                                pseudo_ms2_feature_ids = character(),
                                                step7_pseudo_feature_ids = character()) {
  features <- .fmc_feature_table(object)
  if (anyDuplicated(as.character(features$feature_id))) {
    stop("MFC annotation input feature IDs must be unique.", call. = FALSE)
  }
  feature_id <- as.character(features$feature_id)
  pseudo_from_table <- unique(as.character(pseudo_ms2_feature_ids))
  pseudo_from_table <- pseudo_from_table[!is.na(pseudo_from_table) & nzchar(pseudo_from_table)]
  pseudo_from_step7 <- unique(as.character(step7_pseudo_feature_ids))
  pseudo_from_step7 <- pseudo_from_step7[!is.na(pseudo_from_step7) & nzchar(pseudo_from_step7)]

  pseudo_from_metadata <- character()
  variable_info <- if (methods::is(object, "msi_dataset")) object@variable_info else data.frame()
  if (is.data.frame(variable_info) && nrow(variable_info) &&
      all(c("variable_id", "has_pseudo_ms2") %in% names(variable_info))) {
    flag <- variable_info$has_pseudo_ms2
    flag <- if (is.logical(flag)) {
      flag %in% TRUE
    } else {
      tolower(trimws(as.character(flag))) %in% c("true", "t", "1", "yes", "y")
    }
    pseudo_from_metadata <- as.character(variable_info$variable_id[flag])
  }

  has_table <- feature_id %in% pseudo_from_table
  has_metadata <- feature_id %in% pseudo_from_metadata
  has_step7 <- feature_id %in% pseudo_from_step7
  pseudo_source <- ifelse(
    has_table,
    "pseudo_ms2_table",
    ifelse(
      has_metadata,
      "variable_info_has_pseudo_ms2",
      ifelse(has_step7, "inferred_from_step7_pseudo_ms2_candidate", "none")
    )
  )
  data.frame(
    feature_id = feature_id,
    feature_mz = suppressWarnings(as.numeric(features$feature_mz)),
    annotation_input = TRUE,
    has_pseudo_ms2 = has_table | has_metadata | has_step7,
    pseudo_ms2_flag_source = pseudo_source,
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
}

`.fmc_apply_pure_tissue_feature_leiden_membership` <- 
function(candidates, membership) {
  candidates <- as.data.frame(candidates, stringsAsFactors = FALSE, check.names = FALSE)
  membership <- as.data.frame(membership, stringsAsFactors = FALSE, check.names = FALSE)
  if (!all(c("feature_id", "cluster_id", "cluster_index") %in% names(membership)) ||
      anyDuplicated(as.character(membership$feature_id)) ||
      anyNA(membership$cluster_id) || !is.character(membership$cluster_id)) {
    stop("PURE_TISSUE_WEIGHTED_LEIDEN_CLUSTER_ID_SCHEMA_FAILED", call. = FALSE)
  }
  idx <- match(as.character(candidates$feature_id), as.character(membership$feature_id))
  if (anyNA(idx)) {
    stop("PURE_TISSUE_WEIGHTED_LEIDEN_CANDIDATE_FEATURE_MAPPING_FAILED", call. = FALSE)
  }
  candidates$spatial_cluster_id <- as.character(membership$cluster_id[idx])
  candidates$cluster_index <- as.integer(membership$cluster_index[idx])
  candidates$cluster_id_source <- "weighted_leiden_feature_graph"
  candidates
}

`.fmc_as_df` <- 
function(x) {
  if (is.null(x)) {
    return(data.frame())
  }

  if (inherits(x, "data.table") || inherits(x, "tbl_df") || inherits(x, "tbl")) {
    return(as.data.frame(x, stringsAsFactors = FALSE))
  }

  if (is.data.frame(x)) {
    return(as.data.frame(x, stringsAsFactors = FALSE))
  }

  tryCatch(
    as.data.frame(x, stringsAsFactors = FALSE),
    error = function(e) {
      stop(
        "FMC helper .fmc_as_df() could not convert input to data.frame: ",
        conditionMessage(e),
        call. = FALSE
      )
    }
  )
}

`.fmc_build_pure_tissue_feature_leiden_route` <- 
function(object, config) {
  route_cfg <- .fmc_pure_tissue_feature_leiden_config(config)
  if (is.null(route_cfg)) return(NULL)
  graph_data <- object@feature_graph
  feature_ids <- unique(as.character(graph_data$feature_ids %||% character()))
  if (length(feature_ids) != 919L || anyNA(feature_ids) || any(!nzchar(feature_ids))) {
    stop("PURE_TISSUE_WEIGHTED_LEIDEN_NODE_CONTRACT_FAILED: expected 919 unique selector features.", call. = FALSE)
  }
  pair_stats <- graph_data$pair_stats
  required <- c("feature_i", "feature_j", "SAS", "valid_pair", "pass_joint_filter")
  if (!is.data.frame(pair_stats) || !all(required %in% names(pair_stats))) {
    stop("PURE_TISSUE_WEIGHTED_LEIDEN_PAIR_STATS_UNAVAILABLE", call. = FALSE)
  }
  keep <- !is.na(pair_stats$valid_pair) & pair_stats$valid_pair %in% TRUE &
    !is.na(pair_stats$pass_joint_filter) & pair_stats$pass_joint_filter %in% TRUE &
    is.finite(suppressWarnings(as.numeric(pair_stats$SAS))) &
    suppressWarnings(as.numeric(pair_stats$SAS)) >= route_cfg$edge_sas_threshold
  edge_tbl <- data.frame(
    feature_i = as.character(pair_stats$feature_i[keep]),
    feature_j = as.character(pair_stats$feature_j[keep]),
    weight = as.numeric(pair_stats$SAS[keep]),
    stringsAsFactors = FALSE
  )
  if (nrow(edge_tbl) != 631L || any(!edge_tbl$feature_i %in% feature_ids) ||
      any(!edge_tbl$feature_j %in% feature_ids)) {
    stop("PURE_TISSUE_WEIGHTED_LEIDEN_EDGE_CONTRACT_FAILED: expected 631 in-selector SAS >= 0.75 edges.", call. = FALSE)
  }
  route <- .msi_run_weighted_leiden_stable(
    feature_ids = feature_ids,
    edge_tbl = edge_tbl,
    resolution = route_cfg$resolution,
    restarts = route_cfg$restarts,
    seeds = route_cfg$seeds,
    stability_ari_threshold = route_cfg$stability_ari_threshold,
    prefix = "PT02_LEIDEN"
  )
  membership <- route$membership
  if (!inherits(route$graph, "igraph") || igraph::vcount(route$graph) != 919L ||
      igraph::ecount(route$graph) != 631L || nrow(membership) != 919L ||
      anyDuplicated(membership$feature_id) || anyNA(membership$cluster_id) ||
      !is.character(membership$cluster_id) || !any(membership$cluster_size > 1L)) {
    stop("PURE_TISSUE_WEIGHTED_LEIDEN_MEMBERSHIP_CONTRACT_FAILED", call. = FALSE)
  }
  route$config <- route_cfg
  route
}

`.fmc_chr` <- 
function(x, cols = NULL, default = NA_character_) {
  if (is.data.frame(x) || inherits(x, "data.table") || inherits(x, "tbl_df") || inherits(x, "tbl")) {
    if (is.null(cols)) {
      stop(".fmc_chr() received a data.frame but no candidate columns.", call. = FALSE)
    }
    return(.fmc_value_from_df(x, cols, default = default, mode = "character"))
  }

  if (is.null(x) || length(x) == 0L || all(is.na(x))) {
    return(default)
  }
  value <- as.character(x[[1]])
  if (length(value) == 0L || is.na(value) || !nzchar(value)) {
    return(default)
  }
  value
}

`.fmc_col` <- 
function(df, candidates) {
  if (is.null(df) || !is.data.frame(df)) {
    return(NA_character_)
  }

  nms <- names(df)
  if (is.null(nms) || length(nms) == 0L || length(candidates) == 0L) {
    return(NA_character_)
  }

  idx <- match(tolower(trimws(as.character(candidates))), tolower(trimws(nms)))
  idx <- idx[!is.na(idx)]
  if (length(idx) == 0L) {
    return(NA_character_)
  }
  nms[idx[[1]]]
}

`.fmc_collapse_unique` <- 
function(x, sep = ";") {
  if (is.null(x) || length(x) == 0L) {
    return(NA_character_)
  }
  x <- as.character(x)
  x <- unique(x[!is.na(x) & nzchar(x)])
  if (length(x) == 0L) {
    return(NA_character_)
  }
  paste(x, collapse = sep)
}

`.fmc_config_value` <- 
function(config, key, default = NULL) {
  if (is.null(config)) {
    return(default)
  }

  if (is.list(config) && !is.data.frame(config)) {
    if (!is.null(config[[key]])) {
      return(config[[key]])
    }
    return(default)
  }

  if (is.data.frame(config) || inherits(config, "data.table") ||
      inherits(config, "tbl_df") || inherits(config, "tbl")) {
    cfg <- .fmc_as_df(config)
    if (nrow(cfg) == 0L) {
      return(default)
    }

    key_col <- .fmc_col(cfg, c("parameter", "name", "key", "setting"))
    value_col <- .fmc_col(cfg, c("value", "setting_value", "parameter_value"))
    if (!is.na(key_col) && !is.na(value_col)) {
      hit <- which(tolower(as.character(cfg[[key_col]])) == tolower(as.character(key)))
      if (length(hit) > 0L) {
        return(cfg[[value_col]][hit[[1]]])
      }
    }

    wide_col <- .fmc_col(cfg, key)
    if (!is.na(wide_col)) {
      return(cfg[[wide_col]][[1]])
    }
  }

  default
}

`.fmc_empty` <- 
function(columns = character()) {
  columns <- unique(as.character(columns))
  if (length(columns) == 0L) {
    return(data.frame())
  }

  out <- stats::setNames(
    replicate(length(columns), logical(0), simplify = FALSE),
    columns
  )
  as.data.frame(out, stringsAsFactors = FALSE, check.names = FALSE)
}

`.fmc_evidence_level` <- 
function(types, adducts, rank_groups, score_spatial) {
  types <- unique(as.character(types))
  adducts <- unique(as.character(adducts[!is.na(adducts) & nzchar(adducts)]))
  rank_groups <- unique(as.character(rank_groups[!is.na(rank_groups) & nzchar(rank_groups)]))
  canonical <- any(rank_groups == "canonical") || any(adducts %in% c("[M+H]+", "[M-H]-"))
  has_spatial <- is.finite(score_spatial) && score_spatial > 0
  if (any(types %in% c("pseudo_ms2", "pseudo_ms2_strict"))) return("Level A")
  if (canonical && (("isotope" %in% types) || has_spatial)) return("Level B")
  if (canonical) return("Level C")
  if (any(types %in% c("isotope", "neutral_loss", "ms1_adduct"))) return("Level D")
  "Level E"
}

`.fmc_evidence_level_summary` <- 
function(candidate_table) {
  x <- .fmc_as_df(candidate_table)
  if (!nrow(x) || !"fmc_evidence_level" %in% names(x)) {
    return(.fmc_empty(c("fmc_evidence_level", "n_candidate_pairs", "n_features")))
  }
  levels <- sort(unique(as.character(x$fmc_evidence_level)))
  do.call(rbind, lapply(levels, function(level) {
    sub <- x[as.character(x$fmc_evidence_level) == level, , drop = FALSE]
    data.frame(
      fmc_evidence_level = level,
      n_candidate_pairs = nrow(sub),
      n_features = length(unique(sub$feature_id)),
      stringsAsFactors = FALSE,
      check.names = FALSE
    )
  }))
}

`.fmc_feature_summary_full` <- 
function(object, pseudo, ms1, isotope, neutral_loss,
                                      before_candidates, after_candidates) {
  features <- .fmc_feature_table(object)
  feature_id <- as.character(features$feature_id)
  before_counts <- if (nrow(before_candidates)) table(before_candidates$feature_id) else integer()
  after_counts <- if (nrow(after_candidates)) table(after_candidates$feature_id) else integer()
  out <- data.frame(
    feature_id = feature_id,
    feature_mz = features$feature_mz,
    has_pseudo_ms2_evidence = feature_id %in% as.character(pseudo$feature_id),
    has_ms1_adduct_candidate = feature_id %in% as.character(ms1$feature_id),
    has_isotope_candidate = feature_id %in% as.character(isotope$feature_id),
    has_neutral_loss_candidate = feature_id %in% as.character(neutral_loss$feature_id),
    n_fmc_candidates_before_redundancy = as.integer(before_counts[feature_id]),
    n_fmc_candidates_after_redundancy = as.integer(after_counts[feature_id]),
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
  out$n_fmc_candidates_before_redundancy[is.na(out$n_fmc_candidates_before_redundancy)] <- 0L
  out$n_fmc_candidates_after_redundancy[is.na(out$n_fmc_candidates_after_redundancy)] <- 0L
  out$has_any_fmc_candidate <- out$n_fmc_candidates_after_redundancy > 0L
  top <- after_candidates[after_candidates$rank_within_feature_after_redundancy == 1L, , drop = FALSE]
  if (nrow(top)) {
    top <- top[, c("feature_id", "metabolite_id", "metabolite_name",
                   "score_total_after_redundancy", "fmc_evidence_level"), drop = FALSE]
    names(top) <- c("feature_id", "top_metabolite_id", "top_metabolite_name",
                    "top_score_total_after_redundancy", "top_fmc_evidence_level")
    out <- merge(out, top, by = "feature_id", all.x = TRUE, sort = FALSE)
  }
  out
}

`.fmc_feature_table` <- 
function(object, parse_feature_id_mz = TRUE) {
  variable_info <- NULL

  if (methods::is(object, "msi_dataset")) {
    slot_names <- methods::slotNames(object)
    if ("variable_info" %in% slot_names) {
      variable_info <- methods::slot(object, "variable_info")
    } else if ("variableData" %in% slot_names) {
      variable_info <- methods::slot(object, "variableData")
    }
  }

  if (is.null(variable_info) && is.list(object) && !is.null(object$variable_info)) {
    variable_info <- object$variable_info
  }

  if (is.null(variable_info) && is.data.frame(object)) {
    variable_info <- object
  }

  if (is.null(variable_info)) {
    stop(
      "FMC helper .fmc_feature_table() could not find variable_info in object.",
      call. = FALSE
    )
  }

  variable_info <- .fmc_as_df(variable_info)
  if (nrow(variable_info) == 0L) {
    return(data.frame(
      feature_id = character(),
      feature_mz = numeric(),
      stringsAsFactors = FALSE
    ))
  }

  feature_col <- .fmc_col(variable_info, c("feature_id", "feature", "variable_id", "id", "name"))
  contract_present <- any(c("legacy_feature_mz", "feature_mz_exact", "feature_mz_calibrated") %in% names(variable_info))
  legacy_col <- .fmc_col(variable_info, c("legacy_feature_mz"))
  exact_col <- .fmc_col(variable_info, c("feature_mz_exact"))
  calibrated_col <- .fmc_col(variable_info, c("feature_mz_calibrated"))
  mz_col <- .fmc_col(variable_info, c("feature_mz", "mz", "m_z", "mass_to_charge", "mass"))

  if (!is.na(feature_col)) {
    feature_id <- as.character(variable_info[[feature_col]])
  } else {
    rn <- rownames(variable_info)
    if (!is.null(rn) && length(rn) == nrow(variable_info) &&
        !all(rn %in% as.character(seq_len(nrow(variable_info))))) {
      feature_id <- as.character(rn)
    } else {
      feature_id <- paste0("feature_", seq_len(nrow(variable_info)))
    }
  }

  legacy_feature_mz <- if (!is.na(legacy_col)) suppressWarnings(as.numeric(variable_info[[legacy_col]])) else rep(NA_real_, nrow(variable_info))
  feature_mz_exact <- if (!is.na(exact_col)) suppressWarnings(as.numeric(variable_info[[exact_col]])) else rep(NA_real_, nrow(variable_info))
  feature_mz_calibrated <- if (!is.na(calibrated_col)) suppressWarnings(as.numeric(variable_info[[calibrated_col]])) else rep(NA_real_, nrow(variable_info))
  if (isTRUE(contract_present)) {
    # The contract is row-wise: calibration may be unavailable, in which case
    # exact is the analytical fallback. Legacy m/z is never a matching input
    # once an exact/calibrated contract is present.
    feature_mz <- feature_mz_calibrated
    use_exact <- !is.finite(feature_mz)
    feature_mz[use_exact] <- feature_mz_exact[use_exact]
    if (any(!is.finite(feature_mz))) {
      stop("FMC feature m/z contract has no finite calibrated or exact m/z for every feature.", call. = FALSE)
    }
  } else if (!is.na(mz_col)) {
    feature_mz <- suppressWarnings(as.numeric(variable_info[[mz_col]]))
    legacy_feature_mz <- feature_mz
  } else {
    feature_mz <- if (isTRUE(parse_feature_id_mz)) {
      .fmc_parse_mz_from_feature_id(feature_id)
    } else {
      rep(NA_real_, nrow(variable_info))
    }
    legacy_feature_mz <- feature_mz
  }

  observed_col <- .fmc_col(variable_info, c("observed_mz"))
  observed_mz <- if (!is.na(observed_col)) {
    suppressWarnings(as.numeric(variable_info[[observed_col]]))
  } else if (!is.na(mz_col)) {
    suppressWarnings(as.numeric(variable_info[[mz_col]]))
  } else if (!is.na(legacy_col)) {
    suppressWarnings(as.numeric(variable_info[[legacy_col]]))
  } else {
    rep(NA_real_, nrow(variable_info))
  }
  observed_full_col <- .fmc_col(variable_info, c("observed_mz_full"))
  observed_mz_full <- if (!is.na(observed_full_col)) {
    suppressWarnings(as.numeric(variable_info[[observed_full_col]]))
  } else {
    rep(NA_real_, nrow(variable_info))
  }
  resolved_mz <- .msi_resolve_observed_mz(data.frame(
    observed_mz = observed_mz,
    observed_mz_full = observed_mz_full,
    stringsAsFactors = FALSE,
    check.names = FALSE
  ))
  if (!contract_present) {
    use_full <- is.finite(resolved_mz$observed_mz_full)
    feature_mz[use_full] <- resolved_mz$observed_mz_full[use_full]
  }

  if (all(is.na(feature_mz)) && isTRUE(parse_feature_id_mz)) {
    stop(
      "FMC helper .fmc_feature_table() could not infer feature m/z. ",
      "Expected a column such as feature_mz, mz, m_z, mass_to_charge, ",
      "or feature IDs like m787_8.",
      call. = FALSE
    )
  }

  pure_tissue <- FALSE
  pure_meta <- NULL
  if (methods::is(object, "msi_dataset") && is.list(object@process_info$pure_tissue_mfc)) {
    pure_meta <- object@process_info$pure_tissue_mfc
    pure_tissue <- identical(tolower(as.character(pure_meta$sample_mode %||% "")), "pure_tissue")
  }
  if (isTRUE(pure_tissue)) {
    source <- as.character(pure_meta$ms1_mz_source %||% "")
    if (!identical(source, "cardinal_measured_feature_axis")) {
      stop(
        "pure_tissue MS1 feature m/z provenance is unavailable; expected cardinal_measured_feature_axis.",
        call. = FALSE
      )
    }
    if (any(!is.finite(feature_mz) | feature_mz <= 0)) {
      stop(
        "pure_tissue requires finite measured variable_info$mz for every MS1 feature; no generic m/z fallback is allowed.",
        call. = FALSE
      )
    }
  }

  data.frame(
    feature_id = feature_id,
    legacy_feature_mz = legacy_feature_mz,
    feature_mz_exact = feature_mz_exact,
    feature_mz_calibrated = feature_mz_calibrated,
    feature_mz = feature_mz,
    observed_mz = observed_mz,
    observed_mz_full = observed_mz_full,
    ms1_feature_mz = resolved_mz$mz,
    ms1_mz_source = resolved_mz$source,
    ms1_mz_provenance = ifelse(
      resolved_mz$source == "observed_mz_full",
      "full_precision_feature_centroid",
      ifelse(resolved_mz$source == "observed_mz", "legacy_observed_mz_fallback", "m/z_unavailable")
    ),
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
}

`.fmc_first_nonmissing` <- 
function(x, default = NA) {
  if (is.null(x) || length(x) == 0L) {
    return(default)
  }

  if (is.list(x) && !is.data.frame(x)) {
    x <- unlist(x, recursive = FALSE, use.names = FALSE)
  }

  keep <- !is.na(x)
  if (is.character(x)) {
    keep <- keep & nzchar(x)
  }
  if (!any(keep)) {
    return(default)
  }
  x[which(keep)[[1]]]
}

`.fmc_get_feature_graph_edge` <- 
function(feature_graph, feature_id_a, feature_id_b) {
  if (is.null(feature_graph)) {
    return(NA_real_)
  }

  if (length(feature_id_a) == 0L || length(feature_id_b) == 0L) {
    return(NA_real_)
  }
  feature_id_a <- as.character(feature_id_a[[1]])
  feature_id_b <- as.character(feature_id_b[[1]])
  if (is.na(feature_id_a) || is.na(feature_id_b)) {
    return(NA_real_)
  }
  if (identical(feature_id_a, feature_id_b)) {
    return(1)
  }

  if (is.list(feature_graph) && !is.data.frame(feature_graph)) {
    if (is.data.frame(feature_graph$edges)) {
      feature_graph <- feature_graph$edges
    } else if (is.data.frame(feature_graph$pair_stats)) {
      feature_graph <- feature_graph$pair_stats
    } else if (!is.null(feature_graph$sas_matrix)) {
      feature_graph <- feature_graph$sas_matrix
    }
  }

  if (is.matrix(feature_graph) || methods::is(feature_graph, "Matrix")) {
    rn <- rownames(feature_graph)
    cn <- colnames(feature_graph)
    if (is.null(rn) || is.null(cn) ||
        !feature_id_a %in% rn || !feature_id_b %in% cn) {
      return(NA_real_)
    }
    value <- suppressWarnings(as.numeric(feature_graph[feature_id_a, feature_id_b]))
    if (length(value) == 0L || !is.finite(value)) {
      return(NA_real_)
    }
    return(value)
  }

  edge_table <- tryCatch(.fmc_as_df(feature_graph), error = function(e) data.frame())
  if (nrow(edge_table) == 0L) {
    return(NA_real_)
  }

  col_a <- .fmc_col(edge_table, c("feature_id_a", "from", "feature_a", "source", "feature1"))
  col_b <- .fmc_col(edge_table, c("feature_id_b", "to", "feature_b", "target", "feature2"))
  col_s <- .fmc_col(edge_table, c("SAS", "sas", "edge_weight", "weight", "score"))
  if (is.na(col_a) || is.na(col_b) || is.na(col_s)) {
    return(NA_real_)
  }

  a <- as.character(edge_table[[col_a]])
  b <- as.character(edge_table[[col_b]])
  hit <- (a == feature_id_a & b == feature_id_b) | (a == feature_id_b & b == feature_id_a)
  if (!any(hit)) {
    return(NA_real_)
  }

  value <- suppressWarnings(as.numeric(edge_table[[col_s]][which(hit)[[1]]]))
  if (length(value) == 0L || !is.finite(value)) {
    return(NA_real_)
  }
  value
}

`.fmc_infer_neutral_mass` <- 
function(precursor_mz, adduct, ion_mode = NA_character_,
                                    adduct_rules = default_adduct_rule_table()) {
  precursor_mz <- suppressWarnings(as.numeric(precursor_mz))
  if (!is.finite(precursor_mz) || is.na(adduct) || !nzchar(as.character(adduct))) {
    return(NA_real_)
  }
  rules <- adduct_rules[as.character(adduct_rules$evidence_type) == "ms1_adduct" &
                          as.character(adduct_rules$adduct) == as.character(adduct) &
                          adduct_rules$enabled %in% TRUE, , drop = FALSE]
  if (!is.na(ion_mode) && nzchar(as.character(ion_mode))) {
    rules <- rules[tolower(as.character(rules$ion_mode)) == tolower(as.character(ion_mode)), , drop = FALSE]
  }
  rules <- rules[is.finite(suppressWarnings(as.numeric(rules$mass_shift))) &
                   is.finite(suppressWarnings(as.numeric(rules$charge))) &
                   suppressWarnings(as.numeric(rules$charge)) != 0, , drop = FALSE]
  if (nrow(rules) != 1L) return(NA_real_)
  precursor_mz * abs(as.numeric(rules$charge[[1]])) - as.numeric(rules$mass_shift[[1]])
}

`.fmc_int` <- 
function(x, cols = NULL, default = NA_integer_) {
  if (is.data.frame(x) || inherits(x, "data.table") || inherits(x, "tbl_df") || inherits(x, "tbl")) {
    if (is.null(cols)) {
      stop(".fmc_int() received a data.frame but no candidate columns.", call. = FALSE)
    }
    return(.fmc_value_from_df(x, cols, default = default, mode = "integer"))
  }

  if (is.null(x) || length(x) == 0L || all(is.na(x))) {
    return(default)
  }
  value <- suppressWarnings(as.integer(x[[1]]))
  if (length(value) == 0L || is.na(value)) {
    return(default)
  }
  value
}

`.fmc_ledger_columns` <- 
function(x, branch = NA_character_) {
  x <- .fmc_as_df(x)
  if (!nrow(x)) return(.fmc_rejection_template())
  data.frame(
    source_row_id = as.character(x$source_row_id),
    candidate_key = as.character(x$candidate_key),
    branch = as.character(branch),
    last_reached_stage = as.character(x$last_reached_stage %||% NA_character_),
    rule_fail_reason = as.character(x$rule_fail_reason %||% NA_character_),
    join_fail_reason = as.character(x$join_fail_reason %||% NA_character_),
    removal_reason = as.character(x$removal_reason %||% NA_character_),
    dedup_winner_source_row_id = as.character(x$dedup_winner_source_row_id %||% NA_character_),
    dedup_winner_candidate_key = as.character(x$dedup_winner_candidate_key %||% NA_character_),
    dedup_winner_evidence_source = as.character(x$dedup_winner_evidence_source %||% NA_character_),
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
}

`.fmc_lgl` <- 
function(x, cols = NULL, default = FALSE) {
  if (is.data.frame(x) || inherits(x, "data.table") || inherits(x, "tbl_df") || inherits(x, "tbl")) {
    if (is.null(cols)) {
      stop(".fmc_lgl() received a data.frame but no candidate columns.", call. = FALSE)
    }
    return(.fmc_value_from_df(x, cols, default = default, mode = "logical"))
  }

  default_value <- FALSE
  if (is.logical(default)) {
    default_value <- isTRUE(default[[1]])
  } else if (is.numeric(default)) {
    default_value <- !is.na(default[[1]]) && default[[1]] != 0
  } else if (is.character(default)) {
    y <- tolower(trimws(default[[1]]))
    default_value <- y %in% c("true", "t", "yes", "y", "1")
  }

  if (is.null(x) || length(x) == 0L || all(is.na(x))) {
    return(default_value)
  }

  value <- x[[1]]
  if (is.logical(value)) {
    return(isTRUE(value))
  }
  if (is.numeric(value)) {
    return(!is.na(value) && value != 0)
  }
  if (is.character(value)) {
    y <- tolower(trimws(value))
    if (y %in% c("true", "t", "yes", "y", "1")) return(TRUE)
    if (y %in% c("false", "f", "no", "n", "0")) return(FALSE)
  }
  default_value
}

`.fmc_merge_config` <- 
function(config = NULL, overrides = list()) {
  out <- default_fmc_annotation_config(FALSE)
  if (is.list(config)) {
    for (nm in names(config)) out[[nm]] <- config[[nm]]
  }
  if (length(overrides)) {
    for (nm in names(overrides)) {
      if (!is.null(overrides[[nm]])) out[[nm]] <- overrides[[nm]]
    }
  }
  out$enabled <- .fmc_lgl(out$enabled, FALSE)
  out$enable_ms1_fallback <- .fmc_lgl(out$enable_ms1_fallback, TRUE)
  out$legacy_ms1_fallback <- .fmc_lgl(out$legacy_ms1_fallback, FALSE)
  if (isTRUE(out$legacy_ms1_fallback)) {
    stop("legacy_ms1_fallback is disabled; use load_msi_database_ms1_flat().", call. = FALSE)
  }
  out$ms1_flat_server_mode <- .fmc_lgl(out$ms1_flat_server_mode, FALSE)
  out$ms1_matching_scope <- tolower(trimws(as.character(out$ms1_matching_scope %||%
    if (isTRUE(out$combine_ms1_with_pseudo)) "all_filtered_features" else "pseudo_unresolved_only")))
  if (!out$ms1_matching_scope %in% c("all_filtered_features", "pseudo_unresolved_only", "no_pseudo_spectrum_only")) {
    stop("mfc.ms1_matching_scope must be all_filtered_features, pseudo_unresolved_only, or no_pseudo_spectrum_only.", call. = FALSE)
  }
  out$combine_ms1_with_pseudo <- .fmc_lgl(out$combine_ms1_with_pseudo, FALSE)
  if (identical(out$ms1_matching_scope, "all_filtered_features")) out$combine_ms1_with_pseudo <- TRUE
  out$use_isotope <- .fmc_lgl(out$use_isotope, TRUE)
  out$use_neutral_loss <- .fmc_lgl(out$use_neutral_loss, FALSE)
  out$use_spatial_support <- .fmc_lgl(out$use_spatial_support, TRUE)
  out$count_duplicate_evidence_type_once <- .fmc_lgl(out$count_duplicate_evidence_type_once, TRUE)
  out$ms1_ppm <- suppressWarnings(as.numeric(out$ms1_ppm[[1]]))
  out$ppm_sigma <- suppressWarnings(as.numeric(out$ppm_sigma[[1]]))
  out$max_candidates_per_feature <- suppressWarnings(as.integer(out$max_candidates_per_feature[[1]]))
  out$report_top_k <- suppressWarnings(as.integer((out$report_top_k %||%
    out$max_candidates_per_feature)[[1]]))
  out$max_metabolites_per_feature <- suppressWarnings(as.integer(out$max_metabolites_per_feature[[1]]))
  out$max_features_per_metabolite_soft <- suppressWarnings(as.integer(out$max_features_per_metabolite_soft[[1]]))
  out$max_features_per_metabolite_hard <- suppressWarnings(as.integer(out$max_features_per_metabolite_hard[[1]]))
  if (!is.finite(out$ms1_ppm) || out$ms1_ppm <= 0) out$ms1_ppm <- 10
  if (!is.finite(out$ppm_sigma) || out$ppm_sigma <= 0) out$ppm_sigma <- 5
  if (!is.finite(out$max_candidates_per_feature) || out$max_candidates_per_feature < 1L) out$max_candidates_per_feature <- 100L
  if (!is.finite(out$report_top_k) || out$report_top_k < 1L) out$report_top_k <- 5L
  if (!is.finite(out$max_metabolites_per_feature) || out$max_metabolites_per_feature < 1L) out$max_metabolites_per_feature <- 20L
  if (!is.finite(out$max_features_per_metabolite_soft) || out$max_features_per_metabolite_soft < 1L) out$max_features_per_metabolite_soft <- 30L
  if (!is.finite(out$max_features_per_metabolite_hard) || out$max_features_per_metabolite_hard < 1L) out$max_features_per_metabolite_hard <- 100L
  out
}

`.fmc_metabolite_summary_after` <- 
function(candidate_table) {
  x <- .fmc_as_df(candidate_table)
  if (!nrow(x)) {
    return(.fmc_empty(c("metabolite_id", "metabolite_name", "n_features",
                       "max_score_total_after_redundancy", "best_fmc_evidence_level")))
  }
  groups <- split(seq_len(nrow(x)), as.character(x$metabolite_id), drop = TRUE)
  out <- do.call(rbind, lapply(groups, function(idx) {
    sub <- x[idx, , drop = FALSE]
    best <- which.max(suppressWarnings(as.numeric(sub$score_total_after_redundancy)))
    data.frame(
      metabolite_id = as.character(sub$metabolite_id[[1]]),
      metabolite_name = as.character(sub$metabolite_name[[1]]),
      n_features = length(unique(sub$feature_id)),
      max_score_total_after_redundancy = max(sub$score_total_after_redundancy, na.rm = TRUE),
      best_fmc_evidence_level = as.character(sub$fmc_evidence_level[[best]]),
      stringsAsFactors = FALSE,
      check.names = FALSE
    )
  }))
  rownames(out) <- NULL
  out
}

`.fmc_min_abs` <- 
function(x, default = NA_real_) {
  x <- suppressWarnings(as.numeric(x))
  x <- x[is.finite(x)]
  if (length(x) == 0L) {
    return(default)
  }
  x[which.min(abs(x))]
}

`.fmc_num` <- 
function(x, cols = NULL, default = NA_real_) {
  if (is.data.frame(x) || inherits(x, "data.table") || inherits(x, "tbl_df") || inherits(x, "tbl")) {
    if (is.null(cols)) {
      stop(".fmc_num() received a data.frame but no candidate columns.", call. = FALSE)
    }
    return(.fmc_value_from_df(x, cols, default = default, mode = "numeric"))
  }

  if (is.null(x) || length(x) == 0L || all(is.na(x))) {
    return(default)
  }
  value <- suppressWarnings(as.numeric(x[[1]]))
  if (length(value) == 0L || is.na(value)) {
    return(default)
  }
  value
}

`.fmc_parse_mz_from_feature_id` <- 
function(feature_id) {
  x <- as.character(feature_id)
  out <- rep(NA_real_, length(x))

  hit_decimal <- grepl("^m[0-9]+[.][0-9]+$", x)
  out[hit_decimal] <- suppressWarnings(as.numeric(sub("^m", "", x[hit_decimal])))

  hit_underscore <- grepl("^m[0-9]+_[0-9]+$", x)
  out[hit_underscore] <- suppressWarnings(as.numeric(
    sub("^m([0-9]+)_([0-9]+)$", "\\1.\\2", x[hit_underscore])
  ))

  hit_integer <- grepl("^m[0-9]+$", x)
  out[hit_integer] <- suppressWarnings(as.numeric(sub("^m", "", x[hit_integer])))

  out
}

`.fmc_ppm_error` <- 
function(observed_mz, theoretical_mz) {
  observed_mz <- suppressWarnings(as.numeric(observed_mz))
  theoretical_mz <- suppressWarnings(as.numeric(theoretical_mz))
  out <- (observed_mz - theoretical_mz) / theoretical_mz * 1e6
  out[!is.finite(out)] <- NA_real_
  out
}

`.fmc_prepare_step7_module8_intake` <- 
function(candidates, object, branch = NA_character_) {
  x <- .ann_standardize_mfc_input_candidates(candidates)
  feature_axis <- .mfc2_feature_axis(object = object)
  x <- .mfc2_add_candidate_keys(x, feature_axis = feature_axis)
  if (!nrow(x)) {
    empty <- x
    empty$last_reached_stage <- character()
    empty$retained <- logical()
    empty$rule_fail_reason <- character()
    empty$join_fail_reason <- character()
    empty$removal_reason <- character()
    return(list(
      intake = empty,
      retained = empty,
      rejected = .fmc_rejection_template(),
      candidate_audit = empty,
      required_column_rejections = .fmc_rejection_template(),
      stage_summary = .fmc_stage_summary_template(),
      candidate_hash = .mfc2_candidate_hash(empty)
    ))
  }
  hmdb_missing <- is.na(x$hmdb_key) | !nzchar(as.character(x$hmdb_key))
  missing_match_mode <- is.na(x$match_mode) | !nzchar(as.character(x$match_mode))
  anchor_invalid <- is.na(x$feature_key) | !nzchar(as.character(x$feature_key))
  if (length(feature_axis)) {
    anchor_invalid <- anchor_invalid | !(as.character(x$feature_key) %in% feature_axis)
  }
  retained_flag <- !(hmdb_missing | anchor_invalid | missing_match_mode)
  x$rule_fail_reason <- as.character(x$rule_fail_reason)
  x$rule_fail_reason[hmdb_missing] <- "missing_hmdb_id"
  x$rule_fail_reason[missing_match_mode] <- "missing_match_mode"
  x$rule_fail_reason[retained_flag & (is.na(x$rule_fail_reason) | !nzchar(x$rule_fail_reason))] <- "pass"
  x$join_fail_reason <- NA_character_
  x$join_fail_reason[anchor_invalid] <- "invalid_or_unmapped_anchor_feature_id"
  x$join_fail_reason[missing_match_mode] <- "missing_match_mode"
  x$removal_reason <- NA_character_
  x$removal_reason[hmdb_missing] <- "removed_by_mfc_required_columns"
  x$removal_reason[anchor_invalid] <- "removed_by_mfc_join"
  x$removal_reason[missing_match_mode] <- "missing_match_mode"
  x$last_reached_stage <- ifelse(
    retained_flag, "module8_intake",
    ifelse(
      missing_match_mode, "module8_match_mode_validation",
      ifelse(anchor_invalid, "module8_join", "module8_required_column_validation")
    )
  )
  x$retained <- retained_flag
  retained <- x[retained_flag, , drop = FALSE]
  rejected_rows <- x[!retained_flag, , drop = FALSE]
  rejected <- .fmc_ledger_columns(rejected_rows, branch = branch)
  required <- rejected[rejected$removal_reason %in% "removed_by_mfc_required_columns", , drop = FALSE]
  intake_hash <- .mfc2_candidate_hash(x)
  retained_hash <- .mfc2_candidate_hash(retained)
  stage_summary <- rbind(
    data.frame(
      branch = as.character(branch), stage = "step7_authoritative_input",
      n_input = nrow(x), n_valid_key = sum(!hmdb_missing),
      n_joined = NA_integer_, n_unmatched = NA_integer_,
      n_retained = nrow(x), n_rejected = 0L, n_after_dedup = NA_integer_,
      n_final = NA_integer_, candidate_hash = intake_hash,
      stringsAsFactors = FALSE, check.names = FALSE
    ),
    data.frame(
      branch = as.character(branch), stage = "module8_intake",
      n_input = nrow(x), n_valid_key = sum(!hmdb_missing),
      n_joined = nrow(retained), n_unmatched = sum(anchor_invalid),
      n_retained = nrow(retained), n_rejected = nrow(rejected),
      n_after_dedup = NA_integer_, n_final = NA_integer_,
      candidate_hash = retained_hash,
      stringsAsFactors = FALSE, check.names = FALSE
    )
  )
  list(
    intake = x,
    retained = retained,
    rejected = rejected,
    candidate_audit = x,
    required_column_rejections = required,
    stage_summary = stage_summary,
    candidate_hash = intake_hash,
    retained_hash = retained_hash
  )
}

`.fmc_pure_tissue_candidate_contract` <- 
function(config) {
  contract <- config$pure_tissue_candidate_contract
  required <- c(
    "contract_version", "library_path", "library_sha256", "library_fingerprint_sha256",
    "stage5_cache_schema", "stage5_n_spectra", "stage5_n_unique_spectra",
    "stage5_n_duplicate_spectrum_rows", "stage5_n_peaks", "candidate_key_schema",
    "candidate_count", "candidate_set_sha256"
  )
  if (!is.list(contract) || length(setdiff(required, names(contract)))) {
    stop("PURE_TISSUE_CANDIDATE_CONTRACT_CONFIG_MISSING", call. = FALSE)
  }
  contract$contract_version <- as.character(contract$contract_version[[1L]])
  contract$library_path <- as.character(contract$library_path[[1L]])
  contract$library_sha256 <- tolower(as.character(contract$library_sha256[[1L]]))
  contract$library_fingerprint_sha256 <- tolower(as.character(contract$library_fingerprint_sha256[[1L]]))
  contract$stage5_cache_schema <- as.character(contract$stage5_cache_schema[[1L]])
  contract$candidate_key_schema <- as.character(contract$candidate_key_schema[[1L]])
  contract$candidate_set_sha256 <- tolower(as.character(contract$candidate_set_sha256[[1L]]))
  for (field in c("stage5_n_spectra", "stage5_n_unique_spectra",
                  "stage5_n_duplicate_spectrum_rows", "stage5_n_peaks", "candidate_count")) {
    contract[[field]] <- suppressWarnings(as.integer(contract[[field]][[1L]]))
  }
  hashes <- c(contract$library_sha256, contract$library_fingerprint_sha256,
              contract$candidate_set_sha256)
  if (!identical(contract$contract_version, "pure_tissue_library_candidate_v1") ||
      !identical(contract$candidate_key_schema, "feature_id_hmdb_id_v1") ||
      any(!grepl("^[0-9a-f]{64}$", hashes)) || !nzchar(contract$library_path) ||
      !identical(contract$stage5_cache_schema, "stage5_list_v1") ||
      any(!is.finite(unlist(contract[c("stage5_n_spectra", "stage5_n_unique_spectra",
        "stage5_n_duplicate_spectrum_rows", "stage5_n_peaks", "candidate_count")])) |
        unlist(contract[c("stage5_n_spectra", "stage5_n_unique_spectra",
          "stage5_n_duplicate_spectrum_rows", "stage5_n_peaks", "candidate_count")]) < 0L)) {
    stop("PURE_TISSUE_CANDIDATE_CONTRACT_CONFIG_INVALID", call. = FALSE)
  }
  contract
}

`.fmc_pure_tissue_candidate_key` <- 
function(candidates) {
  candidates <- as.data.frame(candidates, stringsAsFactors = FALSE, check.names = FALSE)
  required <- c("feature_id", "HMDB_ID")
  missing <- setdiff(required, names(candidates))
  if (length(missing)) {
    stop("PURE_TISSUE_CANDIDATE_CONTRACT_KEY_COLUMNS_MISSING: ",
      paste(missing, collapse = ", "), call. = FALSE)
  }
  feature_id <- trimws(as.character(candidates$feature_id))
  hmdb_id <- toupper(trimws(as.character(candidates$HMDB_ID)))
  if (anyNA(feature_id) || anyNA(hmdb_id) || any(!nzchar(feature_id)) || any(!nzchar(hmdb_id))) {
    stop("PURE_TISSUE_CANDIDATE_CONTRACT_KEY_INVALID", call. = FALSE)
  }
  paste(feature_id, hmdb_id, sep = "\r")
}

`.fmc_pure_tissue_candidate_set_sha256` <- 
function(candidates) {
  keys <- .fmc_pure_tissue_candidate_key(candidates)
  if (anyDuplicated(keys)) {
    stop("PURE_TISSUE_CANDIDATE_CONTRACT_DUPLICATE_KEY", call. = FALSE)
  }
  if (!requireNamespace("digest", quietly = TRUE)) {
    stop("PURE_TISSUE_CANDIDATE_CONTRACT_DIGEST_UNAVAILABLE", call. = FALSE)
  }
  digest::digest(paste(sort(keys), collapse = "\n"), algo = "sha256", serialize = FALSE)
}

`.fmc_pure_tissue_feature_leiden_config` <- 
function(config) {
  route <- config$feature_cluster
  if (is.null(route)) return(NULL)
  if (!is.list(route)) {
    stop("pure_tissue mfc.feature_cluster must be a mapping.", call. = FALSE)
  }
  algorithm <- tolower(trimws(as.character(route$algorithm %||% "")))
  threshold <- suppressWarnings(as.numeric(route$edge_sas_threshold %||% NA_real_))
  resolution <- suppressWarnings(as.numeric(route$resolution %||% NA_real_))
  restarts <- suppressWarnings(as.integer(route$restarts %||% NA_integer_))
  seeds <- suppressWarnings(as.integer(unlist(route$seeds %||% integer(), use.names = FALSE)))
  ari_threshold <- suppressWarnings(as.numeric(route$stability_ari_threshold %||% NA_real_))
  if (!identical(algorithm, "weighted_leiden") ||
      !isTRUE(all.equal(threshold, 0.75)) ||
      !isTRUE(all.equal(resolution, 1.0)) ||
      !identical(restarts, 20L) ||
      !identical(seeds, seq_len(20L)) ||
      !isTRUE(all.equal(ari_threshold, 0.90))) {
    stop(
      paste0(
        "pure_tissue weighted-Leiden configuration must be exactly algorithm=weighted_leiden, ",
        "edge_sas_threshold=0.75, resolution=1.0, restarts=20, seeds=1:20, ",
        "stability_ari_threshold=0.90."
      ),
      call. = FALSE
    )
  }
  list(
    algorithm = algorithm,
    edge_sas_threshold = threshold,
    resolution = resolution,
    restarts = restarts,
    seeds = seeds,
    stability_ari_threshold = ari_threshold
  )
}

`.fmc_rank` <- 
function(x, group_col, score_col, decreasing = TRUE) {
  x <- .fmc_as_df(x)
  if (nrow(x) == 0L) {
    return(integer(0))
  }
  if (!group_col %in% names(x)) {
    stop("FMC helper .fmc_rank() missing group column: ", group_col, call. = FALSE)
  }
  if (!score_col %in% names(x)) {
    stop("FMC helper .fmc_rank() missing score column: ", score_col, call. = FALSE)
  }

  group <- as.character(x[[group_col]])
  score <- suppressWarnings(as.numeric(x[[score_col]]))
  score[is.na(score)] <- if (isTRUE(decreasing)) -Inf else Inf

  ranks <- integer(nrow(x))
  for (g in unique(group)) {
    idx <- which(group == g)
    ord <- idx[order(score[idx], decreasing = decreasing, na.last = TRUE)]
    ranks[ord] <- seq_along(ord)
  }
  ranks
}

`.fmc_rejection_template` <- 
function() {
  data.frame(
    source_row_id = character(), candidate_key = character(), branch = character(),
    last_reached_stage = character(), rule_fail_reason = character(),
    join_fail_reason = character(), removal_reason = character(),
    dedup_winner_source_row_id = character(),
    dedup_winner_candidate_key = character(),
    dedup_winner_evidence_source = character(),
    stringsAsFactors = FALSE, check.names = FALSE
  )
}

`.fmc_score_config_with_runtime` <- 
function(score_config, fmc_config) {
  score_config <- .fmc_as_df(score_config)
  extra <- data.frame(
    parameter = c("ms1_ppm", "max_candidates_per_feature", "report_top_k",
                  "max_metabolites_per_feature", "use_isotope",
                  "use_neutral_loss", "use_spatial_support"),
    value = as.character(c(fmc_config$ms1_ppm, fmc_config$max_candidates_per_feature,
                           fmc_config$report_top_k,
                           fmc_config$max_metabolites_per_feature, fmc_config$use_isotope,
                           fmc_config$use_neutral_loss, fmc_config$use_spatial_support)),
    description = "Runtime FMC annotation configuration.",
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
  rbind(score_config, extra)
}

`.fmc_spatial_sas` <- 
function(object, feature_id_a, feature_id_b) {
  if (is.null(object)) {
    return(NA_real_)
  }

  feature_graph <- NULL
  if (methods::is(object, "msi_dataset")) {
    slot_names <- methods::slotNames(object)
    if ("feature_graph" %in% slot_names) {
      feature_graph <- methods::slot(object, "feature_graph")
    }
  }

  if (is.null(feature_graph) && is.list(object) && !is.null(object$feature_graph)) {
    feature_graph <- object$feature_graph
  }

  if (is.null(feature_graph)) {
    return(NA_real_)
  }

  .fmc_get_feature_graph_edge(feature_graph, feature_id_a, feature_id_b)
}

`.fmc_stage_summary_template` <- 
function() {
  data.frame(
    branch = character(), stage = character(), n_input = integer(),
    n_valid_key = integer(), n_joined = integer(), n_unmatched = integer(),
    n_retained = integer(), n_rejected = integer(), n_after_dedup = integer(),
    n_final = integer(), candidate_hash = character(),
    stringsAsFactors = FALSE, check.names = FALSE
  )
}

`.fmc_standardize_evidence` <- 
function(evidence_tables) {
  if (is.null(evidence_tables)) evidence_tables <- list()
  if (is.data.frame(evidence_tables)) evidence_tables <- list(evidence_tables)
  rows <- list()
  evidence_names <- names(evidence_tables)
  if (is.null(evidence_names)) evidence_names <- rep("", length(evidence_tables))
  for (k in seq_along(evidence_tables)) {
    nm <- evidence_names[[k]]
    if (!nzchar(nm)) nm <- paste0("evidence_", k)
    x <- .fmc_as_df(evidence_tables[[k]])
    if (!nrow(x)) next
    if (!"feature_id" %in% names(x)) next
    if (!"metabolite_id" %in% names(x)) {
      if ("candidate_id" %in% names(x)) x$metabolite_id <- as.character(x$candidate_id) else next
    }
    if (!"metabolite_name" %in% names(x)) {
      x$metabolite_name <- if ("candidate_name" %in% names(x)) as.character(x$candidate_name) else NA_character_
    }
    if (!"feature_mz" %in% names(x)) x$feature_mz <- NA_real_
    if (!"evidence_type" %in% names(x)) x$evidence_type <- as.character(nm)
    if (!"evidence_weight" %in% names(x)) x$evidence_weight <- 0
    if (!"adduct" %in% names(x)) x$adduct <- NA_character_
    if (!"rank_group" %in% names(x)) x$rank_group <- NA_character_
    if (!"mz_error_ppm" %in% names(x)) x$mz_error_ppm <- NA_real_
    if (!"spatial_sas_to_parent" %in% names(x)) x$spatial_sas_to_parent <- NA_real_
    if (!"data_source" %in% names(x)) x$data_source <- as.character(nm)
    x$feature_id <- as.character(x$feature_id)
    x$metabolite_id <- as.character(x$metabolite_id)
    x$evidence_weight <- suppressWarnings(as.numeric(x$evidence_weight))
    rows[[length(rows) + 1L]] <- x
  }
  if (!length(rows)) {
    return(.fmc_empty(c("feature_id", "feature_mz", "metabolite_id", "metabolite_name",
                       "evidence_type", "evidence_weight", "adduct", "rank_group",
                       "mz_error_ppm", "spatial_sas_to_parent", "data_source")))
  }
  all_cols <- unique(unlist(lapply(rows, names), use.names = FALSE))
  rows <- lapply(rows, function(x) {
    missing <- setdiff(all_cols, names(x))
    for (col in missing) x[[col]] <- NA
    x[, all_cols, drop = FALSE]
  })
  out <- do.call(rbind, rows)
  out <- out[!is.na(out$feature_id) & nzchar(out$feature_id) &
               !is.na(out$metabolite_id) & nzchar(out$metabolite_id), , drop = FALSE]
  rownames(out) <- NULL
  out
}

`.fmc_validate_pure_tissue_candidate_contract` <- 
function(candidates, config) {
  contract <- .fmc_pure_tissue_candidate_contract(config)
  candidate_hash <- .fmc_pure_tissue_candidate_set_sha256(candidates)
  if (nrow(candidates) != contract$candidate_count ||
      !identical(candidate_hash, contract$candidate_set_sha256)) {
    stop(
      "PURE_TISSUE_WEIGHTED_LEIDEN_CANDIDATE_CONTRACT_FAILED: candidate universe does not match the configured library contract.",
      call. = FALSE
    )
  }
  c(contract, list(observed_candidate_count = nrow(candidates),
    observed_candidate_set_sha256 = candidate_hash))
}

`.fmc_value_from_df` <- 
function(df,
                               cols,
                               default = NA,
                               mode = c("character", "numeric", "integer", "logical")) {
  mode <- match.arg(mode)
  df <- .fmc_as_df(df)
  n <- nrow(df)

  if (mode == "character") {
    default_value <- if (length(default) == 0L) NA_character_ else as.character(default[[1]])
    fallback <- rep(default_value, n)
  } else if (mode == "numeric") {
    default_value <- if (length(default) == 0L) NA_real_ else suppressWarnings(as.numeric(default[[1]]))
    fallback <- rep(default_value, n)
  } else if (mode == "integer") {
    default_value <- if (length(default) == 0L) NA_integer_ else suppressWarnings(as.integer(default[[1]]))
    fallback <- rep(default_value, n)
  } else {
    if (length(default) == 0L || is.na(default[[1]])) {
      default_value <- NA
    } else if (is.logical(default)) {
      default_value <- isTRUE(default[[1]])
    } else if (is.numeric(default)) {
      default_value <- !is.na(default[[1]]) && default[[1]] != 0
    } else {
      default_value <- tolower(trimws(as.character(default[[1]]))) %in% c("true", "t", "yes", "y", "1")
    }
    fallback <- rep(default_value, n)
  }

  if (n == 0L) {
    return(fallback)
  }

  col <- .fmc_col(df, cols)
  if (is.na(col)) {
    return(fallback)
  }

  value <- df[[col]]
  if (mode == "character") {
    value <- as.character(value)
    missing <- is.na(value) | !nzchar(value)
    value[missing] <- fallback[missing]
    return(value)
  }

  if (mode == "numeric") {
    value <- suppressWarnings(as.numeric(value))
    value[is.na(value)] <- fallback[is.na(value)]
    return(value)
  }

  if (mode == "integer") {
    value <- suppressWarnings(as.integer(value))
    value[is.na(value)] <- fallback[is.na(value)]
    return(value)
  }

  out <- fallback
  if (is.logical(value)) {
    out <- value
    out[is.na(out)] <- fallback[is.na(out)]
    return(out)
  }
  if (is.numeric(value)) {
    out <- !is.na(value) & value != 0
    out[is.na(value)] <- fallback[is.na(value)]
    return(out)
  }
  value <- tolower(trimws(as.character(value)))
  out[value %in% c("true", "t", "yes", "y", "1")] <- TRUE
  out[value %in% c("false", "f", "no", "n", "0")] <- FALSE
  out
}

`.mfc_normalize_hmdb_id` <- 
function(x) {
  out <- toupper(trimws(as.character(x)))
  hit <- grepl("^HMDB[0-9]+$", out)
  digits <- sub("^HMDB", "", out[hit])
  short <- nchar(digits) < 7L
  digits[short] <- paste0(strrep("0", 7L - nchar(digits[short])), digits[short])
  out[hit] <- paste0("HMDB", digits)
  out
}

`.mfc_unified_clamp01` <- 
function(x) {
  x <- suppressWarnings(as.numeric(x))
  x[!is.finite(x)] <- 0
  pmin(1, pmax(0, x))
}

`.mfc_unified_evidence_key` <- 
function(x) {
  paste(
    as.character(x$feature_id),
    as.character(x$hmdb_id),
    as.character(x$library_spectrum_id),
    as.character(x$adduct),
    sep = "\r"
  )
}

`.mfc_unified_pseudo_score` <- 
function(x) {
  rmc <- .mfc_unified_clamp01(x$reverse_masked_cosine)
  matched <- .mfc_unified_clamp01(suppressWarnings(as.numeric(x$n_matched_peaks)) / 6)
  lib_fraction <- .mfc_unified_clamp01(x$matched_library_intensity_fraction)
  sas <- .mfc_unified_clamp01(x$S_SAS)
  0.55 * rmc + 0.20 * matched + 0.15 * lib_fraction + 0.10 * sas
}

`.mfc_unified_rejection_reason` <- 
function(x,
                                          min_matched_peaks = 3L,
                                          min_reverse_masked_cosine = 0) {
  n <- nrow(x)
  out <- rep(NA_character_, n)
  n_match <- suppressWarnings(as.integer(x$n_matched_peaks))
  rmc <- suppressWarnings(as.numeric(x$reverse_masked_cosine))
  out[!is.finite(n_match) | n_match < min_matched_peaks] <- "min_matched_peaks"
  fail_rmc <- is.na(out) & (!is.finite(rmc) | rmc <= min_reverse_masked_cosine)
  out[fail_rmc] <- "min_reverse_masked_cosine"
  out[is.na(out)] <- "not_in_passed_pseudo_evidence"
  out
}

`.mfc_unified_value` <- 
function(config, name, default) {
  value <- if (is.list(config)) config[[name]] else NULL
  if (is.null(value) || !length(value)) default else value[[1]]
}

`.mfc2_add_candidate_keys` <- 
function(x, feature_axis = character()) {
  out <- .mfc2_df(x)
  if (!nrow(out)) {
    for (nm in c("source_row_id", "feature_key", "hmdb_key", "adduct_key",
                 "library_spectrum_key", "candidate_key", "dedup_key")) {
      if (!nm %in% names(out)) out[[nm]] <- character()
    }
    return(out)
  }
  if (!"source_row_id" %in% names(out)) out$source_row_id <- NA_character_
  row_index <- seq_len(nrow(out))
  source_row_id <- .mfc2_chr_key(out$source_row_id)
  missing_source <- is.na(source_row_id)
  source_row_id[missing_source] <- paste0("mfc_source_row_", row_index[missing_source])
  out$source_row_id <- source_row_id

  if (!"anchor_feature_id" %in% names(out)) out$anchor_feature_id <- NA_character_
  if (!"feature_id" %in% names(out)) out$feature_id <- out$anchor_feature_id
  feature_key <- .mfc2_resolve_feature_key(out$anchor_feature_id, feature_axis)
  missing_feature_key <- is.na(feature_key)
  feature_id_clean <- .mfc2_chr_key(out$feature_id)
  feature_key[missing_feature_key] <- .mfc2_resolve_feature_key(
    feature_id_clean[missing_feature_key], feature_axis
  )
  if (!length(feature_axis)) {
    feature_key[is.na(feature_key)] <- .mfc2_chr_key(out$anchor_feature_id[is.na(feature_key)])
    feature_key[is.na(feature_key)] <- feature_id_clean[is.na(feature_key)]
  }
  out$feature_key <- feature_key
  mapped <- !is.na(feature_key) & nzchar(feature_key)
  out$feature_id[mapped] <- feature_key[mapped]
  out$anchor_feature_id[mapped] <- feature_key[mapped]

  hmdb_source <- if ("hmdb_id" %in% names(out)) out$hmdb_id else rep(NA_character_, nrow(out))
  missing_hmdb_source <- is.na(.mfc2_chr_key(hmdb_source))
  if (any(missing_hmdb_source) && "metabolite_id" %in% names(out)) {
    metabolite_source <- .mfc2_chr_key(out$metabolite_id)
    use <- missing_hmdb_source & !is.na(metabolite_source)
    hmdb_source[use] <- metabolite_source[use]
    missing_hmdb_source <- is.na(.mfc2_chr_key(hmdb_source))
  }
  if (any(missing_hmdb_source) && "candidate_id" %in% names(out)) {
    candidate_source <- .mfc2_chr_key(out$candidate_id)
    use <- missing_hmdb_source & !is.na(candidate_source)
    hmdb_source[use] <- candidate_source[use]
    missing_hmdb_source <- is.na(.mfc2_chr_key(hmdb_source))
  }
  if (any(missing_hmdb_source) && "compound_id" %in% names(out)) {
    compound_source <- .mfc2_chr_key(out$compound_id)
    use <- missing_hmdb_source & !is.na(compound_source)
    hmdb_source[use] <- compound_source[use]
  }
  out$hmdb_key <- .mfc2_clean_key(hmdb_source)
  if (!"hmdb_id" %in% names(out)) out$hmdb_id <- out$hmdb_key
  if (!"metabolite_id" %in% names(out)) out$metabolite_id <- out$hmdb_key
  out$hmdb_id[is.na(.mfc2_chr_key(out$hmdb_id))] <- out$hmdb_key[is.na(.mfc2_chr_key(out$hmdb_id))]
  out$metabolite_id[is.na(.mfc2_chr_key(out$metabolite_id))] <- out$hmdb_key[is.na(.mfc2_chr_key(out$metabolite_id))]

  if (!"adduct" %in% names(out)) out$adduct <- NA_character_
  if (!"library_spectrum_id" %in% names(out)) out$library_spectrum_id <- NA_character_
  if (!"match_mode" %in% names(out)) out$match_mode <- NA_character_
  out$adduct_key <- .mfc2_key_part(out$adduct, "adduct", out$source_row_id, row_index)
  out$library_spectrum_key <- .mfc2_key_part(
    out$library_spectrum_id, "library_spectrum_id", out$source_row_id, row_index
  )
  match_key <- .mfc2_key_part(out$match_mode, "match_mode", out$source_row_id, row_index)

  existing_candidate_key <- if ("candidate_key" %in% names(out)) {
    .mfc2_chr_key(out$candidate_key)
  } else {
    rep(NA_character_, nrow(out))
  }
  generated_candidate_key <- paste(
    .mfc2_key_part(out$feature_key, "feature_key", out$source_row_id, row_index),
    .mfc2_key_part(out$hmdb_key, "hmdb_key", out$source_row_id, row_index, uppercase = TRUE),
    out$adduct_key,
    out$library_spectrum_key,
    out$source_row_id,
    sep = "|"
  )
  existing_candidate_key[is.na(existing_candidate_key)] <-
    generated_candidate_key[is.na(existing_candidate_key)]
  out$candidate_key <- existing_candidate_key
  out$dedup_key <- paste(
    .mfc2_key_part(out$feature_key, "feature_key", out$source_row_id, row_index),
    .mfc2_key_part(out$hmdb_key, "hmdb_key", out$source_row_id, row_index, uppercase = TRUE),
    out$adduct_key,
    out$library_spectrum_key,
    match_key,
    sep = "|"
  )
  out
}

`.mfc2_candidate_hash` <- 
function(x) {
  x <- .mfc2_df(x)
  if (!nrow(x)) return("empty")
  cols <- intersect(.mfc2_candidate_identity_columns(), names(x))
  if (!length(cols)) return(paste0("rows:", nrow(x)))
  payload <- x[, cols, drop = FALSE]
  payload[] <- lapply(payload, function(col) {
    col <- as.character(col)
    col[is.na(col)] <- "<NA>"
    col
  })
  ord_cols <- intersect(c("source_row_id", "candidate_key", "dedup_key"), names(payload))
  if (length(ord_cols)) {
    ord <- do.call(order, c(payload[ord_cols], list(na.last = TRUE)))
    payload <- payload[ord, , drop = FALSE]
  }
  if (requireNamespace("digest", quietly = TRUE)) {
    digest::digest(payload, algo = "sha256")
  } else {
    paste0("rows:", nrow(payload), "|", paste(utils::capture.output(utils::str(payload)), collapse = "\n"))
  }
}

`.mfc2_candidate_identity_columns` <- 
function() {
  c(
    "source_row_id", "anchor_feature_id", "feature_key", "hmdb_key",
    "adduct_key", "library_spectrum_key", "match_mode", "candidate_key",
    "dedup_key"
  )
}

`.mfc2_chr` <- 
function(x, candidates, default = NA_character_) {
  col <- .mfc2_first_col(x, candidates)
  if (is.na(col)) return(rep(default, nrow(x)))
  out <- as.character(x[[col]])
  out[is.na(out) | !nzchar(out)] <- default
  out
}

`.mfc2_chr_key` <- 
function(x) {
  out <- trimws(as.character(x))
  out[is.na(out) | !nzchar(out)] <- NA_character_
  out
}

`.mfc2_clean_key` <- 
function(x) {
  out <- toupper(trimws(as.character(x)))
  out[is.na(out) | !nzchar(out)] <- NA_character_
  out
}

`.mfc2_cluster_medoid` <- 
function(features, spatial_similarity) {
  features <- unique(as.character(features))
  if (length(features) <= 1L) return(features[[1]])
  sim <- .mfc2_df(spatial_similarity)
  if (!nrow(sim)) return(sort(features)[[1]])
  score <- stats::setNames(rep(0, length(features)), features)
  count <- stats::setNames(rep(0L, length(features)), features)
  for (i in seq_len(nrow(sim))) {
    a <- as.character(sim$feature_id_a[[i]])
    b <- as.character(sim$feature_id_b[[i]])
    value <- suppressWarnings(as.numeric(sim$spatial_similarity[[i]]))
    if (!a %in% features || !b %in% features || !is.finite(value)) next
    score[[a]] <- score[[a]] + value
    score[[b]] <- score[[b]] + value
    count[[a]] <- count[[a]] + 1L
    count[[b]] <- count[[b]] + 1L
  }
  mean_score <- ifelse(count > 0L, score / count, 0)
  sort(names(mean_score)[mean_score == max(mean_score)], method = "radix")[[1]]
}

`.mfc2_cluster_similarity` <- 
function(features, spatial_similarity) {
  if (length(features) < 2L || !is.data.frame(spatial_similarity) || !nrow(spatial_similarity)) return(1)
  key <- paste(pmin(spatial_similarity$feature_id_a, spatial_similarity$feature_id_b),
               pmax(spatial_similarity$feature_id_a, spatial_similarity$feature_id_b), sep = "\r")
  pairs <- utils::combn(sort(features), 2)
  query <- paste(pairs[1, ], pairs[2, ], sep = "\r")
  vals <- spatial_similarity$spatial_similarity[match(query, key)]
  vals <- vals[is.finite(vals)]
  if (length(vals)) stats::median(vals) else 0
}

`.mfc2_component_id` <- 
function(component_id, candidates, analysis_domain = NULL) {
  if (!is.null(component_id) && length(component_id) &&
      !is.na(component_id[[1]]) && nzchar(as.character(component_id[[1]]))) {
    return(list(cluster_label = as.character(component_id[[1]]), component_value = as.character(component_id[[1]])))
  }
  if (identical(tolower(as.character(analysis_domain %||% "")), "all_acquired_pixels")) {
    return(list(cluster_label = "all_acquired_pixels", component_value = NA_character_))
  }
  list(cluster_label = "component1", component_value = "component1")
}

`.mfc2_components` <- 
function(nodes, edges) {
  parent <- seq_along(nodes)
  names(parent) <- nodes
  find <- function(x) {
    i <- match(x, nodes)
    while (parent[[i]] != i) {
      parent[[i]] <<- parent[[parent[[i]]]]
      i <- parent[[i]]
    }
    i
  }
  unite <- function(a, b) {
    ra <- find(a); rb <- find(b)
    if (ra != rb) parent[[rb]] <<- ra
  }
  if (nrow(edges)) {
    for (i in seq_len(nrow(edges))) unite(edges$feature_id_a[[i]], edges$feature_id_b[[i]])
  }
  roots <- vapply(nodes, find, integer(1))
  as.integer(factor(roots, levels = unique(roots)))
}

`.mfc2_df` <- 
function(x) {
  if (is.null(x)) return(data.frame())
  as.data.frame(x, stringsAsFactors = FALSE, check.names = FALSE)
}

`.mfc2_empty_candidate` <- 
function() {
  data.frame(
    component = character(), feature_id = character(),
    anchor_feature_id = character(), mz = numeric(),
    metabolite_id = character(), hmdb_id = character(),
    metabolite_name = character(), compound_name = character(),
    adduct = character(), theoretical_mz = numeric(),
    precursor_mz = numeric(), mz_error_ppm = numeric(),
    match_mode = character(), evidence_source = character(), source = character(),
    lab_id = character(), n_matched_peaks = integer(),
    matched_query_intensity_fraction = numeric(),
    matched_library_intensity_fraction = numeric(),
    library_l2_coverage = numeric(),
    reverse_masked_cosine = numeric(), S_SAS = numeric(),
    S_SAS_score_for_total = numeric(), S_SAS_missing = logical(),
    S_total = numeric(),
    score_ms1 = numeric(), score_ms2 = numeric(), score_spatial = numeric(),
    score_total_pre_mfc = numeric(), confidence_tier = character(),
    mfc_entry_reason = character(),
    feature_mz = numeric(), observed_mz = numeric(),
    observed_mz_full = numeric(), ms1_feature_mz = numeric(),
    ms1_mz_source = character(), ms1_mz_provenance = character(),
    evidence_type = character(), entropy_score = numeric(),
    q_value = numeric(), library_spectrum_id = character(),
    candidate_rank = integer(), raw_score = numeric(),
    evidence_level = character(), raw_evidence_score = numeric(),
    evidence_class_weight = numeric(), final_evidence_score = numeric(),
    source_row_id = character(), feature_key = character(),
    hmdb_key = character(), adduct_key = character(),
    library_spectrum_key = character(), candidate_key = character(),
    dedup_key = character(),
    stringsAsFactors = FALSE, check.names = FALSE
  )
}

`.mfc2_exact_pair_metrics` <- 
function(expression_data, feature_id_a, feature_id_b) {
  empty <- list(AIC = NA_real_, CAR = NA_real_, SAS = NA_real_, n11 = 0L,
                mab = 0L, valid_pair = FALSE)
  ids <- rownames(expression_data)
  if (is.null(ids)) return(empty)
  index <- match(c(feature_id_a, feature_id_b), ids)
  if (anyNA(index)) return(empty)
  a <- suppressWarnings(as.numeric(expression_data[index[[1]], , drop = TRUE]))
  b <- suppressWarnings(as.numeric(expression_data[index[[2]], , drop = TRUE]))
  a[!is.finite(a) | a <= 0] <- 0
  b[!is.finite(b) | b <= 0] <- 0
  active <- a > 0 | b > 0
  n11 <- sum(a > 0 & b > 0)
  mab <- sum(active)
  car <- if (mab > 0L) n11 / mab else NA_real_
  if (mab < 2L || n11 < 1L) {
    return(list(AIC = NA_real_, CAR = car, SAS = NA_real_, n11 = n11,
                mab = mab, valid_pair = FALSE))
  }
  active_a <- a[active]
  active_b <- b[active]
  aic <- if (stats::var(active_a) > 0 && stats::var(active_b) > 0) {
    suppressWarnings(stats::cor(active_a, active_b, method = "spearman"))
  } else {
    NA_real_
  }
  valid <- is.finite(aic) && is.finite(car)
  list(
    AIC = if (valid) as.numeric(aic) else NA_real_,
    CAR = as.numeric(car),
    SAS = if (valid) 0.5 * as.numeric(aic) + 0.5 * as.numeric(car) else NA_real_,
    n11 = as.integer(n11),
    mab = as.integer(mab),
    valid_pair = valid
  )
}

`.mfc2_expression_matrix` <- 
function(expression_data, feature_ids = NULL, pixel_subset = NULL) {
  expr <- if (methods::is(expression_data, "Matrix")) as.matrix(expression_data) else as.matrix(expression_data)
  if (is.null(rownames(expr))) rownames(expr) <- paste0("feature_", seq_len(nrow(expr)))
  if (!is.null(feature_ids)) {
    feature_ids <- as.character(feature_ids)
    idx <- match(feature_ids, rownames(expr))
    if (anyNA(idx)) {
      numeric_idx <- suppressWarnings(as.integer(feature_ids))
      if (all(is.finite(numeric_idx) & numeric_idx >= 1L & numeric_idx <= nrow(expr))) idx <- numeric_idx
    }
    idx <- idx[!is.na(idx)]
    expr <- expr[idx, , drop = FALSE]
  }
  if (!is.null(pixel_subset)) {
    pix <- if (is.logical(pixel_subset)) which(pixel_subset) else suppressWarnings(as.integer(pixel_subset))
    pix <- pix[is.finite(pix) & pix >= 1L & pix <= ncol(expr)]
    expr <- expr[, pix, drop = FALSE]
  }
  expr
}

`.mfc2_feature_axis` <- 
function(object = NULL, feature_metadata = NULL,
                               expression_matrix = NULL) {
  values <- character()
  meta <- .mfc2_df(feature_metadata)
  if (!nrow(meta) && methods::is(object, "msi_dataset")) {
    meta <- .mfc2_df(object@variable_info)
  }
  if (nrow(meta)) {
    for (col in c("feature_id", "variable_id", "id")) {
      if (col %in% names(meta)) values <- c(values, as.character(meta[[col]]))
    }
  }
  expr <- expression_matrix
  if (is.null(expr) && methods::is(object, "msi_dataset")) {
    expr <- object@expression_data
  }
  rn <- tryCatch(rownames(expr), error = function(e) NULL)
  if (!is.null(rn)) values <- c(values, as.character(rn))
  values <- values[!is.na(values) & nzchar(values)]
  unique(values)
}

`.mfc2_first_col` <- 
function(x, candidates) {
  if (!is.data.frame(x)) return(NA_character_)
  hit <- match(tolower(candidates), tolower(names(x)))
  hit <- hit[!is.na(hit)]
  if (length(hit)) names(x)[hit[[1]]] else NA_character_
}

`.mfc2_graph_pair_lookup` <- 
function(feature_graph) {
  if (!is.list(feature_graph)) return(data.frame())
  pair <- feature_graph$pair_stats
  if (!is.data.frame(pair) || !nrow(pair)) pair <- feature_graph$edges
  pair <- .mfc2_df(pair)
  if (!nrow(pair) || !all(c("feature_id_a", "feature_id_b") %in% names(pair))) return(data.frame())
  pair$feature_id_a <- as.character(pair$feature_id_a)
  pair$feature_id_b <- as.character(pair$feature_id_b)
  pair$pair_key <- paste(pmin(pair$feature_id_a, pair$feature_id_b),
                         pmax(pair$feature_id_a, pair$feature_id_b), sep = "\r")
  pair <- pair[!duplicated(pair$pair_key), , drop = FALSE]
  pair
}

`.mfc2_hmdb` <- 
function(...) {
  text <- paste(unlist(list(...), use.names = FALSE), collapse = " ")
  hit <- regmatches(text, regexpr("HMDB[0-9]{5,}", text, perl = TRUE, ignore.case = TRUE))
  if (!length(hit) || is.na(hit) || !nzchar(hit)) NA_character_ else toupper(hit)
}

`.mfc2_hotspot` <- 
function(x, top_fraction = 0.10) {
  x <- suppressWarnings(as.numeric(x))
  x[!is.finite(x)] <- 0
  if (!any(x > 0)) return(rep(FALSE, length(x)))
  cutoff <- suppressWarnings(as.numeric(stats::quantile(x[x > 0], probs = max(0, min(1, 1 - top_fraction)),
                                                        na.rm = TRUE, names = FALSE, type = 8)))
  x >= cutoff & x > 0
}

`.mfc2_int` <- 
function(x, candidates, default = NA_integer_) {
  out <- suppressWarnings(as.integer(.mfc2_num(x, candidates, default)))
  out[!is.finite(out)] <- default
  out
}

`.mfc2_key_part` <- 
function(value, field, source_row_id, row_index,
                           uppercase = FALSE) {
  out <- if (uppercase) .mfc2_clean_key(value) else .mfc2_chr_key(value)
  missing <- is.na(out)
  if (any(missing)) {
    out[missing] <- .mfc2_missing_token(field, source_row_id[missing], row_index[missing])
  }
  out
}

`.mfc2_missing_token` <- 
function(field, source_row_id, row_index) {
  sid <- .mfc2_chr_key(source_row_id)
  sid[is.na(sid)] <- paste0("row", row_index[is.na(sid)])
  paste0("__missing_", field, "__", sid)
}

`.mfc2_num` <- 
function(x, candidates, default = NA_real_) {
  col <- .mfc2_first_col(x, candidates)
  if (is.na(col)) return(rep(default, nrow(x)))
  out <- suppressWarnings(as.numeric(x[[col]]))
  out[!is.finite(out)] <- default
  out
}

`.mfc2_preprocess_expression` <- 
function(expr) {
  expr[!is.finite(expr) | expr < 0] <- 0
  tic <- colSums(expr, na.rm = TRUE)
  finite_tic <- tic[is.finite(tic) & tic > 0]
  scale_to <- if (length(finite_tic)) stats::median(finite_tic) else 1
  keep <- is.finite(tic) & tic > 0
  expr <- expr[, keep, drop = FALSE]
  tic <- tic[keep]
  expr <- sweep(expr, 2, tic, "/") * scale_to
  expr <- log1p(expr)
  for (i in seq_len(nrow(expr))) {
    v <- expr[i, ]
    upper <- suppressWarnings(as.numeric(stats::quantile(v, 0.99, na.rm = TRUE, names = FALSE, type = 8)))
    if (is.finite(upper)) v[v > upper] <- upper
    expr[i, ] <- v
  }
  nonzero_pixels <- colSums(expr > 0, na.rm = TRUE) > 0
  expr[, nonzero_pixels, drop = FALSE]
}

`.mfc2_resolve_feature_key` <- 
function(anchor_feature_id, feature_axis) {
  raw <- trimws(as.character(anchor_feature_id))
  out <- rep(NA_character_, length(raw))
  feature_axis <- unique(as.character(feature_axis))
  exact <- !is.na(raw) & nzchar(raw) & raw %in% feature_axis
  out[exact] <- raw[exact]
  unresolved <- which(is.na(out) & !is.na(raw) & nzchar(raw))
  if (length(unresolved) && length(feature_axis)) {
    numeric_idx <- suppressWarnings(as.integer(raw[unresolved]))
    ok <- is.finite(numeric_idx) & numeric_idx >= 1L & numeric_idx <= length(feature_axis)
    if (any(ok)) out[unresolved[ok]] <- feature_axis[numeric_idx[ok]]
    unresolved <- unresolved[!ok]
  }
  if (length(unresolved) && length(feature_axis)) {
    feature_idx <- suppressWarnings(as.integer(sub("^feature_([0-9]+)$", "\\1", raw[unresolved])))
    ok <- grepl("^feature_[0-9]+$", raw[unresolved]) &
      is.finite(feature_idx) & feature_idx >= 1L & feature_idx <= length(feature_axis)
    if (any(ok)) out[unresolved[ok]] <- feature_axis[feature_idx[ok]]
  }
  out
}

`.mfc2_scale01` <- 
function(x) {
  x <- suppressWarnings(as.numeric(x))
  x[!is.finite(x)] <- NA_real_
  high <- is.finite(x) & x > 1 & x <= 100
  x[high] <- x[high] / 100
  pmin(1, pmax(0, x))
}

`.mfc2_standardize_source` <- 
function(x, source = "annotation") {
  x <- .mfc2_df(x)
  empty <- .mfc2_empty_candidate()
  if (!nrow(x)) return(empty)
  metabolite_id <- .mfc2_chr(x, c("metabolite_id", "hmdb_id", "candidate_id", "compound_id"))
  hmdb_id <- .mfc2_chr(x, c("hmdb_id", "metabolite_id", "candidate_id", "compound_id"))
  hmdb_id <- vapply(seq_along(hmdb_id), function(i) .mfc2_hmdb(hmdb_id[[i]], metabolite_id[[i]]), character(1))
  missing_hmdb <- is.na(hmdb_id) | !nzchar(hmdb_id)
  hmdb_id[missing_hmdb] <- metabolite_id[missing_hmdb]
  out <- data.frame(
    component = .mfc2_chr(x, c("component", "component_id")),
    metabolite_id = metabolite_id,
    hmdb_id = hmdb_id,
    metabolite_name = .mfc2_chr(x, c("metabolite_name", "compound_name", "candidate_name", "name")),
    compound_name = .mfc2_chr(x, c("compound_name", "candidate_name", "metabolite_name", "name")),
    feature_id = .mfc2_chr(x, c("feature_id", "anchor_feature_id", "variable_id")),
    anchor_feature_id = .mfc2_chr(x, c("anchor_feature_id", "feature_id", "variable_id")),
    feature_mz = .mfc2_num(x, c("feature_mz", "anchor_mz", "observed_mz", "mz")),
    observed_mz = .mfc2_num(x, c("observed_mz")),
    observed_mz_full = .mfc2_num(x, c("observed_mz_full")),
    ms1_feature_mz = .mfc2_num(x, c("ms1_feature_mz")),
    ms1_mz_source = .mfc2_chr(x, c("ms1_mz_source")),
    ms1_mz_provenance = .mfc2_chr(x, c("ms1_mz_provenance")),
    mz = .mfc2_num(x, c("mz", "feature_mz", "anchor_mz", "observed_mz")),
    adduct = .mfc2_chr(x, c("adduct", "all_adducts", "precursor_type")),
    evidence_type = .mfc2_chr(x, c("evidence_type", "best_evidence_type", "candidate_source"), source),
    match_mode = .mfc2_chr(x, c("match_mode", "evidence_type", "best_evidence_type"), source),
    evidence_source = .mfc2_chr(x, c("evidence_source", "data_source", "candidate_source"), source),
    source = .mfc2_chr(x, c("source", "database_source", "library_source", "data_source"), source),
    lab_id = .mfc2_chr(x, c("lab_id", "spectrum_id", "library_spectrum_id", "best_spectrum_id", "archive_id")),
    precursor_mz = .mfc2_num(x, c("precursor_mz", "theoretical_mz", "library_precursor_mz")),
    theoretical_mz = .mfc2_num(x, c("theoretical_mz", "precursor_mz", "library_precursor_mz")),
    mz_error_ppm = .mfc2_num(x, c("mz_error_ppm", "precursor_error_ppm", "ppm_error")),
    entropy_score = .mfc2_num(x, c("entropy_score", "entropy_S", "S_FE", "s_fe"), 0),
    n_matched_peaks = .mfc2_int(x, c("n_matched_peaks", "n_reverse_match", "n_match", "n_match_pairs"), 0L),
    matched_query_intensity_fraction = .mfc2_num(x, c("matched_query_intensity_fraction"), NA_real_),
    matched_library_intensity_fraction = .mfc2_num(x, c("matched_library_intensity_fraction", "S_lib_coverage", "u_spec"), NA_real_),
    library_l2_coverage = .mfc2_num(x, c("library_l2_coverage"), NA_real_),
    reverse_masked_cosine = .mfc2_num(x, c("reverse_masked_cosine", "score_ms2", "S_reverse_cos"), NA_real_),
    S_SAS = .mfc2_num(x, c("S_SAS", "S_spatial_reverse", "spatial_sas_to_parent", "score_spatial"), NA_real_),
    S_SAS_score_for_total = .mfc2_num(x, c("S_SAS_score_for_total", "score_spatial", "S_SAS", "S_spatial_reverse"), NA_real_),
    S_SAS_missing = .mfc2_chr(x, c("S_SAS_missing"), NA_character_),
    S_total = .mfc2_num(x, c("S_total", "score_total_pre_mfc", "score_total", "final_score"), NA_real_),
    score_ms1 = .mfc2_num(x, c("score_ms1", "MS1_score", "score_mz_accuracy"), NA_real_),
    score_ms2 = .mfc2_num(x, c("score_ms2", "reverse_masked_cosine", "S_reverse_cos"), NA_real_),
    score_spatial = .mfc2_num(x, c("score_spatial", "S_SAS", "S_spatial_reverse", "spatial_sas_to_parent"), NA_real_),
    score_total_pre_mfc = .mfc2_num(x, c("score_total_pre_mfc", "S_total", "score_total", "final_score"), NA_real_),
    confidence_tier = .mfc2_chr(x, c("confidence_tier", "annotation_confidence")),
    mfc_entry_reason = .mfc2_chr(x, c("mfc_entry_reason", "evidence_failure_reason", "selection_reason"), source),
    q_value = .mfc2_num(x, c("q_value", "annotation_q_value", "qvalue")),
    library_spectrum_id = .mfc2_chr(x, c("library_spectrum_id", "spectrum_id", "best_spectrum_id", "archive_id")),
    candidate_rank = .mfc2_int(x, c("candidate_rank", "rank_within_feature", "rank_total", "rank_entropy", "retrieval_rank")),
    raw_score = .mfc2_num(x, c("raw_score", "final_evidence_score", "score_total", "score_total_after_redundancy", "final_score", "annotation_score"), NA_real_),
    evidence_level = .mfc2_chr(x, c("evidence_level", "confidence_tier", "fmc_evidence_level")),
    raw_evidence_score = .mfc2_num(x, c("raw_evidence_score", "ms1_only_score", "pseudo_ms2_score", "score_total_pre_mfc"), NA_real_),
    evidence_class_weight = .mfc2_num(x, c("evidence_class_weight"), NA_real_),
    final_evidence_score = .mfc2_num(x, c("final_evidence_score", "raw_score", "score_total_pre_mfc"), NA_real_),
    source_row_id = .mfc2_chr(x, c("source_row_id")),
    feature_key = .mfc2_chr(x, c("feature_key")),
    hmdb_key = .mfc2_chr(x, c("hmdb_key")),
    adduct_key = .mfc2_chr(x, c("adduct_key")),
    library_spectrum_key = .mfc2_chr(x, c("library_spectrum_key")),
    candidate_key = .mfc2_chr(x, c("candidate_key")),
    dedup_key = .mfc2_chr(x, c("dedup_key")),
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
  out$match_mode <- ifelse(
    grepl("pseudo", out$match_mode, ignore.case = TRUE),
    "pseudo_ms2",
    ifelse(grepl("ms1|adduct", out$match_mode, ignore.case = TRUE), "ms1_only", out$match_mode)
  )
  pseudo <- out$match_mode %in% "pseudo_ms2"
  out$score_ms2[pseudo] <- .mfc2_scale01(out$reverse_masked_cosine[pseudo])
  out$score_ms2[!pseudo] <- NA_real_
  sas_missing_text <- trimws(as.character(out$S_SAS_missing))
  has_sas_missing_text <- !is.na(sas_missing_text) & nzchar(sas_missing_text)
  out$S_SAS_missing <- ifelse(
    has_sas_missing_text,
    tolower(sas_missing_text) %in% c("true", "t", "1", "yes"),
    !is.finite(out$S_SAS)
  )
  missing_sas_score <- !is.finite(out$S_SAS_score_for_total)
  out$S_SAS_score_for_total[missing_sas_score] <- ifelse(
    is.finite(out$S_SAS[missing_sas_score]),
    .mfc2_scale01(out$S_SAS[missing_sas_score]),
    0
  )
  out$score_spatial[pseudo] <- out$S_SAS_score_for_total[pseudo]
  out$score_spatial[pseudo & !is.finite(out$score_spatial)] <- 0
  ms1_idx <- which(!pseudo & !is.finite(out$score_ms1))
  if (length(ms1_idx)) {
    out$score_ms1[ms1_idx] <- pmax(0, 1 - abs(out$mz_error_ppm[ms1_idx]) / 10)
  }
  pseudo_missing_total <- which(pseudo & !is.finite(out$score_total_pre_mfc))
  if (length(pseudo_missing_total)) {
    out$score_total_pre_mfc[pseudo_missing_total] <-
      0.7 * .mfc2_scale01(out$reverse_masked_cosine[pseudo_missing_total]) +
      0.3 * out$S_SAS_score_for_total[pseudo_missing_total]
  }
  ms1_missing_total <- which(!pseudo & !is.finite(out$score_total_pre_mfc))
  if (length(ms1_missing_total)) {
    out$score_total_pre_mfc[ms1_missing_total] <-
      .mfc2_scale01(out$score_ms1[ms1_missing_total]) * 0.5
  }
  out$S_total[pseudo & !is.finite(out$S_total)] <- out$score_total_pre_mfc[pseudo]
  missing_score <- !is.finite(out$raw_score)
  out$raw_score[missing_score] <- .mfc2_scale01(out$score_total_pre_mfc[missing_score])
  out$raw_score[!is.finite(out$raw_score)] <- 0
  missing_raw_evidence <- !is.finite(out$raw_evidence_score)
  out$raw_evidence_score[missing_raw_evidence] <- out$raw_score[missing_raw_evidence]
  missing_weight <- !is.finite(out$evidence_class_weight)
  out$evidence_class_weight[missing_weight & out$match_mode == "pseudo_ms2"] <- 1.0
  out$evidence_class_weight[missing_weight & out$match_mode == "ms1_only"] <- 0.6
  out$evidence_class_weight[!is.finite(out$evidence_class_weight)] <- 1
  missing_final <- !is.finite(out$final_evidence_score)
  out$final_evidence_score[missing_final] <-
    out$raw_evidence_score[missing_final] * out$evidence_class_weight[missing_final]
  out$final_evidence_score[!is.finite(out$final_evidence_score)] <- out$raw_score[!is.finite(out$final_evidence_score)]
  out$evidence_level[!nzchar(as.character(out$evidence_level)) | is.na(out$evidence_level)] <- ifelse(
    out$match_mode %in% "pseudo_ms2", "pseudo_ms2", "ms1_only"
  )
  out <- .mfc2_add_candidate_keys(out)
  out <- out[!is.na(out$feature_id) & nzchar(out$feature_id) &
               !is.na(out$metabolite_id) & nzchar(out$metabolite_id), , drop = FALSE]
  out <- out[, names(empty), drop = FALSE]
  rownames(out) <- NULL
  out
}

`.ms1_allowed_adducts` <- 
function(allowed_adducts = "auto_by_polarity", polarity = "negative") {
  if (is.null(allowed_adducts) || identical(allowed_adducts, "auto_by_polarity")) {
    polarity <- tolower(as.character(polarity[[1]]))
    if (polarity %in% c("positive", "pos", "+")) {
      return(c("[M+H]+", "[M+Na]+", "[M+K]+", "[M+NH4]+"))
    }
    return(c("[M-H]-", "[M+Cl]-", "[M+FA-H]-"))
  }
  as.character(unlist(allowed_adducts, use.names = FALSE))
}

`.ms1_column` <- 
function(x, candidates, default = NA_character_) {
  col <- .ms1_first_column(x, candidates)
  if (is.na(col)) return(rep(default, nrow(x)))
  out <- x[[col]]
  out
}

`.ms1_first_column` <- 
function(x, candidates) {
  if (!is.data.frame(x)) return(NA_character_)
  normalize <- function(value) gsub("[^a-z0-9]+", "", tolower(as.character(value)))
  hit <- match(normalize(candidates), normalize(names(x)))
  hit <- hit[!is.na(hit)]
  if (length(hit)) names(x)[hit[[1]]] else NA_character_
}

`.ms1_flat_resolve_path` <- 
function(path, project_root = getwd(), must_work = FALSE) {
  if (is.null(path) || length(path) != 1L || is.na(path[[1]]) ||
      !nzchar(trimws(as.character(path[[1]])))) {
    stop("MS1 flat-library path must be one non-empty path.", call. = FALSE)
  }
  value <- path.expand(as.character(path[[1]]))
  if (!grepl("^(/|[A-Za-z]:[/\\\\])", value)) {
    value <- file.path(project_root, value)
  }
  normalizePath(value, winslash = "/", mustWork = must_work)
}

`.ms1_flat_validate_table` <- 
function(table, expected_nrow = NULL) {
  if (!is.data.frame(table) || !identical(class(table), "data.frame")) {
    stop("MS1 flat library must deserialize to a plain data.frame.", call. = FALSE)
  }
  missing_fields <- setdiff(.ms1_flat_required_fields, names(table))
  if (length(missing_fields)) {
    stop(
      "MS1 flat library is missing required field(s): ",
      paste(missing_fields, collapse = ", "),
      call. = FALSE
    )
  }
  if (!is.null(expected_nrow) && !identical(nrow(table), as.integer(expected_nrow))) {
    stop(
      "MS1 flat library row count must be ", as.integer(expected_nrow),
      "; observed ", nrow(table), ".",
      call. = FALSE
    )
  }
  nonempty <- function(x) !is.na(x) & nzchar(trimws(as.character(x)))
  source_id_col <- .ms1_first_column(table, c("source_Lab_IDs", "Lab.ID", "source_lab_ids"))
  if (is.na(source_id_col) || any(!nonempty(table[[source_id_col]]))) {
    stop("MS1 flat library contains missing source Lab.ID provenance.", call. = FALSE)
  }
  if (any(!nonempty(table[["HMDB.ID"]]))) {
    stop("MS1 flat library contains missing HMDB.ID values.", call. = FALSE)
  }
  if (any(!nonempty(table[["Compound.Name"]]))) {
    stop("MS1 flat library contains missing Compound.Name values.", call. = FALSE)
  }
  mz <- suppressWarnings(as.numeric(table[["theoretical_mz"]]))
  if (any(!is.finite(mz) | mz <= 0)) {
    stop("MS1 flat library contains invalid theoretical_mz values.", call. = FALSE)
  }
  invisible(table)
}

`.ms1_normalize_hmdb_id` <- 
function(x) {
  out <- toupper(trimws(as.character(x)))
  hit <- grepl("^HMDB[0-9]+$", out)
  digits <- sub("^HMDB", "", out[hit])
  short <- nchar(digits) < 7L
  digits[short] <- paste0(strrep("0", 7L - nchar(digits[short])), digits[short])
  out[hit] <- paste0("HMDB", digits)
  out
}

`.msi_add_scalar_column_safe` <- 
function(df, name, value) {
  if (!is.data.frame(df)) {
    df <- as.data.frame(df, stringsAsFactors = FALSE, check.names = FALSE)
  }
  if (!is.character(name) || length(name) != 1L || is.na(name) || !nzchar(name)) {
    stop(
      ".msi_add_scalar_column_safe() requires a single non-empty column name.",
      call. = FALSE
    )
  }
  if (is.null(value)) {
    df[[name]] <- rep(NA, nrow(df))
    return(df)
  }
  if (length(value) != 1L) {
    stop(
      ".msi_add_scalar_column_safe() expects a scalar value for column ",
      shQuote(name), "; got length ", length(value), ".",
      call. = FALSE
    )
  }
  df[[name]] <- rep(value, nrow(df))
  df
}

`.msi_adjusted_rand_index` <- 
function(x, y) {
  if (length(x) != length(y) || length(x) < 2L) {
    stop("ARI inputs must have the same length of at least two.", call. = FALSE)
  }
  tab <- table(as.character(x), as.character(y))
  choose2 <- function(z) sum(z * (z - 1) / 2)
  observed <- choose2(tab)
  row_pairs <- choose2(rowSums(tab))
  col_pairs <- choose2(colSums(tab))
  total_pairs <- length(x) * (length(x) - 1) / 2
  expected <- row_pairs * col_pairs / total_pairs
  denominator <- (row_pairs + col_pairs) / 2 - expected
  if (isTRUE(all.equal(denominator, 0))) {
    return(if (isTRUE(all.equal(observed, expected))) 1 else 0)
  }
  (observed - expected) / denominator
}

`.msi_canonicalize_feature_membership` <- 
function(membership,
                                                 prefix = "PT02_LEIDEN") {
  membership <- tibble::as_tibble(membership) |>
    dplyr::transmute(
      feature_id = as.character(feature_id),
      raw_cluster = as.character(cluster_id)
    )
  if (!nrow(membership) || anyNA(membership$feature_id) ||
      any(!nzchar(membership$feature_id)) || anyDuplicated(membership$feature_id) ||
      anyNA(membership$raw_cluster) || any(!nzchar(membership$raw_cluster))) {
    stop("Feature membership must contain one resolved cluster per unique feature.", call. = FALSE)
  }
  ordering <- membership |>
    dplyr::arrange(raw_cluster, feature_id) |>
    dplyr::group_by(raw_cluster) |>
    dplyr::summarise(
      member_key = paste(feature_id, collapse = "|"),
      cluster_size = dplyr::n(),
      .groups = "drop"
    ) |>
    dplyr::arrange(member_key, raw_cluster) |>
    dplyr::mutate(
      cluster_index = as.integer(dplyr::row_number()),
      cluster_id = sprintf("%s_%04d", prefix, cluster_index)
    )
  membership |>
    dplyr::inner_join(
      dplyr::select(ordering, raw_cluster, cluster_id, cluster_index, cluster_size),
      by = "raw_cluster"
    ) |>
    dplyr::transmute(
      feature_id,
      cluster_id = as.character(cluster_id),
      cluster_index = as.integer(cluster_index),
      cluster_size = as.integer(cluster_size),
      is_singleton = cluster_size == 1L
    ) |>
    dplyr::arrange(feature_id)
}

`.msi_missing_summary` <- 
function(context, missing, output_dir = "output", details = NULL) {
  dirs <- .msi_workflow_dirs(output_dir)
  missing <- unique(as.character(missing))
  out <- data.frame(
    context = rep(as.character(context), length(missing)),
    missing_field = missing,
    details = rep(if (is.null(details)) "" else as.character(details), length(missing)),
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
  path <- file.path(dirs$summaries, paste0(gsub("[^A-Za-z0-9_]+", "_", context), "_missing_data.csv"))
  .msi_write_csv(out, path)
  warning(context, " is missing: ", paste(missing, collapse = ", "),
          ". See ", path, call. = FALSE)
  invisible(out)
}

`.msi_resolve_observed_mz` <- 
function(metadata) {
  if (!is.data.frame(metadata)) {
    metadata <- as.data.frame(metadata, stringsAsFactors = FALSE, check.names = FALSE)
  }
  n <- nrow(metadata)
  legacy <- if ("observed_mz" %in% names(metadata)) {
    suppressWarnings(as.numeric(metadata$observed_mz))
  } else {
    rep(NA_real_, n)
  }
  full <- if ("observed_mz_full" %in% names(metadata)) {
    suppressWarnings(as.numeric(metadata$observed_mz_full))
  } else {
    rep(NA_real_, n)
  }
  if (length(legacy) != n || length(full) != n) {
    stop("Observed m/z metadata columns must match the number of rows.", call. = FALSE)
  }
  use_full <- is.finite(full)
  resolved <- legacy
  resolved[use_full] <- full[use_full]
  source <- ifelse(use_full, "observed_mz_full",
                   ifelse(is.finite(legacy), "observed_mz", "unavailable"))
  list(
    mz = resolved,
    source = source,
    observed_mz = legacy,
    observed_mz_full = full
  )
}

`.msi_run_weighted_leiden_stable` <- 
function(feature_ids,
                                            edge_tbl,
                                            resolution = 1,
                                            restarts = 20L,
                                            seeds = seq_len(20L),
                                            stability_ari_threshold = 0.90,
                                            prefix = "PT02_LEIDEN") {
  if (!requireNamespace("igraph", quietly = TRUE) ||
      !exists("cluster_leiden", envir = asNamespace("igraph"), inherits = FALSE)) {
    stop("weighted_leiden requires igraph::cluster_leiden().", call. = FALSE)
  }
  feature_ids <- unique(as.character(feature_ids))
  if (!length(feature_ids) || anyNA(feature_ids) || any(!nzchar(feature_ids))) {
    stop("weighted_leiden requires non-empty unique feature ids.", call. = FALSE)
  }
  if (!is.data.frame(edge_tbl) || !all(c("feature_i", "feature_j", "weight") %in% names(edge_tbl))) {
    stop("weighted_leiden edge_tbl must contain feature_i, feature_j, weight.", call. = FALSE)
  }
  if (!identical(as.integer(restarts), length(seeds)) || restarts < 1L) {
    stop("weighted_leiden restarts must equal the number of supplied seeds.", call. = FALSE)
  }
  if (!is.finite(stability_ari_threshold) || stability_ari_threshold < -1 || stability_ari_threshold > 1) {
    stop("weighted_leiden stability_ari_threshold must be in [-1, 1].", call. = FALSE)
  }
  edge_df <- data.frame(
    from = as.character(edge_tbl$feature_i),
    to = as.character(edge_tbl$feature_j),
    weight = as.numeric(edge_tbl$weight),
    stringsAsFactors = FALSE
  )
  edge_df <- edge_df[
    is.finite(edge_df$weight) & edge_df$from %in% feature_ids & edge_df$to %in% feature_ids &
      edge_df$from != edge_df$to,
    , drop = FALSE
  ]
  edge_df <- edge_df[order(edge_df$from, edge_df$to, -edge_df$weight), , drop = FALSE]
  if (anyDuplicated(paste(edge_df$from, edge_df$to, sep = "\r"))) {
    stop("weighted_leiden edge_tbl contains duplicate undirected edge endpoints.", call. = FALSE)
  }
  graph <- igraph::graph_from_data_frame(
    edge_df,
    directed = FALSE,
    vertices = data.frame(name = feature_ids, stringsAsFactors = FALSE)
  )
  if (igraph::vcount(graph) != length(feature_ids)) {
    stop("weighted_leiden graph vertex count does not match feature ids.", call. = FALSE)
  }
  raw_runs <- lapply(seq_len(as.integer(restarts)), function(i) {
    group_msi_features(
      structure(
        data.frame(feature_i = edge_df$from, feature_j = edge_df$to, weight = edge_df$weight,
                   stringsAsFactors = FALSE),
        feature_ids = feature_ids
      ),
      method = "leiden", resolution = resolution, seed = as.integer(seeds[[i]])
    ) |>
      dplyr::arrange(feature_id)
  })
  canonical_runs <- lapply(raw_runs, .msi_canonicalize_feature_membership, prefix = prefix)
  ari <- diag(1, nrow = length(canonical_runs), ncol = length(canonical_runs))
  if (length(canonical_runs) > 1L) {
    for (i in seq_len(length(canonical_runs) - 1L)) {
      for (j in (i + 1L):length(canonical_runs)) {
        value <- .msi_adjusted_rand_index(canonical_runs[[i]]$cluster_id, canonical_runs[[j]]$cluster_id)
        ari[i, j] <- value
        ari[j, i] <- value
      }
    }
  }
  mean_ari <- rowMeans(ari)
  pairwise_ari <- if (length(canonical_runs) > 1L) ari[upper.tri(ari)] else 1
  mean_pairwise_ari <- mean(pairwise_ari)
  if (mean_pairwise_ari < stability_ari_threshold) {
    stop(sprintf("weighted_leiden stability ARI below threshold %.2f.", stability_ari_threshold), call. = FALSE)
  }
  run_key <- vapply(canonical_runs, function(x) paste(x$cluster_id, collapse = "|"), character(1))
  chosen <- order(-mean_ari, run_key, seq_along(run_key))[[1L]]
  list(
    graph = graph,
    membership = canonical_runs[[chosen]],
    edge_table = tibble::as_tibble(edge_df) |>
      dplyr::transmute(feature_i = from, feature_j = to, weight = as.numeric(weight)),
    diagnostics = list(
      algorithm = "weighted_leiden",
      resolution = as.numeric(resolution),
      restarts = as.integer(restarts),
      seeds = as.integer(seeds),
      stability_ari_threshold = as.numeric(stability_ari_threshold),
      ari = ari,
      min_pairwise_ari = min(pairwise_ari),
      mean_pairwise_ari = mean_pairwise_ari,
      selected_mean_ari = mean_ari[[chosen]],
      selected_restart = as.integer(chosen)
    )
  )
}

`.msi_workflow_dirs` <- 
function(output_dir = "output") {
  dirs <- list(
    root = output_dir,
    figures = file.path(output_dir, "figures"),
    tables = file.path(output_dir, "tables"),
    summaries = file.path(output_dir, "summaries"),
    rds = file.path(output_dir, "rds"),
    logs = file.path(output_dir, "logs"),
    temp = file.path(output_dir, "temp")
  )
  invisible(lapply(dirs, dir.create, recursive = TRUE, showWarnings = FALSE))
  dirs
}

`.msi_write_csv` <- 
function(x, path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  x <- as.data.frame(x, stringsAsFactors = FALSE, check.names = FALSE)
  if (requireNamespace("data.table", quietly = TRUE)) {
    threads <- suppressWarnings(as.integer(getOption("msidataset.n_workers", 1L)))
    if (!is.finite(threads) || threads < 1L) threads <- 1L
    data.table::fwrite(x, path, na = "", nThread = threads)
  } else {
    utils::write.csv(x, path, row.names = FALSE, na = "")
  }
  invisible(path)
}

`%||%` <- 
function(x, y) if (is.null(x) || length(x) == 0L) y else x

`60_mark_chemical_evidence_not_evaluated` <- 
function(reverse_match) {
  scores <- as.data.frame(reverse_match, stringsAsFactors = FALSE, check.names = FALSE)
  if (nrow(scores) == 0L) {
    scores$chem_status <- character()
    scores$chem_role <- character()
    scores$chem_explanation <- character()
    scores$role_delta_mz <- numeric()
    scores$role_ppm_error <- numeric()
    scores$role_evidence_n <- integer()
    return(scores)
  }
  scores$chem_status <- "not_evaluated"
  scores$chem_role <- "unknown"
  scores$chem_explanation <- "chemical evidence was not evaluated in library-conditioned reverse matching"
  scores$role_delta_mz <- NA_real_
  scores$role_ppm_error <- NA_real_
  scores$role_evidence_n <- 0L
  scores
}

`7_annotate_compounds` <- 
function(...,
                                   fmc_annotation = NULL,
                                   fmc_enabled = NULL,
                                   fmc_ms1_ppm = NULL,
                                   fmc_ion_mode = NULL,
                                   fmc_max_candidates_per_feature = NULL,
                                   fmc_max_metabolites_per_feature = NULL,
                                   fmc_max_features_per_metabolite_soft = NULL,
                                   fmc_max_features_per_metabolite_hard = NULL,
                                   fmc_use_isotope = NULL,
                                   fmc_use_neutral_loss = NULL,
                                   fmc_use_spatial_support = NULL,
                                   fmc_metabolite_reference = NULL) {
  warning(
    "7_annotate_compounds is retained for compatibility; use calculate_score for canonical scoring.",
    call. = FALSE
  )
  fmc_args <- list(
    fmc_annotation = fmc_annotation,
    fmc_enabled = fmc_enabled,
    fmc_ms1_ppm = fmc_ms1_ppm,
    fmc_ion_mode = fmc_ion_mode,
    fmc_max_candidates_per_feature = fmc_max_candidates_per_feature,
    fmc_max_metabolites_per_feature = fmc_max_metabolites_per_feature,
    fmc_max_features_per_metabolite_soft = fmc_max_features_per_metabolite_soft,
    fmc_max_features_per_metabolite_hard = fmc_max_features_per_metabolite_hard,
    fmc_use_isotope = fmc_use_isotope,
    fmc_use_neutral_loss = fmc_use_neutral_loss,
    fmc_use_spatial_support = fmc_use_spatial_support,
    fmc_metabolite_reference = fmc_metabolite_reference
  )
  if (isTRUE(.ann_fmc_requested(fmc_args))) {
    warning(
      "FMC is now Module 8. Please call calculate_mfc() after calculate_score().",
      call. = FALSE
    )
  }
  calculate_score(...)
}

`71_calibrate_scores` <- 
function(object, spatial_scores, reverse_match = NULL,
                                  matched_peak_pairs = NULL) {
  empty <- .ann_empty_score_table()
  if (is.null(reverse_match)) reverse_match <- object@annotation_results$reverse_match
  reverse_match <- as.data.frame(reverse_match, stringsAsFactors = FALSE, check.names = FALSE)
  if (nrow(reverse_match) == 0L) {
    return(empty)
  }

  reverse_match <- .ann_add_missing_columns(reverse_match, setdiff(colnames(empty), c(
    "s_spat", "S_FE_cal", "S_reverse_cos_cal", "S_lib_coverage_cal",
    "S_spatial_reverse_cal", "S_precursor_cal", "S_chem_cal",
    "s_fe_cal", "s_arch_cal", "s_rev_cal", "u_spec_cal", "s_chem_cal", "s_spat_cal"
  )))
  reverse_match$entropy_S <- .ann_first_numeric(reverse_match, c("entropy_S", "S_FE", "s_fe"))
  reverse_match$S_FE <- reverse_match$entropy_S
  reverse_match$reverse_masked_cosine <- .ann_first_numeric(
    reverse_match, c("reverse_masked_cosine", "score_ms2", "S_reverse_cos")
  )
  reverse_match$reverse_masked_cosine <- .ann_clip01_zero_na(reverse_match$reverse_masked_cosine)
  reverse_match$S_reverse_cos <- reverse_match$reverse_masked_cosine
  if (all(!is.finite(suppressWarnings(as.numeric(reverse_match$reverse_masked_cosine))))) {
    warning(
      "reverse_masked_cosine is unavailable. Recalculate Step 6 before Step 7.",
      call. = FALSE
    )
  }
  reverse_match$S_lib_coverage <- .ann_first_numeric(reverse_match, c("S_lib_coverage", "u_spec"))
  reverse_match$S_SAS <- .ann_first_numeric(reverse_match, c("S_SAS", "S_spatial_reverse", "s_spat"))
  reverse_match$S_SAS_missing <- !is.finite(suppressWarnings(as.numeric(reverse_match$S_SAS)))
  reverse_match$S_SAS_score_for_total <- ifelse(
    reverse_match$S_SAS_missing,
    0,
    .ann_clip01_zero_na(reverse_match$S_SAS)
  )
  reverse_match$S_spatial_reverse <- reverse_match$S_SAS
  reverse_match$S_chem <- .ann_first_numeric(reverse_match, c("S_chem", "s_chem"))
  sas_from_pairs <- .ann_S_SAS_from_object(object, matched_peak_pairs)
  if (nrow(sas_from_pairs) > 0L) {
    sas_lookup <- stats::setNames(
      suppressWarnings(as.numeric(sas_from_pairs$S_SAS)),
      .ann_split_key(sas_from_pairs, c("cluster_id", "anchor_feature_id", "candidate_key"))
    )
    n_lookup <- stats::setNames(
      suppressWarnings(as.integer(sas_from_pairs$n_SAS_matched_peaks)),
      .ann_split_key(sas_from_pairs, c("cluster_id", "anchor_feature_id", "candidate_key"))
    )
    sum_lookup <- stats::setNames(
      suppressWarnings(as.numeric(sas_from_pairs$sum_matched_library_intensity)),
      .ann_split_key(sas_from_pairs, c("cluster_id", "anchor_feature_id", "candidate_key"))
    )
    status_lookup <- stats::setNames(
      as.character(sas_from_pairs$S_SAS_status),
      .ann_split_key(sas_from_pairs, c("cluster_id", "anchor_feature_id", "candidate_key"))
    )
    rm_key <- .ann_split_key(reverse_match, c("cluster_id", "anchor_feature_id", "candidate_key"))
    sas_value <- suppressWarnings(as.numeric(sas_lookup[rm_key]))
    fill_sas <- is.finite(sas_value)
    reverse_match$S_SAS[fill_sas] <- sas_value[fill_sas]
    reverse_match$S_spatial_reverse[fill_sas] <- sas_value[fill_sas]
    reverse_match$n_SAS_matched_peaks <- suppressWarnings(as.integer(n_lookup[rm_key]))
    reverse_match$sum_matched_library_intensity <- suppressWarnings(as.numeric(sum_lookup[rm_key]))
    reverse_match$S_SAS_status <- as.character(status_lookup[rm_key])
    reverse_match$S_SAS_raw <- suppressWarnings(as.numeric(sas_from_pairs$S_SAS_raw[match(rm_key, names(sas_lookup))]))
  }
  if (all(is.na(suppressWarnings(as.numeric(reverse_match$S_precursor))))) {
    ppm_error <- suppressWarnings(as.numeric(reverse_match$precursor_error_ppm))
    reverse_match$S_precursor <- ifelse(is.finite(ppm_error), 1 / (1 + abs(ppm_error) / 10), 0)
  }
  reverse_match$s_fe <- reverse_match$S_FE
  reverse_match$s_rev <- reverse_match$reverse_masked_cosine
  reverse_match$u_spec <- reverse_match$S_lib_coverage
  reverse_match$s_spat <- reverse_match$S_SAS
  reverse_match$s_chem <- reverse_match$S_chem
  reverse_match$n_matched_peaks <- suppressWarnings(as.integer(.ann_first_numeric(reverse_match, c("n_matched_peaks", "n_reverse_match", "n_match"))))
  reverse_match$n_match_pairs <- suppressWarnings(as.integer(.ann_first_numeric(reverse_match, c("n_match_pairs", "n_reverse_match", "n_match", "n_matched_peaks"))))
  missing_name <- is.na(reverse_match$metabolite_name) | !nzchar(as.character(reverse_match$metabolite_name))
  reverse_match$metabolite_name[missing_name] <- as.character(reverse_match$candidate_name[missing_name])
  spatial_scores <- as.data.frame(spatial_scores, stringsAsFactors = FALSE, check.names = FALSE)
  spatial_scores <- .ann_add_missing_columns(spatial_scores, c(
    "cluster_id", "anchor_feature_id", "candidate_key", "s_spat"
  ))
  if (nrow(spatial_scores) > 0L) {
    spat_lookup <- stats::setNames(
      suppressWarnings(as.numeric(spatial_scores$s_spat)),
      .ann_split_key(spatial_scores, c("cluster_id", "anchor_feature_id", "candidate_key"))
    )
    rm_key <- .ann_split_key(reverse_match, c("cluster_id", "anchor_feature_id", "candidate_key"))
    spat_value <- suppressWarnings(as.numeric(spat_lookup[rm_key]))
    fill_spat <- !is.finite(reverse_match$S_spatial_reverse) & is.finite(spat_value)
    reverse_match$S_spatial_reverse[fill_spat] <- spat_value[fill_spat]
    reverse_match$s_spat[fill_spat] <- spat_value[fill_spat]
  }
  out <- reverse_match
  out <- .ann_add_missing_columns(out, colnames(empty))
  out$s_fe <- out$S_FE
  out$s_rev <- out$reverse_masked_cosine
  out$u_spec <- out$S_lib_coverage
  out$s_spat <- out$S_SAS
  out$s_chem <- out$S_chem
  out$n_matched_peaks <- suppressWarnings(as.integer(.ann_first_numeric(out, c("n_matched_peaks", "n_reverse_match", "n_match"))))
  out$n_match_pairs <- suppressWarnings(as.integer(.ann_first_numeric(out, c("n_match_pairs", "n_reverse_match", "n_match", "n_matched_peaks"))))

  out$S_FE_cal <- .ann_within_anchor_percentile(out, "S_FE")
  out$S_reverse_cos_cal <- .ann_clip01_zero_na(out$reverse_masked_cosine)
  out$S_lib_coverage_cal <- .ann_clip01_zero_na(out$S_lib_coverage)
  out$S_spatial_reverse_cal <- .ann_clip01_zero_na(out$S_spatial_reverse)
  out$S_precursor_cal <- .ann_clip01_zero_na(out$S_precursor)
  out$S_chem_cal <- .ann_clip01_zero_na(out$S_chem)
  out$s_fe_cal <- out$S_FE_cal
  out$s_arch_cal <- .ann_within_anchor_percentile(out, "s_arch")
  out$s_rev_cal <- .ann_clip01_zero_na(out$reverse_masked_cosine)
  out$u_spec_cal <- out$S_lib_coverage_cal
  out$s_chem_cal <- out$S_chem_cal
  out$s_spat_cal <- out$S_spatial_reverse_cal

  out <- out[, colnames(empty), drop = FALSE]
  rownames(out) <- NULL
  out
}

`72_fuse_scores` <- 
function(score_table,
                             weights,
                             reverse_masked_cosine_weight = NULL,
                             S_SAS_weight = NULL,
                             primary_spectral_score = "reverse_masked_cosine") {
  primary_spectral_score <- as.character(primary_spectral_score[[1]])
  if (!identical(primary_spectral_score, "reverse_masked_cosine")) {
    primary_spectral_score <- "reverse_masked_cosine"
  }
  if (is.list(weights) && !is.data.frame(weights)) {
    weights <- unlist(weights, use.names = TRUE)
  }
  legacy_map <- c(
    s_rev = "S_reverse_cos",
    u_spec = "S_lib_coverage",
    s_spat = "S_spatial_reverse",
    s_fe = "S_FE",
    s_chem = "S_chem"
  )
  for (old in names(legacy_map)) {
    if (old %in% names(weights) && !legacy_map[[old]] %in% names(weights)) {
      names(weights)[names(weights) == old] <- legacy_map[[old]]
    }
  }
  if ("s_arch" %in% names(weights) && !"S_precursor" %in% names(weights)) {
    weights <- c(weights, S_precursor = weights[["s_arch"]])
    weights <- weights[names(weights) != "s_arch"]
  }
  required <- c("S_reverse_cos", "S_lib_coverage", "S_spatial_reverse", "S_FE", "S_precursor", "S_chem")
  if (!is.numeric(weights) || any(is.na(weights))) {
    stop("Fusion weights must be numeric and non-missing.", call. = FALSE)
  }
  missing_weights <- setdiff(required, names(weights))
  if (length(missing_weights) > 0L) {
    stop("Fusion weights are missing: ", paste(missing_weights, collapse = ", "), call. = FALSE)
  }
  weights <- weights[required]
  if (any(weights < 0)) {
    stop("Fusion weights must be non-negative.", call. = FALSE)
  }
  weight_sum <- sum(weights)
  if (!is.finite(weight_sum) || weight_sum <= 0) {
    stop("Fusion weights must sum to a positive value.", call. = FALSE)
  }
  weights <- weights / weight_sum

  out <- as.data.frame(score_table, stringsAsFactors = FALSE, check.names = FALSE)
  out <- .ann_add_missing_columns(out, colnames(.ann_empty_score_table()))
  if (nrow(out) == 0L) {
    out$S_total <- numeric()
    out$rank_total <- integer()
    out$s_final <- numeric()
    out$primary_spectral_score <- character()
    out$weight_reverse_masked_cosine <- numeric()
    out$weight_S_reverse_cos <- numeric()
    out$weight_S_lib_coverage <- numeric()
    out$weight_S_spatial_reverse <- numeric()
    out$weight_S_FE <- numeric()
    out$weight_S_precursor <- numeric()
    out$weight_S_chem <- numeric()
    out$weight_s_fe <- numeric()
    out$weight_s_arch <- numeric()
    out$weight_s_rev <- numeric()
    out$weight_u_spec <- numeric()
    out$weight_s_chem <- numeric()
    out$weight_s_spat <- numeric()
    return(out)
  }
  if (!is.null(reverse_masked_cosine_weight) || !is.null(S_SAS_weight)) {
    reverse_masked_cosine_weight <- suppressWarnings(as.numeric(reverse_masked_cosine_weight[[1]]))
    S_SAS_weight <- suppressWarnings(as.numeric(S_SAS_weight[[1]]))
    if (!is.finite(reverse_masked_cosine_weight) || !is.finite(S_SAS_weight) ||
        reverse_masked_cosine_weight < 0 || S_SAS_weight < 0 ||
        reverse_masked_cosine_weight + S_SAS_weight <= 0) {
      stop("reverse_masked_cosine weight and S_SAS_weight must be non-negative and sum to a positive value.", call. = FALSE)
    }
    weight_sum <- reverse_masked_cosine_weight + S_SAS_weight
    reverse_masked_cosine_weight <- reverse_masked_cosine_weight / weight_sum
    S_SAS_weight <- S_SAS_weight / weight_sum
    spectral_score <- .ann_clip01_zero_na(out$reverse_masked_cosine)
    out$S_SAS_missing <- !is.finite(suppressWarnings(as.numeric(out$S_SAS)))
    out$S_SAS_score_for_total <- ifelse(out$S_SAS_missing, 0, .ann_clip01_zero_na(out$S_SAS))
    S_SAS_for_fusion <- out$S_SAS_score_for_total
    out$S_total <- reverse_masked_cosine_weight * spectral_score + S_SAS_weight * S_SAS_for_fusion
    out$s_final <- out$S_total
    out$weight_reverse_masked_cosine <- reverse_masked_cosine_weight
    out$weight_S_SAS <- S_SAS_weight
    out$primary_spectral_score <- primary_spectral_score
    out$score_ms2 <- spectral_score
    out$score_spatial <- S_SAS_for_fusion
    out$score_total_pre_mfc <- out$S_total
    out <- .ann_add_rank_total(out)
    return(out)
  }
  out$s_final <- weights[["S_reverse_cos"]] * .ann_clip01_zero_na(out$S_reverse_cos_cal) +
    weights[["S_lib_coverage"]] * .ann_clip01_zero_na(out$S_lib_coverage_cal) +
    weights[["S_spatial_reverse"]] * .ann_clip01_zero_na(out$S_spatial_reverse_cal) +
    weights[["S_FE"]] * .ann_clip01_zero_na(out$S_FE_cal) +
    weights[["S_precursor"]] * .ann_clip01_zero_na(out$S_precursor_cal) +
    weights[["S_chem"]] * .ann_clip01_zero_na(out$S_chem_cal)
  out$s_final <- pmin(pmax(as.numeric(out$s_final), 0), 1)
  out$S_total <- out$s_final
  out$primary_spectral_score <- primary_spectral_score
  out <- .ann_add_rank_total(out)
  out$weight_S_reverse_cos <- weights[["S_reverse_cos"]]
  out$weight_S_lib_coverage <- weights[["S_lib_coverage"]]
  out$weight_S_spatial_reverse <- weights[["S_spatial_reverse"]]
  out$weight_S_FE <- weights[["S_FE"]]
  out$weight_S_precursor <- weights[["S_precursor"]]
  out$weight_S_chem <- weights[["S_chem"]]
  out$weight_s_fe <- weights[["S_FE"]]
  out$weight_s_arch <- weights[["S_precursor"]]
  out$weight_s_rev <- weights[["S_reverse_cos"]]
  out$weight_u_spec <- weights[["S_lib_coverage"]]
  out$weight_s_chem <- weights[["S_chem"]]
  out$weight_s_spat <- weights[["S_spatial_reverse"]]
  out
}

`73_generate_decoy_library` <- 
function(object,
                                        targets = NULL,
                                        strategy = c("peak_shift", "intensity_permutation", "none"),
                                        decoy_n = 1,
                                        decoy_shift_da = 0.05,
                                        seed = 123) {
  strategy <- match.arg(strategy)
  empty <- .ann_empty_decoy_library()
  if (strategy == "none") {
    return(empty)
  }
  restore_rng <- .msi_rng_restorer()
  on.exit(restore_rng(), add = TRUE)
  set.seed(seed)

  if (is.null(targets)) targets <- object@annotation_results$reverse_match
  targets <- as.data.frame(targets, stringsAsFactors = FALSE, check.names = FALSE)
  if (nrow(targets) == 0L) {
    return(empty)
  }
  decoy_n <- as.integer(decoy_n)
  if (!is.finite(decoy_n) || decoy_n < 1L) {
    decoy_n <- 1L
  }

  spectrum_rows <- list()
  peak_rows <- list()
  for (i in seq_len(nrow(targets))) {
    target <- targets[i, , drop = FALSE]
    ref <- .ann_reference_peaks(target, object)
    if (nrow(ref) == 0L) next
    for (d in seq_len(decoy_n)) {
      decoy_key <- paste("DECOY", target$candidate_key[[1]], d, sep = "::")
      peaks <- ref
      if (strategy == "peak_shift") {
        sign <- if (d %% 2L == 0L) -1 else 1
        peaks$peak_mz <- pmax(abs(peaks$peak_mz + sign * decoy_shift_da * d), .Machine$double.eps)
      } else if (strategy == "intensity_permutation" && nrow(peaks) > 1L) {
        peaks$intensity_sum1 <- sample(peaks$intensity_sum1)
      }
      total <- sum(peaks$intensity_sum1, na.rm = TRUE)
      if (is.finite(total) && total > 0) {
        peaks$intensity_sum1 <- peaks$intensity_sum1 / total
      }
      spectrum_rows[[length(spectrum_rows) + 1L]] <- data.frame(
        cluster_id = target$cluster_id[[1]],
        anchor_feature_id = as.character(target$anchor_feature_id[[1]]),
        anchor_mz = suppressWarnings(as.numeric(target$anchor_mz[[1]])),
        candidate_key = as.character(target$candidate_key[[1]]),
        decoy_candidate_key = decoy_key,
        candidate_id = as.character(target$candidate_id[[1]]),
        spectrum_id = as.character(target$spectrum_id[[1]]),
        candidate_name = as.character(target$candidate_name[[1]]),
        library_source = as.character(target$library_source[[1]]),
        database_source = as.character(target$database_source[[1]]),
        precursor_mz = suppressWarnings(as.numeric(target$precursor_mz[[1]])),
        formula = as.character(target$formula[[1]]),
        inchikey = as.character(target$inchikey[[1]]),
        smiles = as.character(target$smiles[[1]]),
        adduct = as.character(target$adduct[[1]]),
        polarity = as.character(target$polarity[[1]]),
        is_decoy = TRUE,
        stringsAsFactors = FALSE,
        check.names = FALSE
      )
      peak_rows[[length(peak_rows) + 1L]] <- data.frame(
        decoy_candidate_key = decoy_key,
        candidate_key = as.character(target$candidate_key[[1]]),
        spectrum_id = as.character(target$spectrum_id[[1]]),
        peak_mz = peaks$peak_mz,
        intensity_sum1 = peaks$intensity_sum1,
        is_decoy = TRUE,
        stringsAsFactors = FALSE,
        check.names = FALSE
      )
    }
  }

  if (length(spectrum_rows) == 0L || length(peak_rows) == 0L) {
    return(empty)
  }
  list(
    decoy_spectrum_table = do.call(rbind, spectrum_rows),
    decoy_peak_table = do.call(rbind, peak_rows)
  )
}

`74_score_decoys` <- 
function(object,
                              decoy_library,
                              weights,
                              score_table,
                              ppm = 10,
                              tolerance_da = NULL,
                              reverse_masked_cosine_weight = NULL,
                              S_SAS_weight = NULL,
                              primary_spectral_score = "reverse_masked_cosine") {
  primary_spectral_score <- as.character(primary_spectral_score[[1]])
  if (!identical(primary_spectral_score, "reverse_masked_cosine")) {
    primary_spectral_score <- "reverse_masked_cosine"
  }
  empty <- .ann_empty_decoy_scores()
  if (is.null(decoy_library$decoy_spectrum_table) ||
      is.null(decoy_library$decoy_peak_table) ||
      nrow(decoy_library$decoy_spectrum_table) == 0L ||
      nrow(decoy_library$decoy_peak_table) == 0L) {
    return(empty)
  }

  decoy_spectrum <- as.data.frame(decoy_library$decoy_spectrum_table,
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
  decoy_peaks <- as.data.frame(decoy_library$decoy_peak_table,
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
  pseudo <- as.data.frame(object@pseudo_ms2$pseudo_ms2_table,
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
  rows <- list()
  for (i in seq_len(nrow(decoy_spectrum))) {
    decoy <- decoy_spectrum[i, , drop = FALSE]
    query <- pseudo[
      as.character(pseudo$cluster_id) == as.character(decoy$cluster_id[[1]]) &
        as.character(pseudo$anchor_feature_id) == as.character(decoy$anchor_feature_id[[1]]),
      ,
      drop = FALSE
    ]
    ref <- decoy_peaks[decoy_peaks$decoy_candidate_key == decoy$decoy_candidate_key[[1]], , drop = FALSE]
    matches <- .ann_match_query_reference(query, ref, ppm = ppm, tolerance_da = tolerance_da)
    reverse_masked_cosine <- .ann_reverse_cosine(query, ref, matches)
    denom <- sum(suppressWarnings(as.numeric(query$intensity)), na.rm = TRUE)
    matched_query <- unique(matches$query_peak_feature_id)
    u_spec <- if (nrow(matches) > 0L && is.finite(denom) && denom > 0) {
      sum(suppressWarnings(as.numeric(query$intensity[as.character(query$peak_feature_id) %in% matched_query])), na.rm = TRUE) / denom
    } else {
      0
    }
    q <- pmax(suppressWarnings(as.numeric(query$intensity)), 0)
    g <- pmax(suppressWarnings(as.numeric(query$graph_weight)), 0)
    spatial_denom <- sum(q * g, na.rm = TRUE)
    s_spat <- if (nrow(matches) > 0L && is.finite(spatial_denom) && spatial_denom > 0) {
      idx <- as.character(query$peak_feature_id) %in% matched_query
      sum(q[idx] * g[idx], na.rm = TRUE) / spatial_denom
    } else {
      0
    }
    query_sas <- suppressWarnings(as.numeric(query$SAS))
    if (!length(query_sas) || all(!is.finite(query_sas))) {
      query_sas <- suppressWarnings(as.numeric(query$graph_weight))
    }
    sas_lookup <- stats::setNames(query_sas, as.character(query$peak_feature_id))
    matched_sas <- suppressWarnings(as.numeric(sas_lookup[as.character(matches$query_peak_feature_id)]))
    matched_weights <- suppressWarnings(as.numeric(matches$reference_intensity))
    sas_keep <- is.finite(matched_sas) & is.finite(matched_weights) & matched_weights > 0
    S_SAS <- if (any(sas_keep)) {
      stats::weighted.mean(matched_sas[sas_keep], matched_weights[sas_keep])
    } else {
      0
    }
    rows[[length(rows) + 1L]] <- data.frame(
      cluster_id = decoy$cluster_id[[1]],
      anchor_feature_id = as.character(decoy$anchor_feature_id[[1]]),
      anchor_mz = suppressWarnings(as.numeric(decoy$anchor_mz[[1]])),
      candidate_key = as.character(decoy$candidate_key[[1]]),
      decoy_candidate_key = as.character(decoy$decoy_candidate_key[[1]]),
      S_FE = 0,
      reverse_masked_cosine = min(max(reverse_masked_cosine, 0), 1),
      S_SAS = min(max(S_SAS, 0), 1),
      S_SAS_score_for_total = min(max(S_SAS, 0), 1),
      S_SAS_missing = FALSE,
      S_reverse_cos = min(max(reverse_masked_cosine, 0), 1),
      S_lib_coverage = min(max(u_spec, 0), 1),
      S_spatial_reverse = min(max(s_spat, 0), 1),
      S_precursor = 0,
      S_chem = NA_real_,
      n_matched_peaks = nrow(matches),
      query_match_ratio = if (nrow(query)) {
        length(unique(matches$query_peak_feature_id)) / nrow(query)
      } else {
        NA_real_
      },
      library_match_ratio = if (nrow(ref)) nrow(matches) / nrow(ref) else NA_real_,
      s_fe = 0,
      s_arch = 0,
      s_rev = min(max(reverse_masked_cosine, 0), 1),
      u_spec = min(max(u_spec, 0), 1),
      s_chem = NA_real_,
      s_spat = min(max(s_spat, 0), 1),
      is_decoy = TRUE,
      stringsAsFactors = FALSE,
      check.names = FALSE
    )
  }
  if (length(rows) == 0L) {
    return(empty)
  }
  raw <- do.call(rbind, rows)
  raw$S_FE_cal <- 0
  raw$S_reverse_cos_cal <- .ann_clip01_zero_na(raw$S_reverse_cos)
  raw$S_lib_coverage_cal <- .ann_clip01_zero_na(raw$S_lib_coverage)
  raw$S_spatial_reverse_cal <- .ann_clip01_zero_na(raw$S_spatial_reverse)
  raw$S_precursor_cal <- 0
  raw$S_chem_cal <- 0
  raw$s_fe_cal <- 0
  raw$s_arch_cal <- 0
  raw$s_rev_cal <- .ann_clip01_zero_na(raw$s_rev)
  raw$u_spec_cal <- .ann_clip01_zero_na(raw$u_spec)
  raw$s_chem_cal <- 0
  raw$s_spat_cal <- .ann_clip01_zero_na(raw$s_spat)
  fused <- `72_fuse_scores`(
    raw,
    weights = weights,
    reverse_masked_cosine_weight = reverse_masked_cosine_weight,
    S_SAS_weight = S_SAS_weight,
    primary_spectral_score = primary_spectral_score
  )
  out <- .ann_add_missing_columns(fused, colnames(empty))
  out <- out[, colnames(empty), drop = FALSE]
  rownames(out) <- NULL
  out
}

`75_estimate_annotation_q_value` <- 
function(target_scores,
                                             decoy_scores,
                                             method = c("target_decoy", "none")) {
  method <- match.arg(method)
  out <- as.data.frame(target_scores, stringsAsFactors = FALSE, check.names = FALSE)
  if (nrow(out) == 0L) {
    out$annotation_q_value <- numeric()
    out$annotation_fdr_method <- character()
    return(out)
  }
  if (method == "none" || is.null(decoy_scores) || nrow(decoy_scores) == 0L) {
    out$annotation_q_value <- NA_real_
    out$annotation_fdr_method <- if (method == "none") "none" else "none_or_no_decoy"
    return(out)
  }

  decoy_scores <- as.data.frame(decoy_scores, stringsAsFactors = FALSE, check.names = FALSE)
  final_col <- .ann_primary_score_column(out)
  out <- out[order(-suppressWarnings(as.numeric(out[[final_col]]))), , drop = FALSE]
  decoy_final <- suppressWarnings(as.numeric(decoy_scores$s_final))
  target_final <- suppressWarnings(as.numeric(out[[final_col]]))
  fdr_raw <- vapply(target_final, function(threshold) {
    decoy_n <- sum(decoy_final >= threshold, na.rm = TRUE)
    target_n <- sum(target_final >= threshold, na.rm = TRUE)
    (1 + decoy_n) / max(1, target_n)
  }, numeric(1))
  q_value <- rev(cummin(rev(fdr_raw)))
  out$annotation_q_value <- pmin(pmax(q_value, 0), 1)
  out$annotation_fdr_method <- "target_decoy"
  rownames(out) <- NULL
  out
}

`76_filter_annotation_results` <- 
function(annotated,
                                           q_value_cutoff = 0.05,
                                           min_final_score = 0,
                                           S_SAS_threshold = NA_real_,
                                           keep_all = TRUE,
                                           confidence_thresholds = list()) {
  out <- as.data.frame(annotated, stringsAsFactors = FALSE, check.names = FALSE)
  if (nrow(out) == 0L) {
    out$rank_within_anchor <- integer()
    out$pass_annotation_filter <- logical()
    out$annotation_confidence <- character()
    return(out)
  }
  out <- .ann_add_annotation_confidence(out, confidence_thresholds = confidence_thresholds)
  out <- .ann_order_final_annotation(out)
  split_anchor <- .ann_split_key(out, c("cluster_id", "anchor_feature_id"))
  out$rank_within_anchor <- unsplit(lapply(split(seq_len(nrow(out)), split_anchor), seq_along), split_anchor)
  out$rank_total <- suppressWarnings(as.integer(out$rank_within_anchor))
  q <- suppressWarnings(as.numeric(out$annotation_q_value))
  score_col <- .ann_primary_score_column(out)
  score <- suppressWarnings(as.numeric(out[[score_col]]))
  out$qvalue_enabled <- is.finite(q)
  out$passes_score_threshold <- is.finite(score) & score >= min_final_score
  out$passes_q_value <- ifelse(out$qvalue_enabled, q <= q_value_cutoff, TRUE)
  out$passes_q_value[is.na(out$passes_q_value)] <- FALSE
  S_SAS_threshold <- suppressWarnings(as.numeric(S_SAS_threshold[[1]]))
  out$passes_S_SAS_threshold <- if (is.finite(S_SAS_threshold)) {
    suppressWarnings(as.numeric(out$S_SAS)) >= S_SAS_threshold
  } else TRUE
  out$passes_S_SAS_threshold[is.na(out$passes_S_SAS_threshold)] <- FALSE
  out$S_SAS_threshold <- if (is.finite(S_SAS_threshold)) S_SAS_threshold else NA_real_
  out$pass_annotation_filter <- out$passes_score_threshold &
    out$passes_q_value & out$passes_S_SAS_threshold
  out$pass_annotation_filter[is.na(out$pass_annotation_filter)] <- FALSE
  confidence <- as.character(out$annotation_confidence)
  out$pass_review <- out$pass_annotation_filter & confidence %in% "review_only"
  out$pass_primary <- out$pass_annotation_filter & !out$pass_review
  out$annotation_pass_class <- ifelse(
    out$pass_primary, "primary",
    ifelse(out$pass_review, "review", "fail")
  )
  out$rule_fail_reason <- NA_character_
  out$rule_fail_reason[!is.finite(score)] <- "failed_step7_missing_score"
  fail_score <- is.na(out$rule_fail_reason) & !out$passes_score_threshold
  out$rule_fail_reason[fail_score] <- "failed_step7_total_score_rule"
  fail_q <- is.na(out$rule_fail_reason) & out$qvalue_enabled & !out$passes_q_value
  out$rule_fail_reason[fail_q] <- "failed_step7_q_value_rule"
  fail_sas <- is.na(out$rule_fail_reason) & !out$passes_S_SAS_threshold
  out$rule_fail_reason[fail_sas] <- "failed_step7_sas_rule"
  out$rule_fail_reason[out$pass_annotation_filter %in% TRUE] <- "pass"
  if (!isTRUE(keep_all)) {
    out <- out[out$pass_annotation_filter %in% TRUE, , drop = FALSE]
  }
  rownames(out) <- NULL
  out
}

`77_update_annotation_table` <- 
function(object, annotated) {
  annotated <- as.data.frame(annotated, stringsAsFactors = FALSE, check.names = FALSE)
  final_cols <- colnames(.ann_empty_final_annotation())
  annotated <- .ann_add_missing_columns(annotated, final_cols)
  final_annotation <- annotated[, final_cols, drop = FALSE]
  object@annotation_results$final_annotation <- final_annotation
  object@annotation_table <- final_annotation
  object
}

`78_store_process_info` <- 
function(object, parameters) {
  parameters$timestamp <- as.character(Sys.time())
  object@msi_parameters$annotate_compounds <- parameters
  object@process_info$module7_final_annotation <- c(parameters, list(
    n_final_annotations = if (is.data.frame(object@annotation_results$final_annotation)) nrow(object@annotation_results$final_annotation) else NA_integer_
  ))
  object
}

`annotate_features_by_isotopes` <- 
function(object,
                                          fmc_ms1_candidate_table,
                                          ppm = 10,
                                          max_isotope = 2,
                                          isotope_weight = 10,
                                          use_spatial_support = TRUE) {
  features <- .fmc_feature_table(object)
  parents <- .fmc_as_df(fmc_ms1_candidate_table)
  empty <- .fmc_empty(c("feature_id", "parent_feature_id", "feature_mz",
                       "parent_feature_mz", "metabolite_id", "metabolite_name",
                       "isotope_label", "theoretical_mz", "observed_mz",
                       "mz_error_ppm", "evidence_type", "evidence_weight",
                       "spatial_sas_to_parent", "data_source"))
  if (!nrow(features) || !nrow(parents)) return(empty)
  parents <- parents[is.finite(suppressWarnings(as.numeric(parents$feature_mz))) &
                       !is.na(parents$feature_id) & nzchar(as.character(parents$feature_id)), , drop = FALSE]
  if (!nrow(parents)) return(empty)
  ppm <- suppressWarnings(as.numeric(ppm[[1]]))
  max_isotope <- max(1L, min(2L, suppressWarnings(as.integer(max_isotope[[1]]))))
  feat <- features[is.finite(features$feature_mz) & features$feature_mz > 0, , drop = FALSE]
  feat <- feat[order(feat$feature_mz), , drop = FALSE]
  mz_sorted <- feat$feature_mz
  rows <- list()
  tol <- ppm / 1e6
  for (i in seq_len(nrow(parents))) {
    parent_mz <- suppressWarnings(as.numeric(parents$feature_mz[[i]]))
    charge <- if ("charge" %in% names(parents)) {
      suppressWarnings(as.numeric(parents$charge[[i]]))
    } else {
      NA_real_
    }
    if (!length(charge) || !is.finite(charge[[1]]) || charge[[1]] == 0) charge <- 1
    for (iso in seq_len(max_isotope)) {
      spacing <- iso * 1.003355 / abs(charge)
      theoretical <- parent_mz + spacing
      lower <- theoretical / (1 + tol)
      upper <- theoretical / (1 - tol)
      left <- findInterval(lower - 1e-12, mz_sorted) + 1L
      right <- findInterval(upper + 1e-12, mz_sorted)
      if (left <= right && left >= 1L && right <= length(mz_sorted)) {
        hit <- feat[seq.int(left, right), , drop = FALSE]
        hit <- hit[as.character(hit$feature_id) != as.character(parents$feature_id[[i]]), , drop = FALSE]
        if (!nrow(hit)) next
        err <- .fmc_ppm_error(hit$feature_mz, theoretical)
        keep <- abs(err) <= ppm
        if (!any(keep)) next
        hit <- hit[keep, , drop = FALSE]
        err <- err[keep]
        for (j in seq_len(nrow(hit))) {
          rows[[length(rows) + 1L]] <- data.frame(
            feature_id = as.character(hit$feature_id[[j]]),
            parent_feature_id = as.character(parents$feature_id[[i]]),
            feature_mz = as.numeric(hit$feature_mz[[j]]),
            parent_feature_mz = parent_mz,
            metabolite_id = as.character(parents$metabolite_id[[i]]),
            metabolite_name = as.character(parents$metabolite_name[[i]]),
            isotope_label = paste0("M+", iso),
            theoretical_mz = theoretical,
            observed_mz = as.numeric(hit$feature_mz[[j]]),
            mz_error_ppm = err[[j]],
            evidence_type = "isotope",
            evidence_weight = as.numeric(isotope_weight[[1]]),
            spatial_sas_to_parent = if (isTRUE(use_spatial_support)) {
              .fmc_spatial_sas(object, hit$feature_id[[j]], parents$feature_id[[i]])
            } else {
              NA_real_
            },
            data_source = "isotope_spacing_match",
            stringsAsFactors = FALSE,
            check.names = FALSE
          )
        }
      }
    }
  }
  if (!length(rows)) return(empty)
  out <- do.call(rbind, rows)
  out <- out[order(as.character(out$feature_id), as.character(out$metabolite_id),
                   abs(out$mz_error_ppm)), , drop = FALSE]
  rownames(out) <- NULL
  out[, colnames(empty), drop = FALSE]
}

`annotate_features_by_neutral_loss` <- 
function(object,
                                              metabolite_reference,
                                              ppm = 10,
                                              ion_mode,
                                              neutral_loss_rules = default_neutral_loss_rule_table()) {
  features <- .fmc_feature_table(object)
  ref <- .fmc_as_df(metabolite_reference)
  rules <- .fmc_as_df(neutral_loss_rules)
  empty <- .fmc_empty(c("feature_id", "feature_mz", "metabolite_id", "metabolite_name",
                       "neutral_loss", "theoretical_mz", "mz_error_ppm",
                       "evidence_type", "evidence_weight", "data_source"))
  if (!nrow(features) || !nrow(ref) || !nrow(rules) || !"enabled" %in% names(rules)) return(empty)
  rules <- rules[rules$enabled %in% TRUE, , drop = FALSE]
  if (!nrow(rules)) return(empty)
  adducts <- get_adduct_rules_for_mode(ion_mode)
  adducts <- adducts[adducts$evidence_type == "ms1_adduct" & adducts$rank_group %in% "canonical", , drop = FALSE]
  if (!nrow(adducts)) return(empty)
  theo_rows <- list()
  for (a in seq_len(nrow(adducts))) {
    for (r in seq_len(nrow(rules))) {
      theoretical <- (suppressWarnings(as.numeric(ref$neutral_mass)) -
                        as.numeric(rules$mass_loss[[r]]) +
                        as.numeric(adducts$mass_shift[[a]])) / abs(as.numeric(adducts$charge[[a]]))
      theo_rows[[length(theo_rows) + 1L]] <- data.frame(
        metabolite_id = ref$metabolite_id,
        metabolite_name = ref$metabolite_name,
        neutral_loss = as.character(rules$neutral_loss[[r]]),
        theoretical_mz = theoretical,
        evidence_type = "neutral_loss",
        evidence_weight = as.numeric(rules$evidence_weight[[r]]),
        stringsAsFactors = FALSE,
        check.names = FALSE
      )
    }
  }
  theo <- do.call(rbind, theo_rows)
  theo <- theo[is.finite(theo$theoretical_mz) & theo$theoretical_mz > 0, , drop = FALSE]
  if (!nrow(theo)) return(empty)
  theo <- theo[order(theo$theoretical_mz), , drop = FALSE]
  mz_sorted <- theo$theoretical_mz
  tol <- as.numeric(ppm[[1]]) / 1e6
  out_rows <- list()
  for (i in seq_len(nrow(features))) {
    fmz <- suppressWarnings(as.numeric(features$feature_mz[[i]]))
    if (!is.finite(fmz)) next
    left <- findInterval(fmz / (1 + tol) - 1e-12, mz_sorted) + 1L
    right <- findInterval(fmz / (1 - tol) + 1e-12, mz_sorted)
    if (left <= right && left >= 1L && right <= length(mz_sorted)) {
      hit <- theo[seq.int(left, right), , drop = FALSE]
      err <- .fmc_ppm_error(fmz, hit$theoretical_mz)
      keep <- abs(err) <= ppm
      if (!any(keep)) next
      hit <- hit[keep, , drop = FALSE]
      err <- err[keep]
      out_rows[[length(out_rows) + 1L]] <- data.frame(
        feature_id = features$feature_id[[i]],
        feature_mz = fmz,
        hit,
        mz_error_ppm = err,
        data_source = "neutral_loss_match",
        stringsAsFactors = FALSE,
        check.names = FALSE
      )
    }
  }
  if (!length(out_rows)) return(empty)
  out <- do.call(rbind, out_rows)
  out <- out[order(as.character(out$feature_id), -out$evidence_weight,
                   abs(out$mz_error_ppm)), , drop = FALSE]
  rownames(out) <- NULL
  out[, colnames(empty), drop = FALSE]
}

`apply_fmc_redundancy_filter` <- 
function(fmc_candidate_table,
                                        max_metabolites_per_feature = 20,
                                        max_features_per_metabolite_soft = 30,
                                        max_features_per_metabolite_hard = 100) {
  x <- .fmc_as_df(fmc_candidate_table)
  required <- c("feature_id", "metabolite_id", "score_total", "evidence_types")
  missing <- setdiff(required, names(x))
  if (length(missing)) {
    stop("fmc_candidate_table missing required columns: ", paste(missing, collapse = ", "), call. = FALSE)
  }
  if (!nrow(x)) {
    attr(x, "before_redundancy_table") <- x
    attr(x, "fmc_redundancy_summary") <- data.frame(
      n_candidate_pairs_before_redundancy = 0L,
      n_candidate_pairs_after_redundancy = 0L,
      n_promiscuous_metabolites = 0L,
      n_ambiguous_features = 0L,
      stringsAsFactors = FALSE,
      check.names = FALSE
    )
    return(x)
  }
  max_metabolites_per_feature <- suppressWarnings(as.integer(max_metabolites_per_feature[[1]]))
  if (!is.finite(max_metabolites_per_feature) || max_metabolites_per_feature < 1L) {
    max_metabolites_per_feature <- 20L
  }
  feature_counts <- table(x$feature_id)
  metabolite_counts <- tapply(as.character(x$feature_id), as.character(x$metabolite_id),
                              function(z) length(unique(z)))
  x$n_metabolites_per_feature <- as.integer(feature_counts[as.character(x$feature_id)])
  x$n_features_per_metabolite <- as.integer(metabolite_counts[as.character(x$metabolite_id)])
  x$ambiguous_feature <- x$n_metabolites_per_feature > max_metabolites_per_feature
  x$promiscuous_metabolite <- x$n_features_per_metabolite > max_features_per_metabolite_soft
  penalty <- rep(0, nrow(x))
  penalty[x$n_features_per_metabolite > 30 & x$n_features_per_metabolite <= 50] <- -10
  penalty[x$n_features_per_metabolite > 50 & x$n_features_per_metabolite <= 100] <- -20
  penalty[x$n_features_per_metabolite > 100] <- -30
  if (is.finite(max_features_per_metabolite_hard)) {
    penalty[x$n_features_per_metabolite > max_features_per_metabolite_hard] <- pmin(
      penalty[x$n_features_per_metabolite > max_features_per_metabolite_hard],
      -30
    )
  }
  x$redundancy_penalty <- penalty
  x$score_total_after_redundancy <- suppressWarnings(as.numeric(x$score_total)) + penalty
  x$rank_within_feature_after_redundancy <- .fmc_rank(x, "feature_id", "score_total_after_redundancy")
  x$rank_within_metabolite_after_redundancy <- .fmc_rank(x, "metabolite_id", "score_total_after_redundancy")

  before <- x
  has_pseudo_ms2 <- grepl("pseudo_ms2", as.character(x$evidence_types), fixed = TRUE)
  keep <- x$rank_within_feature_after_redundancy <= max_metabolites_per_feature |
    has_pseudo_ms2
  x <- x[keep, , drop = FALSE]
  x$primary_metabolite_feature <- FALSE
  x$supporting_evidence_feature <- FALSE
  metab_groups <- split(seq_len(nrow(x)), as.character(x$metabolite_id), drop = TRUE)
  for (idx in metab_groups) {
    pseudo <- grepl("pseudo_ms2", as.character(x$evidence_types[idx]), fixed = TRUE)
    ord <- order(!pseudo, -x$score_total_after_redundancy[idx], idx)
    primary <- idx[ord[[1]]]
    x$primary_metabolite_feature[primary] <- TRUE
    x$supporting_evidence_feature[idx[idx != primary]] <- TRUE
  }
  x$rank_within_feature_after_redundancy <- .fmc_rank(x, "feature_id", "score_total_after_redundancy")
  x$rank_within_metabolite_after_redundancy <- .fmc_rank(x, "metabolite_id", "score_total_after_redundancy")
  x <- x[order(as.character(x$feature_id), x$rank_within_feature_after_redundancy,
               as.character(x$metabolite_id)), , drop = FALSE]
  rownames(x) <- NULL
  summary <- data.frame(
    n_candidate_pairs_before_redundancy = as.integer(nrow(before)),
    n_candidate_pairs_after_redundancy = as.integer(nrow(x)),
    n_promiscuous_metabolites = as.integer(length(unique(before$metabolite_id[before$promiscuous_metabolite]))),
    n_ambiguous_features = as.integer(length(unique(before$feature_id[before$ambiguous_feature]))),
    max_metabolites_per_feature = as.integer(max_metabolites_per_feature),
    max_features_per_metabolite_soft = as.integer(max_features_per_metabolite_soft),
    max_features_per_metabolite_hard = as.integer(max_features_per_metabolite_hard),
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
  attr(x, "before_redundancy_table") <- before
  attr(x, "fmc_redundancy_summary") <- summary
  x
}

`build_fmc_metabolite_reference` <- 
function(library_object = NULL,
                                           library_table = NULL,
                                           database_name = NULL) {
  source_table <- NULL
  if (is.data.frame(library_table)) {
    source_table <- library_table
  } else if (is.list(library_object)) {
    candidates <- c("metabolite_reference", "ms1_database", "library_spectrum_table",
                    "spectrum_table", "spectra", "metadata")
    for (nm in candidates) {
      if (is.data.frame(library_object[[nm]])) {
        source_table <- library_object[[nm]]
        break
      }
    }
  }
  if (is.null(source_table)) {
    source_table <- data.frame(stringsAsFactors = FALSE, check.names = FALSE)
  }
  tbl <- .fmc_as_df(source_table)
  n_input <- nrow(tbl)
  empty <- .fmc_empty(c("metabolite_id", "metabolite_name", "formula", "neutral_mass",
                       "database", "precursor_mz", "adduct", "ion_mode", "class",
                       "lipid_class", "source_spectrum_id"))
  if (!nrow(tbl)) {
    attr(empty, "fmc_reference_diagnostics") <- data.frame(
      n_input_library_records = 0L,
      n_records_with_neutral_mass = 0L,
      n_records_without_neutral_mass = 0L,
      n_unique_metabolites = 0L,
      database_name = as.character(database_name %||% NA_character_),
      stringsAsFactors = FALSE,
      check.names = FALSE
    )
    return(empty)
  }
  metabolite_id <- .fmc_chr(tbl, c("metabolite_id", "candidate_id", "compound_id",
                                  "hmdb_id", "id", "accession"))
  missing_id <- is.na(metabolite_id) | !nzchar(metabolite_id)
  metabolite_id[missing_id] <- paste0("metabolite_", which(missing_id))
  metabolite_name <- .fmc_chr(tbl, c("metabolite_name", "candidate_name", "compound_name",
                                    "name", "title"), default = NA_character_)
  formula <- .fmc_chr(tbl, c("formula", "molecular_formula"), default = NA_character_)
  neutral_mass <- .fmc_num(tbl, c("neutral_mass", "exact_mass", "monoisotopic_mass",
                                 "monoisotopic_molecular_weight", "exactmass", "mass"),
                           default = NA_real_)
  precursor_mz <- .fmc_num(tbl, c("precursor_mz", "observed_precursor_mz", "parent_mz",
                                 "pepmass"), default = NA_real_)
  adduct <- .fmc_chr(tbl, c("adduct", "precursor_type", "ion_type"), default = NA_character_)
  ion_mode <- tolower(.fmc_chr(tbl, c("ion_mode", "polarity", "mode"), default = NA_character_))
  ion_mode[ion_mode %in% c("pos", "+", "positive ion mode")] <- "positive"
  ion_mode[ion_mode %in% c("neg", "-", "negative ion mode")] <- "negative"
  inferred <- !is.finite(neutral_mass) & is.finite(precursor_mz) & !is.na(adduct) & nzchar(adduct)
  if (any(inferred)) {
    neutral_mass[inferred] <- vapply(which(inferred), function(i) {
      .fmc_infer_neutral_mass(precursor_mz[[i]], adduct[[i]], ion_mode[[i]])
    }, numeric(1))
  }
  database <- .fmc_chr(tbl, c("database", "database_source", "library_source", "source"),
                       default = as.character(database_name %||% NA_character_))
  class_col <- .fmc_chr(tbl, c("class", "superclass", "category"), default = NA_character_)
  lipid_class <- .fmc_chr(tbl, c("lipid_class", "main_class", "subclass"), default = NA_character_)
  spectrum_id <- .fmc_chr(tbl, c("source_spectrum_id", "spectrum_id", "scan", "reference_id"),
                          default = NA_character_)
  out <- data.frame(
    metabolite_id = metabolite_id,
    metabolite_name = metabolite_name,
    formula = formula,
    neutral_mass = neutral_mass,
    database = database,
    precursor_mz = precursor_mz,
    adduct = adduct,
    ion_mode = ion_mode,
    class = class_col,
    lipid_class = lipid_class,
    source_spectrum_id = spectrum_id,
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
  with_mass <- is.finite(out$neutral_mass) & out$neutral_mass > 0
  dropped_without_mass <- sum(!with_mass)
  out <- out[with_mass, , drop = FALSE]
  if (nrow(out)) {
    key <- paste(out$metabolite_id, signif(out$neutral_mass, 10), out$adduct, out$source_spectrum_id, sep = "\r")
    out <- out[!duplicated(key), , drop = FALSE]
  }
  rownames(out) <- NULL
  attr(out, "fmc_reference_diagnostics") <- data.frame(
    n_input_library_records = as.integer(n_input),
    n_records_with_neutral_mass = as.integer(nrow(out)),
    n_records_without_neutral_mass = as.integer(dropped_without_mass),
    n_unique_metabolites = as.integer(length(unique(out$metabolite_id))),
    database_name = as.character(database_name %||% NA_character_),
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
  out
}

`build_mfc_candidate_table` <- 
function(annotation_table = NULL,
                                      ms1_only_match_table = NULL,
                                      pseudo_ms2_match_table = NULL,
                                      feature_metadata = NULL,
                                      expression_matrix = NULL) {
  rows <- list(
    annotation = .mfc2_standardize_source(annotation_table, "annotation"),
    ms1_only = .mfc2_standardize_source(ms1_only_match_table, "ms1_only"),
    pseudo_ms2 = .mfc2_standardize_source(pseudo_ms2_match_table, "pseudo_ms2")
  )
  rows <- rows[vapply(rows, nrow, integer(1)) > 0L]
  if (!length(rows)) return(.mfc2_empty_candidate())
  out <- do.call(rbind, rows)
  feature_metadata <- .mfc2_df(feature_metadata)
  if (nrow(feature_metadata)) {
    fid <- .mfc2_chr(feature_metadata, c("feature_id", "variable_id"))
    fmz <- .mfc2_num(feature_metadata, c("feature_mz", "mz", "observed_mz"))
    lookup <- stats::setNames(fmz, fid)
    missing_mz <- !is.finite(out$feature_mz)
    out$feature_mz[missing_mz] <- lookup[out$feature_id[missing_mz]]
  }
  feature_axis <- .mfc2_feature_axis(feature_metadata = feature_metadata,
                                     expression_matrix = expression_matrix)
  out <- .mfc2_add_candidate_keys(out, feature_axis = feature_axis)
  out$evidence_type <- ifelse(
    grepl("pseudo", out$evidence_type, ignore.case = TRUE) &
      !grepl("strict", out$evidence_type, ignore.case = TRUE),
    "pseudo_ms2",
    ifelse(grepl("ms1", out$evidence_type, ignore.case = TRUE), "ms1_only", out$evidence_type)
  )
  missing_rank <- !is.finite(out$candidate_rank)
  if (any(missing_rank)) {
    out$candidate_rank[missing_rank] <- ave(
      -out$raw_score[missing_rank],
      out$feature_id[missing_rank],
      FUN = function(x) rank(x, ties.method = "first")
    )
  }
  rownames(out) <- NULL
  out
}

`build_mfc_evidence_long` <- 
function(pseudo_ms2_pass_table = NULL,
                                    pseudo_ms2_all_table = NULL,
                                    ms1_match_table = NULL,
                                    ion_mode,
                                    min_matched_peaks = 3L,
                                    min_reverse_masked_cosine = 0,
                                    pseudo_library_name = NA_character_,
                                    ms1_library_name = NA_character_,
                                    sample_mode = "standard_mixture",
                                    pure_tissue_feature_table = NULL) {
  ion_mode <- tolower(as.character(ion_mode[[1]]))
  if (!ion_mode %in% c("positive", "negative")) {
    stop("ion_mode must resolve to positive or negative before MFC evidence construction.", call. = FALSE)
  }
  min_matched_peaks <- as.integer(min_matched_peaks[[1]])
  min_reverse_masked_cosine <- as.numeric(min_reverse_masked_cosine[[1]])

  pseudo_pass <- .mfc2_standardize_source(pseudo_ms2_pass_table, "pseudo_ms2")
  pseudo_all <- .mfc2_standardize_source(pseudo_ms2_all_table, "pseudo_ms2")
  if (!nrow(pseudo_all)) pseudo_all <- pseudo_pass
  pass_key <- unique(.mfc_unified_evidence_key(pseudo_pass))

  pseudo_long <- if (nrow(pseudo_all)) {
    pass <- .mfc_unified_evidence_key(pseudo_all) %in% pass_key
    raw_score <- suppressWarnings(as.numeric(pseudo_all$raw_evidence_score))
    derived <- .mfc_unified_pseudo_score(pseudo_all)
    raw_score[!is.finite(raw_score)] <- derived[!is.finite(raw_score)]
    rejection <- .mfc_unified_rejection_reason(
      pseudo_all,
      min_matched_peaks = min_matched_peaks,
      min_reverse_masked_cosine = min_reverse_masked_cosine
    )
    rejection[pass] <- NA_character_
    data.frame(
      feature_id = as.character(pseudo_all$feature_id),
      source_row_id = as.character(pseudo_all$source_row_id),
      candidate_key = as.character(pseudo_all$candidate_key),
      dedup_key = as.character(pseudo_all$dedup_key),
      feature_mz = suppressWarnings(as.numeric(pseudo_all$feature_mz)),
      HMDB_ID = as.character(pseudo_all$hmdb_id),
      compound_name = as.character(pseudo_all$compound_name),
      adduct = as.character(pseudo_all$adduct),
      ion_mode = ion_mode,
      evidence_source = "pseudo_ms2",
      library_name = ifelse(
        is.na(pseudo_library_name) | !nzchar(as.character(pseudo_library_name)),
        as.character(pseudo_all$source),
        as.character(pseudo_library_name)
      ),
      library_spectrum_id = as.character(pseudo_all$library_spectrum_id),
      pseudo_anchor = as.character(pseudo_all$anchor_feature_id),
      reverse_masked_cosine = suppressWarnings(as.numeric(pseudo_all$reverse_masked_cosine)),
      n_matched_peaks = suppressWarnings(as.integer(pseudo_all$n_matched_peaks)),
      matched_library_intensity_fraction = suppressWarnings(as.numeric(pseudo_all$matched_library_intensity_fraction)),
      S_SAS = suppressWarnings(as.numeric(pseudo_all$S_SAS)),
      step7_q_value = suppressWarnings(as.numeric(pseudo_all$q_value)),
      ms1_ppm_error = NA_real_,
      ms1_ppm_score = NA_real_,
      raw_evidence_score = .mfc_unified_clamp01(raw_score),
      evidence_pass = pass,
      rejection_reason = rejection,
      stringsAsFactors = FALSE,
      check.names = FALSE
    )
  } else {
    data.frame()
  }

  # Preserve the source table separately: `.mfc2_standardize_source()`
  # deliberately projects its historic candidate schema and therefore does
  # not carry pure-tissue provenance columns.  Those columns are only used in
  # the pure-tissue branch below and never alter standard-mixture evidence.
  ms1_source <- as.data.frame(ms1_match_table %||% data.frame(), stringsAsFactors = FALSE,
    check.names = FALSE)
  ms1 <- .mfc2_standardize_source(ms1_match_table, "ms1_only")
  ms1_long <- if (nrow(ms1)) {
    ppm_score <- suppressWarnings(as.numeric(ms1$score_ms1))
    raw_score <- suppressWarnings(as.numeric(ms1$raw_evidence_score))
    ppm_score[!is.finite(ppm_score)] <- raw_score[!is.finite(ppm_score)]
    ppm_score <- .mfc_unified_clamp01(ppm_score)
    raw_score[!is.finite(raw_score)] <- ppm_score[!is.finite(raw_score)]
    data.frame(
      feature_id = as.character(ms1$feature_id),
      source_row_id = as.character(ms1$source_row_id),
      candidate_key = as.character(ms1$candidate_key),
      dedup_key = as.character(ms1$dedup_key),
      feature_mz = suppressWarnings(as.numeric(ms1$feature_mz)),
      HMDB_ID = as.character(ms1$hmdb_id),
      compound_name = as.character(ms1$compound_name),
      adduct = as.character(ms1$adduct),
      ion_mode = ion_mode,
      evidence_source = "ms1",
      library_name = ifelse(
        is.na(ms1_library_name) | !nzchar(as.character(ms1_library_name)),
        as.character(ms1$source),
        as.character(ms1_library_name)
      ),
      library_spectrum_id = as.character(ms1$library_spectrum_id),
      pseudo_anchor = NA_character_,
      reverse_masked_cosine = NA_real_,
      n_matched_peaks = NA_integer_,
      matched_library_intensity_fraction = NA_real_,
      S_SAS = suppressWarnings(as.numeric(ms1$S_SAS)),
      step7_q_value = suppressWarnings(as.numeric(ms1$q_value)),
      ms1_ppm_error = suppressWarnings(as.numeric(ms1$mz_error_ppm)),
      ms1_ppm_score = ppm_score,
      raw_evidence_score = .mfc_unified_clamp01(raw_score),
      evidence_pass = TRUE,
      rejection_reason = NA_character_,
      stringsAsFactors = FALSE,
      check.names = FALSE
    )
  } else {
    data.frame()
  }

  template <- data.frame(
    feature_id = character(), source_row_id = character(),
    candidate_key = character(), dedup_key = character(),
    feature_mz = numeric(), HMDB_ID = character(),
    compound_name = character(), adduct = character(), ion_mode = character(),
    evidence_source = character(), library_name = character(),
    library_spectrum_id = character(), pseudo_anchor = character(),
    reverse_masked_cosine = numeric(), n_matched_peaks = integer(),
    matched_library_intensity_fraction = numeric(), S_SAS = numeric(),
    step7_q_value = numeric(), ms1_ppm_error = numeric(), ms1_ppm_score = numeric(),
    raw_evidence_score = numeric(), evidence_pass = logical(),
    rejection_reason = character(), stringsAsFactors = FALSE, check.names = FALSE
  )
  pure_tissue <- identical(tolower(as.character(sample_mode[[1]])), "pure_tissue")
  if (isTRUE(pure_tissue)) {
    template$ms1_feature_mz <- numeric()
    template$ms1_mz_source <- character()
    template$ms1_mz_provenance <- character()
    template$pseudo_anchor_mz <- numeric()
    template$precursor_mz <- numeric()
    template$precursor_source <- character()
  }
  rows <- list(pseudo_long, ms1_long)
  rows <- rows[vapply(rows, nrow, integer(1)) > 0L]
  if (!length(rows)) return(template)
  out <- do.call(rbind, rows)
  out$HMDB_ID <- .mfc_normalize_hmdb_id(out$HMDB_ID)
  out <- out[nzchar(out$feature_id) & nzchar(out$HMDB_ID), , drop = FALSE]
  out <- out[order(out$feature_id, out$HMDB_ID, out$evidence_source,
                   -out$raw_evidence_score, out$library_spectrum_id), , drop = FALSE]
  if (isTRUE(pure_tissue)) {
    # Pure tissue has no authentic precursor metadata.  The pseudo-MS2 anchor
    # remains explicitly labelled as such; it is not copied into precursor_mz.
    out$ms1_feature_mz <- NA_real_
    out$ms1_mz_source <- NA_character_
    out$ms1_mz_provenance <- NA_character_
    # The measured MS1 axis is a property of every selected feature, not only
    # of features that happened to obtain an MS1-library candidate.  In
    # particular, pseudo-MS2-only evidence must retain the same Cardinal axis
    # provenance without converting its pseudo anchor into a precursor.
    feature_source <- as.data.frame(pure_tissue_feature_table %||% data.frame(),
      stringsAsFactors = FALSE, check.names = FALSE)
    required_feature_source <- c("feature_id", "ms1_feature_mz", "ms1_mz_source", "ms1_mz_provenance")
    if (!nrow(feature_source) || !all(required_feature_source %in% names(feature_source))) {
      stop("HARD_PROVENANCE_BLOCKER_PRECURSOR_OR_MS1_SEMANTICS_INVALID: pure_tissue requires a complete measured MS1 feature table.",
        call. = FALSE)
    }
    source_idx <- match(out$feature_id, as.character(feature_source$feature_id))
    source_rows <- !is.na(source_idx)
    if (!all(source_rows)) {
      stop("HARD_PROVENANCE_BLOCKER_PRECURSOR_OR_MS1_SEMANTICS_INVALID: pure_tissue evidence feature IDs are absent from the measured MS1 feature table.",
        call. = FALSE)
    }
    source_mz <- suppressWarnings(as.numeric(feature_source$ms1_feature_mz[source_idx]))
    source_origin <- as.character(feature_source$ms1_mz_source[source_idx])
    source_provenance <- as.character(feature_source$ms1_mz_provenance[source_idx])
    if (any(!is.finite(source_mz) | source_mz <= 0) ||
        any(source_origin != "cardinal_measured_feature_axis") ||
        any(is.na(source_provenance) | !nzchar(source_provenance))) {
      stop("HARD_PROVENANCE_BLOCKER_PRECURSOR_OR_MS1_SEMANTICS_INVALID: measured MS1 feature provenance is incomplete or not Cardinal-derived.",
        call. = FALSE)
    }
    out$ms1_feature_mz[source_rows] <- source_mz
    out$ms1_mz_source[source_rows] <- source_origin
    out$ms1_mz_provenance[source_rows] <- source_provenance
    out$pseudo_anchor_mz <- ifelse(
      out$evidence_source == "pseudo_ms2", suppressWarnings(as.numeric(out$feature_mz)), NA_real_
    )
    out$precursor_mz <- NA_real_
    out$precursor_source <- "missing_explicit_precursor"
  }
  rownames(out) <- NULL
  out[, names(template), drop = FALSE]
}

`.legacy_calculate_mfc` <- 
function(object,
                          fmc_config = default_fmc_annotation_config(TRUE),
                          metabolite_reference = NULL,
                          library_object = NULL,
                          ion_mode = "negative",
                          ...) {
  if (!methods::is(object, "msi_dataset")) {
    stop("object must be a msi_dataset.", call. = FALSE)
  }
  final_before <- object@annotation_results$final_annotation
  annotation_before <- object@annotation_table
  object <- run_fmc_annotation(
    object = object,
    fmc_config = fmc_config,
    metabolite_reference = metabolite_reference,
    library_object = library_object,
    ion_mode = ion_mode,
    ...
  )
  object@annotation_results$final_annotation <- final_before
  object@annotation_table <- annotation_before
  if (!is.null(object@process_info$module7_fmc_annotation)) {
    object@process_info$module8_mfc <- object@process_info$module7_fmc_annotation
    object@process_info$module8_mfc$module <- "module8_mfc"
    object@process_info$module8_mfc$timestamp <- as.character(Sys.time())
    object@process_info$module7_fmc_annotation <- NULL
  } else if (is.null(object@process_info$module8_mfc)) {
    object@process_info$module8_mfc <- list(
      fmc_enabled = isTRUE((fmc_config %||% list())$enabled),
      module = "module8_mfc",
      timestamp = as.character(Sys.time())
    )
  } else {
    object@process_info$module8_mfc$module <- "module8_mfc"
    object@process_info$module8_mfc$timestamp <- as.character(Sys.time())
  }
  object
}

`calculate_S_SAS` <- 
function(matched_peak_pairs) {
  empty <- data.frame(
    cluster_id = character(),
    anchor_feature_id = character(),
    candidate_key = character(),
    candidate_id = character(),
    spectrum_id = character(),
    S_SAS_raw = numeric(),
    S_SAS = numeric(),
    S_SAS_status = character(),
    n_SAS_matched_peaks = integer(),
    sum_matched_library_intensity = numeric(),
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
  pairs <- as.data.frame(matched_peak_pairs, stringsAsFactors = FALSE, check.names = FALSE)
  if (nrow(pairs) == 0L) {
    return(empty)
  }
  pairs <- .ann_add_missing_columns(pairs, c(
    "cluster_id", "anchor_feature_id", "candidate_key", "candidate_id",
    "spectrum_id", "library_intensity", "SAS", "matched"
  ))
  key <- .ann_split_key(pairs, c("cluster_id", "anchor_feature_id", "candidate_key", "candidate_id", "spectrum_id"))
  rows <- lapply(split(pairs, key), function(df) {
    matched <- df$matched %in% TRUE
    w <- suppressWarnings(as.numeric(df$library_intensity))[matched]
    sas <- suppressWarnings(as.numeric(df$SAS))[matched]
    keep <- is.finite(w) & w > 0 & is.finite(sas)
    sum_w <- sum(w[is.finite(w) & w > 0], na.rm = TRUE)
    S_SAS_raw <- if (!any(matched) || !is.finite(sum_w) || sum_w <= 0 || !any(keep)) {
      NA_real_
    } else {
      sum(w[keep] * sas[keep], na.rm = TRUE) / sum(w[keep], na.rm = TRUE)
    }
    S_SAS_status <- if (is.finite(S_SAS_raw)) "ok" else "missing_spatial_support"
    S_SAS <- if (is.finite(S_SAS_raw)) S_SAS_raw else 0
    data.frame(
      cluster_id = df$cluster_id[[1]],
      anchor_feature_id = as.character(df$anchor_feature_id[[1]]),
      candidate_key = as.character(df$candidate_key[[1]]),
      candidate_id = as.character(df$candidate_id[[1]]),
      spectrum_id = as.character(df$spectrum_id[[1]]),
      S_SAS_raw = as.numeric(S_SAS_raw),
      S_SAS = as.numeric(S_SAS),
      S_SAS_status = S_SAS_status,
      n_SAS_matched_peaks = as.integer(sum(matched, na.rm = TRUE)),
      sum_matched_library_intensity = as.numeric(sum_w),
      stringsAsFactors = FALSE,
      check.names = FALSE
    )
  })
  out <- do.call(rbind, rows)
  out <- out[, colnames(empty), drop = FALSE]
  rownames(out) <- NULL
  out
}

`calculate_score` <- 
function(
  object,
  weights = c(
    S_reverse_cos = 0.35,
    S_lib_coverage = 0.25,
    S_spatial_reverse = 0.20,
    S_FE = 0.10,
    S_precursor = 0.10,
    S_chem = 0.00
  ),
  decoy_strategy = c("peak_shift", "intensity_permutation", "none"),
  fdr_method = c("target_decoy", "none"),
  q_value_cutoff = 0.05,
  min_final_score = 0,
  keep_all = TRUE,
  decoy_n = 1,
  decoy_shift_da = 0.05,
  seed = 123,
  confidence_high_reverse_cos = 0.10,
  confidence_high_reverse_masked_cosine = 0.05,
  confidence_high_lib_coverage = 0.10,
  confidence_high_precursor_ppm = 50,
  confidence_medium_reverse_cos = 0.05,
  confidence_medium_reverse_masked_cosine = 0.01,
  confidence_medium_lib_coverage = 0.05,
  confidence_medium_matched_fraction = 0.02,
  confidence_medium_precursor_ppm = 100,
  reverse_masked_cosine_weight = 0.7,
  S_SAS_weight = 0.3,
  candidate_source_for_step7 = c("reverse_match_filtered", "reverse_match_all"),
  S_SAS_threshold = NA_real_,
  decoy_fragment_tolerance_da = NULL,
  primary_spectral_score = "reverse_masked_cosine",
  output_dir = "output",
  fmc_annotation = NULL,
  fmc_enabled = NULL,
  fmc_ms1_ppm = NULL,
  fmc_ion_mode = NULL,
  fmc_max_candidates_per_feature = NULL,
  fmc_max_metabolites_per_feature = NULL,
  fmc_max_features_per_metabolite_soft = NULL,
  fmc_max_features_per_metabolite_hard = NULL,
  fmc_use_isotope = NULL,
  fmc_use_neutral_loss = NULL,
  fmc_use_spatial_support = NULL,
  fmc_metabolite_reference = NULL,
  ...
) {
  decoy_strategy <- match.arg(decoy_strategy)
  fdr_method <- match.arg(fdr_method)
  candidate_source_for_step7 <- match.arg(candidate_source_for_step7)
  primary_spectral_score <- as.character(primary_spectral_score[[1]])
  if (!identical(primary_spectral_score, "reverse_masked_cosine")) {
    warning(
      "primary_spectral_score is legacy; reverse_masked_cosine is used for new scoring.",
      call. = FALSE
    )
    primary_spectral_score <- "reverse_masked_cosine"
  }

  if (!methods::is(object, "msi_dataset")) {
    stop("object must be a msi_dataset.", call. = FALSE)
  }
  if (is.null(object@pseudo_ms2$pseudo_ms2_table)) {
    .msi_missing_summary("7_annotate_compounds", "pseudo_ms2$pseudo_ms2_table", output_dir)
    return(object)
  }
  if (is.null(object@annotation_results$flash_entropy_hits)) {
    .msi_missing_summary("7_annotate_compounds", "annotation_results$flash_entropy_hits", output_dir)
    return(object)
  }
  candidate_table <- if (identical(candidate_source_for_step7, "reverse_match_all")) {
    object@annotation_results$reverse_match_all
  } else {
    object@annotation_results$reverse_match
  }
  if (is.null(candidate_table)) {
    .msi_missing_summary("7_annotate_compounds", paste0(
      "annotation_results$",
      if (identical(candidate_source_for_step7, "reverse_match_all")) "reverse_match_all" else "reverse_match"
    ), output_dir)
    return(object)
  }
  candidate_table <- as.data.frame(candidate_table, stringsAsFactors = FALSE, check.names = FALSE)
  if (identical(candidate_source_for_step7, "reverse_match_all")) {
    candidate_table <- `60_mark_chemical_evidence_not_evaluated`(candidate_table)
  }
  candidate_table <- .ann_add_missing_columns(candidate_table, c("reverse_masked_cosine", "n_matched_peaks"))
  candidate_table <- .ann_ensure_evidence_pass(candidate_table, object)
  candidate_table$pass_reverse_masked_cosine <-
    suppressWarnings(as.numeric(candidate_table$reverse_masked_cosine)) > 0
  candidate_table$pass_reverse_masked_cosine[is.na(candidate_table$pass_reverse_masked_cosine)] <- FALSE
  candidate_table <- .msi_add_scalar_column_safe(
    candidate_table, "candidate_source_for_step7", candidate_source_for_step7
  )
  candidate_table <- .msi_add_scalar_column_safe(
    candidate_table, "diagnostic_only",
    identical(candidate_source_for_step7, "reverse_match_all")
  )
  annotation_status <- if (nrow(candidate_table) == 0L) {
    if (identical(candidate_source_for_step7, "reverse_match_filtered")) {
      "no_candidates_after_reverse_match_filter"
    } else {
      "no_candidates_available"
    }
  } else {
    "candidates_available"
  }
  candidate_table <- .msi_add_scalar_column_safe(
    candidate_table, "annotation_status", annotation_status
  )
  matched_peak_pairs <- if (identical(candidate_source_for_step7, "reverse_match_all") &&
      is.data.frame(object@annotation_results$reverse_matched_peaks_all)) {
    object@annotation_results$reverse_matched_peaks_all
  } else {
    object@annotation_results$reverse_matched_peaks
  }
  dots <- list(...)
  old_weight_name <- paste0("reverse", "_dp_weight")
  old_threshold_name <- paste0("reverse", "_dp_threshold")
  if (old_weight_name %in% names(dots)) {
    warning(
      "Legacy spectral score weight argument is deprecated; use reverse_masked_cosine_weight.",
      call. = FALSE
    )
    reverse_masked_cosine_weight <- dots[[old_weight_name]]
  }
  if (old_threshold_name %in% names(dots)) {
    warning(
      "Legacy spectral threshold argument is deprecated and ignored; Module 7 uses tiered reverse_masked_cosine rules.",
      call. = FALSE
    )
  }
  ppm <- if ("ppm" %in% names(dots)) dots$ppm else 10

  spatial_scores <- extract_spatial_scores_from_reverse_match(
    object = object,
    reverse_match = candidate_table
  )

  score_table <- `71_calibrate_scores`(
    object = object,
    spatial_scores = spatial_scores,
    reverse_match = candidate_table,
    matched_peak_pairs = matched_peak_pairs
  )

  fused_scores <- `72_fuse_scores`(
    score_table = score_table,
    weights = weights,
    reverse_masked_cosine_weight = reverse_masked_cosine_weight,
    S_SAS_weight = S_SAS_weight,
    primary_spectral_score = primary_spectral_score
  )

  if (fdr_method == "target_decoy" && decoy_strategy != "none") {
    decoy_library <- `73_generate_decoy_library`(
      object = object,
      targets = candidate_table,
      strategy = decoy_strategy,
      decoy_n = decoy_n,
      decoy_shift_da = decoy_shift_da,
      seed = seed
    )

    decoy_scores <- `74_score_decoys`(
      object = object,
      decoy_library = decoy_library,
      weights = weights,
      score_table = fused_scores,
      ppm = ppm,
      tolerance_da = decoy_fragment_tolerance_da,
      reverse_masked_cosine_weight = reverse_masked_cosine_weight,
      S_SAS_weight = S_SAS_weight,
      primary_spectral_score = primary_spectral_score
    )
  } else {
    decoy_library <- .ann_empty_decoy_library()
    decoy_scores <- .ann_empty_decoy_scores()
  }

  annotated <- `75_estimate_annotation_q_value`(
    target_scores = fused_scores,
    decoy_scores = decoy_scores,
    method = fdr_method
  )

  annotated <- `76_filter_annotation_results`(
    annotated = annotated,
    q_value_cutoff = q_value_cutoff,
    min_final_score = min_final_score,
    S_SAS_threshold = S_SAS_threshold,
    keep_all = keep_all,
    confidence_thresholds = list(
      high_reverse_masked_cosine = confidence_high_reverse_masked_cosine,
      high_reverse_cos = confidence_high_reverse_cos,
      high_lib_coverage = confidence_high_lib_coverage,
      high_precursor_ppm = confidence_high_precursor_ppm,
      medium_reverse_masked_cosine = confidence_medium_reverse_masked_cosine,
      medium_reverse_cos = confidence_medium_reverse_cos,
      medium_lib_coverage = confidence_medium_lib_coverage,
      medium_matched_fraction = confidence_medium_matched_fraction,
      medium_precursor_ppm = confidence_medium_precursor_ppm
    )
  )
  if (identical(candidate_source_for_step7, "reverse_match_all") && nrow(annotated)) {
    annotated$pass_annotation_filter <- FALSE
  }

  object <- `77_update_annotation_table`(
    object = object,
    annotated = annotated
  )

  object@annotation_results$final_annotation_summary <- data.frame(
    candidate_source = candidate_source_for_step7,
    n_candidates_input = as.integer(nrow(candidate_table)),
    n_annotations = as.integer(nrow(object@annotation_results$final_annotation)),
    annotation_status = annotation_status,
    diagnostic_only = identical(candidate_source_for_step7, "reverse_match_all"),
    stringsAsFactors = FALSE,
    check.names = FALSE
  )

  object@annotation_results$decoy_library <- decoy_library
  object@annotation_results$decoy_scores <- decoy_scores

  fmc_requested <- .ann_fmc_requested(list(
    fmc_annotation = fmc_annotation,
    fmc_enabled = fmc_enabled,
    fmc_ms1_ppm = fmc_ms1_ppm,
    fmc_ion_mode = fmc_ion_mode,
    fmc_max_candidates_per_feature = fmc_max_candidates_per_feature,
    fmc_max_metabolites_per_feature = fmc_max_metabolites_per_feature,
    fmc_max_features_per_metabolite_soft = fmc_max_features_per_metabolite_soft,
    fmc_max_features_per_metabolite_hard = fmc_max_features_per_metabolite_hard,
    fmc_use_isotope = fmc_use_isotope,
    fmc_use_neutral_loss = fmc_use_neutral_loss,
    fmc_use_spatial_support = fmc_use_spatial_support,
    fmc_metabolite_reference = fmc_metabolite_reference
  ))
  if (isTRUE(fmc_requested)) {
    warning(
      "FMC is now Module 8. Please call calculate_mfc() after calculate_score().",
      call. = FALSE
    )
  }

  mfc_input <- .ann_build_mfc_input_from_pseudo_evidence(object)
  object@annotation_results$mfc_input_candidates <- mfc_input$candidates
  object@annotation_results$mfc_input_summary <- mfc_input$summary
  object@ms2_data$mfc_input_candidates <- mfc_input$candidates
  object@ms2_data$mfc_input_summary <- mfc_input$summary
  object@annotation_results$rdp_candidates <- NULL
  object@annotation_results$rdp_filter_summary <- NULL
  object@ms2_data$rdp_candidates <- NULL
  object@ms2_data$rdp_filter_summary <- NULL

  object <- `78_store_process_info`(
    object = object,
    parameters = list(
      weights = weights,
      decoy_strategy = decoy_strategy,
      fdr_method = fdr_method,
      q_value_cutoff = q_value_cutoff,
      min_final_score = min_final_score,
      keep_all = keep_all,
      decoy_n = decoy_n,
      decoy_shift_da = decoy_shift_da,
      seed = seed,
      ppm = ppm,
      confidence_high_reverse_cos = confidence_high_reverse_cos,
      confidence_high_reverse_masked_cosine = confidence_high_reverse_masked_cosine,
      confidence_high_lib_coverage = confidence_high_lib_coverage,
      confidence_high_precursor_ppm = confidence_high_precursor_ppm,
      confidence_medium_reverse_cos = confidence_medium_reverse_cos,
      confidence_medium_reverse_masked_cosine = confidence_medium_reverse_masked_cosine,
      confidence_medium_lib_coverage = confidence_medium_lib_coverage,
      confidence_medium_matched_fraction = confidence_medium_matched_fraction,
      confidence_medium_precursor_ppm = confidence_medium_precursor_ppm,
      reverse_masked_cosine_weight = reverse_masked_cosine_weight,
      S_SAS_weight = S_SAS_weight,
      candidate_source_for_step7 = candidate_source_for_step7,
      S_SAS_threshold = S_SAS_threshold,
      decoy_fragment_tolerance_da = decoy_fragment_tolerance_da,
      primary_spectral_score = primary_spectral_score,
      diagnostic_only = identical(candidate_source_for_step7, "reverse_match_all"),
      fmc_enabled = FALSE,
      fmc_module = "module8_calculate_mfc",
      S_SAS_missing_policy = "raw_missing_set_to_zero_for_S_total",
      score_formula = "S_total = 0.7 * reverse_masked_cosine + 0.3 * S_SAS_score_for_total",
      primary_reverse_score_column = "reverse_masked_cosine",
      n_mfc_input_candidates = as.integer(nrow(mfc_input$candidates)),
      n_mfc_input_pseudo_ms2 = as.integer(sum(mfc_input$candidates$match_mode == "pseudo_ms2", na.rm = TRUE)),
      n_mfc_input_ms1_only = as.integer(sum(mfc_input$candidates$match_mode %in% c("ms1", "ms1_only"), na.rm = TRUE)),
      mfc_input_source = "pseudo_ms2_candidate_evidence"
    )
  )

  object
}

`cluster_mfc_by_metabolite` <- 
function(mfc_candidate_table,
                                      expression_data,
                                      feature_graph = NULL,
                                      feature_metadata = NULL,
                                      tissue_mask = NULL,
                                      pixel_subset = NULL,
                                      component_id = "component1",
                                      analysis_domain = NULL,
                                      edge_threshold = 0.55) {
  cand <- .mfc2_df(mfc_candidate_table)
  if (!nrow(cand)) {
    out <- data.frame(
      mfc_id = character(), metabolite_id = character(), hmdb_id = character(),
      component_id = character(), cluster_index = integer(), feature_id = character(),
      stringsAsFactors = FALSE, check.names = FALSE
    )
    attr(out, "mfc_edge_table") <- data.frame()
    attr(out, "mfc_spatial_similarity") <- data.frame()
    return(out)
  }
  component_identity <- .mfc2_component_id(component_id, cand, analysis_domain = analysis_domain)
  graph_pairs <- .mfc2_graph_pair_lookup(feature_graph)
  all_edges <- list()
  all_sim <- list()
  memberships <- list()
  m <- 0L
  for (met in unique(as.character(cand$hmdb_id))) {
    sub <- cand[as.character(cand$hmdb_id) == met, , drop = FALSE]
    nodes <- unique(as.character(sub$feature_id))
    if (length(nodes) == 1L) {
      cl <- 1L
      cluster_id <- paste0("MFC_", met, "_", component_identity$cluster_label, "_cluster", cl)
      m <- m + 1L
      memberships[[m]] <- data.frame(
        mfc_id = cluster_id,
        metabolite_id = as.character(sub$metabolite_id[[1]]),
        hmdb_id = met,
        compound_name = as.character(sub$compound_name[[1]]),
        component_id = component_identity$component_value,
        cluster_index = cl,
        feature_id = nodes,
        stringsAsFactors = FALSE,
        check.names = FALSE
      )
      next
    }
    sim <- compute_mfc_spatial_similarity(
      expression_data = expression_data,
      feature_ids = nodes,
      tissue_mask = tissue_mask,
      pixel_subset = pixel_subset
    )
    if (nrow(sim)) {
      sim$hmdb_id <- met
      all_sim[[length(all_sim) + 1L]] <- sim
    }
    edge_rows <- list()
    e <- 0L
    for (i in seq_len(nrow(sim))) {
      fa <- sim$feature_id_a[[i]]
      fb <- sim$feature_id_b[[i]]
      pair <- sub[sub$feature_id %in% c(fa, fb), , drop = FALSE]
      pair_key <- paste(min(fa, fb), max(fa, fb), sep = "\r")
      graph_pair <- graph_pairs[graph_pairs$pair_key %in% pair_key, , drop = FALSE]
      on_demand <- if (!nrow(graph_pair)) {
        .mfc2_exact_pair_metrics(expression_data, fa, fb)
      } else {
        NULL
      }
      pair_metric_source <- if (nrow(graph_pair)) {
        "step3_pair_stats"
      } else {
        "on_demand_active_union_exact"
      }
      sas <- if (nrow(graph_pair) && "SAS" %in% names(graph_pair)) {
        suppressWarnings(as.numeric(graph_pair$SAS[[1]]))
      } else {
        suppressWarnings(as.numeric(on_demand$SAS))
      }
      if (!is.finite(sas)) sas <- 0
      car <- if (nrow(graph_pair) && "CAR" %in% names(graph_pair)) {
        suppressWarnings(as.numeric(graph_pair$CAR[[1]]))
      } else {
        suppressWarnings(as.numeric(on_demand$CAR))
      }
      aic <- if (nrow(graph_pair) && "AIC" %in% names(graph_pair)) {
        suppressWarnings(as.numeric(graph_pair$AIC[[1]]))
      } else {
        suppressWarnings(as.numeric(on_demand$AIC))
      }
      adducts <- unique(as.character(pair$adduct[!is.na(pair$adduct) & nzchar(pair$adduct)]))
      adduct_consistency <- if (length(adducts) == 0L) 0.5 else if (length(adducts) <= 3L) 1 else 0.6
      isotope_compatibility <- if (any(grepl("isotope", pair$evidence_type, ignore.case = TRUE))) 1 else 0.5
      evidence_weight <- if (any(grepl("pseudo", pair$evidence_type, ignore.case = TRUE))) 1 else 0.6
      edge_weight <- 0.40 * sim$spatial_similarity[[i]] +
        0.20 * sas +
        0.15 * adduct_consistency +
        0.05 * isotope_compatibility +
        0.10 * sim$hotspot_jaccard[[i]] +
        0.10 * evidence_weight
      if (is.finite(edge_weight) && edge_weight >= edge_threshold) {
        e <- e + 1L
        edge_rows[[e]] <- data.frame(
          hmdb_id = met,
          feature_id_a = fa,
          feature_id_b = fb,
          SAS = sas,
          CAR = car,
          AIC = aic,
          image_cosine = sim$spatial_cosine[[i]],
          spatial_similarity = sim$spatial_similarity[[i]],
          S_SAS_support = sas,
          adduct_consistency = adduct_consistency,
          isotope_compatibility = isotope_compatibility,
          pair_metric_source = pair_metric_source,
          hotspot_jaccard = sim$hotspot_jaccard[[i]],
          evidence_type_weight = evidence_weight,
          edge_weight = edge_weight,
          stringsAsFactors = FALSE,
          check.names = FALSE
        )
      }
    }
    edges <- if (length(edge_rows)) do.call(rbind, edge_rows) else data.frame(
      hmdb_id = character(), feature_id_a = character(), feature_id_b = character(),
      SAS = numeric(), CAR = numeric(), AIC = numeric(), image_cosine = numeric(),
      spatial_similarity = numeric(), S_SAS_support = numeric(),
      adduct_consistency = numeric(), isotope_compatibility = numeric(), hotspot_jaccard = numeric(),
      pair_metric_source = character(), evidence_type_weight = numeric(), edge_weight = numeric(),
      stringsAsFactors = FALSE, check.names = FALSE
    )
    if (nrow(edges)) all_edges[[length(all_edges) + 1L]] <- edges
    comps <- .mfc2_components(nodes, edges)
    for (cl in sort(unique(comps))) {
      features <- nodes[comps == cl]
      cluster_id <- paste0("MFC_", met, "_", component_identity$cluster_label, "_cluster", cl)
      m <- m + 1L
      memberships[[m]] <- data.frame(
        mfc_id = cluster_id,
        metabolite_id = as.character(sub$metabolite_id[[1]]),
        hmdb_id = met,
        compound_name = as.character(sub$compound_name[[1]]),
        component_id = component_identity$component_value,
        cluster_index = cl,
        feature_id = features,
        stringsAsFactors = FALSE,
        check.names = FALSE
      )
    }
  }
  out <- if (length(memberships)) do.call(rbind, memberships) else data.frame()
  rownames(out) <- NULL
  attr(out, "mfc_edge_table") <- if (length(all_edges)) do.call(rbind, all_edges) else data.frame()
  attr(out, "mfc_spatial_similarity") <- if (length(all_sim)) do.call(rbind, all_sim) else data.frame()
  out
}

`collapse_ms1_candidates_by_hmdb` <- 
function(candidates) {
  x <- as.data.frame(candidates, stringsAsFactors = FALSE, check.names = FALSE)
  if (!nrow(x)) return(x)
  required <- c("feature_id", "HMDB.ID", "abs_ppm_error", "final_evidence_score", "theoretical_mz")
  missing <- setdiff(required, names(x))
  if (length(missing)) {
    stop("MS1 candidate collapse is missing required columns: ", paste(missing, collapse = ", "), call. = FALSE)
  }
  x$HMDB.ID <- .ms1_normalize_hmdb_id(x$HMDB.ID)
  key <- paste(as.character(x$feature_id), x$HMDB.ID, sep = "\r")
  groups <- split(seq_len(nrow(x)), key, drop = TRUE)
  collapse_value <- function(value) {
    value <- trimws(as.character(value))
    value <- value[!is.na(value) & nzchar(value)]
    if (!length(value)) return(NA_character_)
    ord <- order(value)
    value <- value[ord]
    paste(unique(value), collapse = ";")
  }
  rows <- lapply(groups, function(index) {
    group <- x[index, , drop = FALSE]
    ord <- order(group$abs_ppm_error, -group$final_evidence_score, group$theoretical_mz, group$adduct)
    group <- group[ord, , drop = FALSE]
    out <- group[1L, , drop = FALSE]
    out$HMDB.ID <- .ms1_normalize_hmdb_id(out$HMDB.ID)
    out$raw_ms1_match_rows <- as.integer(nrow(group))
    out$matched_adducts <- collapse_value(group$adduct)
    out$matched_theoretical_mz <- collapse_value(sprintf("%.12f", as.numeric(group$theoretical_mz)))
    if ("source_row_count" %in% names(group)) {
      value <- suppressWarnings(as.integer(group$source_row_count))
      out$source_row_count <- as.integer(sum(value[is.finite(value)], na.rm = TRUE))
    }
    if ("source_Lab_IDs" %in% names(group)) out$source_Lab_IDs <- collapse_value(group$source_Lab_IDs)
    if ("source_record_ids" %in% names(group)) out$source_record_ids <- collapse_value(group$source_record_ids)
    out
  })
  out <- do.call(rbind, rows)
  ord <- order(out$feature_id, out$abs_ppm_error, -out$final_evidence_score, out$HMDB.ID)
  out <- out[ord, , drop = FALSE]
  rownames(out) <- NULL
  out
}

`compute_mfc_spatial_similarity` <- 
function(expression_data,
                                           feature_ids = NULL,
                                           tissue_mask = NULL,
                                           pixel_subset = NULL,
                                           coordinates = NULL,
                                           top_fraction = 0.10) {
  subset <- if (!is.null(pixel_subset)) pixel_subset else tissue_mask
  expr <- .mfc2_expression_matrix(expression_data, feature_ids = feature_ids, pixel_subset = subset)
  feature_ids <- rownames(expr)
  expr <- .mfc2_preprocess_expression(expr)
  if (nrow(expr) < 2L || ncol(expr) == 0L) {
    return(data.frame(
      feature_id_a = character(), feature_id_b = character(),
      spatial_cosine = numeric(), spatial_spearman = numeric(),
      hotspot_jaccard = numeric(), spatial_similarity = numeric(),
      stringsAsFactors = FALSE, check.names = FALSE
    ))
  }
  rows <- list()
  k <- 0L
  norms <- sqrt(rowSums(expr^2, na.rm = TRUE))
  hotspots <- lapply(seq_len(nrow(expr)), function(i) .mfc2_hotspot(expr[i, ], top_fraction = top_fraction))
  for (i in seq_len(nrow(expr) - 1L)) {
    for (j in seq.int(i + 1L, nrow(expr))) {
      cosine <- if (norms[[i]] > 0 && norms[[j]] > 0) sum(expr[i, ] * expr[j, ], na.rm = TRUE) / (norms[[i]] * norms[[j]]) else 0
      spear <- suppressWarnings(stats::cor(expr[i, ], expr[j, ], method = "spearman", use = "pairwise.complete.obs"))
      if (!is.finite(spear)) spear <- 0
      inter <- sum(hotspots[[i]] & hotspots[[j]])
      union <- sum(hotspots[[i]] | hotspots[[j]])
      jaccard <- if (union > 0) inter / union else 0
      sim <- 0.50 * max(0, cosine) + 0.30 * max(0, spear) + 0.20 * jaccard
      k <- k + 1L
      rows[[k]] <- data.frame(
        feature_id_a = feature_ids[[i]],
        feature_id_b = feature_ids[[j]],
        spatial_cosine = cosine,
        spatial_spearman = spear,
        hotspot_jaccard = jaccard,
        spatial_similarity = sim,
        stringsAsFactors = FALSE,
        check.names = FALSE
      )
    }
  }
  out <- do.call(rbind, rows)
  rownames(out) <- NULL
  out
}

`deduplicate_mfc_table` <- 
function(mfc_candidate_table,
                                  mfc_cluster_table,
                                  mfc_cluster_score_table = score_mfc_clusters(mfc_candidate_table, mfc_cluster_table),
                                  output_dir = NULL,
                                  score_gap = 0.05) {
  cand <- .mfc2_df(mfc_candidate_table)
  cl <- .mfc2_df(mfc_cluster_table)
  scores <- .mfc2_df(mfc_cluster_score_table)
  if (!nrow(cand)) {
    empty <- list(
      mfc_full_evidence_table = cand,
      mfc_spatial_cluster_table = scores,
      mfc_deduplicated_table = scores,
      mfc_ambiguity_groups = data.frame(),
      mfc_redundancy_log = data.frame()
    )
    return(empty)
  }
  feature_to_mfc <- cl[, c("mfc_id", "hmdb_id", "feature_id"), drop = FALSE]
  full <- merge(cand, feature_to_mfc, by = c("hmdb_id", "feature_id"), all.x = TRUE, sort = FALSE)
  score_cols <- scores[, c("mfc_id", "mfc_score", "mfc_spatial_coherence", "mfc_rank_within_metabolite"), drop = FALSE]
  full <- merge(full, score_cols, by = "mfc_id", all.x = TRUE, sort = FALSE)
  full$mfc_score[!is.finite(full$mfc_score)] <- full$raw_score[!is.finite(full$mfc_score)]
  full$mfc_score[!is.finite(full$mfc_score)] <- 0
  full$rank_metabolites_within_feature <- ave(
    -full$mfc_score,
    full$feature_id,
    FUN = function(x) rank(x, ties.method = "first")
  )
  full$rank_mfc_within_metabolite <- ave(
    -full$mfc_score,
    full$hmdb_id,
    FUN = function(x) rank(x, ties.method = "first")
  )
  duplicate_key <- if ("dedup_key" %in% names(full)) {
    paste(full$mfc_id, full$dedup_key, sep = "\r")
  } else {
    paste(full$mfc_id, full$feature_id, full$evidence_type,
          full$adduct, full$library_spectrum_id, sep = "\r")
  }
  full$is_redundant_evidence <- duplicated(duplicate_key)
  full$redundancy_reason <- ifelse(full$is_redundant_evidence, "duplicate_evidence_row", NA_character_)

  ambiguity_rows <- list()
  group_i <- 0L
  full$ambiguity_group_id <- NA_character_
  full$ambiguity_flag <- FALSE
  for (fid in unique(full$feature_id)) {
    idx <- which(full$feature_id == fid & !full$is_redundant_evidence)
    mets <- unique(full$hmdb_id[idx])
    if (length(mets) < 2L) next
    best <- stats::aggregate(mfc_score ~ hmdb_id, data = full[idx, , drop = FALSE], FUN = max)
    best <- best[order(-best$mfc_score, best$hmdb_id), , drop = FALSE]
    if (nrow(best) >= 2L && is.finite(best$mfc_score[[1]]) &&
        is.finite(best$mfc_score[[2]]) &&
        (best$mfc_score[[1]] - best$mfc_score[[2]]) <= score_gap) {
      group_i <- group_i + 1L
      gid <- sprintf("mfc_ambiguity_%04d", group_i)
      amb_mets <- best$hmdb_id[best$mfc_score >= best$mfc_score[[1]] - score_gap]
      amb_idx <- idx[full$hmdb_id[idx] %in% amb_mets]
      full$ambiguity_group_id[amb_idx] <- gid
      full$ambiguity_flag[amb_idx] <- TRUE
      ambiguity_rows[[group_i]] <- data.frame(
        ambiguity_group_id = gid,
        feature_id = fid,
        hmdb_ids = paste(amb_mets, collapse = ";"),
        n_metabolites = length(amb_mets),
        top_score = best$mfc_score[[1]],
        score_gap = best$mfc_score[[1]] - best$mfc_score[[min(2L, nrow(best))]],
        ambiguity_type = "shared_feature_similar_mfc_score",
        stringsAsFactors = FALSE,
        check.names = FALSE
      )
    }
  }
  dedup <- scores
  if (nrow(dedup)) {
    dedup$rank_mfc_within_metabolite <- ave(-dedup$mfc_score, dedup$hmdb_id, FUN = function(x) rank(x, ties.method = "first"))
    dedup$is_primary_mfc_for_metabolite <- dedup$rank_mfc_within_metabolite == 1L
    feature_map <- strsplit(as.character(dedup$feature_ids), ";", fixed = TRUE)
    dedup$ambiguity_group_id <- vapply(seq_along(feature_map), function(i) {
      features <- feature_map[[i]]
      hit <- full$ambiguity_group_id[full$feature_id %in% features & full$hmdb_id == dedup$hmdb_id[[i]]]
      hit <- hit[!is.na(hit) & nzchar(hit)]
      if (length(hit)) hit[[1]] else NA_character_
    }, character(1))
    dedup$ambiguity_flag <- !is.na(dedup$ambiguity_group_id) & nzchar(dedup$ambiguity_group_id)
  }
  ambiguity <- if (length(ambiguity_rows)) do.call(rbind, ambiguity_rows) else data.frame(
    ambiguity_group_id = character(), feature_id = character(), hmdb_ids = character(),
    n_metabolites = integer(), top_score = numeric(), score_gap = numeric(),
    ambiguity_type = character(), stringsAsFactors = FALSE, check.names = FALSE
  )
  redundancy_cols <- intersect(c(
    "source_row_id", "candidate_key", "dedup_key",
    "mfc_id", "feature_id", "hmdb_id", "evidence_type", "adduct",
    "library_spectrum_id", "redundancy_reason"
  ), names(full))
  redundancy <- full[full$is_redundant_evidence %in% TRUE, redundancy_cols, drop = FALSE]
  if (nrow(redundancy)) {
    redundancy$removal_reason <- "removed_by_mfc_redundancy"
  }
  if (is.null(output_dir)) {
    output_dir <- NA_character_
  }
  if (!is.na(output_dir) && nzchar(output_dir)) {
    dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
    utils::write.csv(full, file.path(output_dir, "mfc_full_evidence_table.csv"), row.names = FALSE, na = "")
    utils::write.csv(scores, file.path(output_dir, "mfc_spatial_cluster_table.csv"), row.names = FALSE, na = "")
    utils::write.csv(dedup, file.path(output_dir, "mfc_deduplicated_table.csv"), row.names = FALSE, na = "")
    utils::write.csv(ambiguity, file.path(output_dir, "mfc_ambiguity_groups.csv"), row.names = FALSE, na = "")
    utils::write.csv(redundancy, file.path(output_dir, "mfc_redundancy_log.csv"), row.names = FALSE, na = "")
  }
  message("MFC spatial deduplication: evidence_rows=", nrow(full),
          ", MFC_before=", nrow(scores), ", MFC_after=", nrow(dedup),
          ", ambiguity_groups=", nrow(ambiguity))
  list(
    mfc_full_evidence_table = full,
    mfc_spatial_cluster_table = scores,
    mfc_deduplicated_table = dedup,
    mfc_ambiguity_groups = ambiguity,
    mfc_redundancy_log = redundancy
  )
}

`default_adduct_rule_table` <- 
function() {
  rows <- list(
    data.frame(ion_mode = "positive", adduct = "pseudo_ms2_strict", mass_shift = NA_real_,
               charge = NA_integer_, evidence_type = "pseudo_ms2_strict",
               evidence_weight = 50, rank_group = "pseudo_ms2",
               description = "Strict reverse-match pseudo-MS2 evidence.", enabled = TRUE),
    data.frame(ion_mode = "negative", adduct = "pseudo_ms2_strict", mass_shift = NA_real_,
               charge = NA_integer_, evidence_type = "pseudo_ms2_strict",
               evidence_weight = 50, rank_group = "pseudo_ms2",
               description = "Strict reverse-match pseudo-MS2 evidence.", enabled = TRUE),
    data.frame(ion_mode = "positive", adduct = c("[M+H]+", "[M+Na]+", "[M+K]+", "[M+NH4]+",
                                                 "[M+H-H2O]+", "[M+ACN+H]+", "[M+2H]2+",
                                                 "[M+H+Na]2+"),
               mass_shift = c(1.007276, 22.989218, 38.963158, 18.033823,
                              1.007276 - 18.010565, 42.033825, 2.014552,
                              23.996494),
               charge = c(1L, 1L, 1L, 1L, 1L, 1L, 2L, 2L),
               evidence_type = "ms1_adduct",
               evidence_weight = c(30, 28, 22, 20, 14, 12, 15, 12),
               rank_group = c("canonical", "secondary", "secondary", "secondary",
                              "in_source_loss", "cluster_adduct", "secondary", "secondary"),
               description = c("Protonated molecule.", "Sodium adduct.", "Potassium adduct.",
                               "Ammonium adduct.", "Water-loss protonated ion.",
                               "Acetonitrile cluster adduct.", "Doubly protonated molecule.",
                               "Mixed proton/sodium double-charge adduct."),
               enabled = TRUE),
    data.frame(ion_mode = "negative", adduct = c("[M-H]-", "[M+FA-H]-", "[M+Ac-H]-",
                                                 "[M+Cl]-", "[M-H-H2O]-", "[M-2H]2-",
                                                 "[M+Br]-"),
               mass_shift = c(-1.007276, 44.998201, 59.013851, 34.969402,
                              -1.007276 - 18.010565, -2.014552, 78.918885),
               charge = c(-1L, -1L, -1L, -1L, -1L, -2L, -1L),
               evidence_type = "ms1_adduct",
               evidence_weight = c(35, 18, 16, 14, 12, 12, 8),
               rank_group = c("canonical", "secondary", "secondary", "secondary",
                              "in_source_loss", "secondary", "secondary"),
               description = c("Deprotonated molecule.", "Formate adduct.", "Acetate adduct.",
                               "Chloride adduct.", "Water-loss deprotonated ion.",
                               "Doubly deprotonated molecule.", "Bromide adduct."),
               enabled = TRUE),
    data.frame(ion_mode = c("positive", "positive", "negative", "negative"),
               adduct = c("M+1 isotope", "M+2 isotope", "M+1 isotope", "M+2 isotope"),
               mass_shift = c(1.003355, 2.006710, 1.003355, 2.006710),
               charge = c(1L, 1L, -1L, -1L),
               evidence_type = "isotope",
               evidence_weight = 10,
               rank_group = "isotope",
               description = c("M+1 isotope spacing.", "M+2 isotope spacing.",
                               "M+1 isotope spacing.", "M+2 isotope spacing."),
               enabled = TRUE)
  )
  out <- do.call(rbind, rows)
  rownames(out) <- NULL
  validate_adduct_rule_table(out)
  out
}

`default_fmc_annotation_config` <- 
function(enabled = FALSE) {
  list(
    enabled = isTRUE(enabled),
    enable_ms1_fallback = TRUE,
    legacy_ms1_fallback = FALSE,
    ms1_matching_scope = "pseudo_unresolved_only",
    combine_ms1_with_pseudo = FALSE,
    ms1_ppm = 10,
    max_candidates_per_feature = 5,
    report_top_k = 5,
    max_metabolites_per_feature = 20,
    max_features_per_metabolite_soft = 30,
    max_features_per_metabolite_hard = 100,
    use_isotope = TRUE,
    use_neutral_loss = FALSE,
    use_spatial_support = TRUE,
    count_duplicate_evidence_type_once = TRUE,
    ms1_library_path = file.path("derived_library", "msi_database_ms1_flat_v2.rds"),
    ms1_flat_server_mode = FALSE,
    ms1_flat_project_root = NULL,
    allowed_adducts = "auto_by_polarity",
    ppm_sigma = 5,
    pseudo_ms2_evidence_class_weight = 1.0,
    ms1_only_evidence_class_weight = 0.6,
    tier3_score_band_min = 0.67,
    tier2_score_band_min = 0.34,
    tier2_score_band_max = 0.66,
    tier1_score_band_max = 0.33,
    pseudo_score_band_min = 0.34,
    ms1_only_score_band_max = 0.33,
    pseudo_raw_weight = 0.70,
    ms1_support_weight = 0.10,
    spatial_weight = 0.20,
    cluster_support_weight = 0.10,
    ambiguity_score_gap = 0.05,
    relative_support_min = 0.20,
    min_priority_score = 0,
    strong_pseudo_raw_min = 0.70,
    spatial_edge_threshold = 0.55
  )
}

`default_fmc_score_config` <- 
function() {
  data.frame(
    parameter = c(
      "ppm_tolerance", "mz_accuracy_max_score", "spatial_weight",
      "count_duplicate_evidence_type_once", "use_neutral_loss",
      "canonical_adducts_positive", "canonical_adducts_negative"
    ),
    value = c("10", "10", "15", "TRUE", "FALSE", "[M+H]+", "[M-H]-"),
    description = c(
      "Default ppm tolerance used for m/z accuracy scoring.",
      "Maximum m/z accuracy bonus added to candidate support score.",
      "Multiplier for optional spatial/SAS support between related evidence features.",
      "Duplicate evidence types for the same feature-metabolite pair count once.",
      "Neutral-loss matching is disabled by default.",
      "Canonical positive adduct labels used for FMC evidence levels.",
      "Canonical negative adduct labels used for FMC evidence levels."
    ),
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
}

`default_neutral_loss_rule_table` <- 
function() {
  data.frame(
    neutral_loss = c("H2O loss", "NH3 loss", "CO2 loss"),
    formula = c("H2O", "NH3", "CO2"),
    mass_loss = c(18.010565, 17.026549, 43.989830),
    evidence_weight = c(12, 10, 10),
    evidence_type = "neutral_loss",
    enabled = FALSE,
    description = c("Water neutral loss.", "Ammonia neutral loss.", "Carbon dioxide neutral loss."),
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
}

`extract_pseudo_ms2_evidence` <- 
function(reverse_match_filtered,
                                        reverse_match_all = NULL,
                                        evidence_weight = 50,
                                        min_matched_peaks_used = NULL,
                                        min_reverse_masked_cosine_used = NULL) {
  df <- .fmc_as_df(reverse_match_filtered)
  empty <- .fmc_empty(c("component", "feature_id", "feature_mz", "metabolite_id", "hmdb_id",
                       "metabolite_name", "candidate_id", "anchor_feature_id", "query_id",
                       "adduct", "theoretical_mz", "mz_error_ppm", "source", "lab_id",
                       "reverse_masked_cosine", "S_SAS", "S_total",
                       "S_SAS_score_for_total", "S_SAS_missing",
                       "n_matched_peaks", "matched_query_intensity_fraction",
                       "matched_library_intensity_fraction", "library_l2_coverage",
                       "match_mode", "score_ms2", "score_spatial",
                       "score_total_pre_mfc", "pseudo_ms2_score", "matched_peak_score",
                       "evidence_level", "raw_evidence_score", "evidence_class_weight",
                       "final_evidence_score", "raw_score",
                       "confidence_tier", "evidence_source",
                       "mfc_entry_reason", "evidence_type", "evidence_weight", "data_source"))
  if (!nrow(df)) return(empty)
  n_matched <- .fmc_int(df, c("n_matched_peaks", "n_reverse_match", "n_match", "n_match_pairs"))
  reverse_masked_cosine <- .fmc_num(df, c("reverse_masked_cosine", "score_ms2", "S_reverse_cos"), 0)
  if ("evidence_pass" %in% names(df)) {
    evidence_pass <- df$evidence_pass
    if (!is.logical(evidence_pass)) {
      evidence_pass <- tolower(trimws(as.character(evidence_pass))) %in% c("true", "t", "1", "yes", "pass")
    }
    pass <- evidence_pass %in% TRUE
  } else {
    threshold_peaks <- if (is.null(min_matched_peaks_used)) NA_integer_ else
      suppressWarnings(as.integer(min_matched_peaks_used[[1]]))
    if (!length(threshold_peaks) || !is.finite(threshold_peaks) || threshold_peaks < 1L) {
      if ("min_matched_peaks_used" %in% names(df) && any(is.finite(suppressWarnings(as.integer(df$min_matched_peaks_used))))) {
        threshold_peaks <- suppressWarnings(as.integer(df$min_matched_peaks_used))
      } else {
        warning(
          "reverse_match_filtered lacks evidence_pass; recomputing pseudo-MS2 FMC evidence from supplied thresholds.",
          call. = FALSE
        )
        threshold_peaks <- 3L
      }
    }
    threshold_rmc <- if (is.null(min_reverse_masked_cosine_used)) NA_real_ else
      suppressWarnings(as.numeric(min_reverse_masked_cosine_used[[1]]))
    if (!length(threshold_rmc) || !is.finite(threshold_rmc) || threshold_rmc < 0) {
      if ("min_reverse_masked_cosine_used" %in% names(df) &&
          any(is.finite(suppressWarnings(as.numeric(df$min_reverse_masked_cosine_used))))) {
        threshold_rmc <- suppressWarnings(as.numeric(df$min_reverse_masked_cosine_used))
      } else {
        threshold_rmc <- 0
      }
    }
    pass <- n_matched >= threshold_peaks & reverse_masked_cosine > threshold_rmc
  }
  df <- df[pass, , drop = FALSE]
  n_matched <- n_matched[pass]
  reverse_masked_cosine <- reverse_masked_cosine[pass]
  if (!nrow(df)) return(empty)
  feature_id <- .fmc_chr(df, c("feature_id", "anchor_feature_id", "query_id"))
  candidate_id <- .fmc_chr(df, c("candidate_id", "metabolite_id", "compound_id", "candidate_key"))
  metabolite_name <- .fmc_chr(df, c("metabolite_name", "candidate_name", "compound_name", "name"))
  S_SAS <- .fmc_num(df, c("S_SAS", "S_spatial_reverse", "s_spat"), NA_real_)
  S_SAS_missing <- !is.finite(S_SAS)
  S_SAS_score <- ifelse(S_SAS_missing, 0, pmin(pmax(S_SAS, 0), 1))
  S_total <- .fmc_num(df, c("S_total", "score_total_pre_mfc"), NA_real_)
  missing_total <- !is.finite(S_total)
  S_total[missing_total] <- 0.7 * pmin(pmax(reverse_masked_cosine[missing_total], 0), 1) +
    0.3 * S_SAS_score[missing_total]
  matched_peak_score <- pmin(n_matched / 6, 1)
  lib_fraction <- .fmc_num(
    df,
    c("matched_library_intensity_fraction", "S_lib_coverage", "u_spec",
      "sum_matched_library_intensity_fraction")
  )
  lib_fraction[!is.finite(lib_fraction)] <- 0
  pseudo_ms2_score <- 0.55 * pmin(pmax(reverse_masked_cosine, 0), 1) +
    0.20 * matched_peak_score +
    0.15 * pmin(pmax(lib_fraction, 0), 1) +
    0.10 * S_SAS_score
  out <- data.frame(
    component = .fmc_chr(df, c("component", "component_id")),
    feature_id = feature_id,
    feature_mz = .fmc_num(df, c("feature_mz", "anchor_mz", "query_mz", "precursor_mz")),
    metabolite_id = candidate_id,
    hmdb_id = candidate_id,
    metabolite_name = metabolite_name,
    candidate_id = candidate_id,
    anchor_feature_id = .fmc_chr(df, c("anchor_feature_id", "feature_id", "query_id")),
    query_id = .fmc_chr(df, c("query_id", "anchor_feature_id", "feature_id")),
    adduct = .fmc_chr(df, c("adduct", "precursor_type", "ion_adduct")),
    theoretical_mz = .fmc_num(df, c("theoretical_mz", "precursor_mz", "library_record_precursor_mz")),
    mz_error_ppm = .fmc_num(df, c("mz_error_ppm", "precursor_error_ppm")),
    source = .fmc_chr(df, c("source", "database_source", "library_source"), "unknown"),
    lab_id = .fmc_chr(df, c("lab_id", "spectrum_id", "library_spectrum_id", "archive_id")),
    reverse_masked_cosine = reverse_masked_cosine,
    S_SAS = S_SAS,
    S_total = S_total,
    S_SAS_score_for_total = S_SAS_score,
    S_SAS_missing = S_SAS_missing,
    n_matched_peaks = n_matched,
    matched_query_intensity_fraction = .fmc_num(df, "matched_query_intensity_fraction"),
    matched_library_intensity_fraction = lib_fraction,
    library_l2_coverage = .fmc_num(df, "library_l2_coverage"),
    match_mode = "pseudo_ms2",
    score_ms2 = reverse_masked_cosine,
    score_spatial = S_SAS_score,
    score_total_pre_mfc = S_total,
    pseudo_ms2_score = pseudo_ms2_score,
    matched_peak_score = matched_peak_score,
    evidence_level = "pseudo_ms2",
    raw_evidence_score = pseudo_ms2_score,
    evidence_class_weight = 1.0,
    final_evidence_score = pseudo_ms2_score,
    raw_score = pseudo_ms2_score,
    confidence_tier = .fmc_chr(df, c("confidence_tier", "annotation_confidence"), "Tier2_pseudoMS2_candidate"),
    evidence_source = "pseudo_ms2_candidate_evidence",
    mfc_entry_reason = "pseudo_ms2_evidence_pass",
    evidence_type = "pseudo_ms2_strict",
    evidence_weight = as.numeric(evidence_weight[[1]]),
    data_source = "pseudo_ms2_reverse_match",
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
  out <- out[!is.na(out$feature_id) & nzchar(out$feature_id) &
               !is.na(out$metabolite_id) & nzchar(out$metabolite_id), , drop = FALSE]
  rownames(out) <- NULL
  out[, colnames(empty), drop = FALSE]
}

`extract_spatial_scores_from_reverse_match` <- 
function(object, reverse_match = NULL) {
  empty <- .ann_empty_spatial_scores()
  if (is.null(reverse_match)) reverse_match <- object@annotation_results$reverse_match
  reverse_match <- as.data.frame(reverse_match, stringsAsFactors = FALSE, check.names = FALSE)
  if (nrow(reverse_match) == 0L) {
    return(empty)
  }
  reverse_match <- .ann_add_missing_columns(reverse_match, c(
    "cluster_id", "anchor_feature_id", "anchor_mz", "candidate_key",
    "S_spatial_reverse", "n_reverse_match"
  ))
  total_counts <- .ann_total_spatial_peak_counts(object)
  out <- data.frame(
    cluster_id = reverse_match$cluster_id,
    anchor_feature_id = as.character(reverse_match$anchor_feature_id),
    anchor_mz = suppressWarnings(as.numeric(reverse_match$anchor_mz)),
    candidate_key = as.character(reverse_match$candidate_key),
    s_spat = suppressWarnings(as.numeric(reverse_match$S_spatial_reverse)),
    explained_spatial_peak_n = suppressWarnings(as.integer(reverse_match$n_reverse_match)),
    total_spatial_peak_n = NA_integer_,
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
  if (nrow(total_counts) > 0L) {
    total_lookup <- stats::setNames(
      total_counts$total_spatial_peak_n,
      .ann_split_key(total_counts, c("cluster_id", "anchor_feature_id"))
    )
    out_key <- .ann_split_key(out, c("cluster_id", "anchor_feature_id"))
    out$total_spatial_peak_n <- suppressWarnings(as.integer(total_lookup[out_key]))
  }
  out <- out[, colnames(empty), drop = FALSE]
  rownames(out) <- NULL
  out
}

`filter_mfc_priority_candidates` <- 
function(scored_candidates, mfc_config = list()) {
  x <- as.data.frame(scored_candidates, stringsAsFactors = FALSE, check.names = FALSE)
  if (!nrow(x)) return(list(retained = x, rejected = x, ambiguity_groups = data.frame()))
  ambiguity_gap <- as.numeric(.mfc_unified_value(mfc_config, "ambiguity_score_gap", 0.05))
  relative_min <- as.numeric(.mfc_unified_value(mfc_config, "relative_support_min", 0.20))
  minimum_score <- as.numeric(.mfc_unified_value(mfc_config, "min_priority_score", 0))
  strong_pseudo_min <- as.numeric(.mfc_unified_value(mfc_config, "strong_pseudo_raw_min", 0.70))
  keep <- is.finite(x$mfc_priority_score) & x$mfc_priority_score >= minimum_score &
    (x$rank_within_feature == 1L |
       x$score_gap <= ambiguity_gap |
       x$relative_support >= relative_min |
       (x$has_pseudo_ms2 & .mfc_unified_clamp01(x$best_pseudo_ms2_score) >= strong_pseudo_min))
  x$retained <- keep
  x$retained_reason <- ifelse(
    !keep, NA_character_,
    ifelse(x$rank_within_feature == 1L, "best_within_feature",
      ifelse(x$score_gap <= ambiguity_gap, "within_ambiguity_score_gap",
        ifelse(x$relative_support >= relative_min, "relative_support",
          "strong_pseudo_ms2_support")))
  )
  x$rejection_reason <- ifelse(
    keep, NA_character_,
    ifelse(!is.finite(x$mfc_priority_score) | x$mfc_priority_score < minimum_score,
      "below_min_priority_score",
      "outside_score_gap_and_relative_support")
  )
  retained <- x[keep, , drop = FALSE]
  rejected <- x[!keep, , drop = FALSE]
  amb <- retained[retained$score_gap <= ambiguity_gap, , drop = FALSE]
  if (nrow(amb)) {
    n_by_feature <- table(amb$feature_id)
    amb <- amb[amb$feature_id %in% names(n_by_feature[n_by_feature > 1L]), , drop = FALSE]
    amb$ambiguity_group <- if (nrow(amb)) paste0("MFC_AMBIG_", match(amb$feature_id, unique(amb$feature_id))) else character()
  }
  list(retained = retained, rejected = rejected, ambiguity_groups = amb)
}

`generate_ms1_only_candidates` <- 
function(object,
                                         ms1_library = .ms1_flat_relative_path,
                                         evidence_pass_features = character(),
                                         ppm = 10,
                                         ppm_sigma = 5,
                                         allowed_adducts = "auto_by_polarity",
                                         polarity = "negative",
                                         max_candidates_per_feature = 100,
                                         report_top_k = NULL) {
  if (!methods::is(object, "msi_dataset")) {
    stop("object must be a msi_dataset.", call. = FALSE)
  }
  lib <- if (is.character(ms1_library) && length(ms1_library) == 1L) {
    load_msi_database_ms1_flat(ms1_library, server_mode = FALSE)
  } else {
    as.data.frame(ms1_library, stringsAsFactors = FALSE, check.names = FALSE)
  }
  .ms1_flat_validate_table(lib)
  # The MS1 query path uses only measured metadata fields. In particular, do
  # not infer a query mass by parsing a feature identifier.
  features <- .fmc_feature_table(object, parse_feature_id_mz = FALSE)
  features <- features[!as.character(features$feature_id) %in% as.character(evidence_pass_features), , drop = FALSE]
  ppm <- suppressWarnings(as.numeric(ppm[[1]]))
  if (!is.finite(ppm) || ppm <= 0) ppm <- 10
  ppm_sigma <- suppressWarnings(as.numeric(ppm_sigma[[1]]))
  if (!is.finite(ppm_sigma) || ppm_sigma <= 0) ppm_sigma <- 5
  report_top_k <- if (is.null(report_top_k)) max_candidates_per_feature else report_top_k
  report_top_k <- suppressWarnings(as.integer(report_top_k[[1]]))
  if (!is.finite(report_top_k) || report_top_k < 1L) report_top_k <- 5L
  allowed <- .ms1_allowed_adducts(allowed_adducts, polarity)
  if (nrow(lib)) {
    adduct <- as.character(.ms1_column(lib, c("Adduct")))
    missing_adduct <- is.na(adduct) | !nzchar(trimws(adduct))
    lib <- lib[missing_adduct | adduct %in% allowed, , drop = FALSE]
    ion <- tolower(as.character(.ms1_column(lib, c("Ion.Mode"))))
    missing_ion <- is.na(ion) | !nzchar(trimws(ion))
    polarity_norm <- tolower(as.character(polarity[[1]]))
    if (polarity_norm %in% c("positive", "pos", "+")) {
      lib <- lib[missing_ion | ion %in% c("positive", "pos", "+"), , drop = FALSE]
    } else {
      lib <- lib[missing_ion | ion %in% c("negative", "neg", "-"), , drop = FALSE]
    }
  }
  empty <- data.frame(
    feature_id = character(), feature_mz = numeric(), observed_mz = numeric(),
    observed_mz_full = numeric(), ms1_feature_mz = numeric(),
    ms1_mz_source = character(), ms1_mz_provenance = character(), HMDB.ID = character(),
    Lab.ID = character(), compound_name = character(), formula = character(),
    adduct = character(), ion_mode = character(), theoretical_mz = numeric(),
    ppm_error = numeric(), abs_ppm_error = numeric(), ppm_score = numeric(),
    adduct_prior = numeric(), isotope_score = numeric(), spatial_quality_score = numeric(),
    ms1_only_score = numeric(), evidence_type = character(), evidence_level = character(),
    raw_evidence_score = numeric(), evidence_class_weight = numeric(),
    final_evidence_score = numeric(), metabolite_id = character(),
    metabolite_name = character(), mz_error_ppm = numeric(),
    score_ms1 = numeric(), raw_score = numeric(), match_mode = character(),
    evidence_source = character(), data_source = character(), evidence_weight = numeric(),
    source_row_count = integer(), source_Lab_IDs = character(),
    source_record_ids = character(), canonical_key = character(),
    raw_ms1_match_rows = integer(), matched_adducts = character(),
    matched_theoretical_mz = character(), rank_within_feature = integer(),
    stringsAsFactors = FALSE, check.names = FALSE
  )
  if (!nrow(features) || !nrow(lib)) return(empty)
  lib <- lib[order(suppressWarnings(as.numeric(lib$theoretical_mz))), , drop = FALSE]
  rows <- vector("list", nrow(features))
  tol <- ppm / 1e6
  resolved_mz <- .msi_resolve_observed_mz(features)
  for (i in seq_len(nrow(features))) {
    fmz <- suppressWarnings(as.numeric(resolved_mz$mz[[i]]))
    if (!is.finite(fmz) || fmz <= 0) next
    lower <- fmz / (1 + tol)
    upper <- fmz / (1 - tol)
    idx <- which(lib$theoretical_mz >= lower & lib$theoretical_mz <= upper)
    if (!length(idx)) next
    hit <- lib[idx, , drop = FALSE]
    err <- (fmz - hit$theoretical_mz) / hit$theoretical_mz * 1e6
    keep <- abs(err) <= ppm
    if (!any(keep)) next
    hit <- hit[keep, , drop = FALSE]
    err <- err[keep]
    ppm_score <- exp(-0.5 * (abs(err) / ppm_sigma)^2)
    hit_adduct <- as.character(.ms1_column(hit, c("Adduct")))
    hit_formula <- as.character(.ms1_column(hit, c("Formula")))
    hit_ion_mode <- as.character(.ms1_column(hit, c("Ion.Mode")))
    adduct_prior <- ifelse(hit_adduct %in% c("[M+H]+", "[M-H]-"), 1.0, 0.7)
    isotope_score <- rep(0, nrow(hit))
    spatial_quality_score <- rep(0, nrow(hit))
    ms1_only_score <- 0.80 * ppm_score + 0.20 * adduct_prior
    rows[[i]] <- data.frame(
      feature_id = as.character(features$feature_id[[i]]),
      feature_mz = fmz,
      observed_mz = as.numeric(resolved_mz$observed_mz[[i]]),
      observed_mz_full = as.numeric(resolved_mz$observed_mz_full[[i]]),
      ms1_feature_mz = fmz,
      ms1_mz_source = as.character(resolved_mz$source[[i]]),
      ms1_mz_provenance = if (identical(resolved_mz$source[[i]], "observed_mz_full")) {
        "full_precision_feature_centroid"
      } else {
        "legacy_observed_mz_fallback"
      },
      HMDB.ID = hit$HMDB.ID,
      Lab.ID = as.character(.ms1_column(hit, c("source_Lab_IDs", "Lab.ID", "source_lab_ids"))),
      compound_name = hit$Compound.Name,
      formula = hit_formula,
      adduct = hit_adduct,
      ion_mode = hit_ion_mode,
      theoretical_mz = hit$theoretical_mz,
      ppm_error = err,
      abs_ppm_error = abs(err),
      ppm_score = ppm_score,
      adduct_prior = adduct_prior,
      isotope_score = isotope_score,
      spatial_quality_score = spatial_quality_score,
      ms1_only_score = ms1_only_score,
      evidence_type = "ms1_only",
      evidence_level = "MS1_only",
      raw_evidence_score = ms1_only_score,
      evidence_class_weight = 0.6,
      final_evidence_score = ms1_only_score * 0.6,
      metabolite_id = hit$HMDB.ID,
      metabolite_name = hit$Compound.Name,
      mz_error_ppm = err,
      score_ms1 = ms1_only_score,
      raw_score = ms1_only_score * 0.6,
      match_mode = "ms1",
      evidence_source = "module8_ms1_only",
      data_source = "module8_ms1_only",
      evidence_weight = ms1_only_score * 0.6,
      source_row_count = suppressWarnings(as.integer(.ms1_column(hit, c("source_row_count"), default = NA_integer_))),
      source_Lab_IDs = as.character(.ms1_column(hit, c("source_Lab_IDs", "Lab.ID"))),
      source_record_ids = as.character(.ms1_column(hit, c("source_record_ids"))),
      canonical_key = as.character(.ms1_column(hit, c("canonical_key"))),
      raw_ms1_match_rows = rep.int(1L, nrow(hit)),
      matched_adducts = hit_adduct,
      matched_theoretical_mz = sprintf("%.12f", as.numeric(hit$theoretical_mz)),
      stringsAsFactors = FALSE,
      check.names = FALSE
    )
  }
  out <- do.call(rbind, rows[!vapply(rows, is.null, logical(1))])
  if (is.null(out) || !nrow(out)) return(empty)
  out <- rank_ms1_candidates_internal(out)
  # Keep the legacy argument/report_top_k in the API for compatibility, but do
  # not apply it here.  Calling code may use ms1_candidate_report_top_k().
  out[, unique(c(colnames(empty), "rank_within_feature")), drop = FALSE]
}

`get_adduct_rules_for_mode` <- 
function(ion_mode,
                                      adduct_rules = default_adduct_rule_table()) {
  ion_mode <- tolower(as.character(ion_mode[[1]]))
  if (!ion_mode %in% c("positive", "negative")) {
    stop("ion_mode must be 'positive' or 'negative'.", call. = FALSE)
  }
  validate_adduct_rule_table(adduct_rules)
  out <- adduct_rules[tolower(as.character(adduct_rules$ion_mode)) == ion_mode &
                        adduct_rules$enabled %in% TRUE, , drop = FALSE]
  rownames(out) <- NULL
  out
}

`group_msi_features` <- 
function(edge_tbl,
                               method = c("louvain", "leiden", "walktrap"),
                               resolution = 1,
                               seed = NULL) {
  if (!requireNamespace("igraph", quietly = TRUE)) {
    stop("Package 'igraph' is required for graph clustering.")
  }
  method <- match.arg(method)
  if (!is.numeric(resolution) || length(resolution) != 1 ||
      !is.finite(resolution) || resolution <= 0) {
    stop("resolution must be a single positive number.")
  }

  if (!is.data.frame(edge_tbl)) {
    stop("edge_tbl must be a data.frame/tibble.")
  }
  required_cols <- c("feature_i", "feature_j", "weight")
  if (!all(required_cols %in% colnames(edge_tbl))) {
    stop("edge_tbl must contain feature_i, feature_j, and weight columns.")
  }

  feature_ids <- attr(edge_tbl, "feature_ids")
  if (is.null(feature_ids) || length(feature_ids) == 0) {
    feature_ids <- unique(c(as.character(edge_tbl$feature_i), as.character(edge_tbl$feature_j)))
  }
  feature_ids <- unique(as.character(feature_ids))

  if (length(feature_ids) == 0) {
    return(tibble::tibble(feature_id = character(), cluster_id = integer()))
  }

  vertices <- data.frame(name = feature_ids, stringsAsFactors = FALSE)
  if (nrow(edge_tbl) > 0) {
    edge_df <- data.frame(
      from = as.character(edge_tbl$feature_i),
      to = as.character(edge_tbl$feature_j),
      weight = as.numeric(edge_tbl$weight),
      stringsAsFactors = FALSE
    )
    edge_df <- edge_df[is.finite(edge_df$weight), , drop = FALSE]
  } else {
    edge_df <- data.frame(from = character(), to = character(), weight = numeric())
  }

  g <- igraph::graph_from_data_frame(
    d = edge_df,
    directed = FALSE,
    vertices = vertices
  )

  if (igraph::ecount(g) == 0) {
    return(tibble::tibble(
      feature_id = feature_ids,
      cluster_id = seq_along(feature_ids)
    ))
  }

  cl <- NULL
  if (method == "louvain") {
    cl <- tryCatch({
      louvain_formals <- names(formals(igraph::cluster_louvain))
      if ("resolution" %in% louvain_formals) {
        igraph::cluster_louvain(
          g,
          weights = igraph::E(g)$weight,
          resolution = as.numeric(resolution)
        )
      } else {
        if (!isTRUE(all.equal(as.numeric(resolution), 1))) {
          warning("Installed igraph::cluster_louvain() does not support resolution; using default behavior.")
        }
        igraph::cluster_louvain(g, weights = igraph::E(g)$weight)
      }
    }, error = function(e) NULL)
    if (is.null(cl)) {
      cl <- igraph::cluster_walktrap(g, weights = igraph::E(g)$weight)
    }
  } else if (method == "leiden") {
    if (!exists("cluster_leiden", envir = asNamespace("igraph"), inherits = FALSE)) {
      stop("Installed igraph does not provide cluster_leiden().", call. = FALSE)
    }
    run_leiden <- function() {
      igraph::cluster_leiden(
        g,
        objective_function = "modularity",
        weights = igraph::E(g)$weight,
        resolution = as.numeric(resolution)
      )
    }
    cl <- if (is.null(seed)) run_leiden() else {
      old_seed <- if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) {
        get(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
      } else {
        NULL
      }
      on.exit({
        if (is.null(old_seed)) {
          if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) {
            rm(".Random.seed", envir = .GlobalEnv)
          }
        } else {
          assign(".Random.seed", old_seed, envir = .GlobalEnv)
        }
      }, add = TRUE)
      set.seed(as.integer(seed[[1L]]))
      run_leiden()
    }
  } else {
    cl <- igraph::cluster_walktrap(g, weights = igraph::E(g)$weight)
  }

  membership <- igraph::membership(cl)
  out <- tibble::tibble(
    feature_id = names(membership),
    cluster_id = as.integer(unname(membership))
  ) |>
    dplyr::arrange(cluster_id, feature_id)

  out
}

`load_msi_database_ms1_flat` <- 
function(
    path = .ms1_flat_relative_path,
    project_root = getwd(),
    server_mode = FALSE,
    expected_nrow = if (isTRUE(server_mode)) 2510438L else NULL) {
  requested_path <- .ms1_flat_resolve_path(path, project_root = project_root, must_work = FALSE)
  if (!identical(tolower(tools::file_ext(requested_path)), "rds")) {
    stop("MS1 flat loader accepts only the published .rds library; legacy RData fallback is disabled.",
         call. = FALSE)
  }
  if (isTRUE(server_mode)) {
    canonical_path <- .ms1_flat_resolve_path(
      .ms1_flat_relative_path, project_root = project_root, must_work = FALSE
    )
    if (!identical(requested_path, canonical_path)) {
      stop(
        "Server mode only permits ", canonical_path,
        "; requested ", requested_path, ".",
        call. = FALSE
      )
    }
  }
  if (!file.exists(requested_path) || dir.exists(requested_path)) {
    stop("MS1 flat RDS does not exist: ", requested_path, call. = FALSE)
  }
  # A plain RDS reader must not initialize a legacy MS1 namespace itself.  A
  # current-source package session can legitimately have ncdf4 present as a
  # transitive namespace while still using only this flat RDS; rejecting that
  # pre-existing namespace made the formal server path impossible to run.
  # Record it for provenance, but fail only if this dependency-free loader
  # caused a forbidden namespace to load.
  prohibited_before <- intersect(.ms1_flat_forbidden_namespaces, loadedNamespaces())
  table <- readRDS(requested_path)
  prohibited_after <- intersect(.ms1_flat_forbidden_namespaces, loadedNamespaces())
  prohibited_loaded_by_loader <- setdiff(prohibited_after, prohibited_before)
  if (isTRUE(server_mode) && length(prohibited_loaded_by_loader)) {
    stop(
      "Flat MS1 RDS loader initialized forbidden legacy namespace(s): ",
      paste(prohibited_loaded_by_loader, collapse = ", "),
      call. = FALSE
    )
  }
  .ms1_flat_validate_table(table, expected_nrow = expected_nrow)
  attr(table, "ms1_flat_loader_id") <- .ms1_flat_loader_id
  attr(table, "ms1_flat_library_path") <- normalizePath(
    requested_path, winslash = "/", mustWork = TRUE
  )
  attr(table, "ms1_flat_server_mode") <- isTRUE(server_mode)
  attr(table, "ms1_flat_preloaded_legacy_namespaces") <- prohibited_before
  table
}

`merge_mfc_feature_hmdb_candidates` <- 
function(evidence_long) {
  evidence <- as.data.frame(evidence_long, stringsAsFactors = FALSE, check.names = FALSE)
  required <- c("feature_id", "HMDB_ID", "evidence_source", "evidence_pass", "raw_evidence_score")
  missing <- setdiff(required, names(evidence))
  if (length(missing)) {
    stop("evidence_long is missing required columns: ", paste(missing, collapse = ", "), call. = FALSE)
  }
  evidence$HMDB_ID <- .mfc_normalize_hmdb_id(evidence$HMDB_ID)
  all_evidence <- evidence
  evidence <- evidence[evidence$evidence_pass %in% TRUE, , drop = FALSE]
  if (!nrow(evidence)) return(data.frame())
  key <- paste(evidence$feature_id, evidence$HMDB_ID, sep = "\r")
  groups <- split(seq_len(nrow(evidence)), key, drop = TRUE)
  rows <- lapply(groups, function(idx) {
    x <- evidence[idx, , drop = FALSE]
    pseudo <- x[x$evidence_source == "pseudo_ms2", , drop = FALSE]
    ms1 <- x[x$evidence_source == "ms1", , drop = FALSE]
    all_group <- all_evidence[
      all_evidence$feature_id == x$feature_id[[1]] &
        all_evidence$HMDB_ID == x$HMDB_ID[[1]],
      , drop = FALSE
    ]
    all_pseudo <- all_group[all_group$evidence_source == "pseudo_ms2", , drop = FALSE]
    best_pseudo <- if (nrow(pseudo)) pseudo[which.max(pseudo$raw_evidence_score), , drop = FALSE] else NULL
    best_ms1 <- if (nrow(ms1)) ms1[which.max(ms1$raw_evidence_score), , drop = FALSE] else NULL
    first <- as.character(x$feature_id[[1]])
    compound <- x$compound_name[!is.na(x$compound_name) & nzchar(x$compound_name)]
    adduct <- if (nrow(ms1)) best_ms1$adduct[[1]] else if (nrow(pseudo)) best_pseudo$adduct[[1]] else NA_character_
    data.frame(
      feature_id = first,
      feature_mz = suppressWarnings(as.numeric(x$feature_mz[[which.max(is.finite(x$feature_mz))]])),
      HMDB_ID = as.character(x$HMDB_ID[[1]]),
      compound_name = if (length(compound)) compound[[1]] else NA_character_,
      best_adduct = as.character(adduct),
      has_pseudo_ms2 = nrow(pseudo) > 0L,
      has_ms1 = nrow(ms1) > 0L,
      evidence_source = if (nrow(pseudo) && nrow(ms1)) {
        "ms1_and_pseudo_ms2"
      } else if (nrow(pseudo)) {
        "pseudo_ms2_only"
      } else {
        "ms1"
      },
      source_row_ids = paste(sort(unique(as.character(x$source_row_id[
        !is.na(x$source_row_id) & nzchar(x$source_row_id)
      ]))), collapse = ";"),
      candidate_keys = paste(sort(unique(as.character(x$candidate_key[
        !is.na(x$candidate_key) & nzchar(x$candidate_key)
      ]))), collapse = ";"),
      dedup_keys = paste(sort(unique(as.character(x$dedup_key[
        !is.na(x$dedup_key) & nzchar(x$dedup_key)
      ]))), collapse = ";"),
      best_source_row_id = if (nrow(pseudo)) {
        as.character(best_pseudo$source_row_id[[1]])
      } else if (nrow(ms1)) {
        as.character(best_ms1$source_row_id[[1]])
      } else {
        NA_character_
      },
      n_library_spectra = length(unique(all_pseudo$library_spectrum_id[
        !is.na(all_pseudo$library_spectrum_id) & nzchar(all_pseudo$library_spectrum_id)
      ])),
      best_pseudo_ms2_spectrum = if (nrow(pseudo)) as.character(best_pseudo$library_spectrum_id[[1]]) else NA_character_,
      best_pseudo_ms2_score = if (nrow(pseudo)) as.numeric(best_pseudo$raw_evidence_score[[1]]) else NA_real_,
      best_ms1_score = if (nrow(ms1)) as.numeric(best_ms1$raw_evidence_score[[1]]) else NA_real_,
      best_reverse_masked_cosine = if (nrow(pseudo)) as.numeric(best_pseudo$reverse_masked_cosine[[1]]) else NA_real_,
      best_n_matched_peaks = if (nrow(pseudo)) as.integer(best_pseudo$n_matched_peaks[[1]]) else NA_integer_,
      best_matched_library_fraction = if (nrow(pseudo)) as.numeric(best_pseudo$matched_library_intensity_fraction[[1]]) else NA_real_,
      best_S_SAS = if (nrow(pseudo)) as.numeric(best_pseudo$S_SAS[[1]]) else NA_real_,
      best_step7_q_value = if (nrow(pseudo)) as.numeric(best_pseudo$step7_q_value[[1]]) else NA_real_,
      best_ms1_ppm_error = if (nrow(ms1)) as.numeric(best_ms1$ms1_ppm_error[[1]]) else NA_real_,
      source_provenance = paste(sort(unique(paste(x$evidence_source, x$library_name, sep = ":"))), collapse = ";"),
      stringsAsFactors = FALSE,
      check.names = FALSE
    )
  })
  out <- do.call(rbind, rows)
  rownames(out) <- NULL
  if (all(c("ms1_feature_mz", "precursor_mz", "precursor_source") %in% names(evidence_long))) {
    out$ms1_feature_mz <- vapply(groups, function(idx) {
      values <- suppressWarnings(as.numeric(evidence$ms1_feature_mz[idx]))
      values <- values[is.finite(values)]
      if (length(values)) values[[1]] else NA_real_
    }, numeric(1))
    out$ms1_mz_source <- vapply(groups, function(idx) {
      values <- as.character(evidence$ms1_mz_source[idx])
      values <- values[!is.na(values) & nzchar(values)]
      if (length(values)) values[[1]] else NA_character_
    }, character(1))
    out$ms1_mz_provenance <- vapply(groups, function(idx) {
      values <- as.character(evidence$ms1_mz_provenance[idx])
      values <- values[!is.na(values) & nzchar(values)]
      if (length(values)) values[[1]] else NA_character_
    }, character(1))
    out$precursor_mz <- NA_real_
    out$precursor_source <- "missing_explicit_precursor"
  }
  out[order(out$feature_id, out$HMDB_ID), , drop = FALSE]
}

`mfc_candidate_report_top_k` <- 
function(candidates, report_top_k = 5L) {
  x <- as.data.frame(candidates, stringsAsFactors = FALSE, check.names = FALSE)
  if (!nrow(x)) return(x)
  k <- suppressWarnings(as.integer(report_top_k[[1]]))
  if (!is.finite(k) || k < 1L) k <- 5L
  if (!"rank_within_feature" %in% names(x)) {
    x <- x[order(x$feature_id, -x$mfc_priority_score, x$HMDB_ID), , drop = FALSE]
    x$rank_within_feature <- ave(seq_len(nrow(x)), x$feature_id, FUN = seq_along)
  }
  x[x$rank_within_feature <= k, , drop = FALSE]
}

`mfc_scoring_parameter_table` <- 
function(mfc_config = list()) {
  defaults <- list(
    tier3_score_band_min = 0.67,
    tier2_score_band_min = 0.34,
    tier2_score_band_max = 0.66,
    tier1_score_band_max = 0.33,
    pseudo_score_band_min = 0.34,
    ms1_only_score_band_max = 0.33,
    pseudo_raw_weight = 0.70,
    ms1_support_weight = 0.10,
    spatial_weight = 0.20,
    cluster_support_weight = 0.10,
    ambiguity_score_gap = 0.05,
    relative_support_min = 0.20,
    max_candidates_per_feature = 5L,
    report_top_k = 5L,
    min_priority_score = 0,
    strong_pseudo_raw_min = 0.70
  )
  value <- vapply(names(defaults), function(nm) {
    as.character(.mfc_unified_value(mfc_config, nm, defaults[[nm]]))
  }, character(1))
  data.frame(
    parameter = names(defaults),
    value = unname(value),
    semantics = c(
      "lower bound for candidates with both valid pseudo-MS2 and MS1 support",
      "lower bound for pseudo-MS2-only candidates",
      "upper bound for pseudo-MS2-only candidates",
      "upper bound for MS1-only candidates",
      "deprecated alias retained for pseudo-MS2 band compatibility",
      "deprecated alias retained for MS1-only band compatibility",
      "weight in pseudo-supported raw score",
      "optional MS1 support weight for pseudo candidates",
      "spatial support weight",
      "cluster support weight",
      "within-feature ambiguity score gap",
      "minimum score relative to best feature candidate",
      "deprecated compatibility alias; reporting cap only",
      "maximum candidates in user-facing report/export",
      "minimum absolute priority score",
      "raw pseudo score retaining strong ambiguous evidence"
    ),
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
}

`ms1_candidate_report_top_k` <- 
function(candidates, report_top_k = 5L) {
  x <- as.data.frame(candidates, stringsAsFactors = FALSE, check.names = FALSE)
  if (!nrow(x)) return(x)
  k <- suppressWarnings(as.integer(report_top_k[[1]]))
  if (!is.finite(k) || k < 1L) k <- 5L
  if (!"rank_within_feature" %in% names(x)) {
    x <- rank_ms1_candidates_internal(x)
  }
  x[x$rank_within_feature <= k, , drop = FALSE]
}

`rank_ms1_candidates_internal` <- 
function(candidates) {
  x <- as.data.frame(candidates, stringsAsFactors = FALSE, check.names = FALSE)
  if (!nrow(x)) return(x)
  x <- collapse_ms1_candidates_by_hmdb(x)
  x <- x[order(x$feature_id, x$abs_ppm_error, -x$final_evidence_score, x$HMDB.ID), , drop = FALSE]
  x$rank_within_feature <- ave(
    seq_len(nrow(x)), x$feature_id,
    FUN = seq_along
  )
  rownames(x) <- NULL
  x
}

`run_fmc_annotation` <- 
function(object,
                               fmc_config = default_fmc_annotation_config(TRUE),
                               metabolite_reference = NULL,
                               library_object = NULL,
                               ion_mode = "negative",
                               adduct_rules = default_adduct_rule_table(),
                               neutral_loss_rules = default_neutral_loss_rule_table(),
                               score_config = default_fmc_score_config(),
                               output_dir = NULL) {
  if (!methods::is(object, "msi_dataset")) {
    stop("object must be a msi_dataset.", call. = FALSE)
  }
  cfg <- .fmc_merge_config(fmc_config)
  if (!isTRUE(cfg$enabled)) return(object)
  ion_mode <- tolower(as.character(ion_mode[[1]]))
  if (!ion_mode %in% c("positive", "negative")) {
    stop("Module 8 requires a resolved ion_mode of positive or negative.", call. = FALSE)
  }
  if (is.null(library_object)) library_object <- object@spectral_library
  ref <- if (is.data.frame(metabolite_reference)) {
    metabolite_reference
  } else {
    build_fmc_metabolite_reference(library_object = library_object)
  }
  reference_diagnostics <- attr(ref, "fmc_reference_diagnostics", exact = TRUE)
  if (!is.data.frame(reference_diagnostics)) {
    reference_diagnostics <- data.frame(
      n_input_library_records = NA_integer_,
      n_records_with_neutral_mass = nrow(ref),
      n_records_without_neutral_mass = NA_integer_,
      n_unique_metabolites = length(unique(ref$metabolite_id)),
      database_name = NA_character_,
      stringsAsFactors = FALSE,
      check.names = FALSE
    )
  }
  reverse_filtered <- if (is.data.frame(object@annotation_results$reverse_match)) {
    object@annotation_results$reverse_match
  } else {
    data.frame()
  }
  reverse_all <- if (is.data.frame(object@annotation_results$reverse_match_all)) {
    object@annotation_results$reverse_match_all
  } else {
    data.frame()
  }
  thresholds <- if (exists(".ann_reverse_evidence_thresholds", mode = "function")) {
    .ann_reverse_evidence_thresholds(object)
  } else {
    list(min_matched_peaks = 3L, min_reverse_masked_cosine = 0)
  }
  legacy_fallback_used <- !is.data.frame(object@annotation_results$mfc_input_candidates)
  if (isTRUE(legacy_fallback_used)) {
    legacy_pseudo <- extract_pseudo_ms2_evidence(
      reverse_filtered, reverse_all, evidence_weight = 50,
      min_matched_peaks_used = thresholds$min_matched_peaks,
      min_reverse_masked_cosine_used = thresholds$min_reverse_masked_cosine
    )
    step7_authoritative_input <- build_mfc_candidate_table(
      pseudo_ms2_match_table = legacy_pseudo,
      ms1_only_match_table = data.frame(),
      feature_metadata = .fmc_feature_table(object),
      expression_matrix = object@expression_data
    )
    step7_authoritative_input <- .ann_standardize_mfc_input_candidates(step7_authoritative_input)
    object@annotation_results$mfc_input_candidates_legacy_rebuilt <- step7_authoritative_input
    object@ms2_data$mfc_input_candidates_legacy_rebuilt <- step7_authoritative_input
    object@annotation_results[["mfc_input_candidates"]] <- FALSE
    object@ms2_data[["mfc_input_candidates"]] <- FALSE
  } else {
    step7_authoritative_input <- object@annotation_results$mfc_input_candidates
  }
  branch_id <- as.character(object@process_info$branch_id %||%
    object@process_info$run_name %||% NA_character_)
  intake_audit <- .fmc_prepare_step7_module8_intake(
    step7_authoritative_input,
    object = object,
    branch = branch_id
  )
  module8_intake_candidates <- intake_audit$intake
  module8_retained_candidates <- intake_audit$retained
  pseudo <- module8_retained_candidates[
    module8_retained_candidates$match_mode %in% "pseudo_ms2",
    ,
    drop = FALSE
  ]
  pseudo_pass_features <- unique(as.character(pseudo$feature_id))
  pseudo_table <- if (is.data.frame(object@pseudo_ms2$pseudo_ms2_table)) {
    object@pseudo_ms2$pseudo_ms2_table
  } else {
    data.frame()
  }
  features_with_pseudo_ms2 <- if (nrow(pseudo_table) && "anchor_feature_id" %in% names(pseudo_table)) {
    unique(as.character(pseudo_table$anchor_feature_id))
  } else if (nrow(reverse_all) && "anchor_feature_id" %in% names(reverse_all)) {
    unique(as.character(reverse_all$anchor_feature_id))
  } else {
    character()
  }
  annotation_input_features <- .fmc_annotation_input_feature_table(
    object = object,
    pseudo_ms2_feature_ids = features_with_pseudo_ms2,
    step7_pseudo_feature_ids = unique(as.character(pseudo$feature_id))
  )
  all_features <- as.character(.fmc_feature_table(object)$feature_id)
  ms1_eligible_features <- switch(
    cfg$ms1_matching_scope,
    all_filtered_features = all_features,
    pseudo_unresolved_only = setdiff(all_features, pseudo_pass_features),
    no_pseudo_spectrum_only = setdiff(all_features, features_with_pseudo_ms2),
    setdiff(all_features, pseudo_pass_features)
  )
  ms1_exclusion_features <- switch(
    cfg$ms1_matching_scope,
    all_filtered_features = character(),
    pseudo_unresolved_only = pseudo_pass_features,
    no_pseudo_spectrum_only = features_with_pseudo_ms2,
    pseudo_pass_features
  )
	  ms1 <- data.frame()
	  ms1_flat_candidate_rows <- 0L
	  ms1_flat_loader_id <- NA_character_
	  ms1_library_path <- cfg$ms1_library_path
	  resolved_ms1_library_path <- NA_character_
	  if (isTRUE(cfg$enable_ms1_fallback)) {
	    ms1_flat <- load_msi_database_ms1_flat(
	      path = ms1_library_path,
	      project_root = as.character(cfg$ms1_flat_project_root %||% getwd()),
	      server_mode = isTRUE(cfg$ms1_flat_server_mode)
	    )
	    ms1_flat_loader_id <- attr(ms1_flat, "ms1_flat_loader_id", exact = TRUE)
	    resolved_ms1_library_path <- attr(ms1_flat, "ms1_flat_library_path", exact = TRUE)
	    ms1 <- generate_ms1_only_candidates(
	      object = object,
	      ms1_library = ms1_flat,
	      evidence_pass_features = ms1_exclusion_features,
      ppm = cfg$ms1_ppm,
      ppm_sigma = cfg$ppm_sigma,
      allowed_adducts = cfg$allowed_adducts,
      polarity = ion_mode,
      max_candidates_per_feature = cfg$max_candidates_per_feature,
      report_top_k = cfg$report_top_k
	    )
	    ms1_flat_candidate_rows <- nrow(ms1)
  }
  if (nrow(ms1)) {
    ms1 <- ms1[as.character(ms1$feature_id) %in% ms1_eligible_features, , drop = FALSE]
  }
  mfc_unified_candidates <- build_mfc_candidate_table(
    pseudo_ms2_match_table = pseudo,
    ms1_only_match_table = ms1,
    feature_metadata = .fmc_feature_table(object),
    expression_matrix = object@expression_data
  )
  isotope <- if (isTRUE(cfg$use_isotope)) {
    annotate_features_by_isotopes(
      object = object,
      fmc_ms1_candidate_table = ms1,
      ppm = cfg$ms1_ppm,
      max_isotope = 2,
      isotope_weight = 10,
      use_spatial_support = cfg$use_spatial_support
    )
  } else {
    .fmc_empty(c("feature_id", "parent_feature_id", "feature_mz", "parent_feature_mz",
                 "metabolite_id", "metabolite_name", "isotope_label", "theoretical_mz",
                 "observed_mz", "mz_error_ppm", "evidence_type", "evidence_weight",
                 "spatial_sas_to_parent", "data_source"))
  }
  neutral <- if (isTRUE(cfg$use_neutral_loss)) {
    annotate_features_by_neutral_loss(
      object = object,
      metabolite_reference = ref,
      ppm = cfg$ms1_ppm,
      ion_mode = ion_mode,
      neutral_loss_rules = neutral_loss_rules
    )
  } else {
    .fmc_empty(c("feature_id", "feature_mz", "metabolite_id", "metabolite_name",
                 "neutral_loss", "theoretical_mz", "mz_error_ppm", "evidence_type",
                 "evidence_weight", "data_source"))
  }
  score_config_runtime <- .fmc_score_config_with_runtime(score_config, cfg)
  scored <- score_fmc_candidates(
    evidence_tables = list(
      pseudo_ms2 = pseudo,
      ms1_adduct = ms1,
      isotope = isotope,
      neutral_loss = neutral
    ),
    object = object,
    feature_graph = object@feature_graph,
    score_config = score_config_runtime
  )
  before <- scored$fmc_candidate_table
  after <- apply_fmc_redundancy_filter(
    before,
    max_metabolites_per_feature = cfg$max_metabolites_per_feature,
    max_features_per_metabolite_soft = cfg$max_features_per_metabolite_soft,
    max_features_per_metabolite_hard = cfg$max_features_per_metabolite_hard
  )
  pure_tissue <- identical(tolower(as.character(cfg$sample_mode %||% "")), "pure_tissue")
  feature_leiden_route <- if (isTRUE(pure_tissue)) {
    .fmc_build_pure_tissue_feature_leiden_route(object, cfg)
  } else {
    NULL
  }
  if (!is.null(feature_leiden_route)) {
    object@feature_graph$graph <- feature_leiden_route$graph
  }
  module8_clusters <- cluster_mfc_by_metabolite(
    mfc_candidate_table = mfc_unified_candidates,
    expression_data = object@expression_data,
    feature_graph = object@feature_graph,
    feature_metadata = .fmc_feature_table(object),
    component_id = if (isTRUE(pure_tissue)) NULL else "component1",
    analysis_domain = if (isTRUE(pure_tissue)) "all_acquired_pixels" else "configured_roi",
    edge_threshold = as.numeric(cfg$spatial_edge_threshold)
  )
  module8_cluster_scores <- score_mfc_clusters(mfc_unified_candidates, module8_clusters)
  module8_dedup <- deduplicate_mfc_table(
    mfc_candidate_table = mfc_unified_candidates,
    mfc_cluster_table = module8_clusters,
    mfc_cluster_score_table = module8_cluster_scores,
    output_dir = NULL
  )
  entropy_process <- object@process_info$module5_entropy_search
  pseudo_library_name <- if (is.list(entropy_process)) {
    as.character(entropy_process$database_name %||% NA_character_)
  } else {
    NA_character_
  }
  mfc_evidence_long <- build_mfc_evidence_long(
    pseudo_ms2_pass_table = pseudo,
    pseudo_ms2_all_table = pseudo,
    ms1_match_table = ms1,
    ion_mode = ion_mode,
    min_matched_peaks = thresholds$min_matched_peaks,
    min_reverse_masked_cosine = thresholds$min_reverse_masked_cosine,
    pseudo_library_name = pseudo_library_name,
    ms1_library_name = as.character(resolved_ms1_library_path[[1]] %||% cfg$ms1_library_path[[1]] %||% NA_character_),
    sample_mode = cfg$sample_mode %||% "standard_mixture",
    pure_tissue_feature_table = if (isTRUE(pure_tissue)) .fmc_feature_table(object) else NULL
  )
  mfc_feature_hmdb_candidates <- merge_mfc_feature_hmdb_candidates(mfc_evidence_long)
  mfc_feature_hmdb_candidates <- score_mfc_priority_candidates(
    feature_hmdb_candidates = mfc_feature_hmdb_candidates,
    mfc_config = cfg,
    cluster_membership = module8_clusters,
    cluster_scores = module8_cluster_scores
  )
  if (!is.null(feature_leiden_route)) {
    candidate_contract <- .fmc_validate_pure_tissue_candidate_contract(
      mfc_feature_hmdb_candidates, cfg
    )
    mfc_feature_hmdb_candidates <- .fmc_apply_pure_tissue_feature_leiden_membership(
      mfc_feature_hmdb_candidates,
      feature_leiden_route$membership
    )
  }
  mfc_priority_filter <- filter_mfc_priority_candidates(mfc_feature_hmdb_candidates, cfg)
  mfc_report_top_k <- mfc_candidate_report_top_k(
    mfc_priority_filter$retained, report_top_k = cfg$report_top_k
  )
  mfc_priority_filter$rejected$removal_reason <- if (nrow(mfc_priority_filter$rejected)) {
    "removed_by_mfc_priority_filter"
  } else {
    character()
  }
  mfc_redundancy_rejections <- module8_dedup$mfc_redundancy_log
  mfc_dedup_rejections <- mfc_redundancy_rejections
  if (nrow(mfc_dedup_rejections)) {
    mfc_dedup_rejections$removal_reason <- "removed_by_mfc_dedup"
  }
  mfc_priority_rejections <- mfc_priority_filter$rejected
  mfc_scoring_parameters <- mfc_scoring_parameter_table(cfg)
  redundancy_summary <- attr(after, "fmc_redundancy_summary", exact = TRUE)
  feature_summary <- .fmc_feature_summary_full(object, pseudo, ms1, isotope, neutral, before, after)
  metabolite_summary <- .fmc_metabolite_summary_after(after)
  evidence_level_summary <- .fmc_evidence_level_summary(after)
  n_input_features <- nrow(.fmc_feature_table(object))
  any_features <- unique(after$feature_id)
  annotation_summary <- data.frame(
    fmc_enabled = TRUE,
    mfc_enabled = TRUE,
    mfc_ms1_fallback_enabled = isTRUE(cfg$enable_ms1_fallback),
    mfc_ms1_fallback_executed = isTRUE(cfg$enable_ms1_fallback) && nrow(ms1) > 0L,
    legacy_ms1_fallback_allowed = FALSE,
    legacy_ms1_fallback_used = FALSE,
    ms1_flat_candidate_rows = as.integer(ms1_flat_candidate_rows),
	  ms1_flat_loader_id = as.character(ms1_flat_loader_id),
    mfc_combined_ms1_pseudo_enabled = isTRUE(cfg$combine_ms1_with_pseudo),
    mfc_ms1_matching_scope = as.character(cfg$ms1_matching_scope),
    fmc_ms1_ppm = as.numeric(cfg$ms1_ppm),
    fmc_ion_mode = ion_mode,
    sample_mode = as.character(cfg$sample_mode %||% "standard_mixture"),
    analysis_domain = as.character(cfg$analysis_domain %||% "configured_roi"),
    n_input_features = as.integer(n_input_features),
    n_features_with_pseudo_ms2_evidence = as.integer(length(unique(pseudo$feature_id))),
    n_features_without_pseudo_ms2 = as.integer(length(setdiff(all_features, features_with_pseudo_ms2))),
    n_features_with_pseudo_ms2_but_no_evidence_pass = as.integer(length(setdiff(features_with_pseudo_ms2, pseudo_pass_features))),
    n_features_with_ms1_adduct_candidates = as.integer(length(unique(ms1$feature_id))),
    n_features_with_isotope_candidates = as.integer(length(unique(isotope$feature_id))),
    n_features_with_any_fmc_candidate = as.integer(length(unique(any_features))),
    n_features_with_no_candidate = as.integer(n_input_features - length(unique(any_features))),
    n_candidate_pairs_before_redundancy = as.integer(nrow(before)),
    n_candidate_pairs_after_redundancy = as.integer(nrow(after)),
    n_mfc_input_candidates = as.integer(nrow(module8_intake_candidates)),
    n_mfc_input_pseudo_ms2 = as.integer(sum(module8_intake_candidates$match_mode == "pseudo_ms2", na.rm = TRUE)),
    n_mfc_input_ms1_only = as.integer(sum(module8_intake_candidates$match_mode %in% c("ms1", "ms1_only"), na.rm = TRUE)),
    n_mfc_unified_candidates = as.integer(nrow(mfc_unified_candidates)),
    step7_candidate_hash = .mfc2_candidate_hash(module8_intake_candidates),
    module8_intake_candidate_hash = .mfc2_candidate_hash(module8_intake_candidates),
    module8_retained_candidate_hash = .mfc2_candidate_hash(module8_retained_candidates),
    legacy_fallback_used = isTRUE(legacy_fallback_used),
    n_mfc_feature_hmdb_candidates = as.integer(nrow(mfc_feature_hmdb_candidates)),
	    n_mfc_retained_candidates = as.integer(nrow(mfc_priority_filter$retained)),
	    n_mfc_report_top_k_candidates = as.integer(nrow(mfc_report_top_k)),
	    report_top_k = as.integer(cfg$report_top_k),
	    n_mfc_rejected_candidates = as.integer(nrow(mfc_priority_filter$rejected)),
    ms1_library_path = as.character(cfg$ms1_library_path[[1]] %||% NA_character_),
	    resolved_ms1_library_path = as.character(resolved_ms1_library_path[[1]] %||% NA_character_),
    ms1_only_definition = paste0("scope=", cfg$ms1_matching_scope),
    score_config_summary = paste(score_config_runtime$parameter, score_config_runtime$value, sep = "=", collapse = ";"),
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
  if (nrow(reference_diagnostics)) {
    annotation_summary$n_reference_metabolites <- reference_diagnostics$n_unique_metabolites[[1]]
    annotation_summary$n_reference_records_with_neutral_mass <- reference_diagnostics$n_records_with_neutral_mass[[1]]
  }
  mfc_result <- list(
    enabled = TRUE,
    mfc_input_candidates = module8_intake_candidates,
    mfc_unified_candidates = mfc_unified_candidates,
    fmc_evidence_table = scored$fmc_evidence_table,
    fmc_candidate_table = after,
    fmc_candidate_table_before_redundancy = before,
    mfc_annotation_input_features = annotation_input_features,
    fmc_feature_summary = feature_summary,
    fmc_metabolite_summary = metabolite_summary,
    fmc_redundancy_summary = redundancy_summary,
    fmc_adduct_rule_table = adduct_rules,
    fmc_score_config = score_config_runtime,
    fmc_annotation_summary = annotation_summary,
    fmc_evidence_level_summary = evidence_level_summary,
    fmc_ms1_candidate_table = ms1,
    module8_mfc_candidate_clusters = module8_cluster_scores,
    module8_mfc_cluster_membership = module8_clusters,
    module8_mfc_representatives = module8_dedup$mfc_deduplicated_table,
    module8_mfc_ambiguity_groups = module8_dedup$mfc_ambiguity_groups,
    module8_mfc_candidates_all = module8_dedup$mfc_full_evidence_table,
    fmc_isotope_candidate_table = isotope,
    fmc_neutral_loss_candidate_table = neutral,
    fmc_metabolite_reference = ref,
    fmc_reference_diagnostics = reference_diagnostics
  )
  mfc_result$step7_to_mfc_stage_summary <- intake_audit$stage_summary
  mfc_result$step7_to_mfc_candidate_audit <- intake_audit$candidate_audit
  mfc_result$step7_to_mfc_retained <- module8_retained_candidates
  mfc_result$step7_to_mfc_rejected <- intake_audit$rejected
  mfc_result$mfc_required_column_rejections <- intake_audit$required_column_rejections
  mfc_result$mfc_redundancy_rejections <- mfc_redundancy_rejections
  mfc_result$mfc_dedup_rejections <- mfc_dedup_rejections
  mfc_result$mfc_priority_rejections <- mfc_priority_rejections
  if (isTRUE(legacy_fallback_used)) {
    mfc_result$legacy_mfc_input_candidates <- module8_intake_candidates
  }
  mfc_result$mfc_evidence_long <- mfc_evidence_long
  mfc_result$mfc_feature_hmdb_candidates <- mfc_feature_hmdb_candidates
  mfc_result$mfc_spatial_clusters <- module8_cluster_scores
  mfc_result$mfc_internal_candidate_pool <- mfc_priority_filter$retained
  mfc_result$mfc_report_top_k_candidates <- mfc_report_top_k
  mfc_result$mfc_filtered_deduplicated <- mfc_priority_filter$retained
  mfc_result$mfc_ambiguity_groups <- mfc_priority_filter$ambiguity_groups
  mfc_result$mfc_rejected_candidates <- mfc_priority_filter$rejected
  mfc_result$mfc_scoring_parameters <- mfc_scoring_parameters
  if (!is.null(feature_leiden_route)) {
    mfc_result$pure_tissue_weighted_leiden_membership <- feature_leiden_route$membership
    mfc_result$pure_tissue_weighted_leiden_edges <- feature_leiden_route$edge_table
    mfc_result$pure_tissue_weighted_leiden_diagnostics <- feature_leiden_route$diagnostics
    mfc_result$pure_tissue_weighted_leiden_config <- feature_leiden_route$config
    mfc_result$pure_tissue_candidate_contract <- candidate_contract
  }
  mfc_result$legacy_field_map <- data.frame(
    legacy_field = c("annotation_results$fmc", "fmc_candidate_table", "module8_mfc_representatives"),
    canonical_field = c("annotation_results$mfc", "mfc_feature_hmdb_candidates", "mfc_filtered_deduplicated"),
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
  object@annotation_results$mfc <- mfc_result
  object@annotation_results$fmc <- mfc_result
  if (!isTRUE(legacy_fallback_used)) {
    object@annotation_results$mfc_input_candidates <- module8_intake_candidates
    object@annotation_results$mfc_input_summary <- annotation_summary
    object@ms2_data$mfc_input_candidates <- module8_intake_candidates
    object@ms2_data$mfc_input_summary <- annotation_summary
  }
  if (!is.null(output_dir) && nzchar(as.character(output_dir[[1]]))) {
    table_dir <- file.path(as.character(output_dir[[1]]), "tables")
    dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)
    utils::write.csv(module8_dedup$mfc_full_evidence_table, file.path(table_dir, "module8_mfc_candidates_all.csv"), row.names = FALSE, na = "")
    utils::write.csv(annotation_input_features, file.path(table_dir, "mfc_annotation_input_features.csv"), row.names = FALSE, na = "")
    utils::write.csv(module8_intake_candidates, file.path(table_dir, "mfc_input_candidates.csv"), row.names = FALSE, na = "")
    utils::write.csv(mfc_unified_candidates, file.path(table_dir, "mfc_unified_candidates.csv"), row.names = FALSE, na = "")
    utils::write.csv(intake_audit$stage_summary, file.path(table_dir, "step7_to_mfc_stage_summary.csv"), row.names = FALSE, na = "")
    utils::write.csv(intake_audit$candidate_audit, file.path(table_dir, "step7_to_mfc_candidate_audit.csv"), row.names = FALSE, na = "")
    utils::write.csv(module8_retained_candidates, file.path(table_dir, "step7_to_mfc_retained.csv"), row.names = FALSE, na = "")
    utils::write.csv(intake_audit$rejected, file.path(table_dir, "step7_to_mfc_rejected.csv"), row.names = FALSE, na = "")
    utils::write.csv(intake_audit$required_column_rejections, file.path(table_dir, "mfc_required_column_rejections.csv"), row.names = FALSE, na = "")
    utils::write.csv(mfc_redundancy_rejections, file.path(table_dir, "mfc_redundancy_rejections.csv"), row.names = FALSE, na = "")
    utils::write.csv(mfc_dedup_rejections, file.path(table_dir, "mfc_dedup_rejections.csv"), row.names = FALSE, na = "")
    utils::write.csv(mfc_priority_rejections, file.path(table_dir, "mfc_priority_rejections.csv"), row.names = FALSE, na = "")
    utils::write.csv(mfc_priority_filter$retained, file.path(table_dir, "mfc_final_candidates.csv"), row.names = FALSE, na = "")
    utils::write.csv(mfc_priority_filter$retained, file.path(table_dir, "mfc_final_candidates_all_ranked.csv"), row.names = FALSE, na = "")
    utils::write.csv(mfc_report_top_k, file.path(table_dir, "mfc_final_candidates_top5_report.csv"), row.names = FALSE, na = "")
    utils::write.csv(module8_cluster_scores, file.path(table_dir, "module8_mfc_candidate_clusters.csv"), row.names = FALSE, na = "")
    utils::write.csv(module8_dedup$mfc_deduplicated_table, file.path(table_dir, "module8_mfc_representatives.csv"), row.names = FALSE, na = "")
    utils::write.csv(module8_dedup$mfc_ambiguity_groups, file.path(table_dir, "module8_mfc_ambiguity_groups.csv"), row.names = FALSE, na = "")
    utils::write.csv(ms1, file.path(table_dir, "module8_ms1_only_candidates.csv"), row.names = FALSE, na = "")
    utils::write.csv(ms1, file.path(table_dir, "module8_ms1_only_candidates_all_ranked.csv"), row.names = FALSE, na = "")
    utils::write.csv(ms1_candidate_report_top_k(ms1, cfg$report_top_k),
      file.path(table_dir, "module8_ms1_only_candidates_top5_report.csv"), row.names = FALSE, na = "")
    utils::write.csv(mfc_evidence_long, file.path(table_dir, "mfc_evidence_long.csv"), row.names = FALSE, na = "")
    utils::write.csv(mfc_feature_hmdb_candidates, file.path(table_dir, "mfc_feature_hmdb_candidates.csv"), row.names = FALSE, na = "")
    utils::write.csv(module8_cluster_scores, file.path(table_dir, "mfc_spatial_clusters.csv"), row.names = FALSE, na = "")
    utils::write.csv(mfc_priority_filter$retained, file.path(table_dir, "mfc_filtered_deduplicated.csv"), row.names = FALSE, na = "")
    utils::write.csv(mfc_priority_filter$ambiguity_groups, file.path(table_dir, "mfc_ambiguity_groups.csv"), row.names = FALSE, na = "")
    utils::write.csv(mfc_priority_filter$rejected, file.path(table_dir, "mfc_rejected_candidates.csv"), row.names = FALSE, na = "")
    utils::write.csv(mfc_scoring_parameters, file.path(table_dir, "mfc_scoring_parameters.csv"), row.names = FALSE, na = "")
    if (!is.null(feature_leiden_route)) {
      utils::write.csv(feature_leiden_route$membership,
        file.path(table_dir, "pure_tissue_weighted_leiden_membership.csv"), row.names = FALSE, na = "")
      utils::write.csv(feature_leiden_route$edge_table,
        file.path(table_dir, "pure_tissue_weighted_leiden_edges.csv"), row.names = FALSE, na = "")
    }
    utils::write.csv(mfc_result$legacy_field_map, file.path(table_dir, "mfc_legacy_field_map.csv"), row.names = FALSE, na = "")
    evidence_summary <- data.frame(
      evidence_type = c("pseudo_ms2", "ms1_only"),
      n_rows = c(sum(mfc_unified_candidates$match_mode == "pseudo_ms2", na.rm = TRUE),
                 sum(mfc_unified_candidates$match_mode == "ms1_only", na.rm = TRUE)),
      evidence_class_weight = c(1.0, 0.6),
      stringsAsFactors = FALSE,
      check.names = FALSE
    )
    utils::write.csv(evidence_summary, file.path(table_dir, "module8_evidence_source_summary.csv"), row.names = FALSE, na = "")
  }
  object@process_info$module8_mfc <- c(
    as.list(annotation_summary[1, , drop = FALSE]),
    list(
      resolved_ion_mode = ion_mode,
      effective_mfc_config = cfg,
      scoring_parameters = mfc_scoring_parameters,
      legacy_fmc_alias_written = TRUE,
      legacy_fallback_used = isTRUE(legacy_fallback_used),
      step7_candidate_hash = .mfc2_candidate_hash(module8_intake_candidates),
      module8_intake_candidate_hash = .mfc2_candidate_hash(module8_intake_candidates),
      module8_retained_candidate_hash = .mfc2_candidate_hash(module8_retained_candidates),
      step7_candidate_row_count = as.integer(nrow(module8_intake_candidates)),
      module8_intake_row_count = as.integer(nrow(module8_intake_candidates)),
      module8_retained_row_count = as.integer(nrow(module8_retained_candidates))
    )
  )
  object
}

`score_fmc_candidates` <- 
function(evidence_tables,
                                 object = NULL,
                                 feature_graph = NULL,
                                 score_config = default_fmc_score_config()) {
  evidence <- .fmc_standardize_evidence(evidence_tables)
  empty_candidate <- .fmc_empty(c("feature_id", "feature_mz", "metabolite_id",
                                 "metabolite_name", "evidence_types",
                                 "best_evidence_type", "all_adducts",
                                 "score_pseudo_ms2", "score_adduct",
                                 "score_isotope", "score_neutral_loss",
                                 "score_mz_accuracy", "score_spatial",
                                 "score_redundancy_penalty", "score_total",
                                 "fmc_evidence_level", "rank_within_feature",
                                 "rank_within_metabolite"))
  if (!nrow(evidence)) {
    return(list(
      fmc_evidence_table = evidence,
      fmc_candidate_table = empty_candidate,
      fmc_feature_summary = .fmc_empty(c("feature_id", "feature_mz", "n_fmc_candidates")),
      fmc_metabolite_summary = .fmc_empty(c("metabolite_id", "metabolite_name", "n_features"))
    ))
  }
  ppm_tolerance <- suppressWarnings(as.numeric(.fmc_config_value(score_config, "ppm_tolerance", 10)))
  if (!is.finite(ppm_tolerance) || ppm_tolerance <= 0) ppm_tolerance <- 10
  mz_accuracy_max <- suppressWarnings(as.numeric(.fmc_config_value(score_config, "mz_accuracy_max_score", 10)))
  if (!is.finite(mz_accuracy_max) || mz_accuracy_max < 0) mz_accuracy_max <- 10
  spatial_weight <- suppressWarnings(as.numeric(.fmc_config_value(score_config, "spatial_weight", 15)))
  if (!is.finite(spatial_weight) || spatial_weight < 0) spatial_weight <- 15

  pair_key <- paste(evidence$feature_id, evidence$metabolite_id, sep = "\r")
  groups <- split(seq_len(nrow(evidence)), pair_key, drop = TRUE)
  rows <- lapply(groups, function(idx) {
    sub <- evidence[idx, , drop = FALSE]
    type_groups <- split(seq_len(nrow(sub)), as.character(sub$evidence_type), drop = TRUE)
    best_by_type <- lapply(type_groups, function(j) {
      w <- suppressWarnings(as.numeric(sub$evidence_weight[j]))
      j[[which.max(ifelse(is.finite(w), w, -Inf))]]
    })
    best_idx <- unlist(best_by_type, use.names = FALSE)
    best <- sub[best_idx, , drop = FALSE]
    score_type <- function(type) {
      hit <- best[as.character(best$evidence_type) == type, , drop = FALSE]
      if (nrow(hit)) max(suppressWarnings(as.numeric(hit$evidence_weight)), na.rm = TRUE) else 0
    }
    adduct_rows <- sub[as.character(sub$evidence_type) == "ms1_adduct", , drop = FALSE]
    score_adduct <- if (nrow(adduct_rows)) {
      max(suppressWarnings(as.numeric(adduct_rows$evidence_weight)), na.rm = TRUE)
    } else {
      0
    }
    best_weight <- suppressWarnings(as.numeric(best$evidence_weight))
    best_evidence <- as.character(best$evidence_type[[which.max(ifelse(is.finite(best_weight), best_weight, -Inf))]])
    mz_error <- .fmc_min_abs(sub$mz_error_ppm)
    score_mz <- if (is.finite(mz_error)) {
      mz_accuracy_max * (1 - abs(mz_error) / ppm_tolerance)
    } else {
      0
    }
    score_mz <- max(0, min(mz_accuracy_max, score_mz))
    spatial_values <- suppressWarnings(as.numeric(sub$spatial_sas_to_parent))
    spatial_values <- spatial_values[is.finite(spatial_values)]
    score_spatial <- if (length(spatial_values)) spatial_weight * max(spatial_values) else 0
    score_pseudo <- max(score_type("pseudo_ms2"), score_type("pseudo_ms2_strict"))
    score_iso <- score_type("isotope")
    score_nl <- score_type("neutral_loss")
    penalty <- 0
    score_total <- score_pseudo + score_adduct + score_iso + score_nl + score_mz + score_spatial - penalty
    data.frame(
      feature_id = as.character(sub$feature_id[[1]]),
      feature_mz = suppressWarnings(as.numeric(.fmc_first_nonmissing(sub$feature_mz, NA_real_))),
      metabolite_id = as.character(sub$metabolite_id[[1]]),
      metabolite_name = as.character(.fmc_first_nonmissing(sub$metabolite_name, NA_character_)),
      evidence_types = .fmc_collapse_unique(best$evidence_type),
      best_evidence_type = best_evidence,
      all_adducts = .fmc_collapse_unique(sub$adduct),
      score_pseudo_ms2 = score_pseudo,
      score_adduct = score_adduct,
      score_isotope = score_iso,
      score_neutral_loss = score_nl,
      score_mz_accuracy = score_mz,
      score_spatial = score_spatial,
      score_redundancy_penalty = penalty,
      score_total = score_total,
      fmc_evidence_level = .fmc_evidence_level(best$evidence_type, sub$adduct, sub$rank_group, score_spatial),
      rank_within_feature = NA_integer_,
      rank_within_metabolite = NA_integer_,
      stringsAsFactors = FALSE,
      check.names = FALSE
    )
  })
  candidate <- do.call(rbind, rows)
  candidate <- candidate[order(as.character(candidate$feature_id), -candidate$score_total,
                               as.character(candidate$metabolite_id)), , drop = FALSE]
  candidate$rank_within_feature <- .fmc_rank(candidate, "feature_id", "score_total")
  candidate$rank_within_metabolite <- .fmc_rank(candidate, "metabolite_id", "score_total")
  rownames(candidate) <- NULL

  features <- if (methods::is(object, "msi_dataset")) .fmc_feature_table(object) else {
    unique(candidate[, c("feature_id", "feature_mz"), drop = FALSE])
  }
  n_by_feature <- table(candidate$feature_id)
  feature_summary <- merge(features, data.frame(
    feature_id = names(n_by_feature),
    n_fmc_candidates = as.integer(n_by_feature),
    stringsAsFactors = FALSE,
    check.names = FALSE
  ), by = "feature_id", all.x = TRUE, sort = FALSE)
  feature_summary$n_fmc_candidates[is.na(feature_summary$n_fmc_candidates)] <- 0L
  top_idx <- which(candidate$rank_within_feature == 1L)
  top <- candidate[top_idx, c("feature_id", "metabolite_id", "metabolite_name",
                             "score_total", "fmc_evidence_level"), drop = FALSE]
  names(top) <- c("feature_id", "top_metabolite_id", "top_metabolite_name",
                  "top_score_total", "top_fmc_evidence_level")
  feature_summary <- merge(feature_summary, top, by = "feature_id", all.x = TRUE, sort = FALSE)

  metabolite_groups <- split(seq_len(nrow(candidate)), as.character(candidate$metabolite_id), drop = TRUE)
  metabolite_summary <- do.call(rbind, lapply(metabolite_groups, function(idx) {
    sub <- candidate[idx, , drop = FALSE]
    data.frame(
      metabolite_id = as.character(sub$metabolite_id[[1]]),
      metabolite_name = as.character(sub$metabolite_name[[1]]),
      n_features = length(unique(sub$feature_id)),
      max_score_total = max(sub$score_total, na.rm = TRUE),
      best_fmc_evidence_level = as.character(sub$fmc_evidence_level[[which.max(sub$score_total)]]),
      stringsAsFactors = FALSE,
      check.names = FALSE
    )
  }))
  rownames(metabolite_summary) <- NULL

  list(
    fmc_evidence_table = evidence,
    fmc_candidate_table = candidate[, colnames(empty_candidate), drop = FALSE],
    fmc_feature_summary = feature_summary,
    fmc_metabolite_summary = metabolite_summary
  )
}

`score_mfc_clusters` <- 
function(mfc_candidate_table,
                               mfc_cluster_table,
                               spatial_similarity = attr(mfc_cluster_table, "mfc_spatial_similarity", exact = TRUE)) {
  cand <- .mfc2_df(mfc_candidate_table)
  cl <- .mfc2_df(mfc_cluster_table)
  if (!nrow(cl)) {
    return(data.frame(
      mfc_id = character(), hmdb_id = character(), metabolite_id = character(),
      compound_name = character(), feature_ids = character(), mfc_spatial_coherence = numeric(),
      cluster_medoid_feature = character(), adduct_diversity = integer(),
      mean_SAS = numeric(), median_SAS = numeric(),
      mfc_best_chem_score = numeric(), mfc_best_S_SAS = numeric(),
      mfc_adduct_consistency = numeric(), mfc_evidence_count = integer(),
      mfc_ms1_only_count = integer(), mfc_pseudo_ms2_count = integer(),
      mfc_score = numeric(), cluster_support_score = numeric(),
      stringsAsFactors = FALSE, check.names = FALSE
    ))
  }
  rows <- lapply(split(seq_len(nrow(cl)), cl$mfc_id, drop = TRUE), function(idx) {
    sub_cl <- cl[idx, , drop = FALSE]
    features <- unique(as.character(sub_cl$feature_id))
    evidence <- cand[cand$hmdb_id %in% sub_cl$hmdb_id & cand$feature_id %in% features, , drop = FALSE]
    chem_values <- .mfc2_scale01(c(
      evidence$entropy_score,
      evidence$reverse_masked_cosine,
      evidence$score_total_pre_mfc,
      evidence$raw_score
    ))
    chem_values <- chem_values[is.finite(chem_values)]
    chem <- if (length(chem_values)) max(chem_values) else 0
    best_sas_values <- .mfc2_scale01(evidence$S_SAS)
    best_sas_values <- best_sas_values[is.finite(best_sas_values)]
    best_sas <- if (length(best_sas_values)) max(best_sas_values) else 0
    adducts <- unique(as.character(evidence$adduct[!is.na(evidence$adduct) & nzchar(evidence$adduct)]))
    adduct_consistency <- if (length(adducts) == 0L) 0.5 else if (length(adducts) <= 3L) 1 else 0.6
    coherence <- .mfc2_cluster_similarity(features, spatial_similarity)
    sas_values <- .mfc2_scale01(evidence$S_SAS)
    sas_values <- sas_values[is.finite(sas_values)]
    mean_sas <- if (length(sas_values)) mean(sas_values) else NA_real_
    median_sas <- if (length(sas_values)) stats::median(sas_values) else NA_real_
    medoid <- .mfc2_cluster_medoid(features, spatial_similarity)
    evidence_count <- nrow(evidence)
    scaled_count <- min(1, log1p(evidence_count) / log1p(10))
    score <- 0.50 * chem + 0.25 * coherence + 0.15 * adduct_consistency + 0.10 * scaled_count
    data.frame(
      mfc_id = as.character(sub_cl$mfc_id[[1]]),
      hmdb_id = as.character(sub_cl$hmdb_id[[1]]),
      metabolite_id = as.character(sub_cl$metabolite_id[[1]]),
      compound_name = as.character(sub_cl$compound_name[[1]]),
      component_id = as.character(sub_cl$component_id[[1]]),
      cluster_index = as.integer(sub_cl$cluster_index[[1]]),
      feature_ids = paste(features, collapse = ";"),
      n_features = length(features),
      cluster_medoid_feature = medoid,
      adduct_diversity = length(adducts),
      mean_SAS = mean_sas,
      median_SAS = median_sas,
      mfc_spatial_coherence = coherence,
      mfc_best_chem_score = chem,
      mfc_best_S_SAS = best_sas,
      mfc_adduct_consistency = adduct_consistency,
      mfc_evidence_count = evidence_count,
      mfc_ms1_only_count = sum(grepl("ms1", evidence$evidence_type, ignore.case = TRUE)),
      mfc_pseudo_ms2_count = sum(grepl("pseudo", evidence$evidence_type, ignore.case = TRUE)),
      mfc_score = pmin(1, pmax(0, score)),
      cluster_support_score = pmin(1, pmax(0, score)),
      stringsAsFactors = FALSE,
      check.names = FALSE
    )
  })
  out <- do.call(rbind, rows)
  out <- out[order(out$hmdb_id, -out$mfc_score, out$mfc_id), , drop = FALSE]
  out$mfc_rank_within_metabolite <- ave(-out$mfc_score, out$hmdb_id, FUN = function(x) rank(x, ties.method = "first"))
  rownames(out) <- NULL
  out
}

`score_mfc_priority_candidates` <- 
function(feature_hmdb_candidates,
                                          mfc_config = list(),
                                          cluster_membership = NULL,
                                          cluster_scores = NULL) {
  x <- as.data.frame(feature_hmdb_candidates, stringsAsFactors = FALSE, check.names = FALSE)
  if (!nrow(x)) return(x)
  defaults <- list(
    tier3_score_band_min = 0.67,
    tier2_score_band_min = 0.34,
    tier2_score_band_max = 0.66,
    tier1_score_band_max = 0.33,
    pseudo_score_band_min = 0.34,
    ms1_only_score_band_max = 0.33,
    pseudo_raw_weight = 0.70,
    ms1_support_weight = 0.10,
    spatial_weight = 0.20,
    cluster_support_weight = 0.10
  )
  cfg <- defaults
  if (is.list(mfc_config)) for (nm in names(mfc_config)) cfg[[nm]] <- mfc_config[[nm]]
  numeric_names <- names(defaults)
  for (nm in numeric_names) {
    cfg[[nm]] <- as.numeric(cfg[[nm]][[1]])
    if (!is.finite(cfg[[nm]]) || cfg[[nm]] < 0) cfg[[nm]] <- defaults[[nm]]
  }
  if (!(0 <= cfg$tier1_score_band_max &&
        cfg$tier1_score_band_max < cfg$tier2_score_band_min &&
        cfg$tier2_score_band_min <= cfg$tier2_score_band_max &&
        cfg$tier2_score_band_max < cfg$tier3_score_band_min &&
        cfg$tier3_score_band_min <= 1)) {
    stop("MFC score bands must satisfy 0 <= tier1_max < tier2_min <= tier2_max < tier3_min <= 1.", call. = FALSE)
  }
  x$evidence_source <- ifelse(
    x$has_pseudo_ms2 & x$has_ms1,
    "ms1_and_pseudo_ms2",
    ifelse(x$has_pseudo_ms2, "pseudo_ms2_only", "ms1_only")
  )

  x$spatial_cluster_id <- NA_character_
  x$cluster_coherence <- 0
  x$cluster_support_score <- 0
  membership <- as.data.frame(cluster_membership, stringsAsFactors = FALSE, check.names = FALSE)
  scores <- as.data.frame(cluster_scores, stringsAsFactors = FALSE, check.names = FALSE)
  if (nrow(membership) && all(c("feature_id", "hmdb_id", "mfc_id") %in% names(membership))) {
    map <- membership[, c("feature_id", "hmdb_id", "mfc_id"), drop = FALSE]
    names(map) <- c("feature_id", "HMDB_ID", "spatial_cluster_id")
    map <- map[!duplicated(paste(map$feature_id, map$HMDB_ID, sep = "\r")), , drop = FALSE]
    x <- merge(x, map, by = c("feature_id", "HMDB_ID"), all.x = TRUE, sort = FALSE,
               suffixes = c("", ".mapped"))
    if ("spatial_cluster_id.mapped" %in% names(x)) {
      x$spatial_cluster_id <- x$spatial_cluster_id.mapped
      x$spatial_cluster_id.mapped <- NULL
    }
  }
  if (nrow(scores) && "mfc_id" %in% names(scores)) {
    coherence <- suppressWarnings(as.numeric(scores$mfc_spatial_coherence))
    support <- suppressWarnings(as.numeric(scores$mfc_score))
    x$cluster_coherence <- .mfc_unified_clamp01(coherence[match(x$spatial_cluster_id, scores$mfc_id)])
    x$cluster_support_score <- .mfc_unified_clamp01(support[match(x$spatial_cluster_id, scores$mfc_id)])
  }
  x$spatial_support <- .mfc_unified_clamp01(x$best_S_SAS)
  pseudo_raw <- .mfc_unified_clamp01(x$best_pseudo_ms2_score)
  ms1_raw <- .mfc_unified_clamp01(x$best_ms1_score)
  pseudo_den <- cfg$pseudo_raw_weight + cfg$spatial_weight + cfg$cluster_support_weight +
    ifelse(x$has_ms1, cfg$ms1_support_weight, 0)
  pseudo_supported <- (
    cfg$pseudo_raw_weight * pseudo_raw +
      cfg$ms1_support_weight * ms1_raw * as.numeric(x$has_ms1) +
      cfg$spatial_weight * x$spatial_support +
      cfg$cluster_support_weight * x$cluster_support_score
  ) / pmax(pseudo_den, .Machine$double.eps)
  ms1_den <- 1 + cfg$spatial_weight + cfg$cluster_support_weight
  ms1_supported <- (
    ms1_raw + cfg$spatial_weight * x$spatial_support +
      cfg$cluster_support_weight * x$cluster_support_score
  ) / ms1_den
  x$pseudo_supported_raw <- ifelse(x$has_pseudo_ms2, .mfc_unified_clamp01(pseudo_supported), NA_real_)
  x$ms1_only_raw <- ifelse(!x$has_pseudo_ms2, .mfc_unified_clamp01(ms1_supported), NA_real_)
  x$mfc_raw_score <- ifelse(
    x$has_pseudo_ms2,
    .mfc_unified_clamp01(pseudo_supported),
    .mfc_unified_clamp01(ms1_supported)
  )
  x$evidence_tier <- ifelse(
    x$has_pseudo_ms2 & x$has_ms1,
    "Tier 1: ms1_and_pseudo_ms2",
    ifelse(x$has_pseudo_ms2, "Tier 2: pseudo_ms2_only", "Tier 3: ms1_only")
  )
  x$mfc_priority_score <- ifelse(
    x$has_pseudo_ms2 & x$has_ms1,
    cfg$tier3_score_band_min + (1 - cfg$tier3_score_band_min) * x$mfc_raw_score,
    ifelse(
      x$has_pseudo_ms2,
      cfg$tier2_score_band_min +
        (cfg$tier2_score_band_max - cfg$tier2_score_band_min) * x$mfc_raw_score,
      cfg$tier1_score_band_max * x$mfc_raw_score
    )
  )
  x$score_semantics <- "priority_score_not_probability"
  x <- x[order(x$feature_id, -x$mfc_priority_score, x$HMDB_ID), , drop = FALSE]
  x$rank_within_feature <- ave(-x$mfc_priority_score, x$feature_id, FUN = function(v) rank(v, ties.method = "first"))
  best <- ave(x$mfc_priority_score, x$feature_id, FUN = max)
  x$score_gap <- best - x$mfc_priority_score
  x$relative_support <- ifelse(best > 0, x$mfc_priority_score / best, 0)
  rownames(x) <- NULL
  x
}

`validate_adduct_rule_table` <- 
function(adduct_rules) {
  required <- c("ion_mode", "adduct", "mass_shift", "charge", "evidence_type",
                "evidence_weight", "rank_group", "description", "enabled")
  if (!is.data.frame(adduct_rules)) {
    stop("adduct_rules must be a data frame.", call. = FALSE)
  }
  missing <- setdiff(required, names(adduct_rules))
  if (length(missing)) {
    stop("adduct_rules missing required columns: ", paste(missing, collapse = ", "), call. = FALSE)
  }
  ion_mode <- tolower(as.character(adduct_rules$ion_mode))
  if (any(!ion_mode %in% c("positive", "negative"))) {
    stop("adduct_rules$ion_mode must be positive or negative.", call. = FALSE)
  }
  weight <- suppressWarnings(as.numeric(adduct_rules$evidence_weight))
  if (any(!is.finite(weight) | weight < 0)) {
    stop("adduct_rules$evidence_weight must be non-negative finite values.", call. = FALSE)
  }
  charge <- suppressWarnings(as.numeric(adduct_rules$charge))
  needs_charge <- as.character(adduct_rules$evidence_type) %in% c("ms1_adduct", "isotope", "neutral_loss")
  if (any(needs_charge & (!is.finite(charge) | charge == 0))) {
    stop("ms1_adduct/isotope/neutral_loss rules require finite non-zero charge.", call. = FALSE)
  }
  invisible(adduct_rules)
}

