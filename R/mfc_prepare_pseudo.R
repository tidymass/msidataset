# Convert the pseudo-MS2 pair table to the scorer query format. Invalid rows
# remain in the returned table with an explicit status for the caller.
.mfc_prepare_pseudo <- function(object, pairs) {
  i <- .msidataset_inputs(object)
  required <- c("precursor_feature_id", "precursor_mz", "fragment_feature_id", "fragment_mz", "sas_d7")
  if (!is.data.frame(pairs) || !all(required %in% names(pairs)))
    stop("Pseudo-MS2 pair columns are missing", call. = FALSE)
  if (any(!is.finite(i$tic) | i$tic <= 0)) stop("TIC must be finite and positive", call. = FALSE)
  values <- i$intensity * (mean(i$tic) / i$tic)
  med <- vapply(seq_len(ncol(values)), function(j) {
    v <- values[, j]
    if (any(!is.finite(v))) return(NA_real_)
    v <- v[v > 0]
    if (length(v)) stats::median(v) else 0
  }, numeric(1))
  idx <- match(pairs$fragment_feature_id, i$feature_id)
  precursor_idx <- match(pairs$precursor_feature_id, i$feature_id)
  out <- pairs
  # Saved pseudo-MS2 pairs may contain rounded labels. Use measured centroids
  # where available and preserve the supplied masses for older inputs.
  full <- .msidataset_feature_mz(object@variable_info)$observed_mz_full
  precursor_full <- is.finite(full[precursor_idx])
  fragment_full <- is.finite(full[idx])
  out$precursor_mz[precursor_full] <- full[precursor_idx[precursor_full]]
  out$fragment_mz[fragment_full] <- full[idx[fragment_full]]
  out$SAS <- out$sas_d7
  out$ms1 <- med[idx]
  max_per <- ave(out$ms1, out$precursor_feature_id, FUN = function(x) {
    if (all(is.na(x))) rep(NA_real_, length(x)) else rep(max(x, na.rm = TRUE), length(x))
  })
  out$intensity <- ifelse(is.na(max_per) | is.na(out$ms1), NA_real_,
                          ifelse(max_per > 0, 100 * out$ms1 / max_per, 0))
  out$peak_id <- paste0("P_", out$precursor_feature_id)
  out$anchor_feature_id <- out$precursor_feature_id
  out$anchor_mz <- out$precursor_mz
  out$fragment_id <- paste0(out$peak_id, "_", out$fragment_feature_id)
  out$cluster_id <- out$peak_id
  out$pseudo_ms2_id <- out$peak_id
  out$peak_feature_id <- out$fragment_feature_id
  out$peak_mz <- out$fragment_mz
  out$graph_weight <- out$SAS
  out$fragment_sas <- out$SAS
  out$pseudo_raw_intensity <- out$intensity
  out$pseudo_intensity <- out$intensity
  out$is_self_peak <- FALSE
  out$peak_role <- "fragment"
  out$preparation_status <- ifelse(is.na(precursor_idx), "missing_precursor_feature",
    ifelse(is.na(idx), "missing_fragment_feature",
      ifelse(!is.finite(out$SAS), "missing_sas", ifelse(!is.finite(out$intensity), "missing_intensity", "ok"))))
  out
}
