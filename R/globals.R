#' @importFrom stats ave
#' @importFrom utils flush.console tail
NULL

# Legacy data-frame evaluation uses these non-standard evaluation symbols.
utils::globalVariables(c('cluster_id', 'cluster_index', 'cluster_size',
  'feature_id', 'from', 'member_key', 'raw_cluster', 'to', 'weight'))
