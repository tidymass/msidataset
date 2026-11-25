#' Run the MSI analysis in cached stages
#'
#' Evidence construction consumes the
#' freshly computed upstream results and external pathway resources. A checkpoint includes
#' hashes of the inputs, the package code and preceding checkpoint. Changing
#' an input invalidates the corresponding chain. A requested `from` stage
#' requires valid earlier checkpoints.
#'
#' @param config Named list with `input` (arguments to [read_msi_data()]),
#'   `labels` (data frame or TSV path), `step6` (RDS path or step6 object),
#'   `vector_library` (directory), and `mfc_output_dir` (writable directory).
#'   If `step6` is absent, `match_spectral_library()` must be available and
#'   `spectral` must contain its configuration. For the evidence stage,
#'   `evidence` provides `output_dir`, `hmdb_pathways`, and `reactome_mapping`.
#'   The vector library's ids.tsv provides taxonomy metadata.
#' @param from First stage to execute; earlier stages must already be cached.
#' @param to Last stage to execute. Use `"build_llm_evidence"` to complete
#'   evidence construction. The default reaches the annotation stage, which
#'   requires calling [annotate_clusters()] separately with an evidence view.
#' @param cache_dir Directory for checkpoint files.
#' @return A list containing the latest results, stage history, and cache keys.
#' @examples
#' \dontrun{run_pipeline(config, to = "retrieve_neighbors", cache_dir = "cache")}
#' @export
run_pipeline <- function(config, from = NULL, to = NULL, cache_dir) {
  stages <- c("read_msi_data", "define_roi", "detect_de_features", "generate_pseudo_ms2",
              "match_spectral_library", "calculate_mfc", "embed_features", "select_k",
              "cluster_features", "select_representatives", "retrieve_neighbors",
              "build_llm_evidence", "annotate_clusters")
  if (!is.list(config) || is.null(names(config))) stop("config must be a named list")
  if (length(cache_dir) != 1L || is.na(cache_dir) || !nzchar(cache_dir)) stop("Invalid cache_dir")
  if (is.null(from)) from <- stages[1L]
  if (is.null(to)) to <- utils::tail(stages, 1L)
  first <- match(from, stages); last <- match(to, stages)
  if (length(first) != 1L || length(last) != 1L || is.na(first) || is.na(last) || first > last)
    stop("Invalid from/to stages")
  dir.create(cache_dir, recursive = TRUE, showWarnings = FALSE)
  if (!dir.exists(cache_dir)) stop("Cannot create cache_dir")
  # Hash every existing path supplied in the configuration, including vector
  # files: metadata alone must not make a stale result appear valid.
  input_config <- config
  input_config$mfc_output_dir <- NULL
  if (!is.null(input_config$spectral)) input_config$spectral$output_dir <- NULL
  if (!is.null(input_config$evidence)) input_config$evidence$output_dir <- NULL
  paths <- unique(unlist(input_config, recursive = TRUE, use.names = FALSE))
  paths <- paths[vapply(paths, function(p) is.character(p) && length(p) == 1L && !is.na(p) && file.exists(p), logical(1))]
  locked_files <- character()
  lock_path <- system.file("config", "locked_d15.yaml", package = "msidataset")
  if (last >= match("match_spectral_library", stages) && length(lock_path) == 1L && nzchar(lock_path) && file.exists(lock_path)) {
    locked <- .msi_locked_resources(lock_path)
    keys <- c("adp_config", "ms1_library", "vector_ids")
    if (is.null(config$step6)) keys <- c(keys, "ms2_library")
    locked_files <- vapply(locked[keys], function(x) if (is.null(x$path)) "" else x$path, character(1))
    locked_files <- c(lock_path, locked_files[file.exists(locked_files)])
  }
  files <- unique(unlist(lapply(paths, function(p) {
    if (dir.exists(p)) list.files(p, recursive = TRUE, full.names = TRUE, all.files = TRUE, no.. = TRUE)
    else p
  }), use.names = FALSE))
  files <- sort(files[file.exists(files) & !dir.exists(files)])
  file_hash <- vapply(files, function(p) digest::digest(file = p, algo = "sha256"), character(1))
  locked_hash <- vapply(locked_files, function(p) digest::digest(file = p, algo = "sha256"), character(1))
  ns <- environment(run_pipeline)
  fns <- mget(ls(ns, all.names = TRUE), envir = ns, inherits = FALSE)
  fns <- fns[vapply(fns, is.function, logical(1))]
  code <- lapply(fns, function(f) list(formals = formals(f), body = body(f)))
  base_key <- digest::digest(list(config = config, files = file_hash, code = code), algo = "sha256")
  state <- list(process_info = list(), cache_keys = list())
  previous_key <- base_key
  for (i in seq_len(last)) {
    step <- stages[i]
    key <- digest::digest(list(previous_key, step,
      if (identical(step, "match_spectral_library")) locked_hash else NULL), algo = "sha256")
    path <- file.path(cache_dir, paste0(sprintf("%02d", i), "_", step, ".rds"))
    checkpoint <- if (file.exists(path)) tryCatch(readRDS(path), error = function(e) NULL) else NULL
    valid <- is.list(checkpoint) && identical(checkpoint$key, key) && is.list(checkpoint$state)
    if (valid) {
      state <- checkpoint$state
      previous_key <- key
      next
    }
    if (i < first) stop("Missing or stale checkpoint before from: ", step)
    started <- as.character(Sys.time())
    state <- switch(step,
      read_msi_data = { state$object <- do.call(read_msi_data, config$input); state },
      define_roi = {
        labels <- config$labels
        if (is.character(labels) && length(labels) == 1L) labels <- data.table::fread(labels, data.table = FALSE, showProgress = FALSE)
        state$object <- define_roi(state$object, labels); state
      },
      detect_de_features = { state$object <- detect_de_features(state$object); state },
      generate_pseudo_ms2 = { state$object <- generate_pseudo_ms2(state$object); state },
      match_spectral_library = {
        if (!is.null(config$step6)) {
          state$step6 <- if (is.character(config$step6) && length(config$step6) == 1L) readRDS(config$step6) else config$step6
        } else {
          if (!exists("match_spectral_library", mode = "function", inherits = TRUE)) stop("step6 is missing and match_spectral_library is unavailable")
          matched <- match_spectral_library(state$object, config$spectral)
          if (methods::is(matched, "mass_dataset")) {
            state$object <- matched
            result <- matched@process_info$match_spectral_library$result
            if (!identical(result$status, "OK") || is.null(result$step6)) stop("Spectral matching did not produce step6")
            state$step6 <- result$step6
          } else {
            state$step6 <- matched
          }
        }
        state
      },
      calculate_mfc = {
        state$object <- calculate_mfc(state$object, state$step6, config$mfc_output_dir)
        mfc <- state$object@process_info$calculate_mfc$result
        if (!identical(mfc$status, "OK")) stop("MFC did not complete: ", mfc$status)
        state$candidates <- mfc$candidates
        state$step6 <- NULL
        state
      },
      embed_features = {
        de <- state$object@process_info$detect_de_features$result
        if (!is.data.frame(de) || !all(c("feature_id", "is_de") %in% names(de))) stop("Unresolved differential-feature results")
        state$unresolved_de <- de[is.na(de$is_de), , drop = FALSE]
        state$embedding <- embed_features(state$candidates, config$vector_library, de$feature_id[which(de$is_de %in% TRUE)]); state
      },
      select_k = { state$selection <- select_k(state$embedding); state },
      cluster_features = { state$clusters <- cluster_features(state$embedding, state$selection); state },
      select_representatives = { state$representatives <- select_representatives(state$embedding, state$clusters$members); state },
      retrieve_neighbors = { state$neighbors <- retrieve_neighbors(state$embedding, state$representatives, config$vector_library); state },
      build_llm_evidence = {
        ev <- config$evidence
        if (is.null(ev) || !all(c("hmdb_pathways", "reactome_mapping", "output_dir") %in% names(ev)))
          stop("evidence config requires hmdb_pathways, reactome_mapping and output_dir")
        de <- state$object@process_info$detect_de_features$result
        variable <- state$object@variable_info
        features <- data.frame(feature_id=as.character(variable$variable_id), mz=variable$mz,
                               log2FC=de$log2FC[match(variable$variable_id,de$feature_id)])
        state$evidence <- build_llm_evidence(candidates=state$candidates, representatives=state$representatives,
          hmdb_pathways=ev$hmdb_pathways, reactome_mapping=ev$reactome_mapping,
          members=state$clusters$members, neighbors=state$neighbors, features=features,
          metadata=file.path(config$vector_library,"ids.tsv"), output_dir=ev$output_dir,
          run_id=if (is.null(ev$run_id)) "d15_resume" else ev$run_id)
        state
      },
      annotate_clusters = stop("Annotation requires a newly constructed evidence view")
    )
    state$process_info[[step]] <- list(function_name = step, input_hash = previous_key,
      cache_key = key, started = started, finished = as.character(Sys.time()), status = "OK")
    state$cache_keys[[step]] <- key
    tmp <- tempfile(pattern = paste0(".", step, "-"), tmpdir = cache_dir)
    saveRDS(list(key = key, state = state), tmp)
    if (!file.rename(tmp, path)) { unlink(tmp); stop("Could not atomically save checkpoint") }
    previous_key <- key
  }
  state
}
