# Internal helpers for entropy retrieval and reverse spectral matching.
# Shared functions already migrated for MFC are provided by zz_legacy_mfc.R.
`.fe_add_missing_columns` <- 
function(df, cols) {
  `50a_add_missing_columns`(df, cols)
}

`.fe_as_data_frame` <- 
function(x) {
  as.data.frame(x, stringsAsFactors = FALSE, check.names = FALSE)
}

`.fe_collect_tables` <- 
function(x) {
  out <- list()
  walk <- function(y) {
    if (is.data.frame(y)) {
      out[[length(out) + 1L]] <<- y
    } else if (is.list(y)) {
      for (item in y) walk(item)
    }
  }
  walk(x)
  out
}

`.fe_compound_key` <- 
function(df) {
  first_nonempty <- function(...) {
    values <- list(...)
    out <- rep(NA_character_, nrow(df))
    for (x in values) {
      x <- as.character(x)
      use <- is.na(out) | !nzchar(out)
      good <- !is.na(x) & nzchar(x)
      out[use & good] <- x[use & good]
    }
    out[is.na(out) | !nzchar(out)] <- paste0("spectrum:", as.character(df$spectrum_id[is.na(out) | !nzchar(out)]))
    out
  }
  name_formula <- paste(as.character(df$candidate_name), as.character(df$formula), sep = "|")
  name_formula[is.na(df$candidate_name) | !nzchar(as.character(df$candidate_name))] <- NA_character_
  first_nonempty(
    ifelse(!is.na(df$inchikey) & nzchar(as.character(df$inchikey)), paste0("inchikey:", df$inchikey), NA_character_),
    ifelse(!is.na(df$candidate_id) & nzchar(as.character(df$candidate_id)), paste0("candidate:", df$candidate_id), NA_character_),
    ifelse(!is.na(df$compound_id) & nzchar(as.character(df$compound_id)), paste0("compound:", df$compound_id), NA_character_),
    ifelse(!is.na(name_formula) & nzchar(name_formula), paste0("name_formula:", name_formula), NA_character_),
    ifelse(!is.na(df$spectrum_id) & nzchar(as.character(df$spectrum_id)), paste0("spectrum:", df$spectrum_id), NA_character_)
  )
}

`.fe_emit_progress` <- 
function(progress_callback = NULL,
                              progress_path = NULL,
                              started_at = Sys.time(),
                              processed_queries = 0L,
                              total_queries = 0L,
                              candidates_so_far = 0L,
                              message = "") {
  elapsed <- as.numeric(difftime(Sys.time(), started_at, units = "secs"))
  row <- data.frame(
    timestamp = format(Sys.time(), "%Y-%m-%d %H:%M:%S"),
    processed_queries = as.integer(processed_queries),
    total_queries = as.integer(total_queries),
    elapsed_seconds = as.numeric(elapsed),
    candidates_so_far = as.integer(candidates_so_far),
    message = as.character(message),
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
  if (!is.null(progress_path) && length(progress_path) && nzchar(as.character(progress_path[[1]]))) {
    path <- as.character(progress_path[[1]])
    dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
    utils::write.table(
      row,
      file = path,
      sep = ",",
      row.names = FALSE,
      col.names = !file.exists(path),
      append = file.exists(path),
      qmethod = "double"
    )
  }
  if (is.function(progress_callback)) {
    try(progress_callback(row), silent = TRUE)
  }
  message(
    "Module 5 progress: processed_queries=", row$processed_queries[[1]],
    "; total_queries=", row$total_queries[[1]],
    "; elapsed_seconds=", round(row$elapsed_seconds[[1]], 2),
    "; candidates_so_far=", row$candidates_so_far[[1]],
    if (nzchar(row$message[[1]])) paste0("; message=", row$message[[1]]) else ""
  )
  try(flush.console(), silent = TRUE)
  invisible(row)
}

`.fe_empty_candidate_hits` <- 
function() {
  data.frame(
    peak_id = character(),
    anchor_feature_id = character(),
    anchor_mz = numeric(),
    compound_key = character(),
    candidate_id = character(),
    compound_id = character(),
    candidate_name = character(),
    inchikey = character(),
    formula = character(),
    smiles = character(),
    database_source = character(),
    best_spectrum_id = character(),
    best_library_source = character(),
    S_FE = numeric(),
    entropy_S = numeric(),
    S_FE_compound = numeric(),
    n_match_fe = integer(),
    matched_query_peak_n = integer(),
    matched_library_peak_n = integer(),
    query_usage_fe = numeric(),
    library_usage_fe = numeric(),
    n_spectra_for_compound = integer(),
    rank_compound_within_anchor = integer(),
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
}

`.fe_empty_entropy_candidates` <- 
function() {
  data.frame(
    query_id = character(),
    pseudo_ms2_id = character(),
    anchor_feature_id = character(),
    anchor_mz = numeric(),
    library_spectrum_id = character(),
    library_id = character(),
    library_source = character(),
    candidate_id = character(),
    spectrum_id = character(),
    metabolite_name = character(),
    candidate_name = character(),
    formula = character(),
    adduct = character(),
    adduct_status = character(),
    library_precursor_mz = numeric(),
    precursor_mz = numeric(),
    precursor_error_da = numeric(),
    precursor_error_ppm = numeric(),
    entropy_score = numeric(),
    entropy_S = numeric(),
    entropy_rank = integer(),
    n_matched_peaks = integer(),
    n_matched_peaks_entropy = integer(),
    matched_peak_summary = character(),
    matched_peak_pairs = character(),
    rank_entropy = integer(),
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
}

`.fe_empty_entropy_filter_summary` <- 
function() {
  data.frame(
    query_id = character(),
    anchor_feature_id = character(),
    n_entropy_candidates_before_top500 = integer(),
    n_entropy_candidates_after_top500 = integer(),
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
}

`.fe_empty_flash_entropy_index` <- 
function() {
  data.frame(
    spectrum_id = character(),
    candidate_id = character(),
    compound_id = character(),
    candidate_name = character(),
    library_source = character(),
    database_source = character(),
    library_fragment_id = character(),
    library_peak_id = character(),
    peak_mz = numeric(),
    library_mz = numeric(),
    library_intensity = numeric(),
    precursor_mz = numeric(),
    formula = character(),
    adduct = character(),
    polarity = character(),
    inchikey = character(),
    smiles = character(),
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
}

`.fe_empty_flash_entropy_index_result` <- 
function() {
  list(peak_index = .fe_empty_flash_entropy_index())
}

`.fe_empty_flash_hits` <- 
function() {
  data.frame(
    peak_id = character(),
    anchor_feature_id = character(),
    query_id = character(),
    anchor_mz = numeric(),
    candidate_id = character(),
    spectrum_id = character(),
    candidate_name = character(),
    library_source = character(),
    precursor_mz = numeric(),
    formula = character(),
    inchikey = character(),
    smiles = character(),
    database_source = character(),
    compound_id = character(),
    adduct = character(),
    polarity = character(),
    entropy_S = numeric(),
    S_FE = numeric(),
    n_match_fe = integer(),
    matched_query_peak_n = integer(),
    matched_library_peak_n = integer(),
    n_FE_match = integer(),
    matched_query_peaks = integer(),
    matched_library_peaks = integer(),
    query_usage_fe = numeric(),
    library_usage_fe = numeric(),
    rank_spectrum_within_anchor = integer(),
    retrieval_rank = integer(),
    rank_entropy = integer(),
    precursor_error_ppm = numeric(),
    s_fe = numeric(),
    cluster_id = character(),
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
}

`.fe_empty_matched_peaks` <- 
function() {
  data.frame(
    peak_id = character(),
    anchor_feature_id = character(),
    anchor_mz = numeric(),
    query_fragment_id = character(),
    query_feature_id = character(),
    query_mz = numeric(),
    query_intensity = numeric(),
    query_SAS = numeric(),
    spectrum_id = character(),
    candidate_id = character(),
    compound_id = character(),
    candidate_name = character(),
    library_source = character(),
    database_source = character(),
    library_fragment_id = character(),
    library_mz = numeric(),
    library_intensity = numeric(),
    ppm_error = numeric(),
    delta_h = numeric(),
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
}

`.fe_empty_ms1_metadata` <- 
function() {
  data.frame(
    candidate_id = character(),
    compound_id = character(),
    candidate_name = character(),
    formula = character(),
    exact_mass = numeric(),
    precursor_mz = numeric(),
    adduct = character(),
    polarity = character(),
    inchikey = character(),
    smiles = character(),
    database_source = character(),
    library_source = character(),
    database_version = character(),
    source_object = character(),
    Kingdom = character(),
    Super_class = character(),
    Class = character(),
    Sub_class = character(),
    CAS.ID = character(),
    PUBCHEM.ID = character(),
    HMDB.ID = character(),
    KEGG.ID = character(),
    CHEBI.ID = character(),
    METLIN.ID = character(),
    LIPIDMAPS.ID = character(),
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
}

`.fe_empty_peak_table` <- 
function(raw_only = FALSE) {
  out <- data.frame(
    spectrum_id = character(),
    candidate_id = character(),
    peak_mz = numeric(),
    intensity_raw = numeric(),
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
  if (!isTRUE(raw_only)) {
    out$intensity_basepeak100 <- numeric()
    out$intensity_sum1 <- numeric()
  }
  out
}

`.fe_empty_spectrum_table` <- 
function() {
  data.frame(
    spectrum_id = character(),
    candidate_id = character(),
    candidate_name = character(),
    library_source = character(),
    compound_id = character(),
    precursor_mz = numeric(),
    precursor_source = character(),
    exact_mass = numeric(),
    formula = character(),
    polarity = character(),
    adduct = character(),
    collision_energy = character(),
    instrument_type = character(),
    predicted = logical(),
    splash = character(),
    inchikey = character(),
    smiles = character(),
    database_source = character(),
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
}

`.fe_entropy_increment` <- 
function(q, r) {
  f <- function(x) ifelse(x > 0, x * log2(x), 0)
  as.numeric(f(q + r) - f(q) - f(r))
}

`.fe_fill_missing_scalar` <- 
function(x, value) {
  if (is.null(x)) {
    return(x)
  }
  if (length(x) == 0L) {
    return(x)
  }
  missing <- is.na(x) | !nzchar(as.character(x))
  if (any(missing)) {
    x[missing] <- value
  }
  x
}

`.fe_fill_spectrum_from_ms1` <- 
function(spectrum, ms1) {
  by_candidate <- merge(
    spectrum,
    ms1,
    by = "candidate_id",
    all.x = TRUE,
    sort = FALSE,
    suffixes = c("", ".ms1")
  )
  if (!any(grepl("\\.ms1$", names(by_candidate)))) {
    by_candidate <- merge(
      spectrum,
      ms1,
      by = "compound_id",
      all.x = TRUE,
      sort = FALSE,
      suffixes = c("", ".ms1")
    )
  }
  for (col in setdiff(colnames(.fe_empty_ms1_metadata()), c("candidate_id"))) {
    ms1_col <- paste0(col, ".ms1")
    if (col %in% names(by_candidate) && ms1_col %in% names(by_candidate)) {
      missing <- is.na(by_candidate[[col]]) | !nzchar(as.character(by_candidate[[col]]))
      by_candidate[[col]][missing] <- by_candidate[[ms1_col]][missing]
    }
  }
  by_candidate[, colnames(.fe_empty_spectrum_table()), drop = FALSE]
}

`.fe_first_nonempty_scalar` <- 
function(...) {
  values <- list(...)
  for (value in values) {
    if (is.null(value) || length(value) == 0L) {
      next
    }
    value <- as.character(value[[1]])
    if (!is.na(value) && nzchar(value)) {
      return(value)
    }
  }
  NA_character_
}

`.fe_lower_bound` <- 
function(x, value) {
  n <- length(x)
  lo <- 1L
  hi <- n + 1L
  while (lo < hi) {
    mid <- floor((lo + hi) / 2)
    if (x[[mid]] < value) {
      lo <- mid + 1L
    } else {
      hi <- mid
    }
  }
  lo
}

`.fe_ms2_from_long_table` <- 
function(df, object_name = "object") {
  df <- .fe_as_data_frame(df)
  df <- .fe_standardize_columns(df, list(
    spectrum_id = c("spectrum_id", "scan", "scan_id", "title", "id", "ID"),
    candidate_id = c("candidate_id", "compound_id", "hmdb_id", "HMDB.ID"),
    compound_id = c("compound_id", "candidate_id", "hmdb_id", "HMDB.ID"),
    candidate_name = c("candidate_name", "name", "compound_name", "metabolite_name", "Name"),
    precursor_mz = c("precursor_mz", "precursorMz", "precursor", "parent_mz"),
    peak_mz = c("peak_mz", "mz", "fragment_mz", "product_mz"),
    intensity_raw = c("intensity_raw", "intensity", "into", "peak_intensity", "abundance"),
    library_source = c("library_source", "source", "database_source"),
    database_source = c("database_source", "source", "library_source"),
    polarity = c("polarity", "ion_mode", "mode"),
    adduct = c("adduct", "adduct_type", "ion_form"),
    collision_energy = c("collision_energy", "ce", "collisionEnergy"),
    instrument_type = c("instrument_type", "instrument"),
    exact_mass = c("exact_mass", "monoisotopic_molecular_weight", "mass", "exactmass"),
    formula = c("formula", "molecular_formula", "Formula"),
    inchikey = c("inchikey", "InChIKey"),
    smiles = c("smiles", "SMILES")
  ))
  if (!"spectrum_id" %in% names(df) || all(is.na(df$spectrum_id) | !nzchar(as.character(df$spectrum_id)))) {
    group_key <- if ("candidate_id" %in% names(df)) as.character(df$candidate_id) else rep(object_name, nrow(df))
    group_key[is.na(group_key) | !nzchar(group_key)] <- object_name
    df$spectrum_id <- paste0("SPECTRUM_", as.integer(factor(group_key, levels = unique(group_key))))
  }
  if (!"candidate_id" %in% names(df) || all(is.na(df$candidate_id) | !nzchar(as.character(df$candidate_id)))) {
    df$candidate_id <- if ("compound_id" %in% names(df)) as.character(df$compound_id) else as.character(df$spectrum_id)
  }

  spectrum_cols <- colnames(.fe_empty_spectrum_table())
  spectrum <- .fe_add_missing_columns(df, spectrum_cols)
  spectrum <- spectrum[!duplicated(spectrum$spectrum_id), spectrum_cols, drop = FALSE]
  peaks <- .fe_add_missing_columns(df, colnames(.fe_empty_peak_table(raw_only = TRUE)))
  peaks <- peaks[, colnames(.fe_empty_peak_table(raw_only = TRUE)), drop = FALSE]
  list(library_spectrum_table = spectrum, library_peak_table = peaks)
}

`.fe_ms2_from_spectrum_list` <- 
function(x, object_name) {
  spectrum_id <- .fe_scalar(x$spectrum_id, NULL)
  if (is.null(spectrum_id)) spectrum_id <- .fe_scalar(x$id, NULL)
  if (is.null(spectrum_id)) spectrum_id <- .fe_scalar(x$title, NULL)
  if (is.null(spectrum_id)) spectrum_id <- paste0("SPECTRUM_", gsub("[^A-Za-z0-9]+", "_", object_name))

  peaks <- NULL
  if (!is.null(x$peaks)) {
    peaks <- .fe_as_data_frame(x$peaks)
    if (ncol(peaks) >= 2L && !all(c("peak_mz", "mz") %in% names(peaks))) {
      names(peaks)[1:2] <- c("peak_mz", "intensity_raw")
    }
  } else {
    peaks <- data.frame(peak_mz = x$mz, intensity_raw = x$intensity, stringsAsFactors = FALSE)
  }
  peaks$spectrum_id <- spectrum_id
  peaks$candidate_id <- .fe_scalar(x$candidate_id, .fe_scalar(x$compound_id, spectrum_id))

  spectrum <- data.frame(
    spectrum_id = spectrum_id,
    candidate_id = peaks$candidate_id[[1]],
    candidate_name = .fe_scalar(x$candidate_name, .fe_scalar(x$name, NA_character_)),
    library_source = .fe_scalar(x$library_source, .fe_scalar(x$source, NA_character_)),
    compound_id = .fe_scalar(x$compound_id, peaks$candidate_id[[1]]),
    precursor_mz = .fe_scalar(x$precursor_mz, NA_real_),
    exact_mass = .fe_scalar(x$exact_mass, NA_real_),
    formula = .fe_scalar(x$formula, NA_character_),
    polarity = .fe_scalar(x$polarity, NA_character_),
    adduct = .fe_scalar(x$adduct, NA_character_),
    collision_energy = .fe_scalar(x$collision_energy, NA_character_),
    instrument_type = .fe_scalar(x$instrument_type, NA_character_),
    predicted = .fe_scalar(x$predicted, NA),
    splash = .fe_scalar(x$splash, NA_character_),
    inchikey = .fe_scalar(x$inchikey, NA_character_),
    smiles = .fe_scalar(x$smiles, NA_character_),
    database_source = .fe_scalar(x$database_source, .fe_scalar(x$source, NA_character_)),
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
  list(library_spectrum_table = spectrum, library_peak_table = peaks)
}

`.fe_peak_index` <- 
function(index) {
  if (is.list(index) && !is.data.frame(index) && is.data.frame(index$peak_index)) {
    return(index$peak_index)
  }
  as.data.frame(index, stringsAsFactors = FALSE, check.names = FALSE)
}

`.fe_rbind_fill` <- 
function(x) {
  x <- x[vapply(x, function(item) is.data.frame(item), logical(1))]
  if (length(x) == 0L) {
    return(data.frame())
  }
  cols <- unique(unlist(lapply(x, names), use.names = FALSE))
  do.call(rbind, lapply(x, function(df) {
    df <- .fe_add_missing_columns(df, cols)
    df[, cols, drop = FALSE]
  }))
}

`.fe_resolve_database_dir` <- 
function(database_dir) {
  if (is.null(database_dir) || dir.exists(database_dir)) {
    return(database_dir)
  }
  fallback <- file.path(dirname(database_dir), "database", basename(database_dir))
  if (dir.exists(fallback)) {
    return(fallback)
  }
  database_dir
}

`.fe_resolve_database_paths` <- 
function(database_root = NULL,
                                       ms1_database_dir = NULL,
                                       ms2_database_dir = NULL) {
  option_root <- getOption("msidataset.database_root", default = NULL)
  env_root <- Sys.getenv("MSIDATASET_DATABASE_ROOT", unset = NA_character_)
  package_root <- system.file("database", package = "msidataset")
  if (!nzchar(package_root)) {
    package_root <- NA_character_
  }

  root <- .fe_first_nonempty_scalar(database_root, option_root, env_root, package_root)
  has_ms1_dir <- !is.na(.fe_first_nonempty_scalar(ms1_database_dir))
  has_ms2_dir <- !is.na(.fe_first_nonempty_scalar(ms2_database_dir))
  if (is.na(root) && !(has_ms1_dir && has_ms2_dir)) {
    stop(
      paste(
        "No spectral database root is configured.",
        "Pass database_root, set options(msidataset.database_root = ...),",
        "or set MSIDATASET_DATABASE_ROOT."
      ),
      call. = FALSE
    )
  }
  root <- if (is.na(root)) NA_character_ else path.expand(as.character(root[[1]]))

  ms1_dir <- if (has_ms1_dir) {
    .fe_first_nonempty_scalar(ms1_database_dir)
  } else {
    c(file.path(root, "ms1_database"), file.path(root, "database", "ms1_database"))
  }
  ms2_dir <- if (has_ms2_dir) {
    .fe_first_nonempty_scalar(ms2_database_dir)
  } else {
    c(file.path(root, "ms2_database"), file.path(root, "database", "ms2_database"))
  }

  ms1_dir <- .fe_resolve_existing_database_dir(ms1_dir, "ms1_database")
  ms2_dir <- .fe_resolve_existing_database_dir(ms2_dir, "ms2_database")
  list(
    database_root = root,
    ms1_database_dir = ms1_dir,
    ms2_database_dir = ms2_dir
  )
}

`.fe_resolve_existing_database_dir` <- 
function(database_dir, label) {
  candidates <- unique(path.expand(as.character(database_dir)))
  candidates <- candidates[!is.na(candidates) & nzchar(candidates)]
  for (candidate in candidates) {
    resolved <- .fe_resolve_database_dir(candidate)
    if (!is.null(resolved) && dir.exists(resolved)) {
      return(resolved)
    }
  }
  stop(
    "Database directory does not exist for ", label, ": ",
    paste(candidates, collapse = " or "),
    call. = FALSE
  )
}

`.fe_scalar` <- 
function(x, default = NA) {
  if (is.null(x) || length(x) == 0L) return(default)
  x[[1]]
}

`.fe_split_key` <- 
function(df, cols) {
  if (nrow(df) == 0L) {
    return(character())
  }
  missing <- setdiff(cols, names(df))
  if (length(missing) > 0L) {
    stop("Cannot build split key; missing columns: ", paste(missing, collapse = ", "), call. = FALSE)
  }
  pieces <- lapply(cols, function(col) {
    value <- as.character(df[[col]])
    value[is.na(value)] <- "<NA>"
    paste0(nchar(value, type = "bytes"), ":", value)
  })
  do.call(paste0, pieces)
}

`.fe_standardize_columns` <- 
function(df, mapping) {
  df <- .fe_as_data_frame(df)
  for (new in names(mapping)) {
    if (new %in% names(df)) next
    old <- mapping[[new]][mapping[[new]] %in% names(df)]
    if (length(old) > 0L) {
      value <- df[[old[[1]]]]
      if (length(old) > 1L) {
        for (candidate in old[-1]) {
          missing <- is.na(value) | !nzchar(as.character(value))
          if (any(missing)) {
            value[missing] <- df[[candidate]][missing]
          }
        }
      }
      df[[new]] <- value
    }
  }
  df
}

`.fe_standardize_library_once` <- 
function(library) {
  if (isTRUE(attr(library, "msidataset_standardized_library", exact = TRUE))) {
    return(library)
  }
  `54_standardize_library_spectra`(library = library)
}

`.fe_standardize_peak_table` <- 
function(df) {
  df <- .fe_as_data_frame(df)
  df <- .fe_standardize_columns(df, list(
    spectrum_id = c("spectrum_id", "scan", "scan_id", "title", "id", "ID"),
    candidate_id = c("candidate_id", "compound_id", "hmdb_id", "HMDB.ID"),
    peak_mz = c("peak_mz", "mz", "MZ", "PeakMZ", "fragment_mz", "product_mz"),
    intensity_raw = c("intensity_raw", "intensity", "Intensity", "IntensityRaw", "into", "peak_intensity", "abundance"),
    intensity_basepeak100 = c("intensity_basepeak100"),
    intensity_sum1 = c("intensity_sum1")
  ))
  if (!"spectrum_id" %in% names(df)) {
    df$spectrum_id <- paste0("SPECTRUM_", seq_len(nrow(df)))
  }
  if (!"candidate_id" %in% names(df)) {
    df$candidate_id <- NA_character_
  }
  .fe_add_missing_columns(df, colnames(.fe_empty_peak_table()))
}

`.fe_standardize_spectrum_table` <- 
function(df) {
  df <- .fe_as_data_frame(df)
  df <- .fe_standardize_columns(df, list(
    spectrum_id = c("spectrum_id", "Lab.ID", "GNPS.ID", "HMDB.ID", "METLIN.ID", "scan", "scan_id", "title", "SpectrumID", "id", "ID"),
    candidate_id = c("candidate_id", "compound_id", "HMDB.ID", "GNPS.ID", "KEGG.ID", "PUBCHEM.ID", "CHEBI.ID", "METLIN.ID", "Lab.ID", "hmdb_id"),
    candidate_name = c("candidate_name", "Compound.name", "name", "compound_name", "metabolite_name", "Name", "IUPAC_name", "Traditional_IUPAC_name"),
    library_source = c("library_source", "source", "database_source"),
    compound_id = c("compound_id", "candidate_id", "HMDB.ID", "GNPS.ID", "KEGG.ID", "PUBCHEM.ID", "CHEBI.ID", "METLIN.ID", "Lab.ID", "hmdb_id"),
    precursor_mz = c("Precursor_mz", "precursor_mz", "precursorMz", "precursor", "parent_mz", "mz", "mz.pos", "mz.neg"),
    exact_mass = c("exact_mass", "monoisotopic_molecular_weight", "mass", "exactmass"),
    formula = c("formula", "molecular_formula", "Formula"),
    polarity = c("polarity", "Polarity", "ion_mode", "mode"),
    adduct = c("adduct", "Adduct", "adduct_type", "ion_form"),
    collision_energy = c("collision_energy", "CE", "ce", "collisionEnergy", "CollisionEnergy"),
    instrument_type = c("instrument_type", "Instrument_type", "Instrument", "instrument"),
    predicted = c("predicted"),
    splash = c("splash"),
    inchikey = c("INCHIKEY.ID", "inchikey", "InChIKey"),
    smiles = c("SMILES.ID", "smiles", "SMILES"),
    database_source = c("database_source", "source", "library_source")
  ))
  df <- .fe_add_missing_columns(df, colnames(.fe_empty_spectrum_table()))
  if (nrow(df) > 0L) {
    missing_sid <- is.na(df$spectrum_id) | !nzchar(as.character(df$spectrum_id))
    df$spectrum_id[missing_sid] <- paste0("SPECTRUM_", which(missing_sid))
    missing_cid <- is.na(df$candidate_id) | !nzchar(as.character(df$candidate_id))
    df$candidate_id[missing_cid] <- as.character(df$compound_id[missing_cid])
    missing_cid <- is.na(df$candidate_id) | !nzchar(as.character(df$candidate_id))
    df$candidate_id[missing_cid] <- as.character(df$spectrum_id[missing_cid])
  }
  df$precursor_mz <- suppressWarnings(as.numeric(df$precursor_mz))
  df$exact_mass <- suppressWarnings(as.numeric(df$exact_mass))
  df[!duplicated(df$spectrum_id), colnames(.fe_empty_spectrum_table()), drop = FALSE]
}

`.fe_upper_bound` <- 
function(x, value) {
  n <- length(x)
  lo <- 1L
  hi <- n + 1L
  while (lo < hi) {
    mid <- floor((lo + hi) / 2)
    if (x[[mid]] <= value) {
      lo <- mid + 1L
    } else {
      hi <- mid
    }
  }
  lo - 1L
}

`.mlc_acquire_lock` <- 
function(cache_path, key_sha256) {
  lock_path <- .mlc_cache_lock_path(cache_path)
  if (dir.exists(lock_path)) .mlc_release_stale_lock(lock_path)
  if (!dir.create(lock_path, recursive = FALSE, showWarnings = FALSE)) {
    stop("A matching-library cache build is already in progress: ", lock_path, call. = FALSE)
  }
  owner <- list(
    pid = Sys.getpid(),
    process_fingerprint = .mlc_pid_fingerprint(Sys.getpid()),
    created_at_utc = format(Sys.time(), tz = "UTC", usetz = TRUE),
    cache_key_sha256 = as.character(key_sha256[[1L]])
  )
  owner_path <- file.path(lock_path, "owner.yml")
  tryCatch({
    yaml::write_yaml(owner, owner_path)
  }, error = function(error) {
    unlink(lock_path)
    stop("Unable to record matching-library cache lock owner: ", conditionMessage(error), call. = FALSE)
  })
  lock_path
}

`.mlc_apply_standard_metadata_overlay` <- 
function(library_metadata,
                                                 fragment_search_index,
                                                 overlay) {
  metadata <- as.data.frame(library_metadata, stringsAsFactors = FALSE, check.names = FALSE)
  fragment_index <- as.data.frame(fragment_search_index, stringsAsFactors = FALSE, check.names = FALSE)
  required_metadata <- c("spectrum_id", "candidate_id", "candidate_name", "compound_id", "precursor_mz", "polarity")
  required_index <- c("spectrum_id", "candidate_id", "candidate_name", "compound_id")
  if (length(setdiff(required_metadata, names(metadata))) ||
      length(setdiff(required_index, names(fragment_index)))) {
    stop("Metadata overlay cannot be applied to an incomplete cache metadata/index schema.", call. = FALSE)
  }
  expected_ids <- as.character(overlay$table$spectrum_id)
  expected_labs <- as.character(overlay$table$Lab.ID)
  row_index <- match(expected_ids, as.character(metadata$spectrum_id))
  if (anyNA(row_index) || anyDuplicated(as.character(metadata$spectrum_id[row_index]))) {
    stop("Metadata overlay target spectrum_id is missing or ambiguous in cache metadata.", call. = FALSE)
  }
  if ("Lab.ID" %in% names(metadata) &&
      any(as.character(metadata$Lab.ID[row_index]) != expected_labs)) {
    stop("Metadata overlay Lab.ID does not match cache metadata.", call. = FALSE)
  }
  old_id <- as.character(metadata$candidate_id[row_index])
  old_name <- as.character(metadata$candidate_name[row_index])
  if (any(old_id != as.character(overlay$table$raw_candidate_id)) ||
      any(old_name != as.character(overlay$table$raw_candidate_name)) ||
      any(abs(suppressWarnings(as.numeric(metadata$precursor_mz[row_index])) - overlay$table$precursor_mz) > 1e-12) ||
      any(tolower(as.character(metadata$polarity[row_index])) != tolower(as.character(overlay$table$polarity)))) {
    stop("Metadata overlay raw identity, precursor, or polarity does not match cache metadata.", call. = FALSE)
  }
  index_target <- as.character(fragment_index$spectrum_id) %in% expected_ids
  if (!any(index_target) ||
      any(as.character(fragment_index$candidate_id[index_target]) !=
          as.character(overlay$table$raw_candidate_id[match(as.character(fragment_index$spectrum_id[index_target]), expected_ids)])) ||
      any(as.character(fragment_index$candidate_name[index_target]) !=
          as.character(overlay$table$raw_candidate_name[match(as.character(fragment_index$spectrum_id[index_target]), expected_ids)])) ) {
    stop("Metadata overlay raw identity does not match the fragment search index.", call. = FALSE)
  }
  metadata_fields <- c("HMDB_ID", "hmdb_id", "compound_name", "metabolite_name", "source_filename", "standard_metadata_overlay_provenance")
  for (field in metadata_fields) if (!field %in% names(metadata)) metadata[[field]] <- NA_character_
  metadata$candidate_id[row_index] <- overlay$table$canonical_HMDB_ID
  metadata$compound_id[row_index] <- overlay$table$canonical_HMDB_ID
  metadata$candidate_name[row_index] <- overlay$table$canonical_compound_name
  metadata$HMDB_ID[row_index] <- overlay$table$canonical_HMDB_ID
  metadata$hmdb_id[row_index] <- overlay$table$canonical_HMDB_ID
  metadata$compound_name[row_index] <- overlay$table$canonical_compound_name
  metadata$metabolite_name[row_index] <- overlay$table$canonical_compound_name
  metadata$source_filename[row_index] <- overlay$table$source_filename
  metadata$standard_metadata_overlay_provenance[row_index] <- overlay$overlay_sha256
  index_overlay_row <- match(as.character(fragment_index$spectrum_id[index_target]), expected_ids)
  fragment_index$candidate_id[index_target] <- overlay$table$canonical_HMDB_ID[index_overlay_row]
  fragment_index$compound_id[index_target] <- overlay$table$canonical_HMDB_ID[index_overlay_row]
  fragment_index$candidate_name[index_target] <- overlay$table$canonical_compound_name[index_overlay_row]
  list(
    library_metadata = metadata,
    fragment_search_index = fragment_index,
    changed_spectrum_ids = expected_ids,
    changed_field_n = as.integer(length(row_index) * 9L),
    unchanged_spectrum_n = as.integer(nrow(metadata) - length(row_index))
  )
}

`.mlc_build_matching_library_cache` <- 
function(source_database_path,
                                              source_object_name,
                                              cache_path,
                                              effective_cleaning_config,
                                              standard_metadata_overlay = NULL) {
  configured_overlay <- .mlc_standard_metadata_overlay_path(effective_cleaning_config)
  overlay_path <- if (!is.null(standard_metadata_overlay) && length(standard_metadata_overlay)) {
    standard_metadata_overlay
  } else {
    configured_overlay
  }
  if (!is.null(overlay_path) && length(overlay_path)) {
    overlay_key <- .mlc_standard_metadata_overlay_key(overlay_path)
    if (is.null(effective_cleaning_config$standard_metadata_overlay)) {
      effective_cleaning_config$standard_metadata_overlay <- overlay_key
    } else if (!identical(effective_cleaning_config$standard_metadata_overlay, overlay_key)) {
      stop("Explicit standard metadata overlay does not match effective cache configuration.", call. = FALSE)
    }
  }
  source_database_path <- normalizePath(as.character(source_database_path[[1L]]), winslash = "/", mustWork = TRUE)
  source_database_sha256 <- .mlc_sha256_file(source_database_path)
  initial_components <- .mlc_key_components(
    source_database_sha256 = source_database_sha256,
    source_object_name = source_object_name,
    effective_cleaning_config = effective_cleaning_config,
    database_schema_digest = "pending-library-read"
  )
  lock_path <- .mlc_acquire_lock(cache_path, .mlc_cache_key(initial_components)$sha256)
  on.exit(.mlc_release_lock(lock_path), add = TRUE)
  if (file.exists(cache_path)) {
    return(.mlc_load_matching_library_cache(
      source_database_path, source_object_name, cache_path, effective_cleaning_config
    ))
  }
  library <- `50_read_spectral_library`(
    spectral_library = source_database_path,
    polarity = effective_cleaning_config$cleaning$polarity,
    min_library_peaks = effective_cleaning_config$cleaning$min_library_peaks
  )
  library <- .fe_standardize_library_once(library)
  database_schema_digest <- .mlc_library_schema_digest(library)
  cleaned_library <- `55_clean_library_peaks`(
    library = library,
    polarity = effective_cleaning_config$cleaning$polarity,
    min_library_peaks = effective_cleaning_config$cleaning$min_library_peaks
  )
  fragment_search_index <- .fe_peak_index(`56_build_flash_entropy_index`(
    library = cleaned_library,
    min_library_peaks = effective_cleaning_config$cleaning$min_library_peaks
  ))
  overlay <- NULL
  overlay_applied <- NULL
  if (!is.null(overlay_path) && length(overlay_path)) {
    standard_ids <- sort(unique(as.character(cleaned_library$library_spectrum_table$spectrum_id)[
      grepl("^STD[0-9]{6}$", as.character(cleaned_library$library_spectrum_table$spectrum_id))
    ]))
    if (length(standard_ids) != 15L) {
      stop("Cold cache build cannot identify exactly 15 standard spectrum_id values for the metadata overlay.", call. = FALSE)
    }
    overlay <- .mlc_read_standard_metadata_overlay(
      overlay_path,
      expected_spectrum_ids = standard_ids,
      expected_lab_ids = standard_ids
    )
    overlay_applied <- .mlc_apply_standard_metadata_overlay(
      cleaned_library$library_spectrum_table, fragment_search_index, overlay
    )
    cleaned_library$library_spectrum_table <- overlay_applied$library_metadata
    fragment_search_index <- overlay_applied$fragment_search_index
  }
  spectrum_id_peak_index <- .rm_build_spectrum_id_peak_index(cleaned_library$library_peak_table)
  components <- .mlc_key_components(
    source_database_sha256 = source_database_sha256,
    source_object_name = source_object_name,
    effective_cleaning_config = effective_cleaning_config,
    database_schema_digest = database_schema_digest
  )
  cache_key <- .mlc_cache_key(components)
  identifier_mappings <- .mlc_library_identifier_mappings(cleaned_library$library_spectrum_table)
  cache <- list(
    cache_schema_version = .mlc_cache_schema_version(),
    cache_key = cache_key,
    source_database_path = source_database_path,
    source_database_sha256 = source_database_sha256,
    source_object_name = as.character(source_object_name),
    effective_cleaning_config = effective_cleaning_config,
    cleaning_code_sha256 = components$cleaning_code_sha256,
    library_metadata = cleaned_library$library_spectrum_table,
    cleaned_peak_table = cleaned_library$library_peak_table,
    spectrum_id_mapping = identifier_mappings$spectrum_id_mapping,
    lab_id_mapping = identifier_mappings$lab_id_mapping,
    fragment_search_index = fragment_search_index,
    spectrum_id_peak_index = spectrum_id_peak_index,
    library_cleaning_fate = .mlc_library_cleaning_fate(
      library, cleaned_library, effective_cleaning_config
    ),
    provenance = list(
      created_at_utc = format(Sys.time(), tz = "UTC", usetz = TRUE),
      builder_pid = Sys.getpid(),
      source_database_loader = "50_read_spectral_library",
      source_database_path_at_build = source_database_path,
      source_database_sha256 = source_database_sha256,
      source_object_name = as.character(source_object_name)
    )
  )
  if (!is.null(overlay)) {
    cache$standard_metadata_overlay <- effective_cleaning_config$standard_metadata_overlay
    cache$parent_cache_lineage <- list(
      parent_cache_path = NA_character_,
      parent_cache_sha256 = NA_character_,
      derivation_mode = "cold_build_with_explicit_standard_metadata_overlay",
      changed_spectrum_n = as.integer(length(overlay_applied$changed_spectrum_ids)),
      changed_field_n = as.integer(overlay_applied$changed_field_n),
      unchanged_spectrum_n = as.integer(overlay_applied$unchanged_spectrum_n)
    )
    cache$provenance$standard_metadata_overlay <- list(
      overlay_path = overlay$path,
      overlay_sha256 = overlay$overlay_sha256,
      overlay_schema_version = overlay$schema_version,
      authoritative_mapping_source_sha256 = overlay$authoritative_mapping_source_sha256,
      changed_spectrum_n = as.integer(length(overlay_applied$changed_spectrum_ids)),
      changed_field_n = as.integer(overlay_applied$changed_field_n)
    )
  }
  cache$spectrum_id_peak_index_identity <- .mlc_index_identity(
    cache$spectrum_id_peak_index, cache$cleaned_peak_table, cache$cache_key, source_database_sha256
  )
  cache$object_digests <- .mlc_object_digests(cache)
  .mlc_validate_cache_object(
    cache, source_database_path, source_database_sha256, source_object_name,
    effective_cleaning_config, require_file_manifest = FALSE
  )
  .mlc_save_cache_atomically(cache, cache_path)
  .mlc_load_matching_library_cache(
    source_database_path, source_object_name, cache_path, effective_cleaning_config
  )
}

`.mlc_cache_helper_code_digest` <- 
function(overlay_enabled = FALSE) {
  if (!isTRUE(overlay_enabled)) return(.mlc_legacy_cache_helper_code_sha256())
  # In a development checkout, hash the helper source itself.  Hashing closures
  # loaded by pkgload is not a stable cache identity across fresh R sessions.
  source_candidates <- unique(c(
    getOption("msidataset.matching_library_cache_source_file", NA_character_),
    file.path(getwd(), "R", "matching_library_cache.R")
  ))
  source_candidates <- source_candidates[!is.na(source_candidates) & nzchar(source_candidates)]
  source_candidates <- source_candidates[file.exists(source_candidates)]
  if (length(source_candidates)) {
    return(.mlc_sha256_file(source_candidates[[1L]]))
  }
  names <- c(
    ".mlc_cache_schema_version", ".mlc_build_matching_library_cache",
    ".mlc_load_matching_library_cache", ".mlc_validate_cache_object",
    ".mlc_validate_prebuilt_spectrum_id_peak_index"
  )
  .mlc_function_code_digest(names)$combined
}

`.mlc_cache_key` <- 
function(components) {
  list(components = components, sha256 = .mlc_sha256_object(components))
}

`.mlc_cache_lock_path` <- 
function(cache_path) {
  paste0(cache_path, ".lock")
}

`.mlc_cache_manifest_path` <- 
function(cache_path) {
  paste0(cache_path, ".manifest.yml")
}

`.mlc_cache_schema_version` <- 
function() {
  "matching-library-cache-schema-0.2"
}

`.mlc_cleaning_code_identity` <- 
function() {
  identity <- .mlc_function_code_digest(c(
    "50_read_spectral_library",
    "54_standardize_library_spectra",
    "55_clean_library_peaks",
    "56_build_flash_entropy_index",
    ".fe_standardize_library_once",
    ".fe_standardize_spectrum_table",
    ".fe_standardize_peak_table",
    ".fe_standardize_columns",
    ".fe_add_missing_columns",
    ".rm_build_spectrum_id_peak_index"
  ))
  # 50/54/55 dispatch to parser helpers in this source unit.  Include its
  # byte-level source identity so an indirect helper change cannot reuse a
  # cache built with different preparation semantics.
  source_path <- file.path(getwd(), "R", "local_database_library.R")
  if (file.exists(source_path)) {
    identity$source_file_sha256 <- .mlc_sha256_file(source_path)
    identity$combined <- .mlc_sha256_object(list(
      functions = identity$functions,
      source_file_sha256 = identity$source_file_sha256
    ))
  }
  identity
}

`.mlc_effective_cache_config` <- 
function(polarity,
                                        min_library_peaks,
                                        ppm,
                                        mz_tol_da,
                                        mz_tol_da_max_mz,
                                        standard_metadata_overlay = NULL,
                                        extra_config = NULL) {
  overlay_path <- standard_metadata_overlay
  if (is.null(overlay_path) && is.list(extra_config) &&
      "standard_metadata_overlay" %in% names(extra_config)) {
    overlay_path <- extra_config$standard_metadata_overlay
    extra_config$standard_metadata_overlay <- NULL
    if (!length(extra_config)) extra_config <- NULL
  }
  config <- list(
    cleaning = list(
      polarity = as.character(polarity[[1L]]),
      min_library_peaks = as.integer(min_library_peaks[[1L]]),
      invalid_peak_rule = "finite_peak_mz_and_finite_positive_intensity_raw",
      minimum_spectrum_rule = "at_least_min_library_peaks_and_positive_finite_total",
      normalization = list(basepeak = "intensity_raw / max(intensity_raw) * 100", sum1 = "intensity_raw / sum(intensity_raw)"),
      duplicate_peak_rule = "preserve_current_55_clean_library_peaks_input_rows",
      row_order_rule = "preserve_filtered_input_order"
    ),
    fragment_search_index = list(
      spectrum_id_field = "spectrum_id",
      ordering = "order(library_mz, spectrum_id)",
      stage5_fragment_ppm = as.numeric(ppm[[1L]]),
      low_mz_abs_tol = if (is.null(mz_tol_da)) NULL else as.numeric(mz_tol_da[[1L]]),
      low_mz_cutoff = as.numeric(mz_tol_da_max_mz[[1L]])
    ),
    stage6_spectrum_index = list(
      spectrum_id_field = "spectrum_id",
      index_semantics = "original_cleaned_peak_table_row_positions"
    )
  )
  if (!is.null(extra_config)) config$caller_supplied <- extra_config
  if (!is.null(overlay_path) && length(overlay_path) &&
      !is.na(as.character(overlay_path[[1L]])) &&
      nzchar(as.character(overlay_path[[1L]]))) {
    config$standard_metadata_overlay <- .mlc_standard_metadata_overlay_key(
      overlay_path
    )
  }
  config
}

`.mlc_fsync_file` <- 
function(path) {
  if (.Platform$OS.type == "windows") return(invisible(TRUE))
  status <- suppressWarnings(system2("sync", c("-f", shQuote(path))))
  if (!identical(status, 0L)) {
    stop("Unable to fsync matching-library cache temporary file: ", path, call. = FALSE)
  }
  invisible(TRUE)
}

`.mlc_function_code_digest` <- 
function(function_names) {
  payload <- lapply(function_names, function(name) {
    fun <- get0(name, mode = "function", inherits = TRUE)
    if (is.null(fun)) {
      stop("Required matching-library cache function is unavailable: ", name, call. = FALSE)
    }
    list(formals = formals(fun), body = body(fun))
  })
  names(payload) <- function_names
  list(
    functions = vapply(payload, .mlc_sha256_object, character(1)),
    combined = .mlc_sha256_object(payload)
  )
}

`.mlc_index_identity` <- 
function(index, peaks, cache_key, source_database_sha256) {
  peaks <- as.data.frame(peaks, stringsAsFactors = FALSE, check.names = FALSE)
  .mlc_validate_index_structure(index, peaks, full = TRUE)
  ids <- as.character(peaks$spectrum_id)
  list(
    cache_schema_version = .mlc_cache_schema_version(),
    cache_key_sha256 = as.character(cache_key$sha256),
    source_database_sha256 = as.character(source_database_sha256),
    cleaned_peak_table_sha256 = .mlc_sha256_object(peaks),
    spectrum_id_field = "spectrum_id",
    n_spectra = as.integer(length(unique(ids[!is.na(ids)]))),
    n_peaks = as.integer(nrow(peaks)),
    n_unique_spectrum_id = as.integer(length(unique(ids[!is.na(ids)]))),
    index_key_sha256 = .mlc_sha256_object(sort(names(index$rows_by_spectrum_id))),
    index_object_sha256 = .mlc_sha256_object(index)
  )
}

`.mlc_key_components` <- 
function(source_database_sha256,
                                source_object_name,
                                effective_cleaning_config,
                                database_schema_digest = NULL) {
  code <- .mlc_cleaning_code_identity()
  overlay <- effective_cleaning_config$standard_metadata_overlay
  overlay_enabled <- is.list(overlay) && isTRUE(overlay$enabled)
  components <- list(
    source_database_sha256 = tolower(as.character(source_database_sha256[[1L]])),
    source_object_name = as.character(source_object_name[[1L]]),
    database_schema_digest = if (is.null(database_schema_digest)) NULL else as.character(database_schema_digest[[1L]]),
    cleaning_code_sha256 = code$combined,
    cleaning_helper_code_sha256 = code$functions,
    cache_helper_code_sha256 = .mlc_cache_helper_code_digest(overlay_enabled),
    effective_cleaning_config_sha256 = .mlc_sha256_object(effective_cleaning_config),
    fragment_index_parameters_sha256 = .mlc_sha256_object(effective_cleaning_config$fragment_search_index),
    normalization_rule_sha256 = .mlc_sha256_object(effective_cleaning_config$cleaning$normalization),
    spectrum_id_field = as.character(effective_cleaning_config$stage6_spectrum_index$spectrum_id_field),
    r_version = R.version$version.string,
    package_versions = .mlc_package_versions(),
    cache_schema_version = .mlc_cache_schema_version()
  )
  if (overlay_enabled) {
    components$standard_metadata_overlay_sha256 <- as.character(overlay$overlay_sha256)
    components$standard_metadata_overlay_schema_version <- as.character(overlay$schema_version)
    components$standard_metadata_overlay_helper_sha256 <- as.character(overlay$helper_code_sha256)
    components$authoritative_mapping_source_sha256 <- as.character(overlay$authoritative_mapping_source_sha256)
  }
  components
}

`.mlc_legacy_cache_helper_code_sha256` <- 
function() {
  # v0.1 cache keys used the byte SHA256 of this source before metadata-overlay
  # support existed.  Retaining this value on the overlay-free branch preserves
  # the approved v0.1 cache contract; overlay-enabled caches use the live file
  # SHA below and therefore cannot be mistaken for a v0.1 cache.
  "f62ce0406af4ad27c0c25abb21233a469ed123b2b7d94f6f0eb809d2fe9c90b7"
}

`.mlc_library_cleaning_fate` <- 
function(raw_library, cleaned_library, effective_cleaning_config) {
  raw_spectrum <- as.data.frame(raw_library$library_spectrum_table, stringsAsFactors = FALSE, check.names = FALSE)
  raw_peaks <- as.data.frame(raw_library$library_peak_table, stringsAsFactors = FALSE, check.names = FALSE)
  cleaned_peaks <- as.data.frame(cleaned_library$library_peak_table, stringsAsFactors = FALSE, check.names = FALSE)
  raw_spectrum <- .fe_standardize_spectrum_table(raw_spectrum)
  raw_peaks <- .fe_standardize_peak_table(raw_peaks)
  all_ids <- unique(c(as.character(raw_spectrum$spectrum_id), as.character(raw_peaks$spectrum_id)))
  all_ids <- all_ids[!is.na(all_ids) & nzchar(all_ids)]
  polarity <- as.character(effective_cleaning_config$cleaning$polarity)
  min_library_peaks <- as.integer(effective_cleaning_config$cleaning$min_library_peaks)
  spectrum_row <- match(all_ids, as.character(raw_spectrum$spectrum_id))
  raw_polarity <- as.character(raw_spectrum$polarity[spectrum_row])
  polarity_eligible <- identical(polarity, "both") | is.na(raw_polarity) | !nzchar(raw_polarity) |
    tolower(raw_polarity) == polarity
  mz <- suppressWarnings(as.numeric(raw_peaks$peak_mz))
  intensity <- suppressWarnings(as.numeric(raw_peaks$intensity_raw))
  valid_peak <- is.finite(mz) & is.finite(intensity) & intensity > 0
  raw_ids <- as.character(raw_peaks$spectrum_id)
  count_ids <- function(values) {
    values <- as.character(values)
    keep <- !is.na(values) & nzchar(values)
    tabulate(match(values[keep], all_ids), nbins = length(all_ids))
  }
  raw_peak_n <- count_ids(raw_ids)
  valid_peak_n <- count_ids(raw_ids[valid_peak])
  cleaned_ids <- as.character(cleaned_peaks$spectrum_id)
  cleaned_peak_n <- count_ids(cleaned_ids)
  retained <- cleaned_peak_n > 0L
  removal_reason <- ifelse(
    retained, "RETAINED",
    ifelse(!polarity_eligible, "POLARITY_EXCLUDED",
      ifelse(valid_peak_n == 0L, "NO_FINITE_POSITIVE_PEAKS",
        ifelse(valid_peak_n < min_library_peaks, "BELOW_MIN_LIBRARY_PEAKS", "REMOVED_BY_55_CLEAN_LIBRARY_PEAKS")
      )
    )
  )
  data.frame(
    spectrum_id = all_ids,
    raw_spectrum_row = as.integer(spectrum_row),
    raw_peak_n = as.integer(raw_peak_n),
    finite_positive_peak_n = as.integer(valid_peak_n),
    cleaned_peak_n = as.integer(cleaned_peak_n),
    retained = as.logical(retained),
    removal_reason = as.character(removal_reason),
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
}

`.mlc_library_identifier_mappings` <- 
function(library_metadata) {
  metadata <- as.data.frame(library_metadata, stringsAsFactors = FALSE, check.names = FALSE)
  required <- c("spectrum_id", "candidate_id")
  if (!all(required %in% names(metadata))) {
    stop("Cannot build cache identifier mappings without spectrum_id and candidate_id metadata.", call. = FALSE)
  }
  spectrum_id <- as.character(metadata$spectrum_id)
  candidate_id <- as.character(metadata$candidate_id)
  lab_id <- spectrum_id
  missing_lab_id <- is.na(lab_id) | !nzchar(lab_id)
  lab_id[missing_lab_id] <- candidate_id[missing_lab_id]
  stable_spectrum_row_id <- as.integer(seq_len(nrow(metadata)))
  list(
    spectrum_id_mapping = data.frame(
      stable_spectrum_row_id = stable_spectrum_row_id,
      spectrum_id = spectrum_id,
      candidate_id = candidate_id,
      stringsAsFactors = FALSE,
      check.names = FALSE
    ),
    lab_id_mapping = data.frame(
      stable_spectrum_row_id = stable_spectrum_row_id,
      spectrum_id = spectrum_id,
      lab_id = lab_id,
      stringsAsFactors = FALSE,
      check.names = FALSE
    )
  )
}

`.mlc_library_schema_digest` <- 
function(library) {
  if (!is.list(library) || !all(c("library_spectrum_table", "library_peak_table") %in% names(library))) {
    stop("Cannot calculate schema digest for an invalid spectral library.", call. = FALSE)
  }
  schema <- lapply(library[c("library_spectrum_table", "library_peak_table")], function(x) {
    x <- as.data.frame(x, stringsAsFactors = FALSE, check.names = FALSE)
    list(names = names(x), classes = vapply(x, function(column) class(column)[[1L]], character(1)))
  })
  .mlc_sha256_object(schema)
}

`.mlc_load_matching_library_cache` <- 
function(source_database_path,
                                             source_object_name,
                                             cache_path,
                                             effective_cleaning_config) {
  cache_path <- normalizePath(as.character(cache_path[[1L]]), winslash = "/", mustWork = TRUE)
  source_database_path <- normalizePath(as.character(source_database_path[[1L]]), winslash = "/", mustWork = TRUE)
  source_database_sha256 <- .mlc_sha256_file(source_database_path)
  cache <- tryCatch(readRDS(cache_path), error = function(error) {
    stop("Matching-library cache is unreadable or truncated: ", conditionMessage(error), call. = FALSE)
  })
  .mlc_validate_cache_object(
    cache = cache,
    source_database_path = source_database_path,
    source_database_sha256 = source_database_sha256,
    source_object_name = source_object_name,
    effective_cleaning_config = effective_cleaning_config,
    cache_path = cache_path,
    require_file_manifest = TRUE
  )
}

`.mlc_object_digests` <- 
function(cache) {
  out <- list(
    library_metadata = .mlc_sha256_object(cache$library_metadata),
    cleaned_peak_table = .mlc_sha256_object(cache$cleaned_peak_table),
    spectrum_id_mapping = .mlc_sha256_object(cache$spectrum_id_mapping),
    lab_id_mapping = .mlc_sha256_object(cache$lab_id_mapping),
    fragment_search_index = .mlc_sha256_object(cache$fragment_search_index),
    spectrum_id_peak_index = .mlc_sha256_object(cache$spectrum_id_peak_index),
    spectrum_id_peak_index_identity = .mlc_sha256_object(cache$spectrum_id_peak_index_identity),
    library_cleaning_fate = .mlc_sha256_object(cache$library_cleaning_fate)
  )
  if ("standard_metadata_overlay" %in% names(cache)) {
    out$standard_metadata_overlay <- .mlc_sha256_object(cache$standard_metadata_overlay)
    out$parent_cache_lineage <- .mlc_sha256_object(cache$parent_cache_lineage)
  }
  out
}

`.mlc_package_versions` <- 
function() {
  packages <- c("digest", "yaml", "methods", "stats", "utils")
  versions <- vapply(packages, function(package) {
    if (!requireNamespace(package, quietly = TRUE)) return(NA_character_)
    as.character(utils::packageVersion(package))
  }, character(1))
  as.list(versions)
}

`.mlc_pid_fingerprint` <- 
function(pid) {
  pid <- suppressWarnings(as.integer(pid[[1L]]))
  if (!is.finite(pid) || pid < 1L || .Platform$OS.type == "windows") return(NA_character_)
  proc_root <- file.path("/proc", as.character(pid))
  stat_path <- file.path(proc_root, "stat")
  command_path <- file.path(proc_root, "cmdline")
  if (!file.exists(stat_path) || !file.exists(command_path)) return(NA_character_)
  stat <- tryCatch(readLines(stat_path, warn = FALSE), error = function(error) character())
  command <- tryCatch(readBin(command_path, what = "raw", n = file.info(command_path)$size), error = function(error) raw())
  if (!length(stat)) return(NA_character_)
  .mlc_sha256_object(list(stat = stat, command = command))
}

`.mlc_pid_is_active` <- 
function(pid) {
  pid <- suppressWarnings(as.integer(pid[[1L]]))
  if (!is.finite(pid) || pid < 1L) return(FALSE)
  if (.Platform$OS.type == "windows") return(NA)
  file.exists(file.path("/proc", as.character(pid)))
}

`.mlc_prebuilt_index_envelope` <- 
function(index, identity, cache_file_sha256 = NA_character_, cache_manifest_sha256 = NA_character_) {
  list(
    index = index,
    identity = identity,
    validated_cache_file_sha256 = as.character(cache_file_sha256),
    validated_cache_manifest_sha256 = as.character(cache_manifest_sha256)
  )
}

`.mlc_read_lock_owner` <- 
function(lock_path) {
  owner_path <- file.path(lock_path, "owner.yml")
  if (!file.exists(owner_path) || !requireNamespace("yaml", quietly = TRUE)) return(NULL)
  tryCatch(yaml::read_yaml(owner_path), error = function(error) NULL)
}

`.mlc_read_standard_metadata_overlay` <- 
function(standard_metadata_overlay,
                                                expected_spectrum_ids = NULL,
                                                expected_lab_ids = NULL) {
  if (is.null(standard_metadata_overlay) || !length(standard_metadata_overlay) ||
      is.na(as.character(standard_metadata_overlay[[1L]])) ||
      !nzchar(as.character(standard_metadata_overlay[[1L]]))) {
    stop("A standard metadata overlay path is required.", call. = FALSE)
  }
  path <- normalizePath(as.character(standard_metadata_overlay[[1L]]), winslash = "/", mustWork = TRUE)
  overlay <- tryCatch(
    utils::read.delim(path, sep = "\t", quote = "", stringsAsFactors = FALSE,
                      check.names = FALSE, fileEncoding = "UTF-8"),
    error = function(error) stop("Standard metadata overlay is unreadable: ", conditionMessage(error), call. = FALSE)
  )
  required <- .mlc_standard_metadata_overlay_required_columns()
  if (!is.data.frame(overlay) || length(setdiff(required, names(overlay)))) {
    stop("Standard metadata overlay has an incomplete schema.", call. = FALSE)
  }
  overlay <- overlay[, required, drop = FALSE]
  overlay[] <- lapply(overlay, function(column) {
    if (is.factor(column)) as.character(column) else column
  })
  overlay$spectrum_id <- trimws(as.character(overlay$spectrum_id))
  overlay$Lab.ID <- trimws(as.character(overlay$Lab.ID))
  overlay$canonical_HMDB_ID <- toupper(trimws(as.character(overlay$canonical_HMDB_ID)))
  overlay$canonical_compound_name <- trimws(as.character(overlay$canonical_compound_name))
  overlay$standard_name <- trimws(as.character(overlay$standard_name))
  overlay$raw_candidate_id <- trimws(as.character(overlay$raw_candidate_id))
  overlay$raw_candidate_name <- trimws(as.character(overlay$raw_candidate_name))
  overlay$source_filename <- trimws(as.character(overlay$source_filename))
  overlay$mapping_source_sha256 <- trimws(as.character(overlay$mapping_source_sha256))
  overlay$approval_status <- trimws(as.character(overlay$approval_status))
  overlay$precursor_mz <- suppressWarnings(as.numeric(overlay$precursor_mz))
  if (nrow(overlay) != 15L) {
    stop("Standard metadata overlay must contain exactly 15 spectra.", call. = FALSE)
  }
  if (any(is.na(overlay$spectrum_id) | !nzchar(overlay$spectrum_id)) ||
      any(is.na(overlay$Lab.ID) | !nzchar(overlay$Lab.ID)) ||
      anyDuplicated(overlay$spectrum_id) || anyDuplicated(overlay$Lab.ID) ||
      any(overlay$spectrum_id != overlay$Lab.ID) ||
      any(!grepl("^STD[0-9]{6}$", overlay$spectrum_id))) {
    stop("Standard metadata overlay requires 15 unique, exact STD spectrum_id/Lab.ID pairs.", call. = FALSE)
  }
  if (any(!grepl("^HMDB[0-9]{7}$", overlay$canonical_HMDB_ID)) ||
      any(is.na(overlay$canonical_compound_name) | !nzchar(overlay$canonical_compound_name)) ||
      any(is.na(overlay$raw_candidate_id) | !nzchar(overlay$raw_candidate_id)) ||
      any(is.na(overlay$raw_candidate_name) | !nzchar(overlay$raw_candidate_name)) ||
      any(is.na(overlay$source_filename) | !nzchar(overlay$source_filename)) ||
      any(!is.finite(overlay$precursor_mz) | overlay$precursor_mz <= 0)) {
    stop("Standard metadata overlay contains missing or invalid identity fields.", call. = FALSE)
  }
  if (any(tolower(overlay$canonical_compound_name) == tolower(overlay$source_filename)) ||
      any(grepl("MSMS[0-9]|NativeID:|File:", overlay$canonical_compound_name, ignore.case = TRUE))) {
    stop("A source filename cannot be used as a canonical compound name in the standard metadata overlay.", call. = FALSE)
  }
  expected_targets <- .mlc_standard_metadata_overlay_expected_targets()
  if (!all(overlay$standard_name %in% names(expected_targets)) ||
      any(unname(expected_targets[overlay$standard_name]) != overlay$canonical_HMDB_ID) ||
      any(overlay$canonical_compound_name != overlay$standard_name) ||
      any(as.integer(table(factor(overlay$standard_name, levels = names(expected_targets)))) != 3L)) {
    stop("Standard metadata overlay does not contain the approved 3x5 target HMDB/name mapping.", call. = FALSE)
  }
  if (any(overlay$standard_name == "ADP" & overlay$canonical_HMDB_ID != "HMDB0001341")) {
    stop("ADP in the standard metadata overlay must be HMDB0001341.", call. = FALSE)
  }
  if (any(tolower(overlay$polarity) != "negative") ||
      any(overlay$approval_status != "APPROVED")) {
    stop("Standard metadata overlay has invalid polarity or non-approved records.", call. = FALSE)
  }
  sha_tokens <- unlist(strsplit(overlay$mapping_source_sha256, "\\|", fixed = FALSE), use.names = FALSE)
  if (!length(sha_tokens) || any(!grepl("^[0-9a-f]{64}$", tolower(sha_tokens)))) {
    stop("Standard metadata overlay mapping_source_sha256 is invalid.", call. = FALSE)
  }
  if (!is.null(expected_spectrum_ids)) {
    expected_spectrum_ids <- sort(unique(as.character(expected_spectrum_ids)))
    if (!identical(sort(overlay$spectrum_id), expected_spectrum_ids)) {
      stop("Standard metadata overlay does not exactly match the permitted 15 spectrum_id values.", call. = FALSE)
    }
  }
  if (!is.null(expected_lab_ids)) {
    expected_lab_ids <- sort(unique(as.character(expected_lab_ids)))
    if (!identical(sort(overlay$Lab.ID), expected_lab_ids)) {
      stop("Standard metadata overlay does not exactly match the permitted 15 Lab.ID values.", call. = FALSE)
    }
  }
  overlay <- overlay[order(overlay$spectrum_id), , drop = FALSE]
  rownames(overlay) <- NULL
  list(
    path = path,
    overlay_sha256 = .mlc_sha256_file(path),
    schema_version = .mlc_standard_metadata_overlay_schema_version(),
    authoritative_mapping_source_sha256 = .mlc_sha256_object(sort(unique(overlay$mapping_source_sha256))),
    table = overlay
  )
}

`.mlc_release_lock` <- 
function(lock_path) {
  if (is.null(lock_path) || !dir.exists(lock_path)) return(invisible(FALSE))
  owner_path <- file.path(lock_path, "owner.yml")
  if (file.exists(owner_path)) unlink(owner_path)
  .mlc_remove_empty_lock_dir(lock_path)
  invisible(TRUE)
}

`.mlc_release_stale_lock` <- 
function(lock_path) {
  owner <- .mlc_read_lock_owner(lock_path)
  active <- if (is.null(owner) || is.null(owner$pid)) NA else .mlc_pid_is_active(owner$pid)
  if (isTRUE(active) && !is.null(owner$process_fingerprint)) {
    active <- identical(as.character(owner$process_fingerprint), as.character(.mlc_pid_fingerprint(owner$pid)))
  }
  if (isTRUE(active) || is.na(active)) {
    stop("Matching-library cache lock is active or unverifiable: ", lock_path, call. = FALSE)
  }
  owner_path <- file.path(lock_path, "owner.yml")
  if (file.exists(owner_path) && unlink(owner_path) != 0L) {
    stop("Unable to remove confirmed-stale matching-library cache lock owner file: ", owner_path, call. = FALSE)
  }
  .mlc_remove_empty_lock_dir(lock_path)
  invisible(TRUE)
}

`.mlc_remove_empty_lock_dir` <- 
function(path) {
  if (!dir.exists(path)) return(invisible(TRUE))
  if (.Platform$OS.type == "windows") {
    if (unlink(path, recursive = FALSE) != 0L) {
      stop("Unable to remove empty matching-library cache lock directory: ", path, call. = FALSE)
    }
  } else {
    status <- suppressWarnings(system2("rmdir", shQuote(path)))
    if (!identical(status, 0L)) {
      stop("Unable to remove empty matching-library cache lock directory: ", path, call. = FALSE)
    }
  }
  invisible(TRUE)
}

`.mlc_resolve_matching_library_cache` <- 
function(source_database_path,
                                                source_object_name,
                                                cache_path,
                                                effective_cleaning_config,
                                                cache_mode = c("error", "build_if_missing")) {
  cache_mode <- match.arg(cache_mode)
  if (is.null(cache_path) || !length(cache_path) || !nzchar(as.character(cache_path[[1L]]))) {
    return(NULL)
  }
  if (is.null(source_database_path) || !length(source_database_path) ||
      !is.character(source_database_path) || length(source_database_path) != 1L) {
    stop("matching_library_cache_path requires one explicit canonical spectral_library file path.", call. = FALSE)
  }
  source_object_name <- .mlc_scalar(
    source_object_name,
    tools::file_path_sans_ext(basename(as.character(source_database_path[[1L]])))
  )
  cache_path <- normalizePath(as.character(cache_path[[1L]]), winslash = "/", mustWork = FALSE)
  if (file.exists(cache_path)) {
    return(.mlc_load_matching_library_cache(
      source_database_path, source_object_name, cache_path, effective_cleaning_config
    ))
  }
  if (!identical(cache_mode, "build_if_missing")) {
    stop("Matching-library cache does not exist and cache_mode is not 'build_if_missing': ", cache_path, call. = FALSE)
  }
  .mlc_build_matching_library_cache(
    source_database_path, source_object_name, cache_path, effective_cleaning_config,
    standard_metadata_overlay = .mlc_standard_metadata_overlay_path(effective_cleaning_config)
  )
}

`.mlc_save_cache_atomically` <- 
function(cache, cache_path) {
  cache_path <- normalizePath(as.character(cache_path[[1L]]), winslash = "/", mustWork = FALSE)
  parent <- dirname(cache_path)
  if (!dir.exists(parent)) {
    stop("Matching-library cache parent directory does not exist: ", parent, call. = FALSE)
  }
  if (file.exists(cache_path)) {
    stop("Refusing to overwrite an existing matching-library cache: ", cache_path, call. = FALSE)
  }
  temporary <- tempfile(pattern = paste0(basename(cache_path), ".tmp-"), tmpdir = parent)
  on.exit(if (file.exists(temporary)) unlink(temporary), add = TRUE)
  saveRDS(cache, temporary, version = 3L, compress = FALSE)
  .mlc_fsync_file(temporary)
  if (!file.rename(temporary, cache_path)) {
    stop("Unable to atomically rename matching-library cache into place: ", cache_path, call. = FALSE)
  }
  cache_file_sha256 <- .mlc_sha256_file(cache_path)
  manifest <- list(
    cache_schema_version = cache$cache_schema_version,
    cache_key_sha256 = cache$cache_key$sha256,
    cache_file_sha256 = cache_file_sha256,
    source_database_sha256 = cache$source_database_sha256,
    source_object_name = cache$source_object_name,
    created_at_utc = cache$provenance$created_at_utc
  )
  manifest_path <- .mlc_cache_manifest_path(cache_path)
  manifest_tmp <- tempfile(pattern = paste0(basename(manifest_path), ".tmp-"), tmpdir = parent)
  on.exit(if (file.exists(manifest_tmp)) unlink(manifest_tmp), add = TRUE)
  yaml::write_yaml(manifest, manifest_tmp)
  .mlc_fsync_file(manifest_tmp)
  if (!file.rename(manifest_tmp, manifest_path)) {
    stop("Unable to atomically write matching-library cache manifest: ", manifest_path, call. = FALSE)
  }
  list(cache_file_sha256 = cache_file_sha256, manifest_path = manifest_path)
}

`.mlc_scalar` <- 
function(x, default = NULL) {
  if (is.null(x) || !length(x)) default else x[[1L]]
}

`.mlc_sha256_file` <- 
function(path) {
  path <- normalizePath(as.character(path[[1L]]), winslash = "/", mustWork = TRUE)
  if (!requireNamespace("digest", quietly = TRUE)) {
    stop("Package 'digest' is required for matching-library cache identity.", call. = FALSE)
  }
  digest::digest(file = path, algo = "sha256")
}

`.mlc_sha256_object` <- 
function(x) {
  if (!requireNamespace("digest", quietly = TRUE)) {
    stop("Package 'digest' is required for matching-library cache identity.", call. = FALSE)
  }
  digest::digest(x, algo = "sha256", serialize = TRUE)
}

`.mlc_standard_metadata_overlay_expected_targets` <- 
function() {
  c(
    "L-Aspartic acid" = "HMDB0000191",
    "L-Glutamic acid" = "HMDB0000148",
    "Citric acid" = "HMDB0000094",
    "AMP" = "HMDB0000045",
    "ADP" = "HMDB0001341"
  )
}

`.mlc_standard_metadata_overlay_key` <- 
function(standard_metadata_overlay) {
  overlay <- .mlc_read_standard_metadata_overlay(standard_metadata_overlay)
  list(
    enabled = TRUE,
    path = overlay$path,
    overlay_sha256 = overlay$overlay_sha256,
    schema_version = overlay$schema_version,
    helper_code_sha256 = .mlc_function_code_digest(c(
      ".mlc_read_standard_metadata_overlay", ".mlc_apply_standard_metadata_overlay"
    ))$combined,
    authoritative_mapping_source_sha256 = overlay$authoritative_mapping_source_sha256
  )
}

`.mlc_standard_metadata_overlay_path` <- 
function(effective_cleaning_config) {
  value <- effective_cleaning_config$standard_metadata_overlay
  if (!is.list(value) || !isTRUE(value$enabled) || is.null(value$path)) return(NULL)
  as.character(value$path[[1L]])
}

`.mlc_standard_metadata_overlay_required_columns` <- 
function() {
  c(
    "spectrum_id", "Lab.ID", "source_filename", "raw_candidate_id",
    "raw_candidate_name", "canonical_HMDB_ID", "canonical_compound_name",
    "standard_name", "replicate_id", "polarity", "precursor_mz",
    "mapping_source", "mapping_source_sha256", "mapping_reason",
    "expected_change", "approval_status"
  )
}

`.mlc_standard_metadata_overlay_schema_version` <- 
function() {
  "standard-metadata-overlay-v0.1"
}

`.mlc_validate_cache_object` <- 
function(cache,
                                       source_database_path,
                                       source_database_sha256,
                                       source_object_name,
                                       effective_cleaning_config,
                                       cache_path = NULL,
                                       require_file_manifest = FALSE) {
  required <- c(
    "cache_schema_version", "cache_key", "source_database_path", "source_database_sha256",
    "source_object_name", "effective_cleaning_config", "cleaning_code_sha256",
    "library_metadata", "cleaned_peak_table", "spectrum_id_mapping", "lab_id_mapping", "fragment_search_index",
    "spectrum_id_peak_index", "spectrum_id_peak_index_identity", "library_cleaning_fate", "provenance", "object_digests"
  )
  if (!is.list(cache) || length(setdiff(required, names(cache)))) {
    stop("Matching-library cache is missing required schema fields.", call. = FALSE)
  }
  if (!identical(as.character(cache$cache_schema_version), .mlc_cache_schema_version())) {
    stop("Matching-library cache schema version does not match this runtime.", call. = FALSE)
  }
  source_database_path <- normalizePath(as.character(source_database_path[[1L]]), winslash = "/", mustWork = TRUE)
  if (!identical(tolower(as.character(cache$source_database_sha256)), tolower(as.character(source_database_sha256)))) {
    stop("Matching-library cache source database SHA256 does not match the requested source database.", call. = FALSE)
  }
  if (!identical(as.character(cache$source_object_name), as.character(source_object_name))) {
    stop("Matching-library cache source object name does not match the requested object.", call. = FALSE)
  }
  expected_components <- .mlc_key_components(
    source_database_sha256 = source_database_sha256,
    source_object_name = source_object_name,
    effective_cleaning_config = effective_cleaning_config,
    database_schema_digest = cache$cache_key$components$database_schema_digest
  )
  expected_key <- .mlc_cache_key(expected_components)
  if (!identical(as.character(cache$cache_key$sha256), as.character(expected_key$sha256))) {
    stop("Matching-library cache key does not match source, code, configuration, R, or package identity.", call. = FALSE)
  }
  if (!identical(.mlc_sha256_object(cache$effective_cleaning_config), .mlc_sha256_object(effective_cleaning_config))) {
    stop("Matching-library cache effective cleaning configuration does not match this call.", call. = FALSE)
  }
  overlay_enabled <- is.list(effective_cleaning_config$standard_metadata_overlay) &&
    isTRUE(effective_cleaning_config$standard_metadata_overlay$enabled)
  if (overlay_enabled &&
      (!all(c("standard_metadata_overlay", "parent_cache_lineage") %in% names(cache)) ||
       !identical(cache$standard_metadata_overlay, effective_cleaning_config$standard_metadata_overlay))) {
    stop("Matching-library cache standard metadata overlay does not match this call.", call. = FALSE)
  }
  if (!overlay_enabled && "standard_metadata_overlay" %in% names(cache)) {
    stop("Matching-library cache requires an explicit standard metadata overlay configuration.", call. = FALSE)
  }
  observed_digests <- .mlc_object_digests(cache)
  if (!identical(unname(observed_digests), unname(cache$object_digests))) {
    stop("Matching-library cache object digest validation failed.", call. = FALSE)
  }
  if (!identical(as.character(cache$cleaning_code_sha256), as.character(expected_components$cleaning_code_sha256))) {
    stop("Matching-library cache cleaning code SHA256 does not match this runtime.", call. = FALSE)
  }
  .mlc_validate_library_identifier_mappings(cache)
  if (overlay_enabled) .mlc_validate_standard_metadata_overlay_applied(cache)
  peaks <- as.data.frame(cache$cleaned_peak_table, stringsAsFactors = FALSE, check.names = FALSE)
  .mlc_validate_index_structure(cache$spectrum_id_peak_index, peaks, full = TRUE)
  identity <- cache$spectrum_id_peak_index_identity
  expected_identity <- .mlc_index_identity(
    cache$spectrum_id_peak_index, peaks, cache$cache_key, source_database_sha256
  )
  if (!identical(identity, expected_identity)) {
    stop("Matching-library cache spectrum_id peak-index identity validation failed.", call. = FALSE)
  }
  cache_file_sha256 <- NA_character_
  cache_manifest_sha256 <- NA_character_
  if (isTRUE(require_file_manifest)) {
    if (is.null(cache_path) || !file.exists(cache_path)) {
      stop("Matching-library cache file is unavailable for file-level validation.", call. = FALSE)
    }
    manifest_path <- .mlc_cache_manifest_path(cache_path)
    if (!file.exists(manifest_path)) {
      stop("Matching-library cache manifest is missing: ", manifest_path, call. = FALSE)
    }
    manifest <- tryCatch(yaml::read_yaml(manifest_path), error = function(error) {
      stop("Matching-library cache manifest is unreadable: ", conditionMessage(error), call. = FALSE)
    })
    cache_file_sha256 <- .mlc_sha256_file(cache_path)
    cache_manifest_sha256 <- .mlc_sha256_file(manifest_path)
    if (!identical(tolower(as.character(manifest$cache_file_sha256)), tolower(cache_file_sha256)) ||
        !identical(as.character(manifest$cache_key_sha256), as.character(cache$cache_key$sha256)) ||
        !identical(tolower(as.character(manifest$source_database_sha256)), tolower(as.character(source_database_sha256)))) {
      stop("Matching-library cache manifest does not validate the cache file identity.", call. = FALSE)
    }
  }
  cached_library <- list(
    library_spectrum_table = cache$library_metadata,
    library_peak_table = cache$cleaned_peak_table
  )
  attr(cached_library, "msidataset_standardized_library") <- TRUE
  list(
    library = cached_library,
    fragment_search_index = cache$fragment_search_index,
    spectrum_id_peak_index = .mlc_prebuilt_index_envelope(
      cache$spectrum_id_peak_index,
      cache$spectrum_id_peak_index_identity,
      cache_file_sha256 = cache_file_sha256,
      cache_manifest_sha256 = cache_manifest_sha256
    ),
    cache_identity = list(
      cache_key_sha256 = cache$cache_key$sha256,
      source_database_path_requested = source_database_path,
      source_database_path_cached = cache$source_database_path,
      source_database_sha256 = cache$source_database_sha256,
      cache_file_sha256 = cache_file_sha256,
      cache_manifest_sha256 = cache_manifest_sha256,
      cache_schema_version = cache$cache_schema_version
    ),
    cache = cache
  )
}

`.mlc_validate_index_structure` <- 
function(index, peaks, full = TRUE) {
  if (!is.list(index) || !isTRUE(index$valid) || !is.list(index$rows_by_spectrum_id)) {
    stop("Prebuilt spectrum_id_peak_index has an invalid schema.", call. = FALSE)
  }
  peaks <- as.data.frame(peaks, stringsAsFactors = FALSE, check.names = FALSE)
  if (!"spectrum_id" %in% names(peaks)) {
    stop("Cannot validate a spectrum index without peak-table spectrum_id.", call. = FALSE)
  }
  ids <- as.character(peaks$spectrum_id)
  valid_rows <- which(!is.na(ids))
  expected_keys <- unique(ids[valid_rows])
  observed_keys <- names(index$rows_by_spectrum_id)
  if (!setequal(observed_keys, expected_keys)) {
    stop("Prebuilt spectrum_id_peak_index key set does not match the cleaned peak table.", call. = FALSE)
  }
  row_counts <- vapply(index$rows_by_spectrum_id, function(rows) length(as.integer(rows)), integer(1))
  if (sum(row_counts) != length(valid_rows)) {
    stop("Prebuilt spectrum_id_peak_index has missing or duplicate row assignments.", call. = FALSE)
  }
  if (!isTRUE(full)) return(invisible(TRUE))
  assigned <- unlist(lapply(index$rows_by_spectrum_id, as.integer), use.names = FALSE)
  if (length(assigned) != length(valid_rows) || any(!is.finite(assigned)) ||
      any(assigned < 1L | assigned > nrow(peaks)) || anyDuplicated(assigned)) {
    stop("Prebuilt spectrum_id_peak_index has invalid, duplicate, or out-of-range row positions.", call. = FALSE)
  }
  if (!identical(sort(assigned), valid_rows)) {
    stop("Prebuilt spectrum_id_peak_index does not cover each valid peak row exactly once.", call. = FALSE)
  }
  for (spectrum_id in observed_keys) {
    rows <- as.integer(index$rows_by_spectrum_id[[spectrum_id]])
    if (!identical(as.character(ids[rows]), rep(as.character(spectrum_id), length(rows)))) {
      stop("Prebuilt spectrum_id_peak_index assigns a row to the wrong spectrum_id.", call. = FALSE)
    }
  }
  invisible(TRUE)
}

`.mlc_validate_library_identifier_mappings` <- 
function(cache) {
  required <- c("spectrum_id_mapping", "lab_id_mapping")
  if (length(setdiff(required, names(cache)))) {
    stop("Matching-library cache is missing stable spectrum_id/Lab.ID mappings.", call. = FALSE)
  }
  expected <- .mlc_library_identifier_mappings(cache$library_metadata)
  if (!identical(cache$spectrum_id_mapping, expected$spectrum_id_mapping) ||
      !identical(cache$lab_id_mapping, expected$lab_id_mapping)) {
    stop("Matching-library cache spectrum_id/Lab.ID mappings do not match library metadata.", call. = FALSE)
  }
  invisible(TRUE)
}

`.mlc_validate_prebuilt_spectrum_id_peak_index` <- 
function(spectrum_id_peak_index,
                                                          peaks,
                                                          matching_library_cache_identity = NULL) {
  if (!is.list(spectrum_id_peak_index) ||
      !all(c("index", "identity", "validated_cache_file_sha256", "validated_cache_manifest_sha256") %in% names(spectrum_id_peak_index))) {
    stop("An explicit spectrum_id_peak_index must be a validated matching-library cache index envelope.", call. = FALSE)
  }
  identity <- spectrum_id_peak_index$identity
  required <- c(
    "cache_schema_version", "cache_key_sha256", "source_database_sha256",
    "cleaned_peak_table_sha256", "spectrum_id_field", "n_spectra", "n_peaks",
    "n_unique_spectrum_id", "index_key_sha256", "index_object_sha256"
  )
  if (!is.list(identity) || length(setdiff(required, names(identity)))) {
    stop("Explicit spectrum_id_peak_index identity is incomplete.", call. = FALSE)
  }
  if (!identical(as.character(identity$cache_schema_version), .mlc_cache_schema_version())) {
    stop("Explicit spectrum_id_peak_index cache schema version is incompatible.", call. = FALSE)
  }
  if (!identical(as.character(identity$spectrum_id_field), "spectrum_id") ||
      !nzchar(as.character(identity$source_database_sha256)) ||
      !nzchar(as.character(spectrum_id_peak_index$validated_cache_file_sha256)) ||
      !nzchar(as.character(spectrum_id_peak_index$validated_cache_manifest_sha256))) {
    stop("Explicit spectrum_id_peak_index lacks validated cache-file/source identity.", call. = FALSE)
  }
  cache_required <- c(
    "cache_key_sha256", "source_database_sha256", "cache_file_sha256",
    "cache_manifest_sha256", "cache_schema_version"
  )
  if (!is.list(matching_library_cache_identity) ||
      length(setdiff(cache_required, names(matching_library_cache_identity)))) {
    stop("Explicit spectrum_id_peak_index requires matching validated cache identity on the Stage 5 library.", call. = FALSE)
  }
  if (!identical(as.character(identity$cache_key_sha256), as.character(matching_library_cache_identity$cache_key_sha256)) ||
      !identical(tolower(as.character(identity$source_database_sha256)), tolower(as.character(matching_library_cache_identity$source_database_sha256))) ||
      !identical(as.character(identity$cache_schema_version), as.character(matching_library_cache_identity$cache_schema_version)) ||
      !identical(tolower(as.character(spectrum_id_peak_index$validated_cache_file_sha256)), tolower(as.character(matching_library_cache_identity$cache_file_sha256)) ) ||
      !identical(tolower(as.character(spectrum_id_peak_index$validated_cache_manifest_sha256)), tolower(as.character(matching_library_cache_identity$cache_manifest_sha256)) )) {
    stop("Explicit spectrum_id_peak_index source database, cache schema, or cache-file identity does not match Stage 5.", call. = FALSE)
  }
  peaks <- as.data.frame(peaks, stringsAsFactors = FALSE, check.names = FALSE)
  index <- spectrum_id_peak_index$index
  .mlc_validate_index_structure(index, peaks, full = TRUE)
  ids <- as.character(peaks$spectrum_id)
  observed <- list(
    cleaned_peak_table_sha256 = .mlc_sha256_object(peaks),
    n_spectra = as.integer(length(unique(ids[!is.na(ids)]))),
    n_peaks = as.integer(nrow(peaks)),
    n_unique_spectrum_id = as.integer(length(unique(ids[!is.na(ids)]))),
    index_key_sha256 = .mlc_sha256_object(sort(names(index$rows_by_spectrum_id))),
    index_object_sha256 = .mlc_sha256_object(index)
  )
  if (!identical(unname(observed), unname(identity[names(observed)]))) {
    stop("Explicit spectrum_id_peak_index does not match the current cleaned peak table.", call. = FALSE)
  }
  index
}

`.mlc_validate_standard_metadata_overlay_applied` <- 
function(cache) {
  overlay_key <- cache$standard_metadata_overlay
  overlay <- .mlc_read_standard_metadata_overlay(overlay_key$path)
  if (!identical(as.character(overlay$overlay_sha256), as.character(overlay_key$overlay_sha256)) ||
      !identical(as.character(overlay$schema_version), as.character(overlay_key$schema_version)) ||
      !identical(as.character(overlay$authoritative_mapping_source_sha256),
                 as.character(overlay_key$authoritative_mapping_source_sha256))) {
    stop("Matching-library cache standard metadata overlay identity is invalid.", call. = FALSE)
  }
  metadata <- as.data.frame(cache$library_metadata, stringsAsFactors = FALSE, check.names = FALSE)
  index <- as.data.frame(cache$fragment_search_index, stringsAsFactors = FALSE, check.names = FALSE)
  target <- match(overlay$table$spectrum_id, as.character(metadata$spectrum_id))
  required <- c("HMDB_ID", "hmdb_id", "compound_name", "metabolite_name", "source_filename")
  if (anyNA(target) || length(setdiff(required, names(metadata))) ||
      any(as.character(metadata$candidate_id[target]) != overlay$table$canonical_HMDB_ID) ||
      any(as.character(metadata$compound_id[target]) != overlay$table$canonical_HMDB_ID) ||
      any(as.character(metadata$HMDB_ID[target]) != overlay$table$canonical_HMDB_ID) ||
      any(as.character(metadata$hmdb_id[target]) != overlay$table$canonical_HMDB_ID) ||
      any(as.character(metadata$candidate_name[target]) != overlay$table$canonical_compound_name) ||
      any(as.character(metadata$compound_name[target]) != overlay$table$canonical_compound_name) ||
      any(as.character(metadata$metabolite_name[target]) != overlay$table$canonical_compound_name) ||
      any(as.character(metadata$source_filename[target]) != overlay$table$source_filename)) {
    stop("Matching-library cache standard metadata overlay was not applied consistently.", call. = FALSE)
  }
  target_index <- as.character(index$spectrum_id) %in% overlay$table$spectrum_id
  ref <- match(as.character(index$spectrum_id[target_index]), overlay$table$spectrum_id)
  if (any(as.character(index$candidate_id[target_index]) != overlay$table$canonical_HMDB_ID[ref]) ||
      any(as.character(index$compound_id[target_index]) != overlay$table$canonical_HMDB_ID[ref]) ||
      any(as.character(index$candidate_name[target_index]) != overlay$table$canonical_compound_name[ref])) {
    stop("Matching-library cache fragment-index overlay metadata was not applied consistently.", call. = FALSE)
  }
  invisible(TRUE)
}

`.msi_adduct_status` <- 
function(adduct) {
  value <- trimws(as.character(adduct))
  ifelse(!is.na(value) & nzchar(value), "explicit", "unknown")
}

`.msi_config_or` <- 
function(x, default = NULL) {
  if (is.null(x) || !length(x)) default else x
}

`.msi_config_scalar` <- 
function(x, default = NULL) {
  x <- .msi_config_or(x, default)
  if (is.null(x) || !length(x)) return(default)
  x[[1]]
}

`.msi_default_tolerance_contract` <- 
function() {
  list(
    contract_version = "tolerance_contract_v1",
    fragment_retrieval = list(
      tolerance_mode = "hybrid_low_mz",
      ppm = 10,
      low_mz_cutoff = 400,
      low_mz_da_floor = 0.01
    ),
    fragment_evidence = list(
      tolerance_mode = "hybrid_low_mz",
      ppm = 10,
      low_mz_cutoff = 400,
      low_mz_da_floor = 0.005
    ),
    precursor_candidate = list(
      mode = "hybrid_feature_floor",
      ppm = 10,
      da_floor_mode = "auto_feature_half_bin"
    )
  )
}

`.msi_effective_precursor_tolerance` <- 
function(observed_mz, contract) {
  mz <- suppressWarnings(as.numeric(observed_mz[[1]]))
  if (!is.finite(mz) || mz <= 0) return(NA_real_)
  node <- contract$precursor_candidate
  ppm <- suppressWarnings(as.numeric(.msi_config_scalar(node$ppm, NA_real_)))
  floor <- suppressWarnings(as.numeric(.msi_config_scalar(node$da_floor, NA_real_)))
  if (!is.finite(ppm) || !is.finite(floor)) return(NA_real_)
  max(mz * ppm / 1e6, floor)
}

`.msi_infer_feature_mz_precision` <- 
function(feature_mz) {
  mz <- suppressWarnings(as.numeric(feature_mz))
  mz <- sort(unique(mz[is.finite(mz) & mz > 0]))
  if (!length(mz)) {
    return(list(method = "unknown", digits = NA_integer_, quantum = NA_real_, half_bin = NA_real_))
  }
  for (digits in 0:8) {
    quantum <- 10^(-digits)
    tolerance <- max(1e-10, quantum * 1e-7)
    if (max(abs(mz - round(mz, digits)), na.rm = TRUE) <= tolerance) {
      return(list(
        method = "decimal_quantized_feature_axis",
        digits = as.integer(digits),
        quantum = quantum,
        half_bin = quantum / 2
      ))
    }
  }
  list(method = "unknown", digits = NA_integer_, quantum = NA_real_, half_bin = NA_real_)
}

`.msi_normalize_tolerance_mode` <- 
function(x, default = "hybrid_low_mz") {
  value <- as.character(.msi_config_scalar(x, default))
  value <- tolower(trimws(value[[1]]))
  if (value %in% c("hybrid_low_mz_fixed_da", "hybrid", "hybrid_low_mz")) {
    return("hybrid_low_mz")
  }
  if (value %in% c("strict_ppm", "ppm")) return("strict_ppm")
  if (value %in% c("bin_aware")) return("bin_aware")
  value
}

`.msi_recursive_merge` <- 
function(default, override) {
  if (is.null(override) || !length(override)) return(default)
  if (!is.list(default) || !is.list(override)) return(override)
  out <- default
  for (nm in names(override)) {
    out[[nm]] <- if (is.list(override[[nm]]) && is.list(out[[nm]])) {
      .msi_recursive_merge(out[[nm]], override[[nm]])
    } else {
      override[[nm]]
    }
  }
  out
}

`.msi_resolve_tolerance_contract` <- 
function(contract = NULL, feature_mz = NULL) {
  if (is.null(contract) || !length(contract)) contract <- list()
  resolved <- .msi_recursive_merge(.msi_default_tolerance_contract(), contract)
  resolved$fragment_retrieval$tolerance_mode <- .msi_normalize_tolerance_mode(
    resolved$fragment_retrieval$tolerance_mode
  )
  resolved$fragment_evidence$tolerance_mode <- .msi_normalize_tolerance_mode(
    resolved$fragment_evidence$tolerance_mode
  )
  .msi_validate_tolerance_contract(resolved)
  precision <- .msi_infer_feature_mz_precision(feature_mz)
  precursor <- resolved$precursor_candidate
  if (identical(as.character(precursor$da_floor_mode), "auto_feature_half_bin")) {
    if (!is.finite(precision$half_bin)) {
      stop(
        "Cannot resolve precursor Da floor: feature m/z quantization is not identifiable.",
        call. = FALSE
      )
    }
    precursor$da_floor <- precision$half_bin
  }
  resolved$precursor_candidate <- precursor
  resolved$feature_mz_precision <- precision
  resolved$resolved <- TRUE
  resolved
}

`.msi_validate_tolerance_contract` <- 
function(contract) {
  required <- c("fragment_retrieval", "fragment_evidence", "precursor_candidate")
  if (!is.list(contract) || !all(required %in% names(contract))) {
    stop("tolerance_contract is missing one or more required policy sections.", call. = FALSE)
  }
  for (section in c("fragment_retrieval", "fragment_evidence")) {
    node <- contract[[section]]
    ppm <- suppressWarnings(as.numeric(.msi_config_scalar(node$ppm, NA_real_)))
    cutoff <- suppressWarnings(as.numeric(.msi_config_scalar(node$low_mz_cutoff, NA_real_)))
    floor <- suppressWarnings(as.numeric(.msi_config_scalar(node$low_mz_da_floor, NA_real_)))
    mode <- as.character(.msi_config_scalar(node$tolerance_mode, ""))
    if (!is.finite(ppm) || ppm <= 0 || !is.finite(cutoff) || cutoff <= 0 ||
        !is.finite(floor) || floor < 0 || !identical(mode, "hybrid_low_mz")) {
      stop("Invalid tolerance_contract.", call. = FALSE)
    }
  }
  precursor <- contract$precursor_candidate
  precursor_ppm <- suppressWarnings(as.numeric(.msi_config_scalar(precursor$ppm, NA_real_)))
  precursor_mode <- as.character(.msi_config_scalar(precursor$mode, ""))
  floor_mode <- as.character(.msi_config_scalar(precursor$da_floor_mode, ""))
  if (!is.finite(precursor_ppm) || precursor_ppm <= 0 ||
      !identical(precursor_mode, "hybrid_feature_floor") ||
      !floor_mode %in% c("auto_feature_half_bin", "explicit")) {
    stop("Invalid tolerance_contract.precursor_candidate.", call. = FALSE)
  }
  if (identical(floor_mode, "explicit")) {
    floor <- suppressWarnings(as.numeric(.msi_config_scalar(precursor$da_floor, NA_real_)))
    if (!is.finite(floor) || floor < 0) stop("Invalid explicit precursor Da floor.", call. = FALSE)
  }
  invisible(TRUE)
}

`.reverse_dp_first_column` <- 
function(x, candidates, required = TRUE) {
  hit <- intersect(candidates, names(x))
  if (length(hit)) return(hit[[1]])
  if (isTRUE(required)) {
    stop("Missing required spectrum column; expected one of: ",
         paste(candidates, collapse = ", "), call. = FALSE)
  }
  NULL
}

`.reverse_dp_normalize_max100` <- 
function(x) {
  x <- suppressWarnings(as.numeric(x))
  x[!is.finite(x) | x < 0] <- 0
  base <- suppressWarnings(max(x, na.rm = TRUE))
  if (!is.finite(base) || base <= 0) return(rep(0, length(x)))
  x / base * 100
}

`.rm_add_missing_columns` <- 
function(df, cols) {
  df <- as.data.frame(df, stringsAsFactors = FALSE, check.names = FALSE)
  for (col in cols) {
    if (!col %in% names(df)) {
      df[[col]] <- rep(NA, nrow(df))
    }
  }
  df
}

`.rm_align_reverse_pairs` <- 
function(pairs) {
  pairs <- as.data.frame(pairs, stringsAsFactors = FALSE, check.names = FALSE)
  if (nrow(pairs) == 0L) {
    return(pairs)
  }
  if (!"library_peak_order" %in% names(pairs) || all(is.na(pairs$library_peak_order))) {
    pairs$library_peak_order <- seq_len(nrow(pairs))
  }
  pairs$library_peak_order <- suppressWarnings(as.integer(pairs$library_peak_order))
  pairs <- pairs[order(pairs$library_peak_order, pairs$library_peak_id), , drop = FALSE]
  pairs$query_intensity[!(pairs$matched %in% TRUE)] <- 0
  pairs
}

`.rm_archive_from_long_table` <- 
function(df) {
  df <- as.data.frame(df, stringsAsFactors = FALSE, check.names = FALSE)
  df <- .rm_standardize_archive_columns(df)
  if (!"spectrum_id" %in% names(df) ||
      all(is.na(df$spectrum_id) | !nzchar(as.character(df$spectrum_id)))) {
    df$spectrum_id <- paste0("ARCHIVE_SPECTRUM_", as.integer(factor(
      if ("candidate_id" %in% names(df)) as.character(df$candidate_id) else seq_len(nrow(df)),
      levels = unique(if ("candidate_id" %in% names(df)) as.character(df$candidate_id) else seq_len(nrow(df)))
    )))
  }
  if (!"candidate_id" %in% names(df)) {
    df$candidate_id <- as.character(df$spectrum_id)
  }
  spectrum <- .rm_add_missing_columns(df, colnames(.rm_empty_archive_spectrum_table()))
  spectrum <- spectrum[!duplicated(spectrum$spectrum_id), colnames(.rm_empty_archive_spectrum_table()), drop = FALSE]
  peaks <- .rm_add_missing_columns(df, colnames(.rm_empty_archive_peak_table()))
  peaks <- peaks[, colnames(.rm_empty_archive_peak_table()), drop = FALSE]
  list(archive_spectrum_table = spectrum, archive_peak_table = peaks)
}

`.rm_bool` <- 
function(x, default = FALSE) {
  if (is.null(x) || length(x) == 0L) return(default)
  if (is.logical(x)) return(isTRUE(x[[1]]))
  value <- tolower(trimws(as.character(x[[1]])))
  if (value %in% c("true", "t", "1", "yes", "y")) return(TRUE)
  if (value %in% c("false", "f", "0", "no", "n")) return(FALSE)
  default
}

`.rm_build_matched_peak_pairs_summary` <- 
function(matched_peak_pairs = NULL,
                                                 matched_peak_pairs_all = NULL,
                                                 reverse_match = NULL,
                                                 reverse_match_all = NULL,
                                                 export_matched_peak_pairs = FALSE,
                                                 export_matched_peak_pairs_all = FALSE,
                                                 store_matched_peak_pairs_in_rds = FALSE,
                                                 max_matched_peak_pairs_export_rows = 100000L,
                                                 compress_large_pair_tables = TRUE,
                                                 skipped_matched_peak_pairs_csv = !isTRUE(export_matched_peak_pairs),
                                                 skipped_matched_peak_pairs_all_csv = !isTRUE(export_matched_peak_pairs_all),
                                                 export_reason = "export_disabled",
                                                 n_matched_peak_pairs_exported = 0L,
                                                 n_matched_peak_pairs_all_exported = 0L,
                                                 matched_peak_pairs_export_truncated = FALSE,
                                                 matched_peak_pairs_all_export_truncated = FALSE,
                                                 matched_peak_pairs_output = NA_character_,
                                                 matched_peak_pairs_all_output = NA_character_,
                                                 timestamp = Sys.time()) {
  pairs <- if (is.data.frame(matched_peak_pairs)) matched_peak_pairs else data.frame()
  pairs_all <- if (is.data.frame(matched_peak_pairs_all)) matched_peak_pairs_all else data.frame()
  rm <- if (is.data.frame(reverse_match)) reverse_match else data.frame()
  rm_all <- if (is.data.frame(reverse_match_all)) reverse_match_all else data.frame()
  candidates_for_counts <- if (nrow(rm_all)) rm_all else if (nrow(rm)) rm else if (nrow(pairs_all)) pairs_all else pairs
  query_key_cols <- if ("query_id" %in% names(candidates_for_counts)) {
    c("query_id", "anchor_feature_id")
  } else {
    c("cluster_id", "anchor_feature_id")
  }
  candidate_key_cols <- c("candidate_key", "cluster_id", "anchor_feature_id")
  if (!all(candidate_key_cols %in% names(candidates_for_counts))) {
    candidate_key_cols <- c("candidate_id", "spectrum_id", "cluster_id", "anchor_feature_id")
  }
  data.frame(
    n_matched_peak_pairs_total = as.integer(nrow(pairs)),
    n_matched_peak_pairs_all_total = as.integer(nrow(pairs_all)),
    n_unique_candidates = .rm_n_unique_key(candidates_for_counts, candidate_key_cols),
    n_unique_queries = .rm_n_unique_key(candidates_for_counts, query_key_cols),
    n_unique_library_spectra = .rm_n_unique_value(candidates_for_counts, c("spectrum_id", "lab_id")),
    export_matched_peak_pairs = isTRUE(export_matched_peak_pairs),
    export_matched_peak_pairs_all = isTRUE(export_matched_peak_pairs_all),
    store_matched_peak_pairs_in_rds = isTRUE(store_matched_peak_pairs_in_rds),
    max_matched_peak_pairs_export_rows = as.integer(max_matched_peak_pairs_export_rows),
    compress_large_pair_tables = isTRUE(compress_large_pair_tables),
    skipped_matched_peak_pairs_csv = isTRUE(skipped_matched_peak_pairs_csv),
    skipped_matched_peak_pairs_all_csv = isTRUE(skipped_matched_peak_pairs_all_csv),
    export_reason = as.character(export_reason[[1]]),
    n_matched_peak_pairs_exported = as.integer(n_matched_peak_pairs_exported),
    n_matched_peak_pairs_all_exported = as.integer(n_matched_peak_pairs_all_exported),
    matched_peak_pairs_export_truncated = isTRUE(matched_peak_pairs_export_truncated),
    matched_peak_pairs_all_export_truncated = isTRUE(matched_peak_pairs_all_export_truncated),
    matched_peak_pairs_output = as.character(matched_peak_pairs_output[[1]]),
    matched_peak_pairs_all_output = as.character(matched_peak_pairs_all_output[[1]]),
    timestamp = format(as.POSIXct(timestamp), "%Y-%m-%d %H:%M:%S"),
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
}

`.rm_build_spectrum_id_peak_index` <- 
function(peaks) {
  if (is.null(peaks)) {
    return(list(rows_by_spectrum_id = list(), valid = FALSE))
  }
  peaks <- as.data.frame(peaks, stringsAsFactors = FALSE, check.names = FALSE)
  required <- c("spectrum_id", "peak_mz", "intensity_sum1", "intensity_raw")
  if (!all(required %in% names(peaks)) ||
      !nrow(peaks) ||
      all(is.na(peaks$intensity_sum1))) {
    return(list(rows_by_spectrum_id = list(), valid = FALSE))
  }
  ids <- as.character(peaks$spectrum_id)
  valid <- !is.na(ids)
  rows_by_spectrum_id <- if (any(valid)) {
    split(which(valid), ids[valid], drop = TRUE)
  } else {
    list()
  }
  list(rows_by_spectrum_id = rows_by_spectrum_id, valid = TRUE)
}

`.rm_emit_progress` <- 
function(progress_callback = NULL,
                              progress_path = NULL,
                              started_at = NULL,
                              processed_candidates = 0L,
                              total_candidates = 0L,
                              reverse_rows_so_far = 0L,
                              pass_rows_so_far = 0L,
                              message = "") {
  started_at <- if (inherits(started_at, "POSIXt")) started_at else Sys.time()
  elapsed <- as.numeric(difftime(Sys.time(), started_at, units = "secs"))
  row <- data.frame(
    timestamp = format(Sys.time(), "%Y-%m-%d %H:%M:%S"),
    processed_candidates = as.integer(processed_candidates),
    total_candidates = as.integer(total_candidates),
    elapsed_seconds = as.numeric(elapsed),
    reverse_rows_so_far = as.integer(reverse_rows_so_far),
    pass_rows_so_far = as.integer(pass_rows_so_far),
    message = as.character(message),
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
  if (!is.null(progress_path) && length(progress_path) && nzchar(as.character(progress_path[[1]]))) {
    path <- as.character(progress_path[[1]])
    dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
    utils::write.table(
      row,
      file = path,
      sep = ",",
      row.names = FALSE,
      col.names = !file.exists(path),
      append = file.exists(path),
      qmethod = "double"
    )
  }
  if (is.function(progress_callback)) {
    try(progress_callback(row), silent = TRUE)
  }
  message(
    "Module 6 progress: processed_candidates=", row$processed_candidates[[1]],
    "; total_candidates=", row$total_candidates[[1]],
    "; elapsed_seconds=", round(row$elapsed_seconds[[1]], 2),
    "; reverse_rows_so_far=", row$reverse_rows_so_far[[1]],
    "; pass_rows_so_far=", row$pass_rows_so_far[[1]],
    if (nzchar(row$message[[1]])) paste0("; message=", row$message[[1]]) else ""
  )
  try(flush.console(), silent = TRUE)
  invisible(row)
}

`.rm_empty_archive_hits` <- 
function() {
  data.frame(
    cluster_id = character(),
    anchor_feature_id = character(),
    anchor_mz = numeric(),
    delta = numeric(),
    expected_precursor = numeric(),
    archive_id = character(),
    spectrum_id = character(),
    candidate_id = character(),
    candidate_name = character(),
    library_source = character(),
    precursor_mz = numeric(),
    formula = character(),
    inchikey = character(),
    smiles = character(),
    database_source = character(),
    adduct = character(),
    adduct_status = character(),
    polarity = character(),
    s_arch = numeric(),
    n_match_arch = integer(),
    matched_query_peak_n_arch = integer(),
    matched_archive_peak_n = integer(),
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
}

`.rm_empty_archive_peak_table` <- 
function() {
  data.frame(
    archive_id = character(),
    spectrum_id = character(),
    candidate_id = character(),
    peak_mz = numeric(),
    intensity_raw = numeric(),
    intensity_sum1 = numeric(),
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
}

`.rm_empty_archive_spectrum_table` <- 
function() {
  data.frame(
    archive_id = character(),
    spectrum_id = character(),
    candidate_id = character(),
    candidate_name = character(),
    library_source = character(),
    precursor_mz = numeric(),
    formula = character(),
    inchikey = character(),
    smiles = character(),
    database_source = character(),
    adduct = character(),
    polarity = character(),
    source_file = character(),
    scan = character(),
    metadata_json = character(),
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
}

`.rm_empty_merged_candidates` <- 
function() {
  data.frame(
    cluster_id = character(),
    anchor_feature_id = character(),
    anchor_mz = numeric(),
    query_id = character(),
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
    adduct_status = character(),
    polarity = character(),
    entropy_S = numeric(),
    S_FE = numeric(),
    s_fe = numeric(),
    retrieval_rank = integer(),
    n_FE_match = integer(),
    matched_query_peaks = integer(),
    matched_library_peaks = integer(),
    precursor_error_ppm = numeric(),
    s_arch = numeric(),
    matched_query_peak_n_arch = integer(),
    matched_archive_peak_n = integer(),
    candidate_source = character(),
    candidate_retrieved = logical(),
    precursor_gate_pass = logical(),
    fragment_retrieval_hits = integer(),
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
}

`.rm_empty_query_spectrum` <- 
function() {
  data.frame(
    cluster_id = character(),
    anchor_feature_id = character(),
    anchor_mz = numeric(),
    query_peak_id = character(),
    query_mz = numeric(),
    query_intensity = numeric(),
    query_intensity_raw = numeric(),
    SAS = numeric(),
    graph_support = numeric(),
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
}

`.rm_empty_rdp_filter_summary` <- 
function() {
  data.frame(
    query_id = character(),
    anchor_feature_id = character(),
    n_candidates_after_entropy = integer(),
    n_pseudo_ms2_candidate_evidence = integer(),
    min_matched_peaks_evidence = integer(),
    min_reverse_masked_cosine_evidence = numeric(),
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
}

`.rm_empty_reverse_matched_peaks` <- 
function() {
  data.frame(
    library_peak_order = integer(),
    library_peak_id = character(),
    library_mz = numeric(),
    library_fragment_mz = numeric(),
    library_intensity = numeric(),
    library_fragment_intensity = numeric(),
    library_intensity_raw = numeric(),
    library_intensity_norm = numeric(),
    library_relative_intensity = numeric(),
    query_peak_id = character(),
    query_fragment_feature_id = character(),
    query_mz = numeric(),
    query_fragment_mz = numeric(),
    query_intensity = numeric(),
    query_intensity_raw = numeric(),
    query_intensity_norm = numeric(),
    query_relative_intensity = numeric(),
    mz_error = numeric(),
    mz_error_ppm = numeric(),
    SAS = numeric(),
    query_sas = numeric(),
    contribution_to_dot = numeric(),
    contribution_to_S_SAS = numeric(),
    query_full_l2 = numeric(),
    query_matched_l2 = numeric(),
    library_all_l2 = numeric(),
    query_total_intensity = numeric(),
    matched_query_intensity = numeric(),
    matched_library_intensity = numeric(),
    total_library_intensity = numeric(),
    matched_query_intensity_fraction = numeric(),
    matched_library_intensity_fraction = numeric(),
    library_l2_coverage = numeric(),
    reverse_masked_cosine = numeric(),
    match_id = character(),
    n_query_peaks_total = integer(),
    n_library_peaks_total = integer(),
    matched = logical(),
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
}

`.rm_empty_reverse_matched_peaks_full` <- 
function() {
  cbind(
    data.frame(
      cluster_id = character(),
      anchor_feature_id = character(),
      query_id = character(),
      candidate_key = character(),
      candidate_id = character(),
      spectrum_id = character(),
      compound_name = character(),
      candidate_name = character(),
      stringsAsFactors = FALSE,
      check.names = FALSE
    ),
    .rm_empty_reverse_matched_peaks()
  )
}

`.rm_empty_reverse_scores` <- 
function(base_only = FALSE) {
  base <- data.frame(
    cluster_id = character(),
    anchor_feature_id = character(),
    anchor_mz = numeric(),
    candidate_key = character(),
    candidate_id = character(),
    spectrum_id = character(),
    archive_id = character(),
    reverse_masked_cosine = numeric(),
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
  if (isTRUE(base_only)) {
    return(base)
  }
  data.frame(
    component = character(),
    cluster_id = character(),
    anchor_feature_id = character(),
    feature_id = character(),
    anchor_mz = numeric(),
    candidate_key = character(),
    candidate_id = character(),
    spectrum_id = character(),
    archive_id = character(),
    query_id = character(),
    query_spectrum_id = character(),
    target_hmdb_id = character(),
    metabolite_name = character(),
    candidate_name = character(),
    library_source = character(),
    database_source = character(),
    source = character(),
    lab_id = character(),
    instrument = character(),
    precursor_mz = numeric(),
    formula = character(),
    inchikey = character(),
    smiles = character(),
    adduct = character(),
    adduct_status = character(),
    polarity = character(),
    candidate_source = character(),
    candidate_retrieved = logical(),
    precursor_gate_pass = logical(),
    fragment_retrieval_hits = integer(),
    entropy_S = numeric(),
    S_FE = numeric(),
    retrieval_rank = integer(),
    rank_entropy = integer(),
    precursor_error_ppm = numeric(),
    S_precursor = numeric(),
    reverse_masked_cosine = numeric(),
    mass_error_ppm_summary = character(),
    mass_error_ppm_median = numeric(),
    mass_error_ppm_max = numeric(),
    matched_peak_pairs_id = character(),
    query_unmatched_peak_policy = character(),
    reference_unmatched_peak_policy = character(),
    S_lib_coverage = numeric(),
    library_l2_coverage = numeric(),
    matched_library_intensity_fraction = numeric(),
    matched_query_intensity_fraction = numeric(),
    matched_peak_fraction = numeric(),
    n_matched_peaks = integer(),
    n_reverse_match = integer(),
    S_SAS = numeric(),
    S_SAS_score_for_total = numeric(),
    S_SAS_missing = logical(),
    S_total = numeric(),
    rank_total = integer(),
    rank_within_feature = integer(),
    rank_margin_to_next = numeric(),
    S_spatial_reverse = numeric(),
    n_SAS_matched_peaks = integer(),
    sum_matched_library_intensity = numeric(),
    matched_library_intensity = numeric(),
    matched_query_intensity = numeric(),
    total_library_intensity = numeric(),
    n_pseudo_peaks = integer(),
    n_library_peaks = integer(),
    reverse_match_status = character(),
    s_fe = numeric(),
    s_arch = numeric(),
    s_rev = numeric(),
    n_match = integer(),
    u_spec = numeric(),
    s_spat = numeric(),
    matched_query_peak_n = integer(),
    matched_reference_peak_n = integer(),
    n_match_pairs = integer(),
    unique_query_matched_peaks = integer(),
    unique_library_matched_peaks = integer(),
    n_query_peaks_total = integer(),
    n_library_peaks_total = integer(),
    query_match_ratio = numeric(),
    library_match_ratio = numeric(),
    query_peak_reuse_n = integer(),
    query_peak_reuse_fraction = numeric(),
    has_query_peak_reuse = logical(),
    match_mode = character(),
    score_ms2 = numeric(),
    score_spatial = numeric(),
    score_total_pre_mfc = numeric(),
    pass_n_matched_peaks_ge3 = logical(),
    pass_reverse_masked_cosine = logical(),
    evidence_pass = logical(),
    fragment_evidence_hits = integer(),
    evidence_eligible = logical(),
    evidence_eligibility_reason = character(),
    min_matched_peaks_used = integer(),
    min_reverse_masked_cosine_used = numeric(),
    evidence_failure_reason = character(),
    evidence_filter_reason = character(),
    confidence_tier = character(),
    chem_status = character(),
    chem_role = character(),
    chem_explanation = character(),
    role_delta_mz = numeric(),
    role_ppm_error = numeric(),
    role_evidence_n = integer(),
    spectral_shape_support = character(),
    passes_spectral_shape_filter = logical(),
    pass_min_matched_peaks = logical(),
    pass_reverse_match = logical(),
    reverse_filter_mode = I(character()),
    threshold_validation_status = character(),
    selection_reason = character(),
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
}

`.rm_merge_candidate_rows` <- 
function(df) {
  out <- df[1, , drop = FALSE]
  for (col in names(df)) {
    if (col == "candidate_source") next
    values <- df[[col]]
    if (is.numeric(values) || is.integer(values)) {
      if (all(is.na(values))) {
        out[[col]] <- values[[1]]
      } else if (col %in% c("retrieval_rank")) {
        out[[col]] <- min(values, na.rm = TRUE)
      } else {
        out[[col]] <- max(values, na.rm = TRUE)
      }
    } else {
      idx <- which(!is.na(values) & nzchar(as.character(values)))
      out[[col]] <- if (length(idx) > 0L) values[[idx[[1]]]] else values[[1]]
    }
  }
  sources <- unique(as.character(df$candidate_source))
  out$candidate_source <- if (all(c("flash", "archive") %in% sources)) "both" else sources[[1]]
  out
}

`.rm_n_unique_key` <- 
function(df, cols) {
  if (!is.data.frame(df) || nrow(df) == 0L) return(0L)
  cols <- cols[cols %in% names(df)]
  if (!length(cols)) return(0L)
  key <- do.call(paste, c(lapply(cols, function(col) {
    value <- as.character(df[[col]])
    value[is.na(value)] <- ""
    value
  }), sep = "\r"))
  as.integer(length(unique(key[nzchar(key)])))
}

`.rm_n_unique_value` <- 
function(df, cols) {
  if (!is.data.frame(df) || nrow(df) == 0L) return(0L)
  cols <- cols[cols %in% names(df)]
  if (!length(cols)) return(0L)
  values <- unlist(lapply(cols, function(col) as.character(df[[col]])), use.names = FALSE)
  values <- values[!is.na(values) & nzchar(values)]
  as.integer(length(unique(values)))
}

`.rm_normalize_library_spectrum` <- 
function(x) {
  if (nrow(x) == 0L) return(x)
  total <- sum(as.numeric(x$library_intensity), na.rm = TRUE)
  if (is.finite(total) && total > 0) {
    x$library_intensity <- as.numeric(x$library_intensity) / total
  }
  rownames(x) <- NULL
  x
}

`.rm_read_archive_library` <- 
function(archive_library) {
  if (is.character(archive_library) && length(archive_library) == 1L) {
    path <- archive_library[[1]]
    if (!file.exists(path)) {
      stop("archive_library file does not exist: ", path, call. = FALSE)
    }
    ext <- tolower(tools::file_ext(path))
    if (ext == "rds") {
      archive_library <- readRDS(path)
    } else if (ext %in% c("rdata", "rda")) {
      env <- new.env(parent = emptyenv())
      loaded_names <- load(path, envir = env)
      archive_library <- lapply(loaded_names, function(nm) get(nm, envir = env))
      names(archive_library) <- loaded_names
    } else if (ext == "csv") {
      archive_library <- utils::read.csv(path, stringsAsFactors = FALSE, check.names = FALSE)
    } else if (ext %in% c("tsv", "txt")) {
      archive_library <- utils::read.delim(path, stringsAsFactors = FALSE, check.names = FALSE)
    } else {
      stop("Unsupported archive_library file type: ", ext, call. = FALSE)
    }
  }

  if (is.data.frame(archive_library)) {
    return(.rm_archive_from_long_table(archive_library))
  }

  if (is.list(archive_library) &&
      all(c("archive_spectrum_table", "archive_peak_table") %in% names(archive_library))) {
    return(list(
      archive_spectrum_table = archive_library$archive_spectrum_table,
      archive_peak_table = archive_library$archive_peak_table
    ))
  }

  if (is.list(archive_library) &&
      all(c("library_spectrum_table", "library_peak_table") %in% names(archive_library))) {
    spectrum <- archive_library$library_spectrum_table
    peaks <- archive_library$library_peak_table
    return(list(
      archive_spectrum_table = spectrum,
      archive_peak_table = peaks
    ))
  }

  if (is.list(archive_library)) {
    for (item in archive_library) {
      if (is.data.frame(item) ||
          (is.list(item) && any(c(
            "archive_spectrum_table", "archive_peak_table",
            "library_spectrum_table", "library_peak_table"
          ) %in% names(item)))) {
        return(.rm_read_archive_library(item))
      }
    }
  }

  stop("Unsupported archive_library input.", call. = FALSE)
}

`.rm_spectrum_id_index_rows` <- 
function(index, spectrum_id) {
  if (is.null(index) || !isTRUE(index$valid) ||
      length(spectrum_id) != 1L || is.na(spectrum_id)) {
    return(integer())
  }
  rows <- index$rows_by_spectrum_id[[as.character(spectrum_id)]]
  if (is.null(rows)) integer() else as.integer(rows)
}

`.rm_split_key` <- 
function(df, cols) {
  df <- as.data.frame(df, stringsAsFactors = FALSE, check.names = FALSE)
  if (nrow(df) == 0L) {
    return(character())
  }
  missing <- setdiff(cols, names(df))
  if (length(missing) > 0L) {
    stop("Cannot build reverse-match key; missing columns: ", paste(missing, collapse = ", "), call. = FALSE)
  }
  pieces <- lapply(cols, function(col) {
    value <- as.character(df[[col]])
    value[is.na(value)] <- "<NA>"
    paste0(nchar(value, type = "bytes"), ":", value)
  })
  do.call(paste0, pieces)
}

`.rm_standardize_archive_columns` <- 
function(df) {
  mapping <- list(
    archive_id = c("archive_id", "spectrum_id", "scan", "scan_id", "id", "ID"),
    spectrum_id = c("spectrum_id", "scan", "scan_id", "title", "id", "ID"),
    candidate_id = c("candidate_id", "compound_id", "hmdb_id", "HMDB.ID", "id", "ID"),
    candidate_name = c("candidate_name", "name", "compound_name", "Compound.name", "Name"),
    library_source = c("library_source", "source", "database_source"),
    precursor_mz = c("precursor_mz", "Precursor_mz", "precursorMz", "precursor", "parent_mz", "mz"),
    formula = c("formula", "Formula", "molecular_formula"),
    inchikey = c("inchikey", "INCHIKEY.ID", "InChIKey"),
    smiles = c("smiles", "SMILES", "SMILES.ID"),
    database_source = c("database_source", "source", "library_source"),
    adduct = c("adduct", "Adduct", "adduct_type"),
    polarity = c("polarity", "Polarity", "ion_mode", "mode"),
    source_file = c("source_file", "file"),
    scan = c("scan", "scan_id"),
    metadata_json = c("metadata_json"),
    peak_mz = c("peak_mz", "mz", "MZ", "fragment_mz", "product_mz"),
    intensity_raw = c("intensity_raw", "intensity", "Intensity", "into", "peak_intensity", "abundance"),
    intensity_sum1 = c("intensity_sum1")
  )
  for (new in names(mapping)) {
    if (new %in% names(df)) next
    old <- mapping[[new]][mapping[[new]] %in% names(df)]
    if (length(old) > 0L) {
      value <- df[[old[[1]]]]
      if (length(old) > 1L) {
        for (candidate in old[-1]) {
          missing <- is.na(value) | !nzchar(as.character(value))
          if (any(missing)) {
            value[missing] <- df[[candidate]][missing]
          }
        }
      }
      df[[new]] <- value
    }
  }
  df
}

`.rm_standardize_archive_peak_table` <- 
function(df) {
  df <- as.data.frame(df, stringsAsFactors = FALSE, check.names = FALSE)
  df <- .rm_standardize_archive_columns(df)
  df <- .rm_add_missing_columns(df, colnames(.rm_empty_archive_peak_table()))
  df[, colnames(.rm_empty_archive_peak_table()), drop = FALSE]
}

`.rm_standardize_archive_spectrum_table` <- 
function(df) {
  df <- as.data.frame(df, stringsAsFactors = FALSE, check.names = FALSE)
  df <- .rm_standardize_archive_columns(df)
  df <- .rm_add_missing_columns(df, colnames(.rm_empty_archive_spectrum_table()))
  if (nrow(df) > 0L) {
    missing_sid <- is.na(df$spectrum_id) | !nzchar(as.character(df$spectrum_id))
    df$spectrum_id[missing_sid] <- paste0("ARCHIVE_SPECTRUM_", which(missing_sid))
    missing_cid <- is.na(df$candidate_id) | !nzchar(as.character(df$candidate_id))
    df$candidate_id[missing_cid] <- as.character(df$spectrum_id[missing_cid])
    missing_archive <- is.na(df$archive_id) | !nzchar(as.character(df$archive_id))
    df$archive_id[missing_archive] <- as.character(df$spectrum_id[missing_archive])
  }
  df$precursor_mz <- suppressWarnings(as.numeric(df$precursor_mz))
  df[!duplicated(df$spectrum_id), colnames(.rm_empty_archive_spectrum_table()), drop = FALSE]
}

`5_calculate_flash_entropy` <- 
function(object,
                                        spectral_library = NULL,
                                        database_root = NULL,
                                        ms1_database_dir = NULL,
                                        ms2_database_dir = NULL,
                                        ppm = 10,
                                        mz_tol_da = NULL,
                                        mz_tol_da_max_mz = Inf,
                                        top_n_spectra = 500,
                                        top_n_compounds = 100,
                                        max_spectra_per_compound = 3,
                                        min_match_fe = 1,
                                        min_s_fe = 0,
                                        top_n_entropy_candidates = 500,
                                        precursor_filter = FALSE,
                                        precursor_tol = 0.05,
  polarity = c("both", "positive", "negative"),
                                        min_library_peaks = 2,
                                        max_library_spectra = Inf,
                                        top_n = NULL,
                                        min_score = NULL,
                                        progress_every = NULL,
                                        progress_path = NULL,
                                        progress_callback = NULL,
                                        progress_started_at = NULL,
                                        output_dir = "output",
                                        matching_library_cache_path = NULL,
                                        matching_library_cache_mode = c("error", "build_if_missing"),
                                        matching_library_cache_source_object_name = NULL,
                                        matching_library_cache_config = NULL,
                                        tolerance_contract = NULL,
                                        ...) {
  polarity <- match.arg(polarity)
  if (!methods::is(object, "msi_dataset")) {
    stop("object must be a msi_dataset.", call. = FALSE)
  }
  if (!is.list(object@pseudo_ms2) ||
      is.null(object@pseudo_ms2$pseudo_ms2_table) ||
      !is.data.frame(object@pseudo_ms2$pseudo_ms2_table)) {
    .msi_missing_summary("5_calculate_flash_entropy", "pseudo_ms2$pseudo_ms2_table", output_dir)
    return(object)
  }
  if (!is.null(tolerance_contract)) {
    if (!exists(".msi_resolve_tolerance_contract", mode = "function")) {
      stop("Shared tolerance contract helpers are unavailable.", call. = FALSE)
    }
    tolerance_contract <- .msi_resolve_tolerance_contract(
      tolerance_contract,
      feature_mz = object@variable_info$mz
    )
    retrieval_policy <- tolerance_contract$fragment_retrieval
    ppm <- as.numeric(retrieval_policy$ppm)
    mz_tol_da <- as.numeric(retrieval_policy$low_mz_da_floor)
    mz_tol_da_max_mz <- as.numeric(retrieval_policy$low_mz_cutoff)
  }
  dots <- list(...)
  if (!is.null(top_n)) {
    top_n_spectra <- top_n
  } else if (!is.null(dots$top_n)) {
    top_n_spectra <- dots$top_n
  }
  if (!is.null(min_score)) {
    min_s_fe <- min_score
  } else if (!is.null(dots$min_score)) {
    min_s_fe <- dots$min_score
  }
  top_n_entropy_candidates <- as.integer(top_n_entropy_candidates[[1]])
  if (!is.finite(top_n_entropy_candidates) || top_n_entropy_candidates < 1L) {
    stop("top_n_entropy_candidates must be a positive integer.", call. = FALSE)
  }
  if (is.finite(suppressWarnings(as.numeric(top_n_compounds[[1]]))) &&
      top_n_entropy_candidates > as.integer(top_n_compounds[[1]])) {
    message(
      "flash_entropy.top_n_candidates exceeds top_n_compounds; compound-level tables ",
      "are capped by top_n_compounds, while entropy_candidates remain spectrum-level."
    )
  }

  matching_library_cache <- NULL
  if (!is.null(matching_library_cache_path) && length(matching_library_cache_path) &&
      nzchar(as.character(matching_library_cache_path[[1L]]))) {
    matching_library_cache_mode <- match.arg(matching_library_cache_mode)
    effective_cache_config <- .mlc_effective_cache_config(
      polarity = polarity,
      min_library_peaks = min_library_peaks,
      ppm = ppm,
      mz_tol_da = mz_tol_da,
      mz_tol_da_max_mz = mz_tol_da_max_mz,
      extra_config = matching_library_cache_config
    )
    matching_library_cache <- .mlc_resolve_matching_library_cache(
      source_database_path = spectral_library,
      source_object_name = matching_library_cache_source_object_name,
      cache_path = matching_library_cache_path,
      effective_cleaning_config = effective_cache_config,
      cache_mode = matching_library_cache_mode
    )
  }

  if (!is.null(matching_library_cache)) {
    library <- matching_library_cache$library
  } else if (is.null(spectral_library) &&
      is.list(object@spectral_library) &&
      all(c("library_spectrum_table", "library_peak_table") %in% names(object@spectral_library))) {
    library <- list(
      library_spectrum_table = object@spectral_library$library_spectrum_table,
      library_peak_table = object@spectral_library$library_peak_table
    )
    if (!is.null(object@spectral_library$flash_entropy_index)) {
      library$flash_entropy_index <- object@spectral_library$flash_entropy_index
    }
    if (isTRUE(attr(object@spectral_library, "msidataset_standardized_library", exact = TRUE))) {
      attr(library, "msidataset_standardized_library") <- TRUE
    }
    for (attr_name in c("database_root", "ms1_database_dir", "ms2_database_dir")) {
      attr_value <- attr(object@spectral_library, attr_name, exact = TRUE)
      if (!is.null(attr_value)) {
        attr(library, attr_name) <- attr_value
      }
    }
  } else {
    library <- `50_read_spectral_library`(
      spectral_library = spectral_library,
      database_root = database_root,
      ms1_database_dir = ms1_database_dir,
      ms2_database_dir = ms2_database_dir,
      polarity = polarity,
      min_library_peaks = min_library_peaks,
      max_library_spectra = max_library_spectra,
      ...
    )
  }

  library <- .fe_standardize_library_once(library)
  resolved_database_root <- attr(library, "database_root", exact = TRUE)
  if (is.null(resolved_database_root) || length(resolved_database_root) == 0L) {
    resolved_database_root <- database_root
  }
  resolved_ms1_database_dir <- attr(library, "ms1_database_dir", exact = TRUE)
  if (is.null(resolved_ms1_database_dir) || length(resolved_ms1_database_dir) == 0L) {
    resolved_ms1_database_dir <- ms1_database_dir
  }
  resolved_ms2_database_dir <- attr(library, "ms2_database_dir", exact = TRUE)
  if (is.null(resolved_ms2_database_dir) || length(resolved_ms2_database_dir) == 0L) {
    resolved_ms2_database_dir <- ms2_database_dir
  }

  if (is.null(matching_library_cache)) {
    library <- `55_clean_library_peaks`(
      library = library,
      polarity = polarity,
      min_library_peaks = min_library_peaks
    )
    index <- `56_build_flash_entropy_index`(
      library = library,
      min_library_peaks = min_library_peaks
    )
  } else {
    index <- list(peak_index = matching_library_cache$fragment_search_index)
  }
  peak_index <- .fe_peak_index(index)

  matched_peaks <- `57_match_and_score_flash_entropy`(
    pseudo_ms2 = object@pseudo_ms2$pseudo_ms2_table,
    index = peak_index,
    ppm = ppm,
    mz_tol_da = mz_tol_da,
    mz_tol_da_max_mz = mz_tol_da_max_mz,
    precursor_filter = precursor_filter,
    precursor_tol = precursor_tol,
    precursor_tolerance_contract = tolerance_contract,
    progress_every = progress_every,
    progress_path = progress_path,
    progress_callback = progress_callback,
    progress_started_at = progress_started_at
  )

  spectrum_hits <- `58_aggregate_entropy_score`(
    matched_peaks = matched_peaks,
    pseudo_ms2 = object@pseudo_ms2$pseudo_ms2_table,
    index = peak_index,
    top_n_spectra = top_n_spectra,
    min_match_fe = min_match_fe,
    min_s_fe = min_s_fe
  )

  candidate_hits <- `58b_collapse_flash_entropy_candidates`(
    spectrum_hits = spectrum_hits,
    top_n_compounds = top_n_compounds,
    max_spectra_per_compound = max_spectra_per_compound
  )

  entropy_stage <- `58c_build_entropy_candidates`(
    spectrum_hits = spectrum_hits,
    matched_peaks = matched_peaks,
    top_n_entropy_candidates = top_n_entropy_candidates
  )

  object <- `59_store_flash_entropy_hits`(
    object = object,
    library = library,
    index = peak_index,
    matched_peaks = matched_peaks,
    spectrum_hits = spectrum_hits,
    candidate_hits = candidate_hits,
    entropy_candidates = entropy_stage$entropy_candidates,
    entropy_filter_summary = entropy_stage$entropy_filter_summary,
    spectrum_id_peak_index = if (is.null(matching_library_cache)) NULL else matching_library_cache$spectrum_id_peak_index,
    matching_library_cache_identity = if (is.null(matching_library_cache)) NULL else matching_library_cache$cache_identity,
    parameters = list(
      database_root = resolved_database_root,
      ms1_database_dir = resolved_ms1_database_dir,
      ms2_database_dir = resolved_ms2_database_dir,
      ppm = ppm,
      mz_tol_da = mz_tol_da,
      mz_tol_da_max_mz = mz_tol_da_max_mz,
      top_n_spectra = top_n_spectra,
      top_n_compounds = top_n_compounds,
      top_n_entropy_candidates = top_n_entropy_candidates,
      max_spectra_per_compound = max_spectra_per_compound,
      min_match_fe = min_match_fe,
      min_s_fe = min_s_fe,
      precursor_filter = precursor_filter,
      precursor_tol = precursor_tol,
      tolerance_contract = tolerance_contract,
      fragment_retrieval_tolerance = if (is.null(tolerance_contract)) {
        NA_character_
      } else {
        paste0(
          tolerance_contract$fragment_retrieval$tolerance_mode,
          ": max(", tolerance_contract$fragment_retrieval$ppm,
          " ppm, ", tolerance_contract$fragment_retrieval$low_mz_da_floor,
          " Da) below ", tolerance_contract$fragment_retrieval$low_mz_cutoff, " m/z"
        )
      },
      precursor_candidate_tolerance = if (is.null(tolerance_contract)) {
        NA_character_
      } else {
        paste0(
          tolerance_contract$precursor_candidate$mode, ": max(",
          tolerance_contract$precursor_candidate$ppm,
          " ppm, ", tolerance_contract$precursor_candidate$da_floor, " Da)"
        )
      },
      polarity = polarity,
      min_library_peaks = min_library_peaks,
      max_library_spectra = max_library_spectra,
      n_query_peaks = nrow(object@pseudo_ms2$pseudo_ms2_table),
      n_library_index_peaks = nrow(peak_index),
      n_matched_peaks = nrow(matched_peaks),
      n_spectrum_hits = nrow(spectrum_hits),
      n_raw_spectrum_hits = nrow(spectrum_hits),
      n_after_top_n_spectra = nrow(spectrum_hits),
      n_after_compound_collapse = nrow(candidate_hits),
      n_after_top_n_compounds = nrow(candidate_hits),
      n_final_entropy_candidates = nrow(entropy_stage$entropy_candidates),
      n_candidate_hits = nrow(candidate_hits),
      n_entropy_candidates = nrow(entropy_stage$entropy_candidates)
    )
  )

  object
}

`50_read_spectral_library` <- 
function(spectral_library = NULL,
                                       database_root = NULL,
                                       ms1_database_dir = NULL,
                                       ms2_database_dir = NULL,
                                       polarity = c("both", "positive", "negative"),
                                       min_library_peaks = 2,
                                       max_library_spectra = Inf,
                                       ...) {
  polarity <- match.arg(polarity)

  if (is.null(spectral_library)) {
    database_paths <- .fe_resolve_database_paths(
      database_root = database_root,
      ms1_database_dir = ms1_database_dir,
      ms2_database_dir = ms2_database_dir
    )
    db <- `51_load_local_database_files`(
      ms1_database_dir = database_paths$ms1_database_dir,
      ms2_database_dir = database_paths$ms2_database_dir
    )
    ms1_metadata <- `52_extract_ms1_metadata`(db$ms1_objects)
    ms2_library <- `53_extract_ms2_spectra`(
      db$ms2_objects,
      max_library_spectra = max_library_spectra
    )
    library <- `54_standardize_library_spectra`(
      library = ms2_library,
      ms1_metadata = ms1_metadata
    )
    attr(library, "msidataset_standardized_library") <- TRUE
    attr(library, "database_root") <- database_paths$database_root
    attr(library, "ms1_database_dir") <- database_paths$ms1_database_dir
    attr(library, "ms2_database_dir") <- database_paths$ms2_database_dir
    return(library)
  }

  if (is.character(spectral_library) && length(spectral_library) == 1L) {
    path <- spectral_library[[1]]
    if (!file.exists(path)) {
      stop("spectral_library file does not exist: ", path, call. = FALSE)
    }
    ext <- tolower(tools::file_ext(path))
    if (ext == "rds") {
      spectral_library <- readRDS(path)
    } else if (ext %in% c("rdata", "rda")) {
      env <- new.env(parent = emptyenv())
      loaded_names <- load(path, envir = env)
      spectral_library <- lapply(loaded_names, function(nm) get(nm, envir = env))
      names(spectral_library) <- loaded_names
    } else if (ext == "csv") {
      spectral_library <- utils::read.csv(path, stringsAsFactors = FALSE, check.names = FALSE)
    } else if (ext %in% c("tsv", "txt")) {
      spectral_library <- utils::read.delim(path, stringsAsFactors = FALSE, check.names = FALSE)
    } else {
      stop("Unsupported spectral_library file type: ", ext, ". Supported: .rds, .RData, .rda, .csv, .tsv.", call. = FALSE)
    }
  }

  if (is.list(spectral_library) &&
      all(c("library_spectrum_table", "library_peak_table") %in% names(spectral_library))) {
    library <- list(
      library_spectrum_table = spectral_library$library_spectrum_table,
      library_peak_table = spectral_library$library_peak_table
    )
    if (!is.null(spectral_library$flash_entropy_index)) {
      library$flash_entropy_index <- spectral_library$flash_entropy_index
    }
    if (isTRUE(attr(spectral_library, "msidataset_standardized_library", exact = TRUE))) {
      attr(library, "msidataset_standardized_library") <- TRUE
    }
    if (isTRUE(attr(spectral_library, "msidataset_hmdb_processed_long", exact = TRUE))) {
      attr(library, "msidataset_hmdb_processed_long") <- TRUE
    }
    return(library)
  }

  if (is.data.frame(spectral_library)) {
    return(`53_extract_ms2_spectra`(
      list(spectral_library = spectral_library),
      max_library_spectra = max_library_spectra
    ))
  }

  if (is.list(spectral_library)) {
    return(`53_extract_ms2_spectra`(
      spectral_library,
      max_library_spectra = max_library_spectra
    ))
  }

  stop("Unsupported spectral_library input.", call. = FALSE)
}

`50a_add_missing_columns` <- 
function(x, required_cols) {
  x <- as.data.frame(x, stringsAsFactors = FALSE, check.names = FALSE)

  for (col in required_cols) {
    if (!col %in% names(x)) {
      x[[col]] <- rep(NA, nrow(x))
    }
  }

  x
}

`50b_database_field` <- 
function(x, field) {
  if (isS4(x)) {
    return(attributes(x)[[field]])
  }
  if (is.list(x)) {
    return(x[[field]])
  }
  NULL
}

`50b_database_names` <- 
function(x) {
  if (isS4(x)) {
    return(names(attributes(x)))
  }
  if (is.list(x)) {
    return(names(x))
  }
  character()
}

`50b_is_database_class` <- 
function(x) {
  (isS4(x) || is.list(x)) &&
    all(c("database.info", "spectra.info", "spectra.data") %in%
          `50b_database_names`(x))
}

`50c_extract_database_info` <- 
function(x, object_name = NA_character_) {
  out <- list(
    database_source = NA_character_,
    library_source = NA_character_,
    database_version = NA_character_,
    source_object = object_name
  )

  if (!`50b_is_database_class`(x)) {
    return(out)
  }

  info <- `50b_database_field`(x, "database.info")

  if (is.list(info)) {
    if ("Source" %in% names(info)) {
      out$database_source <- as.character(info$Source[[1]])
      out$library_source <- as.character(info$Source[[1]])
    }

    if ("Version" %in% names(info)) {
      out$database_version <- as.character(info$Version[[1]])
    }
  }

  out
}

`51_load_local_database_files` <- 
function(ms1_database_dir, ms2_database_dir) {
  load_dir <- function(database_dir) {
    database_dir <- .fe_resolve_database_dir(database_dir)
    if (is.null(database_dir) || !dir.exists(database_dir)) {
      stop("Database directory does not exist: ", database_dir, call. = FALSE)
    }

    files <- list.files(
      database_dir,
      pattern = "\\.(RData|rda)$",
      full.names = TRUE,
      ignore.case = TRUE
    )

    if (length(files) == 0) {
      stop("No .RData or .rda files found in: ", database_dir, call. = FALSE)
    }

    objects <- list()

    for (file in files) {
      env <- new.env(parent = emptyenv())
      loaded_names <- load(file, envir = env)

      for (nm in loaded_names) {
        key <- paste(basename(file), nm, sep = "::")
        objects[[key]] <- get(nm, envir = env)
      }
    }

    objects
  }

  list(
    ms1_objects = load_dir(ms1_database_dir),
    ms2_objects = load_dir(ms2_database_dir)
  )
}

`52_extract_ms1_metadata` <- 
function(ms1_objects) {
  tables <- list()
  if (is.null(names(ms1_objects))) {
    names(ms1_objects) <- paste0("MS1_OBJECT_", seq_along(ms1_objects))
  }
  for (nm in names(ms1_objects)) {
    x <- ms1_objects[[nm]]
    if (`50b_is_database_class`(x)) {
      tbl <- `52b_extract_ms1_from_database_class`(x, object_name = nm)
      if (nrow(tbl) > 0L) {
        tables[[length(tables) + 1L]] <- tbl
      }
    } else {
      collected <- .fe_collect_tables(x)
      tables <- c(tables, collected)
    }
  }
  if (length(tables) == 0L) {
    return(.fe_empty_ms1_metadata())
  }

  out <- lapply(seq_along(tables), function(i) {
    df <- .fe_as_data_frame(tables[[i]])
    df <- .fe_standardize_columns(df, list(
      candidate_id = c("HMDB.ID", "KEGG.ID", "PUBCHEM.ID", "CHEBI.ID", "METLIN.ID", "Lab.ID", "candidate_id", "compound_id", "hmdb_id", "id", "ID"),
      compound_id = c("HMDB.ID", "KEGG.ID", "PUBCHEM.ID", "CHEBI.ID", "METLIN.ID", "Lab.ID", "compound_id", "candidate_id", "hmdb_id", "id", "ID"),
      candidate_name = c("Compound.name", "candidate_name", "compound_name", "name", "Name", "metabolite_name", "IUPAC_name", "Traditional_IUPAC_name"),
      formula = c("formula", "molecular_formula", "Formula"),
      exact_mass = c("monisotopic_molecular_weight", "exact_mass", "ExactMass", "mz", "MZ", "mass", "exactmass"),
      precursor_mz = c("mz", "mz.pos", "mz.neg", "precursor_mz", "precursorMz", "Precursor_mz", "precursor"),
      adduct = c("adduct", "Adduct", "adduct_type", "ion_form"),
      polarity = c("polarity", "Polarity", "ion_mode", "mode"),
      inchikey = c("INCHIKEY.ID", "InChIKey", "inchikey"),
      smiles = c("SMILES.ID", "SMILES", "smiles"),
      database_source = c("database_source", "Source", "Submitter", "source", "library_source"),
      library_source = c("library_source", "database_source", "Submitter", "source")
    ))
    .fe_add_missing_columns(df, colnames(.fe_empty_ms1_metadata()))
  })
  ms1_cols <- colnames(.fe_empty_ms1_metadata())
  out <- do.call(rbind, lapply(out, function(df) df[, ms1_cols, drop = FALSE]))

  missing_candidate <- is.na(out$candidate_id) | !nzchar(as.character(out$candidate_id))
  if (any(missing_candidate)) {
    compound <- as.character(out$compound_id)
    use_compound <- missing_candidate & !is.na(compound) & nzchar(compound)
    out$candidate_id[use_compound] <- compound[use_compound]
    missing_candidate <- is.na(out$candidate_id) | !nzchar(as.character(out$candidate_id))
    if (any(missing_candidate)) {
      source <- as.character(out$database_source)
      source[is.na(source) | !nzchar(source)] <- as.character(out$source_object[is.na(source) | !nzchar(source)])
      source[is.na(source) | !nzchar(source)] <- "MS1"
      out$candidate_id[missing_candidate] <- paste0(source[missing_candidate], "_", which(missing_candidate))
    }
  }
  missing_compound <- is.na(out$compound_id) | !nzchar(as.character(out$compound_id))
  if (any(missing_compound)) {
    out$compound_id[missing_compound] <- as.character(out$candidate_id[missing_compound])
  }

  out$candidate_id <- as.character(out$candidate_id)
  out$compound_id <- as.character(out$compound_id)
  out$exact_mass <- suppressWarnings(as.numeric(out$exact_mass))
  out$precursor_mz <- suppressWarnings(as.numeric(out$precursor_mz))
  out <- out[!duplicated(out$candidate_id), colnames(.fe_empty_ms1_metadata()), drop = FALSE]
  rownames(out) <- NULL
  out
}

`52b_extract_ms1_from_database_class` <- 
function(x, object_name = "databaseClass") {
  if (!`50b_is_database_class`(x)) {
    return(data.frame())
  }

  spectra_info <- `50b_database_field`(x, "spectra.info")

  if (!is.data.frame(spectra_info)) {
    return(data.frame())
  }

  db_info <- `50c_extract_database_info`(x, object_name = object_name)

  out <- as.data.frame(spectra_info, stringsAsFactors = FALSE, check.names = FALSE)

  out$database_source <- db_info$database_source
  out$library_source <- db_info$library_source
  out$database_version <- db_info$database_version
  out$source_object <- db_info$source_object

  out
}

`53_extract_ms2_spectra` <- 
function(ms2_objects, max_library_spectra = Inf) {
  spectrum_tables <- list()
  peak_tables <- list()

  add_standard_library <- function(x) {
    spectrum_tables[[length(spectrum_tables) + 1L]] <<- x$library_spectrum_table
    peak_tables[[length(peak_tables) + 1L]] <<- x$library_peak_table
  }

  process_object <- function(x, object_name = "object") {
    if (`50b_is_database_class`(x)) {
      remaining <- max_library_spectra
      if (is.finite(max_library_spectra) && length(peak_tables) > 0L) {
        parsed <- length(unique(unlist(lapply(peak_tables, function(tbl) as.character(tbl$spectrum_id)))))
        remaining <- max(max_library_spectra - parsed, 0)
      }
      extracted <- `53a_extract_ms2_from_database_class`(
        x,
        object_name = object_name,
        max_library_spectra = remaining
      )
      spectrum_tables[[length(spectrum_tables) + 1L]] <<- extracted$library_spectrum_table
      peak_tables[[length(peak_tables) + 1L]] <<- extracted$library_peak_table
      return(invisible(TRUE))
    }
    if (is.data.frame(x)) {
      parsed <- .fe_ms2_from_long_table(x, object_name = object_name)
      spectrum_tables[[length(spectrum_tables) + 1L]] <<- parsed$library_spectrum_table
      peak_tables[[length(peak_tables) + 1L]] <<- parsed$library_peak_table
      return(invisible(TRUE))
    }
    if (is.list(x) && all(c("library_spectrum_table", "library_peak_table") %in% names(x))) {
      add_standard_library(x)
      return(invisible(TRUE))
    }
    if (is.list(x) && all(c("spectrum_table", "peak_table") %in% names(x))) {
      spectrum_tables[[length(spectrum_tables) + 1L]] <<- x$spectrum_table
      peak_tables[[length(peak_tables) + 1L]] <<- x$peak_table
      return(invisible(TRUE))
    }
    if (is.list(x)) {
      if (!is.null(x$mz) && !is.null(x$intensity)) {
        parsed <- .fe_ms2_from_spectrum_list(x, object_name)
        spectrum_tables[[length(spectrum_tables) + 1L]] <<- parsed$library_spectrum_table
        peak_tables[[length(peak_tables) + 1L]] <<- parsed$library_peak_table
        return(invisible(TRUE))
      }
      if (!is.null(x$peaks)) {
        parsed <- .fe_ms2_from_spectrum_list(x, object_name)
        spectrum_tables[[length(spectrum_tables) + 1L]] <<- parsed$library_spectrum_table
        peak_tables[[length(peak_tables) + 1L]] <<- parsed$library_peak_table
        return(invisible(TRUE))
      }
      child_names <- names(x)
      if (is.null(child_names)) {
        child_names <- paste0("ITEM_", seq_along(x))
      }
      for (i in seq_along(x)) {
        process_object(x[[i]], object_name = paste(object_name, child_names[[i]], sep = "::"))
      }
    }
    invisible(TRUE)
  }

  if (is.null(names(ms2_objects))) {
    names(ms2_objects) <- paste0("MS2_OBJECT_", seq_along(ms2_objects))
  }
  for (nm in names(ms2_objects)) {
    process_object(ms2_objects[[nm]], object_name = nm)
  }

  spectrum_table <- if (length(spectrum_tables) > 0L) {
    do.call(rbind, lapply(spectrum_tables, .fe_standardize_spectrum_table))
  } else {
    .fe_empty_spectrum_table()
  }
  peak_table <- if (length(peak_tables) > 0L) {
    do.call(rbind, lapply(peak_tables, .fe_standardize_peak_table))
  } else {
    .fe_empty_peak_table(raw_only = TRUE)
  }

  if (nrow(peak_table) == 0L) {
    stop(
      "No MS2 peaks could be extracted from ms2_database. databaseClass@spectra.data was recognized but no peak_mz/intensity table could be parsed.",
      call. = FALSE
    )
  }

  list(
    library_spectrum_table = spectrum_table,
    library_peak_table = peak_table
  )
}

`53a_extract_ms2_from_database_class` <- 
function(x,
                                                  object_name = "databaseClass",
                                                  max_library_spectra = Inf) {
  if (!`50b_is_database_class`(x)) {
    return(list(
      library_spectrum_table = data.frame(),
      library_peak_table = data.frame()
    ))
  }

  db_info <- `50c_extract_database_info`(x, object_name = object_name)

  spectra_info <- `50b_database_field`(x, "spectra.info")
  spectra_data <- `50b_database_field`(x, "spectra.data")

  if (!is.data.frame(spectra_info)) {
    spectrum_table <- data.frame()
  } else {
    spectrum_table <- as.data.frame(spectra_info, stringsAsFactors = FALSE, check.names = FALSE)
  }

  spectrum_table <- .fe_standardize_columns(spectrum_table, list(
    spectrum_id = c("Lab.ID", "GNPS.ID", "HMDB.ID", "METLIN.ID", "candidate_id", "spectrum_id", "SpectrumID", "id", "ID"),
    candidate_id = c("HMDB.ID", "GNPS.ID", "KEGG.ID", "PUBCHEM.ID", "CHEBI.ID", "METLIN.ID", "Lab.ID", "candidate_id", "compound_id"),
    compound_id = c("HMDB.ID", "GNPS.ID", "KEGG.ID", "PUBCHEM.ID", "CHEBI.ID", "METLIN.ID", "Lab.ID", "compound_id"),
    candidate_name = c("Compound.name", "candidate_name", "compound_name", "Name", "name", "IUPAC_name", "Traditional_IUPAC_name"),
    precursor_mz = c("Precursor_mz", "precursor_mz"),
    exact_mass = c("monisotopic_molecular_weight", "exact_mass", "mz"),
    formula = c("Formula", "formula"),
    polarity = c("Polarity", "polarity"),
    adduct = c("Adduct", "adduct"),
    collision_energy = c("CE", "collision_energy", "CollisionEnergy"),
    instrument_type = c("Instrument_type", "instrument_type", "Instrument", "instrument"),
    inchikey = c("INCHIKEY.ID", "InChIKey", "inchikey"),
    smiles = c("SMILES.ID", "SMILES", "smiles"),
    library_source = c("library_source", "Submitter", "Source"),
    database_source = c("database_source", "Submitter", "Source")
  ))

  precursor_value <- suppressWarnings(as.numeric(spectrum_table$precursor_mz))
  spectrum_table$precursor_source <- ifelse(
    is.finite(precursor_value) & precursor_value > 0,
    "EXPLICIT_LIBRARY_PRECURSOR",
    "MISSING_LIBRARY_PRECURSOR"
  )

  if (nrow(spectrum_table) > 0L) {
    spectrum_table <- .fe_add_missing_columns(
      spectrum_table,
      c("library_source", "database_source")
    )
    spectrum_table$library_source <- .fe_fill_missing_scalar(
      spectrum_table$library_source,
      db_info$library_source
    )
    spectrum_table$database_source <- .fe_fill_missing_scalar(
      spectrum_table$database_source,
      db_info$database_source
    )
    spectrum_table$database_version <- db_info$database_version
    spectrum_table$source_object <- db_info$source_object
    if (!"spectrum_id" %in% names(spectrum_table)) {
      spectrum_table$spectrum_id <- rep(NA_character_, nrow(spectrum_table))
    }
    missing_sid <- is.na(spectrum_table$spectrum_id) |
      !nzchar(as.character(spectrum_table$spectrum_id))
    source <- as.character(spectrum_table$database_source)
    source[is.na(source) | !nzchar(source)] <- "DATABASE"
    spectrum_table$spectrum_id[missing_sid] <- paste0(
      source[missing_sid],
      "_SPECTRUM_",
      which(missing_sid)
    )
  }

  spectrum_table <- .fe_add_missing_columns(
    spectrum_table,
    c(
      colnames(.fe_empty_spectrum_table()),
      "Lab.ID", "GNPS.ID", "HMDB.ID", "KEGG.ID", "PUBCHEM.ID",
      "CHEBI.ID", "METLIN.ID", "Ontology", "Kingdom", "Super_class",
      "Class", "Sub_class", "database_version", "source_object"
    )
  )

  spectrum_rows <- list()
  peak_rows <- list()
  parsed_spectra <- 0L
  max_library_spectra <- suppressWarnings(as.numeric(max_library_spectra))
  if (length(max_library_spectra) == 0L || is.na(max_library_spectra)) {
    max_library_spectra <- Inf
  }

  if (is.list(spectra_data) && length(spectra_data) > 0L) {
    group_names <- names(spectra_data)
    if (is.null(group_names)) {
      group_names <- paste0("Spectra.", seq_along(spectra_data))
    }
    metadata_index <- `53b_build_spectrum_metadata_index`(spectrum_table)

    stop_parsing <- FALSE
    for (group_i in seq_along(spectra_data)) {
      group <- spectra_data[[group_i]]
      group_name <- group_names[[group_i]]
      group_polarity <- NA_character_
      if (grepl("positive", group_name, ignore.case = TRUE)) {
        group_polarity <- "positive"
      } else if (grepl("negative", group_name, ignore.case = TRUE)) {
        group_polarity <- "negative"
      }
      if (!is.list(group) || length(group) == 0L) {
        next
      }
      spectrum_ids <- names(group)
      if (is.null(spectrum_ids)) {
        spectrum_ids <- paste0(group_name, "_", seq_along(group))
      }
      for (spectrum_i in seq_along(group)) {
        if (is.finite(max_library_spectra) && parsed_spectra >= max_library_spectra) {
          stop_parsing <- TRUE
          break
        }
        spectrum_id <- spectrum_ids[[spectrum_i]]
        peak_table <- `53c_extract_peak_frame`(group[[spectrum_i]])
        if (nrow(peak_table) == 0L) {
          next
        }
        metadata <- `53b_match_spectrum_metadata`(spectrum_id, spectrum_table, metadata_index)
        metadata <- .fe_add_missing_columns(metadata, names(spectrum_table))
        metadata$library_source <- .fe_fill_missing_scalar(
          metadata$library_source,
          db_info$library_source
        )
        metadata$database_source <- .fe_fill_missing_scalar(
          metadata$database_source,
          db_info$database_source
        )
        if ("database_version" %in% names(metadata)) {
          metadata$database_version <- .fe_fill_missing_scalar(
            metadata$database_version,
            db_info$database_version
          )
        }
        if ("source_object" %in% names(metadata)) {
          metadata$source_object <- .fe_fill_missing_scalar(
            metadata$source_object,
            db_info$source_object
          )
        }
        missing_pol <- is.na(metadata$polarity) | !nzchar(as.character(metadata$polarity))
        if (any(missing_pol) && !is.na(group_polarity)) {
          metadata$polarity[missing_pol] <- group_polarity
        }
        missing_sid <- is.na(metadata$spectrum_id) |
          !nzchar(as.character(metadata$spectrum_id))
        metadata$spectrum_id[missing_sid] <- as.character(spectrum_id)
        missing_cid <- is.na(metadata$candidate_id) |
          !nzchar(as.character(metadata$candidate_id))
        metadata$candidate_id[missing_cid] <- metadata$spectrum_id[missing_cid]

        peak_table$spectrum_id <- as.character(metadata$spectrum_id[[1]])
        peak_table$candidate_id <- as.character(metadata$candidate_id[[1]])

        spectrum_rows[[length(spectrum_rows) + 1L]] <- metadata
        peak_rows[[length(peak_rows) + 1L]] <- peak_table[, c("spectrum_id", "candidate_id", "peak_mz", "intensity_raw"), drop = FALSE]
        parsed_spectra <- parsed_spectra + 1L
      }
      if (isTRUE(stop_parsing)) {
        break
      }
    }
  }

  if (length(spectrum_rows) > 0L) {
    spectrum_out <- .fe_rbind_fill(spectrum_rows)
  } else {
    spectrum_out <- spectrum_table[0, , drop = FALSE]
  }
  if (length(peak_rows) > 0L) {
    peak_out <- .fe_rbind_fill(peak_rows)
  } else {
    peak_out <- .fe_empty_peak_table(raw_only = TRUE)
  }

  list(
    library_spectrum_table = spectrum_out,
    library_peak_table = peak_out
  )
}

`53b_build_spectrum_metadata_index` <- 
function(spectrum_table) {
  spectrum_table <- as.data.frame(spectrum_table, stringsAsFactors = FALSE, check.names = FALSE)
  index <- new.env(parent = emptyenv(), hash = TRUE)
  if (!nrow(spectrum_table)) {
    return(index)
  }
  for (col in c("spectrum_id", "Lab.ID", "GNPS.ID", "HMDB.ID", "candidate_id")) {
    if (!col %in% names(spectrum_table)) {
      next
    }
    values <- as.character(spectrum_table[[col]])
    keep <- !is.na(values) & nzchar(values)
    if (!any(keep)) {
      next
    }
    rows <- which(keep)
    values <- values[keep]
    first <- !duplicated(values)
    rows <- rows[first]
    values <- values[first]
    keys <- paste(col, values, sep = "\r")
    for (i in seq_along(keys)) {
      index[[keys[[i]]]] <- rows[[i]]
    }
  }
  index
}

`53b_match_spectrum_metadata` <- 
function(spectrum_id, spectrum_table, metadata_index = NULL) {
  spectrum_table <- as.data.frame(spectrum_table, stringsAsFactors = FALSE, check.names = FALSE)
  if (nrow(spectrum_table) > 0L) {
    for (col in c("spectrum_id", "Lab.ID", "GNPS.ID", "HMDB.ID", "candidate_id")) {
      if (col %in% names(spectrum_table)) {
        if (!is.null(metadata_index)) {
          key <- paste(col, as.character(spectrum_id), sep = "\r")
          idx <- metadata_index[[key]]
          if (!is.null(idx)) {
            return(spectrum_table[idx[[1]], , drop = FALSE])
          }
          next
        }
        idx <- which(as.character(spectrum_table[[col]]) == as.character(spectrum_id))
        if (length(idx) > 0L) {
          return(spectrum_table[idx[[1]], , drop = FALSE])
        }
      }
    }
  }

  out <- .fe_add_missing_columns(data.frame(), names(spectrum_table))
  out[1, "spectrum_id"] <- as.character(spectrum_id)
  if ("candidate_id" %in% names(out)) {
    out[1, "candidate_id"] <- as.character(spectrum_id)
  }
  out
}

`53c_extract_peak_frame` <- 
function(peak_object) {
  empty <- data.frame(
    peak_mz = numeric(),
    intensity_raw = numeric(),
    stringsAsFactors = FALSE,
    check.names = FALSE
  )

  extract_rectangular <- function(x) {
    df <- as.data.frame(x, stringsAsFactors = FALSE, check.names = FALSE)
    if (ncol(df) < 2L) {
      return(empty)
    }
    mz_names <- c("peak_mz", "mz", "MZ", "PeakMZ", "fragment_mz", "product_mz")
    int_names <- c("intensity_raw", "intensity", "Intensity", "IntensityRaw", "peak_intensity", "abundance")
    mz_col <- mz_names[mz_names %in% names(df)]
    int_col <- int_names[int_names %in% names(df)]
    if (length(mz_col) == 0L || length(int_col) == 0L) {
      mz_col <- names(df)[[1]]
      int_col <- names(df)[[2]]
    } else {
      mz_col <- mz_col[[1]]
      int_col <- int_col[[1]]
    }
    data.frame(
      peak_mz = suppressWarnings(as.numeric(df[[mz_col]])),
      intensity_raw = suppressWarnings(as.numeric(df[[int_col]])),
      stringsAsFactors = FALSE,
      check.names = FALSE
    )
  }

  if (is.data.frame(peak_object) || is.matrix(peak_object)) {
    return(extract_rectangular(peak_object))
  }

  if (is.list(peak_object)) {
    if (!is.null(peak_object$mz) && !is.null(peak_object$intensity)) {
      return(data.frame(
        peak_mz = suppressWarnings(as.numeric(peak_object$mz)),
        intensity_raw = suppressWarnings(as.numeric(peak_object$intensity)),
        stringsAsFactors = FALSE,
        check.names = FALSE
      ))
    }
    if (!is.null(peak_object$MZ) && !is.null(peak_object$Intensity)) {
      return(data.frame(
        peak_mz = suppressWarnings(as.numeric(peak_object$MZ)),
        intensity_raw = suppressWarnings(as.numeric(peak_object$Intensity)),
        stringsAsFactors = FALSE,
        check.names = FALSE
      ))
    }
    pieces <- lapply(peak_object, `53c_extract_peak_frame`)
    pieces <- pieces[vapply(pieces, nrow, integer(1)) > 0L]
    if (length(pieces) > 0L) {
      return(do.call(rbind, pieces))
    }
  }

  empty
}

`54_standardize_library_spectra` <- 
function(library, ms1_metadata = NULL) {
  if (!is.list(library) || !all(c("library_spectrum_table", "library_peak_table") %in% names(library))) {
    stop("library must contain library_spectrum_table and library_peak_table.", call. = FALSE)
  }

  spectrum <- .fe_standardize_spectrum_table(library$library_spectrum_table)
  peaks <- .fe_standardize_peak_table(library$library_peak_table)

  if (!is.null(ms1_metadata) && nrow(ms1_metadata) > 0L) {
    ms1 <- `52_extract_ms1_metadata`(list(ms1_metadata))
    spectrum <- .fe_fill_spectrum_from_ms1(spectrum, ms1)
  }

  if (!"candidate_id" %in% names(peaks) || all(is.na(peaks$candidate_id))) {
    candidate_map <- stats::setNames(as.character(spectrum$candidate_id), as.character(spectrum$spectrum_id))
    peaks$candidate_id <- unname(candidate_map[as.character(peaks$spectrum_id)])
  } else {
    missing_candidate <- is.na(peaks$candidate_id) | !nzchar(as.character(peaks$candidate_id))
    if (any(missing_candidate)) {
      candidate_map <- stats::setNames(as.character(spectrum$candidate_id), as.character(spectrum$spectrum_id))
      peaks$candidate_id[missing_candidate] <- unname(candidate_map[as.character(peaks$spectrum_id[missing_candidate])])
    }
  }

  peaks <- peaks[!duplicated(.fe_split_key(peaks, c("spectrum_id", "peak_mz", "intensity_raw"))), , drop = FALSE]
  spectrum <- spectrum[, colnames(.fe_empty_spectrum_table()), drop = FALSE]
  peaks <- .fe_add_missing_columns(peaks, colnames(.fe_empty_peak_table()))
  peaks <- peaks[, colnames(.fe_empty_peak_table()), drop = FALSE]
  rownames(spectrum) <- NULL
  rownames(peaks) <- NULL
  out <- list(library_spectrum_table = spectrum, library_peak_table = peaks)
  attr(out, "msidataset_standardized_library") <- TRUE
  out
}

`55_clean_library_peaks` <- 
function(library,
                                     polarity = c("both", "positive", "negative"),
                                     min_library_peaks = 2) {
  polarity <- match.arg(polarity)
  min_library_peaks <- as.integer(min_library_peaks)
  if (!is.finite(min_library_peaks) || min_library_peaks < 1L) {
    stop("min_library_peaks must be a positive integer.", call. = FALSE)
  }

  spectrum <- .fe_standardize_spectrum_table(library$library_spectrum_table)
  peaks <- .fe_standardize_peak_table(library$library_peak_table)

  if (!identical(polarity, "both") && "polarity" %in% colnames(spectrum)) {
    keep_spec <- is.na(spectrum$polarity) |
      !nzchar(as.character(spectrum$polarity)) |
      tolower(as.character(spectrum$polarity)) == polarity
    spectrum <- spectrum[keep_spec, , drop = FALSE]
    peaks <- peaks[as.character(peaks$spectrum_id) %in% as.character(spectrum$spectrum_id), , drop = FALSE]
  }

  peaks$peak_mz <- suppressWarnings(as.numeric(peaks$peak_mz))
  peaks$intensity_raw <- suppressWarnings(as.numeric(peaks$intensity_raw))
  peaks <- peaks[
    is.finite(peaks$peak_mz) &
      is.finite(peaks$intensity_raw) &
      peaks$intensity_raw > 0,
    ,
    drop = FALSE
  ]

  if (nrow(peaks) == 0L) {
    stop(
      "No valid library peaks remain after cleaning. Check ms2_database parsing, especially databaseClass@spectra.data.",
      call. = FALSE
    )
  }

  split_rows <- split(seq_len(nrow(peaks)), as.character(peaks$spectrum_id))
  peaks$intensity_basepeak100 <- NA_real_
  peaks$intensity_sum1 <- NA_real_
  valid_ids <- character(length(split_rows))
  valid_count <- 0L
  for (sid in names(split_rows)) {
    idx <- split_rows[[sid]]
    ints <- as.numeric(peaks$intensity_raw[idx])
    if (length(idx) < min_library_peaks) {
      next
    }
    base <- max(ints, na.rm = TRUE)
    total <- sum(ints, na.rm = TRUE)
    if (!is.finite(base) || base <= 0 || !is.finite(total) || total <= 0) {
      next
    }
    peaks$intensity_basepeak100[idx] <- ints / base * 100
    peaks$intensity_sum1[idx] <- ints / total
    valid_count <- valid_count + 1L
    valid_ids[[valid_count]] <- sid
  }
  valid_ids <- valid_ids[seq_len(valid_count)]

  peaks <- peaks[as.character(peaks$spectrum_id) %in% valid_ids, , drop = FALSE]
  spectrum <- spectrum[as.character(spectrum$spectrum_id) %in% valid_ids, , drop = FALSE]
  if (nrow(peaks) == 0L) {
    stop(
      "No valid library peaks remain after cleaning. Check ms2_database parsing, especially databaseClass@spectra.data.",
      call. = FALSE
    )
  }
  rownames(peaks) <- NULL
  rownames(spectrum) <- NULL
  out <- list(library_spectrum_table = spectrum, library_peak_table = peaks)
  if (is.list(library) && !is.null(library$flash_entropy_index)) {
    fast_index <- as.data.frame(library$flash_entropy_index, stringsAsFactors = FALSE, check.names = FALSE)
    if ("spectrum_id" %in% names(fast_index)) {
      fast_index <- fast_index[as.character(fast_index$spectrum_id) %in% as.character(spectrum$spectrum_id), , drop = FALSE]
    }
    out$flash_entropy_index <- fast_index
  }
  if (isTRUE(attr(library, "msidataset_standardized_library", exact = TRUE))) {
    attr(out, "msidataset_standardized_library") <- TRUE
  }
  if (isTRUE(attr(library, "msidataset_hmdb_processed_long", exact = TRUE))) {
    attr(out, "msidataset_hmdb_processed_long") <- TRUE
  }
  out
}

`56_build_flash_entropy_index` <- 
function(library, min_library_peaks = 2) {
  if (is.list(library) && !is.null(library$flash_entropy_index)) {
    fast_index <- as.data.frame(library$flash_entropy_index, stringsAsFactors = FALSE, check.names = FALSE)
    required <- colnames(.fe_empty_flash_entropy_index())
    if (all(required %in% names(fast_index))) {
      fast_index <- fast_index[, required, drop = FALSE]
      fast_index$library_mz <- suppressWarnings(as.numeric(fast_index$library_mz))
      fast_index$library_intensity <- suppressWarnings(as.numeric(fast_index$library_intensity))
      fast_index <- fast_index[
        is.finite(fast_index$library_mz) &
          is.finite(fast_index$library_intensity) &
          fast_index$library_intensity > 0,
        ,
        drop = FALSE
      ]
      fast_index <- fast_index[order(fast_index$library_mz, fast_index$spectrum_id), , drop = FALSE]
      rownames(fast_index) <- NULL
      message("Using prebuilt Flash Entropy peak index with ", nrow(fast_index), " peak rows.")
      return(list(peak_index = fast_index))
    }
  }
  peak_table <- as.data.frame(library$library_peak_table, stringsAsFactors = FALSE, check.names = FALSE)
  spectrum_table <- as.data.frame(library$library_spectrum_table, stringsAsFactors = FALSE, check.names = FALSE)
  peak_table <- .fe_standardize_peak_table(peak_table)
  spectrum_table <- .fe_standardize_spectrum_table(spectrum_table)
  if (!"intensity_sum1" %in% names(peak_table)) {
    peak_table$intensity_sum1 <- NA_real_
  }
  peak_table$peak_mz <- suppressWarnings(as.numeric(peak_table$peak_mz))
  peak_table$intensity_raw <- suppressWarnings(as.numeric(peak_table$intensity_raw))
  peak_table$intensity_sum1 <- suppressWarnings(as.numeric(peak_table$intensity_sum1))
  peak_table <- peak_table[
    is.finite(peak_table$peak_mz) &
      is.finite(peak_table$intensity_raw) &
      peak_table$intensity_raw > 0,
    ,
    drop = FALSE
  ]
  if (nrow(peak_table) == 0L) {
    return(.fe_empty_flash_entropy_index_result())
  }

  split_idx <- split(seq_len(nrow(peak_table)), as.character(peak_table$spectrum_id))
  valid <- logical(nrow(peak_table))
  for (sid in names(split_idx)) {
    idx <- split_idx[[sid]]
    if (length(idx) < as.integer(min_library_peaks)) {
      next
    }
    total <- sum(peak_table$intensity_raw[idx], na.rm = TRUE)
    if (!is.finite(total) || total <= 0) {
      next
    }
    peak_table$intensity_sum1[idx] <- peak_table$intensity_raw[idx] / total
    valid[idx] <- TRUE
  }
  peak_table <- peak_table[valid, , drop = FALSE]
  if (nrow(peak_table) == 0L) {
    return(.fe_empty_flash_entropy_index_result())
  }

  meta_cols <- c(
    "spectrum_id", "candidate_id", "compound_id", "candidate_name",
    "library_source", "database_source", "precursor_mz", "formula",
    "adduct", "polarity", "inchikey", "smiles"
  )
  spectrum_table <- .fe_add_missing_columns(spectrum_table, meta_cols)
  peak_table <- merge(
    peak_table,
    spectrum_table[, meta_cols, drop = FALSE],
    by = c("spectrum_id", "candidate_id"),
    all.x = TRUE,
    sort = FALSE
  )
  peak_table$library_mz <- suppressWarnings(as.numeric(peak_table$peak_mz))
  peak_table$library_intensity <- suppressWarnings(as.numeric(peak_table$intensity_sum1))
  peak_table <- peak_table[
    is.finite(peak_table$library_mz) &
      is.finite(peak_table$library_intensity) &
      peak_table$library_intensity > 0,
    ,
    drop = FALSE
  ]
  peak_table <- peak_table[order(peak_table$library_mz, peak_table$spectrum_id), , drop = FALSE]
  peak_table$library_fragment_id <- paste0(as.character(peak_table$spectrum_id), "_", seq_len(nrow(peak_table)))
  peak_table$library_peak_id <- peak_table$library_fragment_id
  index <- peak_table[, colnames(.fe_empty_flash_entropy_index()), drop = FALSE]
  rownames(index) <- NULL
  list(peak_index = index)
}

`57_match_and_score_flash_entropy` <- 
function(pseudo_ms2,
                                               index,
                                               ppm = 10,
                                               mz_tol_da = NULL,
                                               mz_tol_da_max_mz = Inf,
                                               precursor_filter = FALSE,
                                               precursor_tol = 0.05,
                                               progress_every = NULL,
                                               progress_path = NULL,
                                               progress_callback = NULL,
                                               progress_started_at = NULL,
                                               precursor_tolerance_contract = NULL) {
  empty <- .fe_empty_matched_peaks()
  progress_started_at <- if (inherits(progress_started_at, "POSIXt")) progress_started_at else Sys.time()
  query <- as.data.frame(pseudo_ms2, stringsAsFactors = FALSE, check.names = FALSE)
  peak_index <- .fe_peak_index(index)

  q_required <- c("peak_id", "anchor_feature_id", "anchor_mz", "fragment_id", "fragment_feature_id", "fragment_mz", "SAS", "intensity")
  missing_q <- setdiff(q_required, names(query))
  if (length(missing_q) > 0L) {
    stop("pseudo_ms2 missing required columns: ", paste(missing_q, collapse = ", "), call. = FALSE)
  }
  if (nrow(query) == 0L || nrow(peak_index) == 0L) {
    return(empty)
  }
  query_key <- .fe_split_key(query, c("peak_id", "anchor_feature_id"))
  total_queries <- length(unique(query_key))
  progress_every <- if (is.null(progress_every) || length(progress_every) == 0L) {
    NA_integer_
  } else {
    suppressWarnings(as.integer(progress_every[[1]]))
  }
  progress_enabled <- is.finite(progress_every) && progress_every > 0L
  processed_queries <- 0L
  last_progress <- 0L
  seen_queries <- new.env(parent = emptyenv(), hash = TRUE)

  query$fragment_mz <- suppressWarnings(as.numeric(query$fragment_mz))
  query$intensity <- suppressWarnings(as.numeric(query$intensity))
  query$anchor_mz <- suppressWarnings(as.numeric(query$anchor_mz))
  query$SAS <- suppressWarnings(as.numeric(query$SAS))
  peak_index$library_mz <- suppressWarnings(as.numeric(peak_index$library_mz))
  peak_index$library_intensity <- suppressWarnings(as.numeric(peak_index$library_intensity))
  peak_index$precursor_mz <- suppressWarnings(as.numeric(peak_index$precursor_mz))
  valid_library <- is.finite(peak_index$library_mz) &
    peak_index$library_mz > 0 &
    is.finite(peak_index$library_intensity) &
    peak_index$library_intensity > 0
  if (!any(valid_library)) {
    warning(
      "Flash Entropy library index has no peaks with positive finite library_mz and library_intensity.",
      call. = FALSE
    )
    return(empty)
  }
  invalid_library_n <- sum(!valid_library)
  peak_index <- peak_index[valid_library, , drop = FALSE]
  peak_index <- peak_index[order(peak_index$library_mz, peak_index$spectrum_id), , drop = FALSE]
  rownames(peak_index) <- NULL
  library_mz <- as.numeric(peak_index$library_mz)
  if (invalid_library_n > 0L) {
    warning(
      "Dropped ", invalid_library_n,
      " Flash Entropy library index peak(s) with non-positive or non-finite m/z/intensity.",
      call. = FALSE
    )
  }

  rows <- vector("list", 0L)
  row_i <- 0L
  ppm <- as.numeric(ppm[[1]])
  if (!is.finite(ppm) || ppm < 0) {
    stop("ppm must be a non-negative finite numeric value.", call. = FALSE)
  }
  mz_tol <- if (!is.null(mz_tol_da)) as.numeric(mz_tol_da[[1]]) else NA_real_
  if (!is.null(mz_tol_da) && (!is.finite(mz_tol) || mz_tol < 0)) {
    stop("mz_tol_da must be a non-negative finite numeric value.", call. = FALSE)
  }
  mz_tol_da_max_mz <- as.numeric(mz_tol_da_max_mz[[1]])
  if (!is.finite(mz_tol_da_max_mz)) {
    mz_tol_da_max_mz <- Inf
  }
  if (!is.infinite(mz_tol_da_max_mz) && mz_tol_da_max_mz <= 0) {
    stop("mz_tol_da_max_mz must be a positive numeric value or Inf.", call. = FALSE)
  }
  for (i in seq_len(nrow(query))) {
    q_key <- as.character(query_key[[i]])
    if (!exists(q_key, envir = seen_queries, inherits = FALSE)) {
      assign(q_key, TRUE, envir = seen_queries)
      processed_queries <- processed_queries + 1L
      if (isTRUE(progress_enabled) &&
          (processed_queries %% progress_every == 0L || processed_queries == total_queries)) {
        .fe_emit_progress(
          progress_callback = progress_callback,
          progress_path = progress_path,
          started_at = progress_started_at,
          processed_queries = processed_queries,
          total_queries = total_queries,
          candidates_so_far = row_i,
          message = "matching query peaks"
        )
        last_progress <- processed_queries
      }
    }
    q_mz <- as.numeric(query$fragment_mz[[i]])
    q_int <- as.numeric(query$intensity[[i]])
    if (!is.finite(q_mz) || q_mz <= 0 || !is.finite(q_int) || q_int <= 0) {
      next
    }
	    ppm_da <- q_mz * ppm / 1e6
	    use_da_tol <- !is.null(mz_tol_da) && q_mz < mz_tol_da_max_mz
	    effective_tolerance_da <- if (isTRUE(use_da_tol)) {
	      pmax(ppm_da, mz_tol)
	    } else {
	      ppm_da
	    }
	    lower <- q_mz - effective_tolerance_da
	    upper <- q_mz + effective_tolerance_da
    start <- .fe_lower_bound(library_mz, lower)
    end <- .fe_upper_bound(library_mz, upper)
    hit_idx <- if (start <= end) seq.int(start, end) else integer()
    if (length(hit_idx) == 0L) {
      next
    }
    if (isTRUE(precursor_filter)) {
      precursor <- as.numeric(peak_index$precursor_mz[hit_idx])
      precursor_tolerance <- if (!is.null(precursor_tolerance_contract)) {
        .msi_effective_precursor_tolerance(
          query$anchor_mz[[i]], precursor_tolerance_contract
        )
      } else {
        as.numeric(precursor_tol)
      }
      keep_precursor <- is.finite(precursor) & is.finite(precursor_tolerance) &
        abs(precursor - as.numeric(query$anchor_mz[[i]])) <= precursor_tolerance
      hit_idx <- hit_idx[keep_precursor]
      if (length(hit_idx) == 0L) {
        next
      }
    }
    for (j in hit_idx) {
      lib_int <- as.numeric(peak_index$library_intensity[[j]])
      if (!is.finite(lib_int) || lib_int <= 0) {
        next
      }
      row_i <- row_i + 1L
      rows[[row_i]] <- data.frame(
        peak_id = as.character(query$peak_id[[i]]),
        anchor_feature_id = as.character(query$anchor_feature_id[[i]]),
        anchor_mz = as.numeric(query$anchor_mz[[i]]),
        query_fragment_id = as.character(query$fragment_id[[i]]),
        query_feature_id = as.character(query$fragment_feature_id[[i]]),
        query_mz = q_mz,
        query_intensity = q_int,
        query_SAS = as.numeric(query$SAS[[i]]),
        spectrum_id = as.character(peak_index$spectrum_id[[j]]),
        candidate_id = as.character(peak_index$candidate_id[[j]]),
        compound_id = as.character(peak_index$compound_id[[j]]),
        candidate_name = as.character(peak_index$candidate_name[[j]]),
        library_source = as.character(peak_index$library_source[[j]]),
        database_source = as.character(peak_index$database_source[[j]]),
        library_fragment_id = as.character(peak_index$library_fragment_id[[j]]),
        library_mz = as.numeric(peak_index$library_mz[[j]]),
        library_intensity = lib_int,
        ppm_error = as.numeric((q_mz - peak_index$library_mz[[j]]) / peak_index$library_mz[[j]] * 1e6),
        delta_h = .fe_entropy_increment(q_int, lib_int),
        stringsAsFactors = FALSE,
        check.names = FALSE
      )
    }
  }
  if (isTRUE(progress_enabled) && last_progress < total_queries) {
    .fe_emit_progress(
      progress_callback = progress_callback,
      progress_path = progress_path,
      started_at = progress_started_at,
      processed_queries = total_queries,
      total_queries = total_queries,
      candidates_so_far = row_i,
      message = "matching complete"
    )
  }

  if (length(rows) == 0L) {
    return(empty)
  }
  out <- do.call(rbind, rows)
  out$delta_h <- pmax(as.numeric(out$delta_h), 0)
  rownames(out) <- NULL
  out
}

`58_aggregate_entropy_score` <- 
function(matched_peaks,
                                         pseudo_ms2,
                                         index,
                                         top_n_spectra = 200,
                                         min_match_fe = 1,
                                         min_s_fe = 0) {
  empty <- .fe_empty_flash_hits()
  matched_peaks <- as.data.frame(matched_peaks, stringsAsFactors = FALSE, check.names = FALSE)
  if (nrow(matched_peaks) == 0L) {
    return(empty)
  }
  pseudo_ms2 <- as.data.frame(pseudo_ms2, stringsAsFactors = FALSE, check.names = FALSE)
  index <- .fe_peak_index(index)
  query_total <- stats::aggregate(
    intensity ~ peak_id + anchor_feature_id,
    data = pseudo_ms2,
    FUN = function(x) sum(as.numeric(x), na.rm = TRUE)
  )
  names(query_total)[names(query_total) == "intensity"] <- "query_total_intensity"
  library_total <- stats::aggregate(
    library_intensity ~ spectrum_id,
    data = index,
    FUN = function(x) sum(as.numeric(x), na.rm = TRUE)
  )
  names(library_total)[names(library_total) == "library_intensity"] <- "library_total_intensity"

  key_cols <- c("peak_id", "anchor_feature_id", "anchor_mz", "candidate_id", "spectrum_id")
  split_key <- .fe_split_key(matched_peaks, key_cols)
  pieces <- split(matched_peaks, split_key, drop = TRUE)
  out <- lapply(pieces, function(df) {
    query_unique <- df[!duplicated(df$query_fragment_id), c("query_fragment_id", "query_intensity"), drop = FALSE]
    lib_unique <- df[!duplicated(df$library_fragment_id), c("library_fragment_id", "library_intensity"), drop = FALSE]
    data.frame(
      peak_id = as.character(df$peak_id[[1]]),
      anchor_feature_id = as.character(df$anchor_feature_id[[1]]),
      anchor_mz = as.numeric(df$anchor_mz[[1]]),
      candidate_id = as.character(df$candidate_id[[1]]),
      spectrum_id = as.character(df$spectrum_id[[1]]),
      entropy_S = sum(as.numeric(df$delta_h), na.rm = TRUE),
      n_match_fe = nrow(df),
      matched_query_peak_n = length(unique(as.character(df$query_fragment_id))),
      matched_library_peak_n = length(unique(as.character(df$library_fragment_id))),
      matched_query_intensity_sum = sum(as.numeric(query_unique$query_intensity), na.rm = TRUE),
      matched_library_intensity_sum = sum(as.numeric(lib_unique$library_intensity), na.rm = TRUE),
      stringsAsFactors = FALSE,
      check.names = FALSE
    )
  })
  out <- do.call(rbind, out)
  out <- merge(out, query_total, by = c("peak_id", "anchor_feature_id"), all.x = TRUE, sort = FALSE)
  out <- merge(out, library_total, by = "spectrum_id", all.x = TRUE, sort = FALSE)
  out$query_usage_fe <- out$matched_query_intensity_sum / out$query_total_intensity
  out$library_usage_fe <- out$matched_library_intensity_sum / out$library_total_intensity
  out$query_usage_fe[!is.finite(out$query_usage_fe)] <- 0
  out$library_usage_fe[!is.finite(out$library_usage_fe)] <- 0
  out <- out[
    is.finite(out$entropy_S) &
      out$entropy_S >= as.numeric(min_s_fe) &
      out$n_match_fe >= as.integer(min_match_fe),
    ,
    drop = FALSE
  ]
  if (nrow(out) == 0L) {
    return(empty)
  }

  meta <- index[!duplicated(.fe_split_key(index, c("candidate_id", "spectrum_id"))), c(
    "candidate_id", "spectrum_id", "candidate_name", "library_source",
    "database_source", "compound_id", "precursor_mz", "formula",
    "inchikey", "smiles", "adduct", "polarity"
  ), drop = FALSE]
  out <- merge(out, meta, by = c("candidate_id", "spectrum_id"), all.x = TRUE, sort = FALSE)
  out$query_id <- as.character(out$peak_id)
  out$precursor_error_ppm <- (as.numeric(out$anchor_mz) - as.numeric(out$precursor_mz)) /
    as.numeric(out$precursor_mz) * 1e6
  out$precursor_error_ppm[!is.finite(out$precursor_error_ppm)] <- NA_real_
  out$S_FE <- out$entropy_S
  out <- out[!duplicated(.fe_split_key(out, key_cols)), , drop = FALSE]
  out <- out[order(out$peak_id, out$anchor_feature_id, -out$entropy_S, out$candidate_id, out$spectrum_id), , drop = FALSE]
  top_n_spectra <- as.integer(top_n_spectra)
  if (!is.finite(top_n_spectra) || top_n_spectra < 1L) {
    top_n_spectra <- nrow(out)
  }
  split_anchor <- .fe_split_key(out, c("peak_id", "anchor_feature_id"))
  keep <- unlist(lapply(split(seq_len(nrow(out)), split_anchor), function(idx) idx[seq_len(min(length(idx), top_n_spectra))]), use.names = FALSE)
  out <- out[sort(keep), , drop = FALSE]
  out <- do.call(rbind, lapply(split(out, .fe_split_key(out, c("peak_id", "anchor_feature_id"))), function(df) {
    df <- df[order(-df$entropy_S, df$candidate_id, df$spectrum_id), , drop = FALSE]
    df$rank_spectrum_within_anchor <- seq_len(nrow(df))
    df$retrieval_rank <- df$rank_spectrum_within_anchor
    df$rank_entropy <- df$rank_spectrum_within_anchor
    df
  }))
  out$s_fe <- out$entropy_S
  out$n_FE_match <- out$n_match_fe
  out$matched_query_peaks <- out$matched_query_peak_n
  out$matched_library_peaks <- out$matched_library_peak_n
  out$cluster_id <- out$peak_id
  out <- out[, colnames(empty), drop = FALSE]
  rownames(out) <- NULL
  out
}

`58b_collapse_flash_entropy_candidates` <- 
function(spectrum_hits,
                                                    top_n_compounds = 100,
                                                    max_spectra_per_compound = 3) {
  empty <- .fe_empty_candidate_hits()
  spectrum_hits <- as.data.frame(spectrum_hits, stringsAsFactors = FALSE, check.names = FALSE)
  if (nrow(spectrum_hits) == 0L) {
    return(empty)
  }
  spectrum_hits$compound_key <- .fe_compound_key(spectrum_hits)
  spectrum_hits <- spectrum_hits[order(
    spectrum_hits$peak_id,
    spectrum_hits$anchor_feature_id,
    spectrum_hits$compound_key,
    -spectrum_hits$entropy_S,
    spectrum_hits$spectrum_id
  ), , drop = FALSE]
  max_spectra_per_compound <- as.integer(max_spectra_per_compound)
  if (!is.finite(max_spectra_per_compound) || max_spectra_per_compound < 1L) {
    max_spectra_per_compound <- nrow(spectrum_hits)
  }
  split_compound <- .fe_split_key(spectrum_hits, c("peak_id", "anchor_feature_id", "compound_key"))
  limited <- do.call(rbind, lapply(split(spectrum_hits, split_compound), function(df) {
    df[seq_len(min(nrow(df), max_spectra_per_compound)), , drop = FALSE]
  }))
  split_limited <- .fe_split_key(limited, c("peak_id", "anchor_feature_id", "compound_key"))

  out <- lapply(split(limited, split_limited), function(df) {
    best <- df[order(-df$entropy_S, df$spectrum_id), , drop = FALSE][1, , drop = FALSE]
    data.frame(
      peak_id = best$peak_id[[1]],
      anchor_feature_id = best$anchor_feature_id[[1]],
      anchor_mz = as.numeric(best$anchor_mz[[1]]),
      compound_key = best$compound_key[[1]],
      candidate_id = best$candidate_id[[1]],
      compound_id = best$compound_id[[1]],
      candidate_name = best$candidate_name[[1]],
      inchikey = best$inchikey[[1]],
      formula = best$formula[[1]],
      smiles = best$smiles[[1]],
      database_source = best$database_source[[1]],
      best_spectrum_id = best$spectrum_id[[1]],
      best_library_source = best$library_source[[1]],
      entropy_S = max(as.numeric(df$entropy_S), na.rm = TRUE),
      S_FE = max(as.numeric(df$entropy_S), na.rm = TRUE),
      S_FE_compound = max(as.numeric(df$entropy_S), na.rm = TRUE),
      n_match_fe = best$n_match_fe[[1]],
      matched_query_peak_n = best$matched_query_peak_n[[1]],
      matched_library_peak_n = best$matched_library_peak_n[[1]],
      query_usage_fe = best$query_usage_fe[[1]],
      library_usage_fe = best$library_usage_fe[[1]],
      n_spectra_for_compound = length(unique(as.character(df$spectrum_id))),
      stringsAsFactors = FALSE,
      check.names = FALSE
    )
  })
  out <- do.call(rbind, out)
  out <- out[order(out$peak_id, out$anchor_feature_id, -out$entropy_S, out$compound_key), , drop = FALSE]
  top_n_compounds <- as.integer(top_n_compounds)
  if (!is.finite(top_n_compounds) || top_n_compounds < 1L) {
    top_n_compounds <- nrow(out)
  }
  out <- do.call(rbind, lapply(split(out, .fe_split_key(out, c("peak_id", "anchor_feature_id"))), function(df) {
    df <- df[order(-df$entropy_S, df$compound_key), , drop = FALSE]
    df <- df[seq_len(min(nrow(df), top_n_compounds)), , drop = FALSE]
    df$rank_compound_within_anchor <- seq_len(nrow(df))
    df
  }))
  out <- out[, colnames(empty), drop = FALSE]
  rownames(out) <- NULL
  out
}

`58c_build_entropy_candidates` <- 
function(spectrum_hits,
                                           matched_peaks,
                                           top_n_entropy_candidates = 500) {
  entropy_candidates <- .fe_empty_entropy_candidates()
  entropy_filter_summary <- .fe_empty_entropy_filter_summary()
  spectrum_hits <- as.data.frame(spectrum_hits, stringsAsFactors = FALSE, check.names = FALSE)
  if (nrow(spectrum_hits) == 0L) {
    return(list(
      entropy_candidates = entropy_candidates,
      entropy_filter_summary = entropy_filter_summary
    ))
  }
  top_n_entropy_candidates <- as.integer(top_n_entropy_candidates[[1]])
  if (!is.finite(top_n_entropy_candidates) || top_n_entropy_candidates < 1L) {
    top_n_entropy_candidates <- 500L
  }
  spectrum_hits <- .fe_add_missing_columns(spectrum_hits, c(
    "query_id", "peak_id", "anchor_feature_id", "anchor_mz", "candidate_id",
    "spectrum_id", "candidate_name", "precursor_mz", "entropy_S",
    "n_match_fe", "rank_entropy", "library_source", "database_source",
    "formula", "adduct", "adduct_status"
  ))
  if (all(is.na(spectrum_hits$query_id))) {
    spectrum_hits$query_id <- as.character(spectrum_hits$peak_id)
  }
  spectrum_hits$entropy_S <- suppressWarnings(as.numeric(spectrum_hits$entropy_S))
  if (all(is.na(spectrum_hits$entropy_S)) && "S_FE" %in% names(spectrum_hits)) {
    spectrum_hits$entropy_S <- suppressWarnings(as.numeric(spectrum_hits$S_FE))
  }
  matched_peaks <- as.data.frame(matched_peaks, stringsAsFactors = FALSE, check.names = FALSE)
  matched_lookup <- character()
  if (nrow(matched_peaks) > 0L) {
    matched_peaks <- .fe_add_missing_columns(matched_peaks, c(
      "peak_id", "anchor_feature_id", "candidate_id", "spectrum_id",
      "query_fragment_id", "query_mz", "library_mz"
    ))
    matched_key <- .fe_split_key(matched_peaks, c("peak_id", "anchor_feature_id", "candidate_id", "spectrum_id"))
    pair_text <- paste0(
      as.character(matched_peaks$query_fragment_id), ":",
      signif(suppressWarnings(as.numeric(matched_peaks$query_mz)), 8), "->",
      signif(suppressWarnings(as.numeric(matched_peaks$library_mz)), 8)
    )
    matched_lookup <- vapply(split(pair_text, matched_key), function(x) {
      paste(unique(x), collapse = ";")
    }, character(1))
  }

  split_query <- .fe_split_key(spectrum_hits, c("query_id", "anchor_feature_id"))
  pieces <- lapply(split(spectrum_hits, split_query), function(df) {
    df <- df[order(-df$entropy_S, df$candidate_id, df$spectrum_id), , drop = FALSE]
    df$rank_entropy <- seq_len(nrow(df))
    df[df$rank_entropy <= top_n_entropy_candidates, , drop = FALSE]
  })
  kept <- do.call(rbind, pieces)
  kept_key <- .fe_split_key(kept, c("peak_id", "anchor_feature_id", "candidate_id", "spectrum_id"))
  matched_peak_pairs <- unname(matched_lookup[kept_key])
  matched_peak_pairs[is.na(matched_peak_pairs)] <- ""
  entropy_candidates <- data.frame(
    query_id = as.character(kept$query_id),
    pseudo_ms2_id = as.character(kept$query_id),
    anchor_feature_id = as.character(kept$anchor_feature_id),
    anchor_mz = suppressWarnings(as.numeric(kept$anchor_mz)),
    library_spectrum_id = as.character(kept$spectrum_id),
    library_id = as.character(kept$candidate_id),
    library_source = as.character(kept$library_source),
    candidate_id = as.character(kept$candidate_id),
    spectrum_id = as.character(kept$spectrum_id),
    metabolite_name = as.character(kept$candidate_name),
    candidate_name = as.character(kept$candidate_name),
    formula = as.character(kept$formula),
    adduct = as.character(kept$adduct),
    adduct_status = .msi_adduct_status(kept$adduct),
    library_precursor_mz = suppressWarnings(as.numeric(kept$precursor_mz)),
    precursor_mz = suppressWarnings(as.numeric(kept$precursor_mz)),
    precursor_error_da = suppressWarnings(as.numeric(kept$anchor_mz) - as.numeric(kept$precursor_mz)),
    precursor_error_ppm = suppressWarnings((as.numeric(kept$anchor_mz) - as.numeric(kept$precursor_mz)) /
      as.numeric(kept$precursor_mz) * 1e6),
    entropy_score = suppressWarnings(as.numeric(kept$entropy_S)),
    entropy_S = suppressWarnings(as.numeric(kept$entropy_S)),
    entropy_rank = suppressWarnings(as.integer(kept$rank_entropy)),
    n_matched_peaks = suppressWarnings(as.integer(kept$n_match_fe)),
    n_matched_peaks_entropy = suppressWarnings(as.integer(kept$n_match_fe)),
    matched_peak_summary = matched_peak_pairs,
    matched_peak_pairs = matched_peak_pairs,
    rank_entropy = suppressWarnings(as.integer(kept$rank_entropy)),
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
  entropy_candidates <- .fe_add_missing_columns(entropy_candidates, colnames(.fe_empty_entropy_candidates()))
  entropy_candidates <- entropy_candidates[, colnames(.fe_empty_entropy_candidates()), drop = FALSE]

  before <- stats::aggregate(
    candidate_id ~ query_id + anchor_feature_id,
    data = spectrum_hits,
    FUN = length
  )
  names(before)[names(before) == "candidate_id"] <- "n_entropy_candidates_before_top500"
  after <- stats::aggregate(
    candidate_id ~ query_id + anchor_feature_id,
    data = entropy_candidates,
    FUN = length
  )
  names(after)[names(after) == "candidate_id"] <- "n_entropy_candidates_after_top500"
  first_anchor <- spectrum_hits[!duplicated(.fe_split_key(spectrum_hits, c("query_id", "anchor_feature_id"))), c("query_id", "anchor_feature_id"), drop = FALSE]
  entropy_filter_summary <- merge(first_anchor, before, by = c("query_id", "anchor_feature_id"), all.x = TRUE, sort = FALSE)
  entropy_filter_summary <- merge(entropy_filter_summary, after, by = c("query_id", "anchor_feature_id"), all.x = TRUE, sort = FALSE)
  entropy_filter_summary$n_entropy_candidates_before_top500[is.na(entropy_filter_summary$n_entropy_candidates_before_top500)] <- 0L
  entropy_filter_summary$n_entropy_candidates_after_top500[is.na(entropy_filter_summary$n_entropy_candidates_after_top500)] <- 0L
  entropy_filter_summary <- entropy_filter_summary[, colnames(.fe_empty_entropy_filter_summary()), drop = FALSE]
  rownames(entropy_candidates) <- NULL
  rownames(entropy_filter_summary) <- NULL
  list(
    entropy_candidates = entropy_candidates,
    entropy_filter_summary = entropy_filter_summary
  )
}

`59_store_flash_entropy_hits` <- 
function(object,
                                          library,
                                          index,
                                          matched_peaks,
                                          spectrum_hits,
                                          candidate_hits,
                                          entropy_candidates,
                                          entropy_filter_summary,
                                          spectrum_id_peak_index = NULL,
                                          matching_library_cache_identity = NULL,
                                          parameters) {
  spectral_library <- list(
    library_spectrum_table = library$library_spectrum_table,
    library_peak_table = library$library_peak_table,
    flash_entropy_index = index
  )
  if (!is.null(spectrum_id_peak_index)) {
    spectral_library$spectrum_id_peak_index <- spectrum_id_peak_index
  }
  if (!is.null(matching_library_cache_identity)) {
    spectral_library$matching_library_cache_identity <- matching_library_cache_identity
  }
  object@spectral_library <- spectral_library
  attr(object@spectral_library, "msidataset_standardized_library") <- TRUE
  attr(object@spectral_library, "database_root") <- parameters$database_root
  attr(object@spectral_library, "ms1_database_dir") <- parameters$ms1_database_dir
  attr(object@spectral_library, "ms2_database_dir") <- parameters$ms2_database_dir
  object@annotation_results$flash_entropy_matched_peaks <- matched_peaks
  object@annotation_results$flash_entropy_spectrum_hits <- spectrum_hits
  object@annotation_results$flash_entropy_candidate_hits <- candidate_hits
  object@annotation_results$flash_entropy_hits <- spectrum_hits
  object@annotation_results$entropy_candidates <- entropy_candidates
  object@annotation_results$entropy_filter_summary <- entropy_filter_summary
  object@ms2_data$entropy_candidates <- entropy_candidates
  object@ms2_data$entropy_filter_summary <- entropy_filter_summary
  parameters$note <- paste(
    "Entropy S is used for high-recall candidate retrieval.",
    "Module 5 keeps the top Entropy S candidates per pseudo-MS2 query."
  )
  object@process_info$module5_entropy_search <- parameters
  object@msi_parameters$calculate_flash_entropy <- parameters
  object
}

`6_reverse_match` <- 
function(object,
                              archive_library = NULL,
                              archive_index = NULL,
                              ppm = 10,
                              fragment_tolerance_ppm = ppm,
                              fragment_tolerance_da = NULL,
                              fragment_tolerance_mode = "max",
                              low_mz_cutoff = 400,
                              fixed_da_low_mz = 0.005,
                              mz_bin_width = 0.01,
                              bin_margin = 0.5,
                              precursor_tol = 0.05,
                              delta_grid = c(0),
                              archive_top_n = 200,
                              use_flash_hits = TRUE,
                              use_archive_hits = !is.null(archive_library) || !is.null(archive_index),
                              min_candidate_peaks = 1,
                              min_matched_peaks = 3L,
                              min_matched_peaks_evidence = min_matched_peaks,
                              min_reverse_masked_cosine_evidence = 0,
                              matched_library_intensity_fraction_min = 0.30,
                              reverse_match_mode = c("strict", "permissive"),
                              reverse_filter_mode = "reverse_masked_cosine_evidence",
                              moderate_min_lib_coverage = 0.6,
                              moderate_min_matched_peaks = 3L,
                              moderate_max_precursor_ppm = 50,
                              workers = 1L,
                              progress_every = NULL,
                              progress_path = NULL,
                              progress_started_at = NULL,
                              progress_callback = NULL,
                              export_matched_peak_pairs = FALSE,
                              export_matched_peak_pairs_all = FALSE,
                              max_matched_peak_pairs_export_rows = 100000L,
                              compress_large_pair_tables = TRUE,
                              store_matched_peak_pairs_in_rds = TRUE,
                              output_dir = "output",
                              tolerance_contract = NULL,
                              ...,
                              spectrum_id_peak_index = NULL) {
  reverse_match_mode <- match.arg(reverse_match_mode)
  export_matched_peak_pairs <- .rm_bool(export_matched_peak_pairs, FALSE)
  export_matched_peak_pairs_all <- .rm_bool(export_matched_peak_pairs_all, FALSE)
  compress_large_pair_tables <- .rm_bool(compress_large_pair_tables, TRUE)
  store_matched_peak_pairs_in_rds <- .rm_bool(store_matched_peak_pairs_in_rds, TRUE)
  max_matched_peak_pairs_export_rows <- suppressWarnings(as.integer(max_matched_peak_pairs_export_rows[[1]]))
  if (!is.finite(max_matched_peak_pairs_export_rows) || max_matched_peak_pairs_export_rows < 1L) {
    max_matched_peak_pairs_export_rows <- 100000L
  }
  legacy_args <- list(...)
  legacy_arg_names <- names(legacy_args)
  old_arg_names <- c(
    paste0("reverse", "_dp_threshold"),
    paste0("reverse", "_dp_candidate_threshold"),
    paste0("reverse", "_dp_strong_threshold"),
    paste0("reverse", "_cos_candidate_threshold"),
    paste0("reverse", "_cos_strong_threshold")
  )
  if (length(intersect(legacy_arg_names, old_arg_names))) {
    warning(
      "Deprecated reverse-match threshold arguments are ignored; Module 6 uses reverse_masked_cosine evidence.",
      call. = FALSE
    )
  }
  reverse_filter_mode <- as.character(reverse_filter_mode[[1]])
  if (reverse_filter_mode %in% c(paste0("reverse", "_dp"), paste0("reverse", "_dp_tiered"), "spectral_shape")) {
    warning(
      "Deprecated reverse-match filter modes are ignored; using reverse_masked_cosine_evidence.",
      call. = FALSE
    )
    reverse_filter_mode <- "reverse_masked_cosine_evidence"
  }
  if (!identical(reverse_filter_mode, "reverse_masked_cosine_evidence")) {
    stop("reverse_filter_mode must be 'reverse_masked_cosine_evidence'.", call. = FALSE)
  }
  if (!is.null(fragment_tolerance_da)) {
    warning(
      "fragment_tolerance_da is deprecated for the main reverse masked cosine path; ",
      "using it as fixed_da_low_mz for backward compatibility. Use ",
      "reverse_matching.low_mz_abs_tol / fixed_da_low_mz instead.",
      call. = FALSE
    )
    fixed_da_low_mz <- fragment_tolerance_da
  }
  if (!missing(mz_bin_width) || !missing(bin_margin)) {
    warning(
      "mz_bin_width and bin_margin no longer control main reverse masked cosine matching; ",
      "they are only used by archive_search compatibility paths.",
      call. = FALSE
    )
  }
  if (!methods::is(object, "msi_dataset")) {
    stop("object must be a msi_dataset.", call. = FALSE)
  }
  if (!is.null(tolerance_contract)) {
    if (!exists(".msi_resolve_tolerance_contract", mode = "function")) {
      stop("Shared tolerance contract helpers are unavailable.", call. = FALSE)
    }
    tolerance_contract <- .msi_resolve_tolerance_contract(
      tolerance_contract,
      feature_mz = object@variable_info$mz
    )
    evidence_policy <- tolerance_contract$fragment_evidence
    fragment_tolerance_ppm <- as.numeric(evidence_policy$ppm)
    fragment_tolerance_mode <- as.character(evidence_policy$tolerance_mode)
    low_mz_cutoff <- as.numeric(evidence_policy$low_mz_cutoff)
    fixed_da_low_mz <- as.numeric(evidence_policy$low_mz_da_floor)
  }
  if (is.null(object@pseudo_ms2$pseudo_ms2_table)) {
    .msi_missing_summary("6_reverse_match", "pseudo_ms2$pseudo_ms2_table", output_dir)
    return(object)
  }
  if (isTRUE(use_flash_hits) && is.null(object@annotation_results$flash_entropy_hits)) {
    .msi_missing_summary("6_reverse_match", "annotation_results$flash_entropy_hits", output_dir)
    return(object)
  }
  if (isTRUE(use_flash_hits) &&
      (is.null(object@spectral_library$library_spectrum_table) ||
       is.null(object@spectral_library$library_peak_table))) {
    .msi_missing_summary("6_reverse_match", c("spectral_library$library_spectrum_table", "spectral_library$library_peak_table"), output_dir)
    return(object)
  }

  if (isTRUE(use_archive_hits)) {
    if (is.null(archive_index)) {
      archive_index <- `60_build_archive_index`(
        archive_library = archive_library,
        ppm = fragment_tolerance_ppm,
        min_candidate_peaks = min_candidate_peaks,
        ...
      )
    } else if (!all(c("archive_spectrum_table", "archive_peak_table", "precursor_index", "fragment_index") %in% names(archive_index))) {
      archive_index <- `60_build_archive_index`(
        archive_library = archive_index,
        ppm = fragment_tolerance_ppm,
        min_candidate_peaks = min_candidate_peaks,
        ...
      )
    }
    object@archive_index <- archive_index
    archive_hits <- `60_search_archive_candidates`(
      object = object,
      archive_index = archive_index,
      fragment_tolerance_ppm = fragment_tolerance_ppm,
      fragment_tolerance_da = fragment_tolerance_da,
      fragment_tolerance_mode = fragment_tolerance_mode,
      low_mz_cutoff = low_mz_cutoff,
      fixed_da_low_mz = fixed_da_low_mz,
      mz_bin_width = mz_bin_width,
      bin_margin = bin_margin,
      precursor_tol = precursor_tol,
      tolerance_contract = tolerance_contract,
      delta_grid = delta_grid,
      archive_top_n = archive_top_n
    )
  } else {
    archive_index <- list()
    object@archive_index <- archive_index
    archive_hits <- .rm_empty_archive_hits()
  }

  candidates <- `60_prepare_reverse_candidates`(
    flash_hits = if (isTRUE(use_flash_hits)) `60_entropy_candidates_for_reverse`(object) else NULL,
    archive_hits = archive_hits,
    use_flash_hits = use_flash_hits,
    use_archive_hits = use_archive_hits
  )

  spectrum_index_started_at <- Sys.time()
  if (is.null(spectrum_id_peak_index)) {
    spectrum_id_peak_index <- .rm_build_spectrum_id_peak_index(
      object@spectral_library$library_peak_table
    )
    spectrum_index_build_count <- 1L
    spectrum_index_origin <- "built_internally"
    spectrum_index_validation_seconds <- 0
  } else {
    spectrum_id_peak_index <- .mlc_validate_prebuilt_spectrum_id_peak_index(
      spectrum_id_peak_index,
      object@spectral_library$library_peak_table,
      object@spectral_library$matching_library_cache_identity
    )
    spectrum_index_build_count <- 0L
    spectrum_index_origin <- "prebuilt_matching_library_cache"
    spectrum_index_validation_seconds <- as.numeric(difftime(Sys.time(), spectrum_index_started_at, units = "secs"))
  }
  spectrum_index_build_seconds <- if (isTRUE(spectrum_index_build_count == 1L)) {
    as.numeric(difftime(Sys.time(), spectrum_index_started_at, units = "secs"))
  } else {
    0
  }
  spectrum_index_telemetry <- new.env(parent = emptyenv())
  spectrum_index_telemetry$direct_accesses <- 0L
  spectrum_index_telemetry$fallback_linear_scans <- 0L

  matched <- `60_reverse_match_library_peaks`(
    object = object,
    candidates = candidates,
    fragment_tolerance_ppm = fragment_tolerance_ppm,
    fragment_tolerance_da = fragment_tolerance_da,
    fragment_tolerance_mode = fragment_tolerance_mode,
    low_mz_cutoff = low_mz_cutoff,
    fixed_da_low_mz = fixed_da_low_mz,
    mz_bin_width = mz_bin_width,
    bin_margin = bin_margin,
    workers = workers,
    progress_every = progress_every,
    progress_path = progress_path,
    progress_started_at = progress_started_at,
    progress_callback = progress_callback,
    spectrum_id_peak_index = spectrum_id_peak_index,
    spectrum_index_telemetry = spectrum_index_telemetry
  )

	  reverse_match_all <- `60_summarise_reverse_match`(
	    candidates = candidates,
	    matched_peak_pairs = matched,
	    min_matched_peaks_evidence = min_matched_peaks_evidence,
	    min_reverse_masked_cosine_evidence = min_reverse_masked_cosine_evidence,
	    progress_every = progress_every,
    progress_path = progress_path,
    progress_started_at = progress_started_at,
    progress_callback = progress_callback
  )

  rdp_filter <- `60_apply_reverse_masked_cosine_filter`(
    reverse_match = reverse_match_all,
    matched_peak_pairs = matched,
    min_matched_peaks_evidence = min_matched_peaks_evidence,
    min_reverse_masked_cosine_evidence = min_reverse_masked_cosine_evidence
  )
  if (is.data.frame(rdp_filter$all_candidates)) {
    reverse_match_all <- rdp_filter$all_candidates
  }

  reverse_match <- `60_mark_chemical_evidence_not_evaluated`(
    reverse_match = rdp_filter$pseudo_ms2_candidate_evidence
  )
  matched_peak_pairs_summary <- .rm_build_matched_peak_pairs_summary(
    matched_peak_pairs = rdp_filter$matched_peak_pairs,
    matched_peak_pairs_all = matched,
    reverse_match = reverse_match,
    reverse_match_all = reverse_match_all,
    export_matched_peak_pairs = export_matched_peak_pairs,
    export_matched_peak_pairs_all = export_matched_peak_pairs_all,
    store_matched_peak_pairs_in_rds = store_matched_peak_pairs_in_rds,
    max_matched_peak_pairs_export_rows = max_matched_peak_pairs_export_rows,
    compress_large_pair_tables = compress_large_pair_tables,
    skipped_matched_peak_pairs_csv = !isTRUE(export_matched_peak_pairs),
    skipped_matched_peak_pairs_all_csv = !isTRUE(export_matched_peak_pairs_all),
    export_reason = if (!isTRUE(export_matched_peak_pairs) && !isTRUE(export_matched_peak_pairs_all)) {
      "export_disabled"
    } else {
      "pending_caller_export"
    }
  )
  .rm_emit_progress(
    progress_callback = progress_callback,
    progress_path = progress_path,
    started_at = progress_started_at,
    processed_candidates = nrow(candidates),
    total_candidates = nrow(candidates),
    reverse_rows_so_far = nrow(reverse_match_all),
    pass_rows_so_far = nrow(reverse_match),
    message = "reverse scoring completed"
  )

  object <- `60_store_reverse_match_results`(
    object = object,
    archive_index = archive_index,
    archive_hits = archive_hits,
    merged_candidates = candidates,
    reverse_matched_peaks_all = matched,
    reverse_matched_peaks = rdp_filter$matched_peak_pairs,
    reverse_match = reverse_match,
    reverse_match_all = reverse_match_all,
    pseudo_ms2_candidate_evidence_summary = rdp_filter$pseudo_ms2_candidate_evidence_summary,
    matched_peak_pairs_summary = matched_peak_pairs_summary,
    store_matched_peak_pairs_in_rds = store_matched_peak_pairs_in_rds,
    parameters = list(
      fragment_ppm = fragment_tolerance_ppm,
      tolerance_mode = "hybrid_low_mz",
      low_mz_cutoff = low_mz_cutoff,
      low_mz_abs_tol = fixed_da_low_mz,
      fixed_da_low_mz = fixed_da_low_mz,
      archive_search = list(
        mz_bin_width = if (isTRUE(use_archive_hits)) mz_bin_width else NA_real_,
        bin_margin = if (isTRUE(use_archive_hits)) bin_margin else NA_real_,
        fragment_tolerance_da = if (isTRUE(use_archive_hits)) fragment_tolerance_da else NA_real_
      ),
      precursor_tol = precursor_tol,
      delta_grid = delta_grid,
      archive_top_n = archive_top_n,
      use_flash_hits = use_flash_hits,
      use_archive_hits = use_archive_hits,
      min_candidate_peaks = min_candidate_peaks,
      min_matched_peaks_evidence = max(1L, as.integer(min_matched_peaks_evidence[[1]])),
      min_reverse_masked_cosine_evidence = as.numeric(min_reverse_masked_cosine_evidence[[1]]),
      matched_library_intensity_fraction_role = "diagnostic_only",
      reverse_match_mode = reverse_match_mode,
      reverse_filter_mode = "reverse_masked_cosine_evidence",
      moderate_min_lib_coverage = moderate_min_lib_coverage,
      moderate_min_matched_peaks = moderate_min_matched_peaks,
      moderate_max_precursor_ppm = moderate_max_precursor_ppm,
      tolerance_contract = tolerance_contract,
      fragment_evidence_tolerance = if (is.null(tolerance_contract)) {
        NA_character_
      } else {
        paste0(
          tolerance_contract$fragment_evidence$tolerance_mode,
          ": max(", tolerance_contract$fragment_evidence$ppm,
          " ppm, ", tolerance_contract$fragment_evidence$low_mz_da_floor,
          " Da) below ", tolerance_contract$fragment_evidence$low_mz_cutoff, " m/z"
        )
      },
      precursor_candidate_tolerance = if (is.null(tolerance_contract)) {
        NA_character_
      } else {
        paste0(
          tolerance_contract$precursor_candidate$mode, ": max(",
          tolerance_contract$precursor_candidate$ppm,
          " ppm, ", tolerance_contract$precursor_candidate$da_floor, " Da)"
        )
      },
      workers = workers,
      spectrum_id_peak_index_build_seconds = spectrum_index_build_seconds,
      spectrum_id_peak_index_validation_seconds = spectrum_index_validation_seconds,
      spectrum_id_peak_index_build_count = as.integer(spectrum_index_build_count),
      spectrum_id_peak_index_origin = spectrum_index_origin,
      spectrum_id_peak_index_bytes = as.numeric(utils::object.size(spectrum_id_peak_index)),
      spectrum_id_peak_index_direct_accesses = as.integer(spectrum_index_telemetry$direct_accesses),
      spectrum_id_peak_index_fallback_linear_scans = as.integer(spectrum_index_telemetry$fallback_linear_scans),
      export_matched_peak_pairs = export_matched_peak_pairs,
      export_matched_peak_pairs_all = export_matched_peak_pairs_all,
      max_matched_peak_pairs_export_rows = max_matched_peak_pairs_export_rows,
      compress_large_pair_tables = compress_large_pair_tables,
      store_matched_peak_pairs_in_rds = store_matched_peak_pairs_in_rds,
      note = paste(
        "Entropy S is used for high-recall candidate retrieval.",
        "Module 6 uses reverse_masked_cosine as the primary library-space reverse-match score.",
        "matched_library_intensity_fraction is retained as a diagnostic field only."
      )
    )
  )

  object
}

`60_apply_reverse_masked_cosine_filter` <- 
function(reverse_match,
                                                    matched_peak_pairs,
                                                    min_matched_peaks_evidence = 3L,
                                                    min_reverse_masked_cosine_evidence = 0) {
  empty_scores <- .rm_empty_reverse_scores()
  empty_pairs <- .rm_empty_reverse_matched_peaks_full()
  reverse_match <- as.data.frame(reverse_match, stringsAsFactors = FALSE, check.names = FALSE)
  matched_peak_pairs <- as.data.frame(matched_peak_pairs, stringsAsFactors = FALSE, check.names = FALSE)
  min_matched_peaks_evidence <- suppressWarnings(as.integer(min_matched_peaks_evidence[[1]]))
  if (!is.finite(min_matched_peaks_evidence) || min_matched_peaks_evidence < 1L) {
    min_matched_peaks_evidence <- 3L
  }
  min_reverse_masked_cosine_evidence <- suppressWarnings(as.numeric(min_reverse_masked_cosine_evidence[[1]]))
  if (!is.finite(min_reverse_masked_cosine_evidence) || min_reverse_masked_cosine_evidence < 0) {
    min_reverse_masked_cosine_evidence <- 0
  }
  if (nrow(reverse_match) == 0L) {
    return(list(
      all_candidates = empty_scores,
      pseudo_ms2_candidate_evidence = empty_scores,
      pseudo_ms2_candidate_evidence_summary = .rm_empty_rdp_filter_summary(),
      matched_peak_pairs = empty_pairs
    ))
  }
  reverse_match <- .rm_add_missing_columns(reverse_match, colnames(empty_scores))
  missing_query <- is.na(reverse_match$query_id) |
    !nzchar(as.character(reverse_match$query_id))
  reverse_match$query_id[missing_query] <-
    as.character(reverse_match$cluster_id[missing_query])
  missing_candidate <- is.na(reverse_match$candidate_id) |
    !nzchar(as.character(reverse_match$candidate_id))
  reverse_match$candidate_id[missing_candidate] <-
    as.character(reverse_match$candidate_key[missing_candidate])

  n_matched <- suppressWarnings(as.integer(reverse_match$n_matched_peaks))
  n_matched[!is.finite(n_matched)] <- 0L
  rmc <- suppressWarnings(as.numeric(reverse_match$reverse_masked_cosine))
  rmc[!is.finite(rmc)] <- 0
  sas <- suppressWarnings(as.numeric(reverse_match$S_SAS))
  sas_score <- ifelse(is.finite(sas), pmin(pmax(sas, 0), 1), 0)
  s_total <- 0.7 * pmin(pmax(rmc, 0), 1) + 0.3 * sas_score
  reverse_match$reverse_masked_cosine <- pmin(pmax(rmc, 0), 1)
  reverse_match$S_SAS_score_for_total <- sas_score
  reverse_match$S_SAS_missing <- !is.finite(sas)
  reverse_match$S_total <- pmin(pmax(s_total, 0), 1)
  reverse_match$score_ms2 <- reverse_match$reverse_masked_cosine
  reverse_match$score_spatial <- reverse_match$S_SAS_score_for_total
  reverse_match$score_total_pre_mfc <- reverse_match$S_total
  reverse_match$match_mode <- "pseudo_ms2"
  reverse_match$pass_n_matched_peaks_ge3 <- n_matched >= min_matched_peaks_evidence
  reverse_match$pass_reverse_masked_cosine <- reverse_match$reverse_masked_cosine > min_reverse_masked_cosine_evidence
  reverse_match$evidence_pass <- reverse_match$pass_n_matched_peaks_ge3 &
    reverse_match$pass_reverse_masked_cosine
  reverse_match$min_matched_peaks_used <- as.integer(min_matched_peaks_evidence)
  reverse_match$min_reverse_masked_cosine_used <- as.numeric(min_reverse_masked_cosine_evidence)
  reverse_match$pass_min_matched_peaks <- reverse_match$pass_n_matched_peaks_ge3
  reverse_match$pass_reverse_match <- reverse_match$evidence_pass
  reverse_match$passes_spectral_shape_filter <- reverse_match$evidence_pass
  reverse_match$spectral_shape_support <- ifelse(reverse_match$evidence_pass, "evidence_pass", "evidence_fail")
  reverse_match$reverse_filter_mode <- "reverse_masked_cosine_evidence"
  reverse_match$threshold_validation_status <- "reverse_masked_cosine_primary"
  reverse_match$selection_reason <- ifelse(
    reverse_match$evidence_pass,
    "n_matched_peaks_and_reverse_masked_cosine_evidence_pass",
    "below_reverse_masked_cosine_evidence_rule"
  )
  reverse_match$evidence_failure_reason <- ifelse(
    reverse_match$evidence_pass,
    "pass",
    ifelse(
      as.character(reverse_match$reverse_match_status) == "no_library_spectrum",
      "missing library spectrum",
      ifelse(
        n_matched < min_matched_peaks_evidence,
        paste0("n_matched_peaks < ", min_matched_peaks_evidence),
        "reverse_masked_cosine <= threshold"
      )
    )
  )
  reverse_match$evidence_filter_reason <- reverse_match$evidence_failure_reason

  rank_one <- function(df) {
    df <- df[order(
      -(df$evidence_pass %in% TRUE),
      -suppressWarnings(as.numeric(df$S_total)),
      -suppressWarnings(as.numeric(df$reverse_masked_cosine)),
      -suppressWarnings(as.numeric(df$n_matched_peaks)),
      df$candidate_id,
      df$spectrum_id,
      na.last = TRUE
    ), , drop = FALSE]
    df$rank_within_feature <- seq_len(nrow(df))
    score <- suppressWarnings(as.numeric(df$S_total))
    score[!is.finite(score)] <- NA_real_
    next_score <- c(score[-1L], NA_real_)
    df$rank_margin_to_next <- ifelse(is.finite(score) & is.finite(next_score), score - next_score, NA_real_)
    df
  }
  split_query <- .rm_split_key(reverse_match, c("query_id", "anchor_feature_id"))
  reverse_match <- do.call(rbind, lapply(split(reverse_match, split_query), rank_one))
  rownames(reverse_match) <- NULL

  sas_for_tier <- suppressWarnings(as.numeric(reverse_match$S_SAS))
  rank_for_tier <- suppressWarnings(as.integer(reverse_match$rank_within_feature))
  n_for_tier <- suppressWarnings(as.integer(reverse_match$n_matched_peaks))
  rmc_for_tier <- suppressWarnings(as.numeric(reverse_match$reverse_masked_cosine))
  st_for_tier <- suppressWarnings(as.numeric(reverse_match$S_total))
  tier1 <- n_for_tier >= 3L & rmc_for_tier >= 0.05 & sas_for_tier >= 0.60 &
    st_for_tier >= 0.20 & (is.na(rank_for_tier) | rank_for_tier <= 3L)
  tier2 <- !tier1 & n_for_tier >= 3L & rmc_for_tier > 0
  tier3 <- !tier1 & !tier2 & (
    (n_for_tier == 2L & rmc_for_tier >= 0.05) |
      (n_for_tier == 1L & rmc_for_tier >= 0.20 & sas_for_tier >= 0.80)
  )
  reverse_match$confidence_tier <- ifelse(
    tier1, "Tier1_pseudoMS2_spatial",
    ifelse(tier2, "Tier2_pseudoMS2_candidate",
      ifelse(tier3, "Tier3_review_only", "No_pseudoMS2_evidence"))
  )

  selected <- reverse_match[reverse_match$evidence_pass %in% TRUE, , drop = FALSE]
  before <- stats::aggregate(candidate_id ~ query_id + anchor_feature_id, data = reverse_match, FUN = length)
  names(before)[names(before) == "candidate_id"] <- "n_candidates_after_entropy"
  after <- if (nrow(selected)) {
    out <- stats::aggregate(candidate_id ~ query_id + anchor_feature_id, data = selected, FUN = length)
    names(out)[names(out) == "candidate_id"] <- "n_pseudo_ms2_candidate_evidence"
    out
  } else {
    data.frame(query_id = character(), anchor_feature_id = character(),
      n_pseudo_ms2_candidate_evidence = integer(), stringsAsFactors = FALSE,
      check.names = FALSE)
  }
  first_query <- reverse_match[!duplicated(.rm_split_key(reverse_match,
    c("query_id", "anchor_feature_id"))), c("query_id", "anchor_feature_id"), drop = FALSE]
  summary <- merge(first_query, before, by = c("query_id", "anchor_feature_id"), all.x = TRUE, sort = FALSE)
  summary <- merge(summary, after, by = c("query_id", "anchor_feature_id"), all.x = TRUE, sort = FALSE)
  summary$n_candidates_after_entropy[is.na(summary$n_candidates_after_entropy)] <- 0L
  summary$n_pseudo_ms2_candidate_evidence[is.na(summary$n_pseudo_ms2_candidate_evidence)] <- 0L
  summary$min_matched_peaks_evidence <- min_matched_peaks_evidence
  summary$min_reverse_masked_cosine_evidence <- min_reverse_masked_cosine_evidence
  summary <- summary[, colnames(.rm_empty_rdp_filter_summary()), drop = FALSE]

  matched_out <- empty_pairs
  if (nrow(selected) > 0L && nrow(matched_peak_pairs) > 0L) {
    matched_peak_pairs <- .rm_add_missing_columns(matched_peak_pairs, colnames(empty_pairs))
    keep_key <- .rm_split_key(selected, c("candidate_key", "cluster_id", "anchor_feature_id"))
    pair_key <- .rm_split_key(matched_peak_pairs, c("candidate_key", "cluster_id", "anchor_feature_id"))
    matched_out <- matched_peak_pairs[pair_key %in% keep_key, , drop = FALSE]
    matched_out <- matched_out[, colnames(empty_pairs), drop = FALSE]
  }
  selected <- .rm_add_missing_columns(selected, colnames(empty_scores))
  reverse_match <- .rm_add_missing_columns(reverse_match, colnames(empty_scores))
  selected <- selected[, colnames(empty_scores), drop = FALSE]
  reverse_match <- reverse_match[, colnames(empty_scores), drop = FALSE]
  rownames(selected) <- NULL
  rownames(summary) <- NULL
  rownames(matched_out) <- NULL
  list(
    all_candidates = reverse_match,
    pseudo_ms2_candidate_evidence = selected,
    pseudo_ms2_candidate_evidence_summary = summary,
    matched_peak_pairs = matched_out
  )
}

`60_build_archive_index` <- 
function(archive_library,
                                     ppm = 10,
                                     min_candidate_peaks = 1,
                                     ...) {
  if (is.null(archive_library)) {
    stop("archive_library must be provided when archive_index is NULL.", call. = FALSE)
  }

  archive_library <- .rm_read_archive_library(archive_library)
  archive_spectrum_table <- .rm_standardize_archive_spectrum_table(
    archive_library$archive_spectrum_table
  )
  archive_peak_table <- .rm_standardize_archive_peak_table(
    archive_library$archive_peak_table
  )

  archive_peak_table$peak_mz <- suppressWarnings(as.numeric(archive_peak_table$peak_mz))
  archive_peak_table$intensity_raw <- suppressWarnings(as.numeric(archive_peak_table$intensity_raw))
  archive_peak_table <- archive_peak_table[
    is.finite(archive_peak_table$peak_mz) &
      is.finite(archive_peak_table$intensity_raw) &
      archive_peak_table$intensity_raw > 0,
    ,
    drop = FALSE
  ]

  if (nrow(archive_peak_table) == 0L) {
    stop("No valid archive peaks remain after cleaning.", call. = FALSE)
  }

  min_candidate_peaks <- as.integer(min_candidate_peaks)
  if (!is.finite(min_candidate_peaks) || min_candidate_peaks < 1L) {
    min_candidate_peaks <- 1L
  }

  archive_peak_table$intensity_sum1 <- NA_real_
  valid_spectra <- character()
  split_rows <- split(seq_len(nrow(archive_peak_table)), as.character(archive_peak_table$spectrum_id))
  for (sid in names(split_rows)) {
    idx <- split_rows[[sid]]
    ints <- as.numeric(archive_peak_table$intensity_raw[idx])
    if (length(idx) < min_candidate_peaks) {
      next
    }
    total <- sum(ints, na.rm = TRUE)
    if (!is.finite(total) || total <= 0) {
      next
    }
    archive_peak_table$intensity_sum1[idx] <- ints / total
    valid_spectra <- c(valid_spectra, sid)
  }

  archive_peak_table <- archive_peak_table[
    as.character(archive_peak_table$spectrum_id) %in% valid_spectra,
    ,
    drop = FALSE
  ]
  archive_spectrum_table <- archive_spectrum_table[
    as.character(archive_spectrum_table$spectrum_id) %in% valid_spectra,
    ,
    drop = FALSE
  ]

  if (nrow(archive_peak_table) == 0L) {
    stop("No archive spectra meet min_candidate_peaks.", call. = FALSE)
  }

  archive_map <- stats::setNames(
    as.character(archive_spectrum_table$archive_id),
    as.character(archive_spectrum_table$spectrum_id)
  )
  archive_peak_table$archive_id <- unname(archive_map[as.character(archive_peak_table$spectrum_id)])
  missing_archive <- is.na(archive_peak_table$archive_id) |
    !nzchar(as.character(archive_peak_table$archive_id))
  archive_peak_table$archive_id[missing_archive] <- as.character(archive_peak_table$spectrum_id[missing_archive])

  precursor_index <- archive_spectrum_table[, c(
    "archive_id", "spectrum_id", "candidate_id", "precursor_mz"
  ), drop = FALSE]
  precursor_index$precursor_mz <- suppressWarnings(as.numeric(precursor_index$precursor_mz))
  precursor_index <- precursor_index[is.finite(precursor_index$precursor_mz), , drop = FALSE]

  fragment_index <- archive_peak_table[, c(
    "archive_id", "spectrum_id", "candidate_id", "peak_mz", "intensity_sum1"
  ), drop = FALSE]
  fragment_index$peak_mz <- suppressWarnings(as.numeric(fragment_index$peak_mz))
  fragment_index$intensity_sum1 <- suppressWarnings(as.numeric(fragment_index$intensity_sum1))
  fragment_index <- fragment_index[
    is.finite(fragment_index$peak_mz) &
      is.finite(fragment_index$intensity_sum1) &
      fragment_index$intensity_sum1 > 0,
    ,
    drop = FALSE
  ]
  fragment_index <- fragment_index[order(fragment_index$peak_mz, fragment_index$spectrum_id), , drop = FALSE]
  fragment_index$archive_peak_id <- seq_len(nrow(fragment_index))

  archive_spectrum_table <- archive_spectrum_table[, colnames(.rm_empty_archive_spectrum_table()), drop = FALSE]
  archive_peak_table <- archive_peak_table[, colnames(.rm_empty_archive_peak_table()), drop = FALSE]
  rownames(archive_spectrum_table) <- NULL
  rownames(archive_peak_table) <- NULL
  rownames(precursor_index) <- NULL
  rownames(fragment_index) <- NULL

  list(
    archive_spectrum_table = archive_spectrum_table,
    archive_peak_table = archive_peak_table,
    precursor_index = precursor_index,
    fragment_index = fragment_index
  )
}

`60_calculate_library_coverage` <- 
function(library_intensity, matched) {
  l <- suppressWarnings(as.numeric(library_intensity))
  matched <- matched %in% TRUE
  if (length(l) == 0L) {
    return(list(S_lib_coverage = 0, matched_peak_fraction = 0, n_reverse_match = 0L))
  }
  l[!is.finite(l)] <- 0
  total <- sum(l, na.rm = TRUE)
  coverage <- if (is.finite(total) && total > 0) sum(l[matched], na.rm = TRUE) / total else 0
  list(
    S_lib_coverage = pmin(pmax(as.numeric(coverage), 0), 1),
    matched_peak_fraction = sum(matched, na.rm = TRUE) / length(l),
    n_reverse_match = as.integer(sum(matched, na.rm = TRUE))
  )
}

`60_calculate_reverse_cosine` <- 
function(library_intensity, matched_query_intensity) {
  l <- suppressWarnings(as.numeric(library_intensity))
  q <- suppressWarnings(as.numeric(matched_query_intensity))
  if (length(l) == 0L) {
    return(NA_real_)
  }
  if (length(q) != length(l)) {
    stop("library_intensity and matched_query_intensity must have the same length.", call. = FALSE)
  }
  l[!is.finite(l)] <- 0
  q[!is.finite(q)] <- 0
  if (all(q == 0)) {
    return(0)
  }
  denom <- sqrt(sum(l^2, na.rm = TRUE) * sum(q^2, na.rm = TRUE))
  if (!is.finite(denom) || denom <= 0) {
    return(0)
  }
  pmin(pmax(sum(l * q, na.rm = TRUE) / denom, 0), 1)
}

`60_calculate_spatial_reverse_score` <- 
function(SAS, library_intensity, matched) {
  matched <- matched %in% TRUE
  if (!any(matched)) {
    return(NA_real_)
  }
  sas <- suppressWarnings(as.numeric(SAS))[matched]
  w <- suppressWarnings(as.numeric(library_intensity))[matched]
  if (length(sas) == 0L || all(is.na(sas))) {
    return(NA_real_)
  }
  keep <- is.finite(sas) & is.finite(w) & w > 0
  if (!any(keep)) {
    return(NA_real_)
  }
  as.numeric(stats::weighted.mean(sas[keep], w = w[keep]))
}

`60_entropy_candidates_for_reverse` <- 
function(object) {
  entropy_candidates <- NULL
  if (is.data.frame(object@ms2_data$entropy_candidates)) {
    entropy_candidates <- object@ms2_data$entropy_candidates
  } else if (is.data.frame(object@annotation_results$entropy_candidates)) {
    entropy_candidates <- object@annotation_results$entropy_candidates
  }
  flash_hits <- object@annotation_results$flash_entropy_hits
  if (is.null(entropy_candidates) || nrow(entropy_candidates) == 0L || is.null(flash_hits)) {
    return(flash_hits)
  }
  entropy_candidates <- as.data.frame(entropy_candidates, stringsAsFactors = FALSE, check.names = FALSE)
  flash_hits <- as.data.frame(flash_hits, stringsAsFactors = FALSE, check.names = FALSE)
  entropy_candidates <- .rm_add_missing_columns(entropy_candidates, c("query_id", "anchor_feature_id", "candidate_id", "spectrum_id", "rank_entropy"))
  flash_hits <- .rm_add_missing_columns(flash_hits, c("query_id", "anchor_feature_id", "candidate_id", "spectrum_id", "rank_entropy", "entropy_S"))
  keep_key <- .rm_split_key(entropy_candidates, c("query_id", "anchor_feature_id", "candidate_id", "spectrum_id"))
  flash_key <- .rm_split_key(flash_hits, c("query_id", "anchor_feature_id", "candidate_id", "spectrum_id"))
  keep <- flash_key %in% keep_key
  out <- flash_hits[keep, , drop = FALSE]
  if (nrow(out) == 0L) {
    return(entropy_candidates)
  }
  rank_lookup <- stats::setNames(suppressWarnings(as.integer(entropy_candidates$rank_entropy)), keep_key)
  out$rank_entropy <- suppressWarnings(as.integer(rank_lookup[flash_key[keep]]))
  out
}

`60_get_library_spectrum` <- 
function(candidate,
                                      object,
                                      spectrum_id_peak_index = NULL,
                                      spectrum_index_telemetry = NULL) {
  empty <- data.frame(
    library_peak_id = character(),
    library_mz = numeric(),
    library_fragment_mz = numeric(),
    library_intensity = numeric(),
    library_intensity_raw = numeric(),
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
  candidate <- as.data.frame(candidate, stringsAsFactors = FALSE, check.names = FALSE)
  source <- as.character(candidate$candidate_source[[1]])

  if (source %in% c("flash", "both")) {
    peaks <- object@spectral_library$library_peak_table
    if (!is.null(peaks)) {
      peaks <- as.data.frame(peaks, stringsAsFactors = FALSE, check.names = FALSE)
      peaks <- .rm_add_missing_columns(peaks, c("spectrum_id", "candidate_id", "peak_mz", "intensity_sum1", "intensity_raw"))
      if (!"intensity_sum1" %in% names(peaks) || all(is.na(peaks$intensity_sum1))) {
        peaks$intensity_sum1 <- suppressWarnings(as.numeric(peaks$intensity_raw)) /
          ave(suppressWarnings(as.numeric(peaks$intensity_raw)), peaks$spectrum_id, FUN = function(x) sum(x, na.rm = TRUE))
      }
      spectrum_id <- as.character(candidate$spectrum_id[[1]])
      ref <- peaks[0, , drop = FALSE]
      indexed_rows <- if (!is.null(spectrum_id_peak_index)) {
        .rm_spectrum_id_index_rows(spectrum_id_peak_index, spectrum_id)
      } else {
        integer()
      }
      if (length(indexed_rows)) {
        ref <- peaks[indexed_rows, , drop = FALSE]
        if (!is.null(spectrum_index_telemetry)) {
          spectrum_index_telemetry$direct_accesses <- spectrum_index_telemetry$direct_accesses + 1L
        }
      } else {
        if (!is.null(spectrum_index_telemetry)) {
          spectrum_index_telemetry$fallback_linear_scans <- spectrum_index_telemetry$fallback_linear_scans + 1L
        }
        ref <- peaks[as.character(peaks$spectrum_id) == spectrum_id, , drop = FALSE]
      }
      if (nrow(ref) == 0L) {
        ref <- peaks[as.character(peaks$candidate_id) == as.character(candidate$candidate_id[[1]]), , drop = FALSE]
      }
      if (nrow(ref) > 0L) {
        out <- data.frame(
          library_peak_id = if ("library_peak_id" %in% names(ref)) as.character(ref$library_peak_id) else paste0(as.character(candidate$spectrum_id[[1]]), "_", seq_len(nrow(ref))),
          library_mz = suppressWarnings(as.numeric(ref$peak_mz)),
          library_intensity = suppressWarnings(as.numeric(ref$intensity_sum1)),
          library_intensity_raw = suppressWarnings(as.numeric(ref$intensity_raw)),
          stringsAsFactors = FALSE,
          check.names = FALSE
        )
        out <- out[is.finite(out$library_mz) & is.finite(out$library_intensity) & out$library_intensity > 0, , drop = FALSE]
        return(.rm_normalize_library_spectrum(out))
      }
    }
  }

  if (source %in% c("archive", "both") && length(object@archive_index) > 0L) {
    peaks <- object@archive_index$archive_peak_table
    if (!is.null(peaks)) {
      peaks <- as.data.frame(peaks, stringsAsFactors = FALSE, check.names = FALSE)
      peaks <- .rm_add_missing_columns(peaks, c("archive_id", "spectrum_id", "candidate_id", "peak_mz", "intensity_sum1", "intensity_raw"))
      ref <- peaks[as.character(peaks$archive_id) == as.character(candidate$archive_id[[1]]), , drop = FALSE]
      if (nrow(ref) == 0L) {
        ref <- peaks[as.character(peaks$spectrum_id) == as.character(candidate$spectrum_id[[1]]), , drop = FALSE]
      }
      if (nrow(ref) > 0L) {
        out <- data.frame(
          library_peak_id = paste0(as.character(candidate$spectrum_id[[1]]), "_", seq_len(nrow(ref))),
          library_mz = suppressWarnings(as.numeric(ref$peak_mz)),
          library_intensity = suppressWarnings(as.numeric(ref$intensity_sum1)),
          library_intensity_raw = suppressWarnings(as.numeric(ref$intensity_raw)),
          stringsAsFactors = FALSE,
          check.names = FALSE
        )
        out <- out[is.finite(out$library_mz) & is.finite(out$library_intensity) & out$library_intensity > 0, , drop = FALSE]
        return(.rm_normalize_library_spectrum(out))
      }
    }
  }

  empty
}

`60_match_library_to_query` <- 
function(library_peaks,
                                        query_peaks,
                                        tolerance_ppm = 10,
                                        tolerance_da = NULL,
                                        tolerance_mode = "max",
                                        low_mz_cutoff = 400,
                                        fixed_da_low_mz = 0.005,
                                        mz_bin_width = 0.01,
                                        bin_margin = 0.5) {
  empty <- .rm_empty_reverse_matched_peaks()
  lib <- as.data.frame(library_peaks, stringsAsFactors = FALSE, check.names = FALSE)
  query <- as.data.frame(query_peaks, stringsAsFactors = FALSE, check.names = FALSE)
  if (nrow(lib) == 0L) {
    return(empty)
  }
  if (!is.null(tolerance_da)) {
    warning(
      "tolerance_da is deprecated for 60_match_library_to_query(); using it as ",
      "fixed_da_low_mz for backward compatibility. The canonical effective ",
      "tolerance is max(library_mz * tolerance_ppm / 1e6, fixed_da_low_mz) ",
      "only below low_mz_cutoff.",
      call. = FALSE
    )
    fixed_da_low_mz <- tolerance_da
  }
  lib <- .rm_add_missing_columns(lib, c(
    "library_peak_id", "library_mz", "library_intensity", "library_intensity_raw"
  ))
  query <- .rm_add_missing_columns(query, colnames(.rm_empty_query_spectrum()))
  lib$library_intensity_raw <- suppressWarnings(as.numeric(lib$library_intensity_raw))
  use_norm_as_raw <- !is.finite(lib$library_intensity_raw)
  lib$library_intensity_raw[use_norm_as_raw] <- suppressWarnings(as.numeric(lib$library_intensity[use_norm_as_raw]))
  query$query_intensity_raw <- suppressWarnings(as.numeric(query$query_intensity_raw))
  use_query_norm_as_raw <- !is.finite(query$query_intensity_raw)
  query$query_intensity_raw[use_query_norm_as_raw] <- suppressWarnings(as.numeric(query$query_intensity[use_query_norm_as_raw]))
  score <- compute_reverse_masked_cosine(
    library_peaks = lib,
    query_peaks = query,
    tolerance_ppm = tolerance_ppm,
    tolerance_da_lowmz = fixed_da_low_mz,
    lowmz_threshold = low_mz_cutoff,
    one_to_one = TRUE
  )
  matched <- score$matched_peaks
  out <- data.frame(
    library_peak_order = as.integer(matched$library_peak_order),
    library_peak_id = as.character(matched$library_peak_id),
    library_mz = as.numeric(matched$library_mz),
    library_fragment_mz = as.numeric(matched$library_mz),
    library_intensity = as.numeric(matched$library_intensity_norm),
    library_fragment_intensity = as.numeric(matched$library_intensity_norm),
    library_intensity_raw = as.numeric(matched$library_intensity_raw),
    library_intensity_norm = as.numeric(matched$library_intensity_norm),
    library_relative_intensity = as.numeric(matched$library_intensity_norm),
    query_peak_id = as.character(matched$query_peak_id),
    query_fragment_feature_id = as.character(matched$query_peak_id),
    query_mz = as.numeric(matched$query_mz),
    query_fragment_mz = as.numeric(matched$query_mz),
    query_intensity = as.numeric(matched$query_intensity_norm),
    query_intensity_raw = as.numeric(matched$query_intensity_raw),
    query_intensity_norm = as.numeric(matched$query_intensity_norm),
    query_relative_intensity = as.numeric(matched$query_intensity_norm),
    mz_error = as.numeric(matched$mz_error_da),
    mz_error_ppm = as.numeric(matched$mz_error_ppm),
    SAS = as.numeric(matched$SAS),
    query_sas = as.numeric(matched$SAS),
    contribution_to_dot = as.numeric(matched$contribution_to_dot),
    contribution_to_S_SAS = as.numeric(matched$contribution_to_S_SAS),
    query_full_l2 = rep(score$query_full_l2, nrow(matched)),
    query_matched_l2 = rep(score$query_masked_l2, nrow(matched)),
    library_all_l2 = rep(score$library_all_l2, nrow(matched)),
    query_total_intensity = rep(score$total_query_intensity, nrow(matched)),
    matched_query_intensity = rep(score$matched_query_intensity, nrow(matched)),
    matched_library_intensity = rep(score$matched_library_intensity, nrow(matched)),
    total_library_intensity = rep(score$total_library_intensity, nrow(matched)),
    matched_query_intensity_fraction = rep(score$matched_query_intensity_fraction, nrow(matched)),
    matched_library_intensity_fraction = rep(score$matched_library_intensity_fraction, nrow(matched)),
    library_l2_coverage = rep(score$library_l2_coverage, nrow(matched)),
    reverse_masked_cosine = rep(score$reverse_masked_cosine, nrow(matched)),
    match_id = NA_character_,
    n_query_peaks_total = as.integer(nrow(score$query_peaks)),
    n_library_peaks_total = as.integer(nrow(matched)),
    matched = matched$matched %in% TRUE,
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
  matched_idx <- which(out$matched %in% TRUE)
  if (length(matched_idx) > 0L) {
    out$match_id[matched_idx] <- paste0("M", seq_along(matched_idx))
  }
  out <- out[, colnames(empty), drop = FALSE]
  rownames(out) <- NULL
  out
}

`60_prepare_reverse_candidates` <- 
function(flash_hits = NULL,
                                            archive_hits = NULL,
                                            use_flash_hits = TRUE,
                                            use_archive_hits = TRUE) {
  empty <- .rm_empty_merged_candidates()
  pieces <- list()

  if (isTRUE(use_flash_hits) && !is.null(flash_hits) && nrow(flash_hits) > 0L) {
    flash <- as.data.frame(flash_hits, stringsAsFactors = FALSE, check.names = FALSE)
    flash <- .rm_add_missing_columns(flash, colnames(empty))
    flash$entropy_S <- suppressWarnings(as.numeric(flash$entropy_S))
    flash$S_FE <- suppressWarnings(as.numeric(flash$S_FE))
    fill_entropy <- !is.finite(flash$entropy_S) & is.finite(flash$S_FE)
    flash$entropy_S[fill_entropy] <- flash$S_FE[fill_entropy]
    fill_s_fe <- !is.finite(flash$S_FE) & is.finite(flash$entropy_S)
    flash$S_FE[fill_s_fe] <- flash$entropy_S[fill_s_fe]
    flash$s_fe <- flash$entropy_S
    if ("rank_spectrum_within_anchor" %in% names(flash)) flash$retrieval_rank <- suppressWarnings(as.integer(flash$rank_spectrum_within_anchor))
    if ("rank_entropy" %in% names(flash) && any(!is.na(flash$rank_entropy))) flash$retrieval_rank <- suppressWarnings(as.integer(flash$rank_entropy))
    if ("n_match_fe" %in% names(flash)) flash$n_FE_match <- suppressWarnings(as.integer(flash$n_match_fe))
    if ("matched_query_peak_n" %in% names(flash)) flash$matched_query_peaks <- suppressWarnings(as.integer(flash$matched_query_peak_n))
    if ("matched_library_peak_n" %in% names(flash)) flash$matched_library_peaks <- suppressWarnings(as.integer(flash$matched_library_peak_n))
    if (!"query_id" %in% names(flash) || all(is.na(flash$query_id))) flash$query_id <- as.character(flash$peak_id)
    flash$archive_id <- NA_character_
    flash$s_arch <- NA_real_
    flash$candidate_source <- "flash"
    pieces[[length(pieces) + 1L]] <- flash[, colnames(empty), drop = FALSE]
  }

  if (isTRUE(use_archive_hits) && !is.null(archive_hits) && nrow(archive_hits) > 0L) {
    archive <- as.data.frame(archive_hits, stringsAsFactors = FALSE, check.names = FALSE)
    archive <- .rm_add_missing_columns(archive, colnames(empty))
    archive$entropy_S <- 0
    archive$S_FE <- 0
    archive$s_fe <- 0
    archive$retrieval_rank <- NA_integer_
    archive$n_FE_match <- NA_integer_
    archive$matched_query_peaks <- archive$matched_query_peak_n_arch
    archive$matched_library_peaks <- archive$matched_archive_peak_n
    archive$query_id <- as.character(archive$cluster_id)
    archive$candidate_source <- "archive"
    pieces[[length(pieces) + 1L]] <- archive[, colnames(empty), drop = FALSE]
  }

  if (length(pieces) == 0L) {
    return(empty)
  }

  out <- do.call(rbind, pieces)
  out$library_source[is.na(out$library_source) | !nzchar(as.character(out$library_source))] <- "unknown"
  if (!"adduct_status" %in% names(out)) out$adduct_status <- NA_character_
  out$adduct_status <- ifelse(
    !is.na(trimws(as.character(out$adduct))) & nzchar(trimws(as.character(out$adduct))),
    "explicit", "unknown"
  )
  out$candidate_retrieved <- TRUE
  out$precursor_gate_pass <- TRUE
  out$fragment_retrieval_hits <- suppressWarnings(as.integer(out$n_FE_match))
  candidate_id <- as.character(out$candidate_id)
  spectrum_id <- as.character(out$spectrum_id)
  candidate_id[is.na(candidate_id) | !nzchar(candidate_id)] <- "NA_candidate"
  spectrum_id[is.na(spectrum_id) | !nzchar(spectrum_id)] <- "NA_spectrum"
  out$candidate_key <- paste(out$library_source, candidate_id, spectrum_id, sep = "::")
  split_key <- .rm_split_key(out, c("cluster_id", "anchor_feature_id", "candidate_key"))
  out <- do.call(rbind, lapply(split(out, split_key), .rm_merge_candidate_rows))
  out <- out[order(out$cluster_id, out$anchor_feature_id, out$retrieval_rank, -out$entropy_S, out$candidate_key), , drop = FALSE]
  out <- out[, colnames(empty), drop = FALSE]
  rownames(out) <- NULL
  out
}

`60_reverse_match_library_peaks` <- 
function(object,
                                          candidates,
                                          fragment_tolerance_ppm = 10,
                                          fragment_tolerance_da = NULL,
                                          fragment_tolerance_mode = "max",
                                          low_mz_cutoff = 400,
                                          fixed_da_low_mz = 0.005,
                                          mz_bin_width = 0.01,
                                          bin_margin = 0.5,
                                          workers = 1L,
                                          progress_every = NULL,
                                          progress_path = NULL,
                                          progress_started_at = NULL,
                                          progress_callback = NULL,
                                          spectrum_id_peak_index = NULL,
                                          spectrum_index_telemetry = NULL) {
  empty <- .rm_empty_reverse_matched_peaks_full()
  candidates <- as.data.frame(candidates, stringsAsFactors = FALSE, check.names = FALSE)
  if (nrow(candidates) == 0L) {
    return(empty)
  }
  pseudo <- `60_standardize_reverse_query_spectrum`(object@pseudo_ms2$pseudo_ms2_table)
  pseudo_split <- split(pseudo, .rm_split_key(pseudo, c("cluster_id", "anchor_feature_id")), drop = TRUE)
  if (is.null(spectrum_id_peak_index)) {
    spectrum_id_peak_index <- .rm_build_spectrum_id_peak_index(object@spectral_library$library_peak_table)
  }
  match_one <- function(i) {
    candidate <- candidates[i, , drop = FALSE]
    query_key <- .rm_split_key(candidate, c("cluster_id", "anchor_feature_id"))
    query <- pseudo_split[[query_key]]
    if (is.null(query)) {
      query <- pseudo[0, , drop = FALSE]
    }
    lib <- `60_get_library_spectrum`(
      candidate,
      object,
      spectrum_id_peak_index = spectrum_id_peak_index,
      spectrum_index_telemetry = spectrum_index_telemetry
    )
    if (nrow(lib) == 0L) return(NULL)
    matched <- `60_match_library_to_query`(
      library_peaks = lib,
      query_peaks = query,
      tolerance_ppm = fragment_tolerance_ppm,
      tolerance_da = fragment_tolerance_da,
      tolerance_mode = fragment_tolerance_mode,
      low_mz_cutoff = low_mz_cutoff,
      fixed_da_low_mz = fixed_da_low_mz,
      mz_bin_width = mz_bin_width,
      bin_margin = bin_margin
    )
    if (nrow(matched) == 0L) return(NULL)
    matched$anchor_feature_id <- as.character(candidate$anchor_feature_id[[1]])
    matched$query_id <- as.character(candidate$query_id[[1]])
    matched$candidate_id <- as.character(candidate$candidate_id[[1]])
    matched$spectrum_id <- as.character(candidate$spectrum_id[[1]])
    matched$candidate_key <- as.character(candidate$candidate_key[[1]])
    matched$compound_name <- as.character(candidate$candidate_name[[1]])
    matched$candidate_name <- as.character(candidate$candidate_name[[1]])
    matched$cluster_id <- candidate$cluster_id[[1]]
    matched[, colnames(empty), drop = FALSE]
  }
  workers <- max(1L, as.integer(workers[[1]]))
  progress_started_at <- if (inherits(progress_started_at, "POSIXt")) progress_started_at else Sys.time()
  progress_every <- if (is.null(progress_every) || length(progress_every) == 0L) {
    nrow(candidates)
  } else {
    suppressWarnings(as.integer(progress_every[[1]]))
  }
  if (!is.finite(progress_every) || progress_every < 1L) {
    progress_every <- nrow(candidates)
  }
  candidate_idx <- seq_len(nrow(candidates))
  chunks <- split(candidate_idx, ceiling(seq_along(candidate_idx) / progress_every), drop = TRUE)
  rows <- vector("list", 0L)
  reverse_rows_so_far <- 0L
  processed_candidates <- 0L
  for (chunk in chunks) {
    chunk_rows <- if (workers > 1L && .Platform$OS.type != "windows" && length(chunk) > 1L) {
      parallel::mclapply(chunk, match_one, mc.cores = min(workers, length(chunk)))
    } else {
      lapply(chunk, match_one)
    }
    chunk_rows <- Filter(Negate(is.null), chunk_rows)
    if (length(chunk_rows)) {
      rows <- c(rows, chunk_rows)
      reverse_rows_so_far <- reverse_rows_so_far + sum(vapply(chunk_rows, nrow, integer(1)))
    }
    processed_candidates <- processed_candidates + length(chunk)
    .rm_emit_progress(
      progress_callback = progress_callback,
      progress_path = progress_path,
      started_at = progress_started_at,
      processed_candidates = processed_candidates,
      total_candidates = nrow(candidates),
      reverse_rows_so_far = reverse_rows_so_far,
      pass_rows_so_far = 0L,
      message = "matching library peaks"
    )
  }
  rows <- Filter(Negate(is.null), rows)
  if (length(rows) == 0L) {
    return(empty)
  }
  out <- do.call(rbind, rows)
  out <- out[, colnames(empty), drop = FALSE]
  rownames(out) <- NULL
  out
}

`60_search_archive_candidates` <- 
function(object,
                                           archive_index,
                                           fragment_tolerance_ppm = 10,
                                           fragment_tolerance_da = NULL,
                                           fragment_tolerance_mode = "max",
                                           low_mz_cutoff = 400,
                                           fixed_da_low_mz = 0.005,
                                           mz_bin_width = 0.01,
                                           bin_margin = 0.5,
                                           precursor_tol = 0.05,
                                           delta_grid = c(0),
                                           archive_top_n = 200,
                                           tolerance_contract = NULL) {
  empty <- .rm_empty_archive_hits()
  if (is.null(archive_index) || length(archive_index) == 0L) {
    return(empty)
  }

  pseudo <- `60_standardize_reverse_query_spectrum`(object@pseudo_ms2$pseudo_ms2_table)
  precursor_index <- as.data.frame(archive_index$precursor_index, stringsAsFactors = FALSE, check.names = FALSE)
  fragment_index <- as.data.frame(archive_index$fragment_index, stringsAsFactors = FALSE, check.names = FALSE)
  spectrum <- as.data.frame(archive_index$archive_spectrum_table, stringsAsFactors = FALSE, check.names = FALSE)
  if (nrow(pseudo) == 0L || nrow(precursor_index) == 0L || nrow(fragment_index) == 0L) {
    return(empty)
  }

  precursor_index$precursor_mz <- suppressWarnings(as.numeric(precursor_index$precursor_mz))
  fragment_index$peak_mz <- suppressWarnings(as.numeric(fragment_index$peak_mz))
  fragment_index$intensity_sum1 <- suppressWarnings(as.numeric(fragment_index$intensity_sum1))

  anchors <- unique(pseudo[, c("cluster_id", "anchor_feature_id", "anchor_mz"), drop = FALSE])
  rows <- list()
  for (i in seq_len(nrow(anchors))) {
    anchor <- anchors[i, , drop = FALSE]
    query <- pseudo[
      as.character(pseudo$cluster_id) == as.character(anchor$cluster_id) &
        as.character(pseudo$anchor_feature_id) == as.character(anchor$anchor_feature_id),
      ,
      drop = FALSE
    ]
    query <- query[is.finite(query$query_mz) & is.finite(query$query_intensity) & query$query_intensity > 0, , drop = FALSE]
    if (nrow(query) == 0L) next

    for (delta in delta_grid) {
      expected_precursor <- as.numeric(anchor$anchor_mz) + as.numeric(delta)
      precursor_tolerance <- if (!is.null(tolerance_contract)) {
        .msi_effective_precursor_tolerance(expected_precursor, tolerance_contract)
      } else {
        as.numeric(precursor_tol)
      }
      keep <- is.finite(precursor_index$precursor_mz) &
        is.finite(precursor_tolerance) &
        abs(precursor_index$precursor_mz - expected_precursor) <= precursor_tolerance
      candidates <- precursor_index[keep, , drop = FALSE]
      if (nrow(candidates) == 0L) next

      for (j in seq_len(nrow(candidates))) {
        ref <- fragment_index[
          as.character(fragment_index$spectrum_id) == as.character(candidates$spectrum_id[[j]]),
          ,
          drop = FALSE
        ]
        if (nrow(ref) == 0L) next
        ref <- data.frame(
          library_peak_id = seq_len(nrow(ref)),
          library_mz = ref$peak_mz,
          library_intensity = ref$intensity_sum1,
          stringsAsFactors = FALSE,
          check.names = FALSE
        )
        matches <- `60_match_library_to_query`(
          library_peaks = ref,
          query_peaks = query,
          tolerance_ppm = fragment_tolerance_ppm,
          tolerance_da = fragment_tolerance_da,
          tolerance_mode = fragment_tolerance_mode,
          low_mz_cutoff = low_mz_cutoff,
          fixed_da_low_mz = fixed_da_low_mz,
          mz_bin_width = mz_bin_width,
          bin_margin = bin_margin
        )
        matched <- matches[matches$matched %in% TRUE, , drop = FALSE]
        if (nrow(matched) == 0L) next
        rows[[length(rows) + 1L]] <- data.frame(
          cluster_id = anchor$cluster_id[[1]],
          anchor_feature_id = as.character(anchor$anchor_feature_id[[1]]),
          anchor_mz = as.numeric(anchor$anchor_mz[[1]]),
          delta = as.numeric(delta),
          expected_precursor = expected_precursor,
          archive_id = as.character(candidates$archive_id[[j]]),
          spectrum_id = as.character(candidates$spectrum_id[[j]]),
          candidate_id = as.character(candidates$candidate_id[[j]]),
          s_arch = `60_calculate_reverse_cosine`(ref$library_intensity, matches$query_intensity),
          n_match_arch = nrow(matched),
          matched_query_peak_n_arch = length(unique(matched$query_peak_id)),
          matched_archive_peak_n = length(unique(matched$library_peak_id)),
          stringsAsFactors = FALSE,
          check.names = FALSE
        )
      }
    }
  }

  if (length(rows) == 0L) {
    return(empty)
  }

  out <- do.call(rbind, rows)
  meta_cols <- c(
    "archive_id", "spectrum_id", "candidate_id", "candidate_name", "library_source",
    "precursor_mz", "formula", "inchikey", "smiles", "database_source", "adduct", "polarity"
  )
  spectrum <- .rm_add_missing_columns(spectrum, meta_cols)
  out <- merge(out, spectrum[, meta_cols, drop = FALSE],
    by = c("archive_id", "spectrum_id", "candidate_id"),
    all.x = TRUE,
    sort = FALSE
  )
  out <- out[order(out$cluster_id, out$anchor_feature_id, -out$s_arch), , drop = FALSE]
  archive_top_n <- as.integer(archive_top_n)
  if (!is.finite(archive_top_n) || archive_top_n < 1L) {
    archive_top_n <- nrow(out)
  }
  split_anchor <- .rm_split_key(out, c("cluster_id", "anchor_feature_id"))
  keep <- unlist(lapply(split(seq_len(nrow(out)), split_anchor), function(idx) {
    idx[seq_len(min(length(idx), archive_top_n))]
  }), use.names = FALSE)
  out <- out[sort(keep), , drop = FALSE]
  out <- .rm_add_missing_columns(out, colnames(empty))
  out <- out[, colnames(empty), drop = FALSE]
  rownames(out) <- NULL
  out
}

`60_standardize_reverse_query_spectrum` <- 
function(query) {
  empty <- .rm_empty_query_spectrum()
  query <- as.data.frame(query, stringsAsFactors = FALSE, check.names = FALSE)
  if (nrow(query) == 0L) {
    return(empty)
  }
  query <- .rm_add_missing_columns(query, c(
    "cluster_id", "anchor_feature_id", "anchor_mz", "peak_feature_id",
    "fragment_feature_id", "peak_mz", "fragment_mz", "intensity", "SAS",
    "graph_weight", "peak_id", "fragment_id", "pseudo_raw_intensity",
    "is_self_peak", "peak_role"
  ))
  query <- query[!(query$is_self_peak %in% TRUE) &
    !(tolower(as.character(query$peak_role)) %in% "anchor_self"), , drop = FALSE]
  query_mz <- ifelse(is.finite(suppressWarnings(as.numeric(query$peak_mz))),
    suppressWarnings(as.numeric(query$peak_mz)),
    suppressWarnings(as.numeric(query$fragment_mz))
  )
  query_peak_id <- as.character(query$peak_feature_id)
  missing_id <- is.na(query_peak_id) | !nzchar(query_peak_id)
  query_peak_id[missing_id] <- as.character(query$fragment_feature_id[missing_id])
  missing_id <- is.na(query_peak_id) | !nzchar(query_peak_id)
  query_peak_id[missing_id] <- as.character(query$fragment_id[missing_id])
  missing_id <- is.na(query_peak_id) | !nzchar(query_peak_id)
  query_peak_id[missing_id] <- paste0("query_peak_", which(missing_id))
  sas <- suppressWarnings(as.numeric(query$SAS))
  graph <- suppressWarnings(as.numeric(query$graph_weight))
  sas[!is.finite(sas) & is.finite(graph)] <- graph[!is.finite(sas) & is.finite(graph)]

  out <- data.frame(
    cluster_id = query$cluster_id,
    anchor_feature_id = as.character(query$anchor_feature_id),
    anchor_mz = suppressWarnings(as.numeric(query$anchor_mz)),
    query_peak_id = query_peak_id,
    query_mz = query_mz,
    query_intensity = suppressWarnings(as.numeric(query$intensity)),
    query_intensity_raw = ifelse(
      is.finite(suppressWarnings(as.numeric(query$pseudo_raw_intensity))),
      suppressWarnings(as.numeric(query$pseudo_raw_intensity)),
      suppressWarnings(as.numeric(query$intensity))
    ),
    SAS = sas,
    graph_support = graph,
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
  out <- out[is.finite(out$query_mz) & is.finite(out$query_intensity) & out$query_intensity > 0, , drop = FALSE]
  out <- out[, colnames(empty), drop = FALSE]
  rownames(out) <- NULL
  out
}

`60_store_reverse_match_results` <- 
function(object,
                                             archive_index,
                                             archive_hits,
                                             merged_candidates,
                                             reverse_matched_peaks_all,
                                             reverse_matched_peaks,
                                             reverse_match,
                                             reverse_match_all,
                                             pseudo_ms2_candidate_evidence_summary = NULL,
                                             parameters = list(),
                                             matched_peak_pairs_summary = NULL,
                                             store_matched_peak_pairs_in_rds = TRUE,
                                             ...) {
  legacy_args <- list(...)
  if (is.null(pseudo_ms2_candidate_evidence_summary) &&
      is.data.frame(legacy_args$rdp_filter_summary)) {
    pseudo_ms2_candidate_evidence_summary <- legacy_args$rdp_filter_summary
  }
  if (is.null(pseudo_ms2_candidate_evidence_summary)) {
    pseudo_ms2_candidate_evidence_summary <- .rm_empty_rdp_filter_summary()
  }
  store_matched_peak_pairs_in_rds <- .rm_bool(store_matched_peak_pairs_in_rds, TRUE)
  if (is.null(matched_peak_pairs_summary) || !is.data.frame(matched_peak_pairs_summary)) {
    matched_peak_pairs_summary <- .rm_build_matched_peak_pairs_summary(
      matched_peak_pairs = reverse_matched_peaks,
      matched_peak_pairs_all = reverse_matched_peaks_all,
      reverse_match = reverse_match,
      reverse_match_all = reverse_match_all,
      store_matched_peak_pairs_in_rds = store_matched_peak_pairs_in_rds
    )
  }
  object@archive_index <- archive_index
  object@annotation_results$archive_hits <- archive_hits
  object@annotation_results$merged_candidates <- merged_candidates
  object@annotation_results$reverse_match_all <- reverse_match_all
  object@annotation_results$reverse_match <- reverse_match
  object@annotation_results$pseudo_ms2_candidate_evidence_all <- reverse_match_all
  object@annotation_results$pseudo_ms2_candidate_evidence <- reverse_match
  object@annotation_results$pseudo_ms2_candidate_evidence_summary <- pseudo_ms2_candidate_evidence_summary
  object@annotation_results$matched_peak_pairs_summary <- matched_peak_pairs_summary
  object@annotation_results$rdp_candidates <- NULL
  object@annotation_results$rdp_filter_summary <- NULL
  object@ms2_data$pseudo_ms2_candidate_evidence_all <- reverse_match_all
  object@ms2_data$pseudo_ms2_candidate_evidence <- reverse_match
  object@ms2_data$pseudo_ms2_candidate_evidence_summary <- pseudo_ms2_candidate_evidence_summary
  object@ms2_data$matched_peak_pairs_summary <- matched_peak_pairs_summary
  object@ms2_data$rdp_candidates <- NULL
  object@ms2_data$rdp_filter_summary <- NULL
  if (isTRUE(store_matched_peak_pairs_in_rds)) {
    object@annotation_results$reverse_matched_peaks_all <- reverse_matched_peaks_all
    object@annotation_results$reverse_matched_peaks <- reverse_matched_peaks
    object@annotation_results$matched_peak_pairs_all <- reverse_matched_peaks_all
    object@annotation_results$matched_peak_pairs <- reverse_matched_peaks
    object@ms2_data$reverse_matched_peaks_all <- reverse_matched_peaks_all
    object@ms2_data$reverse_matched_peaks <- reverse_matched_peaks
    object@ms2_data$matched_peak_pairs_all <- reverse_matched_peaks_all
    object@ms2_data$matched_peak_pairs <- reverse_matched_peaks
  } else {
    object <- `60_strip_matched_peak_pairs_from_object`(object, matched_peak_pairs_summary)
  }
  if (is.list(parameters)) {
    parameters$store_matched_peak_pairs_in_rds <- store_matched_peak_pairs_in_rds
  }
  object@msi_parameters$reverse_match <- parameters
  object@process_info$module6_reverse_match <- parameters
  object
}

`60_strip_matched_peak_pairs_from_object` <- 
function(object, matched_peak_pairs_summary = NULL) {
  if (!methods::is(object, "msi_dataset")) {
    return(object)
  }
  if (is.data.frame(matched_peak_pairs_summary)) {
    object@annotation_results$matched_peak_pairs_summary <- matched_peak_pairs_summary
    object@ms2_data$matched_peak_pairs_summary <- matched_peak_pairs_summary
  }
  for (name in c("reverse_matched_peaks", "reverse_matched_peaks_all",
                 "matched_peak_pairs", "matched_peak_pairs_all")) {
    object@annotation_results[[name]] <- NULL
    object@ms2_data[[name]] <- NULL
  }
  if (is.list(object@msi_parameters$reverse_match)) {
    object@msi_parameters$reverse_match$store_matched_peak_pairs_in_rds <- FALSE
  }
  if (is.list(object@process_info$module6_reverse_match)) {
    object@process_info$module6_reverse_match$store_matched_peak_pairs_in_rds <- FALSE
  }
  object
}

`60_summarise_reverse_match` <- 
function(candidates,
	                                         matched_peak_pairs,
	                                         min_matched_peaks_evidence = 3L,
	                                         min_reverse_masked_cosine_evidence = 0,
	                                         progress_every = NULL,
                                         progress_path = NULL,
                                         progress_started_at = NULL,
                                         progress_callback = NULL) {
  empty <- .rm_empty_reverse_scores()
  candidates <- as.data.frame(candidates, stringsAsFactors = FALSE, check.names = FALSE)
	  matched_peak_pairs <- as.data.frame(matched_peak_pairs, stringsAsFactors = FALSE, check.names = FALSE)
	  min_matched_peaks_evidence <- suppressWarnings(as.integer(min_matched_peaks_evidence[[1]]))
	  if (!is.finite(min_matched_peaks_evidence) || min_matched_peaks_evidence < 1L) {
	    min_matched_peaks_evidence <- 3L
	  }
	  min_reverse_masked_cosine_evidence <- suppressWarnings(as.numeric(min_reverse_masked_cosine_evidence[[1]]))
	  if (!is.finite(min_reverse_masked_cosine_evidence) || min_reverse_masked_cosine_evidence < 0) {
	    min_reverse_masked_cosine_evidence <- 0
	  }
	  if (nrow(candidates) == 0L) {
    return(empty)
  }
  progress_started_at <- if (inherits(progress_started_at, "POSIXt")) progress_started_at else Sys.time()
  progress_every <- if (is.null(progress_every) || length(progress_every) == 0L) {
    nrow(candidates)
  } else {
    suppressWarnings(as.integer(progress_every[[1]]))
  }
  if (!is.finite(progress_every) || progress_every < 1L) {
    progress_every <- nrow(candidates)
  }
  matched_split <- list()
  if (nrow(matched_peak_pairs)) {
    matched_pair_keys <- .rm_split_key(matched_peak_pairs, c("candidate_key", "cluster_id", "anchor_feature_id"))
    matched_split <- split(matched_peak_pairs, matched_pair_keys, drop = TRUE)
  }

  rows <- vector("list", nrow(candidates))
  for (i in seq_len(nrow(candidates))) {
    candidate <- candidates[i, , drop = FALSE]
    pair_key <- .rm_split_key(candidate, c("candidate_key", "cluster_id", "anchor_feature_id"))
    pairs <- matched_split[[pair_key]]
    if (is.null(pairs)) {
      pairs <- matched_peak_pairs[0, , drop = FALSE]
    }
    pairs <- .rm_align_reverse_pairs(pairs)
    status <- "matched"
    if (nrow(pairs) == 0L) {
      status <- "no_library_spectrum"
    } else if (all(!(pairs$matched %in% TRUE))) {
      status <- "no_matched_peaks"
    }
    reverse_masked_cosine <- if (nrow(pairs) > 0L && "reverse_masked_cosine" %in% names(pairs)) {
      value <- suppressWarnings(as.numeric(pairs$reverse_masked_cosine))
      value <- value[is.finite(value)]
      if (length(value)) value[[1]] else `60_calculate_reverse_cosine`(pairs$library_intensity, pairs$query_intensity)
    } else if (nrow(pairs) > 0L) {
      `60_calculate_reverse_cosine`(pairs$library_intensity, pairs$query_intensity)
    } else {
      0
    }
    reverse_masked_cosine <- pmin(pmax(as.numeric(reverse_masked_cosine), 0), 1)
    coverage <- if (nrow(pairs) > 0L) {
      `60_calculate_library_coverage`(pairs$library_intensity, pairs$matched)
    } else {
      list(S_lib_coverage = 0, matched_peak_fraction = 0, n_reverse_match = 0L)
    }
    library_l2_coverage <- if (nrow(pairs) > 0L && "library_l2_coverage" %in% names(pairs)) {
      value <- suppressWarnings(as.numeric(pairs$library_l2_coverage))
      value <- value[is.finite(value)]
      if (length(value)) value[[1]] else {
        lib <- suppressWarnings(as.numeric(pairs$library_intensity))
        matched_for_l2 <- pairs$matched %in% TRUE
        denom <- sum(lib^2, na.rm = TRUE)
        if (is.finite(denom) && denom > 0) sqrt(sum(lib[matched_for_l2]^2, na.rm = TRUE) / denom) else 0
      }
    } else if (nrow(pairs) > 0L) {
      lib <- suppressWarnings(as.numeric(pairs$library_intensity))
      matched_for_l2 <- pairs$matched %in% TRUE
      denom <- sum(lib^2, na.rm = TRUE)
      if (is.finite(denom) && denom > 0) sqrt(sum(lib[matched_for_l2]^2, na.rm = TRUE) / denom) else 0
    } else {
      0
    }
    spatial <- if (nrow(pairs) > 0L) {
      `60_calculate_spatial_reverse_score`(pairs$SAS, pairs$library_intensity, pairs$matched)
    } else {
      NA_real_
    }
    S_SAS_score_for_total <- if (is.finite(spatial)) pmin(pmax(as.numeric(spatial), 0), 1) else 0
    S_SAS_missing <- !is.finite(spatial)
    S_total <- pmin(pmax(0.7 * reverse_masked_cosine + 0.3 * S_SAS_score_for_total, 0), 1)
    matched_bool <- if (nrow(pairs) > 0L) pairs$matched %in% TRUE else logical()
    n_unique_query <- if (nrow(pairs) > 0L) {
      length(unique(as.character(pairs$query_peak_id[matched_bool & !is.na(pairs$query_peak_id)])))
    } else {
      0L
    }
    n_unique_library <- if (nrow(pairs) > 0L) {
      length(unique(as.character(pairs$library_peak_id[matched_bool & !is.na(pairs$library_peak_id)])))
    } else {
      0L
    }
    n_query_total <- if (nrow(pairs) > 0L && "n_query_peaks_total" %in% names(pairs)) {
      value <- suppressWarnings(as.integer(pairs$n_query_peaks_total))
      value <- value[is.finite(value)]
      if (length(value)) max(value) else 0L
    } else {
      0L
    }
    if (!is.finite(n_query_total)) n_query_total <- 0L
    n_library_total <- if (nrow(pairs) > 0L && "n_library_peaks_total" %in% names(pairs)) {
      value <- suppressWarnings(as.integer(pairs$n_library_peaks_total))
      value <- value[is.finite(value)]
      if (length(value)) max(value) else nrow(pairs)
    } else if (nrow(pairs) > 0L) {
      nrow(pairs)
    } else {
      0L
    }
    if (!is.finite(n_library_total)) n_library_total <- 0L
    matched_pair_n <- if (nrow(pairs) > 0L) sum(matched_bool, na.rm = TRUE) else 0L
    query_peak_reuse_n <- max(0L, as.integer(matched_pair_n - n_unique_query))
    query_peak_reuse_fraction <- if (matched_pair_n > 0L) {
      query_peak_reuse_n / matched_pair_n
    } else {
      0
    }
    precursor_error <- suppressWarnings(as.numeric(candidate$precursor_error_ppm[[1]]))
    S_precursor <- if (is.finite(precursor_error)) 1 / (1 + abs(precursor_error) / 10) else 0
    matched_ppm <- suppressWarnings(as.numeric(pairs$mz_error_ppm[matched_bool]))
    mass_error_ppm_median <- if (length(matched_ppm) && any(is.finite(matched_ppm))) {
      stats::median(matched_ppm, na.rm = TRUE)
    } else {
      NA_real_
    }
    mass_error_ppm_max <- if (length(matched_ppm) && any(is.finite(matched_ppm))) {
      max(abs(matched_ppm), na.rm = TRUE)
    } else {
      NA_real_
    }
    matched_pair_ids <- as.character(pairs$match_id[matched_bool])
    matched_pair_ids <- matched_pair_ids[!is.na(matched_pair_ids) & nzchar(matched_pair_ids)]
    mass_error_ppm_summary <- if (length(matched_ppm) && any(is.finite(matched_ppm))) {
      paste0(
        "median=", signif(stats::median(matched_ppm, na.rm = TRUE), 6),
        ";max_abs=", signif(max(abs(matched_ppm), na.rm = TRUE), 6)
      )
    } else {
      NA_character_
    }
    matched_query_intensity <- if (nrow(pairs) > 0L && "matched_query_intensity" %in% names(pairs)) {
      value <- suppressWarnings(as.numeric(pairs$matched_query_intensity))
      value <- value[is.finite(value)]
      if (length(value)) value[[1]] else sum(as.numeric(pairs$query_intensity[matched_bool]), na.rm = TRUE)
    } else {
      sum(as.numeric(pairs$query_intensity[matched_bool]), na.rm = TRUE)
    }
    matched_library_intensity <- if (nrow(pairs) > 0L && "matched_library_intensity" %in% names(pairs)) {
      value <- suppressWarnings(as.numeric(pairs$matched_library_intensity))
      value <- value[is.finite(value)]
      if (length(value)) value[[1]] else sum(as.numeric(pairs$library_intensity[matched_bool]), na.rm = TRUE)
    } else {
      sum(as.numeric(pairs$library_intensity[matched_bool]), na.rm = TRUE)
    }
    matched_query_intensity_fraction <- if (nrow(pairs) > 0L &&
        "matched_query_intensity_fraction" %in% names(pairs)) {
      value <- suppressWarnings(as.numeric(pairs$matched_query_intensity_fraction))
      value <- value[is.finite(value)]
      if (length(value)) max(value) else 0
    } else {
      0
    }
    matched_library_intensity_fraction <- if (nrow(pairs) > 0L &&
        "matched_library_intensity_fraction" %in% names(pairs)) {
      value <- suppressWarnings(as.numeric(pairs$matched_library_intensity_fraction))
      value <- value[is.finite(value)]
      if (length(value)) max(value) else as.numeric(coverage$S_lib_coverage)
    } else {
      as.numeric(coverage$S_lib_coverage)
    }
    target_hmdb_id <- regmatches(
      paste(candidate$candidate_id[[1]], candidate$candidate_name[[1]], candidate$metabolite_name[[1]]),
      regexpr("HMDB[0-9]{5,}", paste(candidate$candidate_id[[1]], candidate$candidate_name[[1]], candidate$metabolite_name[[1]]),
              perl = TRUE, ignore.case = TRUE)
    )
    target_hmdb_id <- if (length(target_hmdb_id) && !is.na(target_hmdb_id) && nzchar(target_hmdb_id)) {
      toupper(target_hmdb_id)
    } else {
      NA_character_
    }
    source_value <- as.character(candidate$database_source[[1]])
    if (is.na(source_value) || !nzchar(source_value)) source_value <- as.character(candidate$library_source[[1]])
    lab_id_value <- as.character(candidate$spectrum_id[[1]])
    if (is.na(lab_id_value) || !nzchar(lab_id_value)) lab_id_value <- as.character(candidate$candidate_id[[1]])
    rows[[i]] <- data.frame(
      component = if ("component" %in% names(candidate)) as.character(candidate$component[[1]]) else NA_character_,
      cluster_id = candidate$cluster_id[[1]],
      anchor_feature_id = as.character(candidate$anchor_feature_id[[1]]),
      feature_id = as.character(candidate$anchor_feature_id[[1]]),
      anchor_mz = as.numeric(candidate$anchor_mz[[1]]),
      candidate_key = as.character(candidate$candidate_key[[1]]),
      query_id = as.character(candidate$query_id[[1]]),
      query_spectrum_id = as.character(candidate$query_id[[1]]),
      candidate_id = as.character(candidate$candidate_id[[1]]),
      spectrum_id = as.character(candidate$spectrum_id[[1]]),
      archive_id = as.character(candidate$archive_id[[1]]),
      target_hmdb_id = target_hmdb_id,
      metabolite_name = as.character(candidate$candidate_name[[1]]),
      candidate_name = as.character(candidate$candidate_name[[1]]),
      library_source = as.character(candidate$library_source[[1]]),
      database_source = as.character(candidate$database_source[[1]]),
      source = source_value,
      lab_id = lab_id_value,
      instrument = if ("instrument" %in% names(candidate)) as.character(candidate$instrument[[1]]) else NA_character_,
      precursor_mz = as.numeric(candidate$precursor_mz[[1]]),
      formula = as.character(candidate$formula[[1]]),
      inchikey = as.character(candidate$inchikey[[1]]),
      smiles = as.character(candidate$smiles[[1]]),
      adduct = as.character(candidate$adduct[[1]]),
      adduct_status = if ("adduct_status" %in% names(candidate)) {
        as.character(candidate$adduct_status[[1]])
      } else {
        .msi_adduct_status(candidate$adduct[[1]])
      },
      polarity = as.character(candidate$polarity[[1]]),
      candidate_source = as.character(candidate$candidate_source[[1]]),
      candidate_retrieved = if ("candidate_retrieved" %in% names(candidate)) {
        isTRUE(candidate$candidate_retrieved[[1]])
      } else {
        TRUE
      },
      precursor_gate_pass = if ("precursor_gate_pass" %in% names(candidate)) {
        isTRUE(candidate$precursor_gate_pass[[1]])
      } else {
        TRUE
      },
      fragment_retrieval_hits = if ("fragment_retrieval_hits" %in% names(candidate)) {
        as.integer(candidate$fragment_retrieval_hits[[1]])
      } else if ("n_FE_match" %in% names(candidate)) {
        as.integer(candidate$n_FE_match[[1]])
      } else {
        NA_integer_
      },
      entropy_S = as.numeric(candidate$entropy_S[[1]]),
      S_FE = as.numeric(candidate$entropy_S[[1]]),
      retrieval_rank = as.integer(candidate$retrieval_rank[[1]]),
      rank_entropy = as.integer(candidate$retrieval_rank[[1]]),
      precursor_error_ppm = precursor_error,
      S_precursor = as.numeric(S_precursor),
      reverse_masked_cosine = as.numeric(reverse_masked_cosine),
      mass_error_ppm_summary = mass_error_ppm_summary,
      mass_error_ppm_median = as.numeric(mass_error_ppm_median),
      mass_error_ppm_max = as.numeric(mass_error_ppm_max),
      matched_peak_pairs_id = paste(matched_pair_ids, collapse = ";"),
      query_unmatched_peak_policy = "ignored",
      reference_unmatched_peak_policy = "masked_query_zero_in_library_space",
      S_lib_coverage = as.numeric(coverage$S_lib_coverage),
      library_l2_coverage = as.numeric(library_l2_coverage),
      matched_library_intensity_fraction = as.numeric(matched_library_intensity_fraction),
      matched_query_intensity_fraction = as.numeric(matched_query_intensity_fraction),
      matched_peak_fraction = as.numeric(coverage$matched_peak_fraction),
      n_matched_peaks = as.integer(coverage$n_reverse_match),
      n_reverse_match = as.integer(coverage$n_reverse_match),
      S_SAS = as.numeric(spatial),
      S_SAS_score_for_total = as.numeric(S_SAS_score_for_total),
      S_SAS_missing = isTRUE(S_SAS_missing),
      S_total = as.numeric(S_total),
      score_ms2 = as.numeric(reverse_masked_cosine),
      score_spatial = as.numeric(S_SAS_score_for_total),
      score_total_pre_mfc = as.numeric(S_total),
      S_spatial_reverse = as.numeric(spatial),
      n_SAS_matched_peaks = as.integer(coverage$n_reverse_match),
      sum_matched_library_intensity = as.numeric(matched_library_intensity),
      matched_library_intensity = as.numeric(matched_library_intensity),
      matched_query_intensity = as.numeric(matched_query_intensity),
      total_library_intensity = sum(as.numeric(pairs$library_intensity), na.rm = TRUE),
      n_pseudo_peaks = as.integer(n_query_total),
      n_library_peaks = as.integer(n_library_total),
      reverse_match_status = status,
      s_fe = as.numeric(candidate$entropy_S[[1]]),
      s_rev = as.numeric(reverse_masked_cosine),
      n_match = as.integer(coverage$n_reverse_match),
      u_spec = as.numeric(coverage$S_lib_coverage),
      s_spat = as.numeric(spatial),
      s_arch = as.numeric(candidate$s_arch[[1]]),
      matched_query_peak_n = if (nrow(pairs) > 0L) length(unique(pairs$query_peak_id[matched_bool])) else 0L,
      matched_reference_peak_n = as.integer(coverage$n_reverse_match),
      n_match_pairs = as.integer(coverage$n_reverse_match),
      unique_query_matched_peaks = as.integer(n_unique_query),
      unique_library_matched_peaks = as.integer(n_unique_library),
      n_query_peaks_total = as.integer(n_query_total),
      n_library_peaks_total = as.integer(n_library_total),
      query_match_ratio = if (n_query_total > 0L) n_unique_query / n_query_total else NA_real_,
      library_match_ratio = if (n_library_total > 0L) n_unique_library / n_library_total else NA_real_,
      query_peak_reuse_n = as.integer(query_peak_reuse_n),
      query_peak_reuse_fraction = as.numeric(query_peak_reuse_fraction),
      has_query_peak_reuse = query_peak_reuse_n > 0L,
      match_mode = "pseudo_ms2",
	      pass_n_matched_peaks_ge3 = as.integer(coverage$n_reverse_match) >= min_matched_peaks_evidence,
      pass_reverse_masked_cosine = is.finite(reverse_masked_cosine) &
	        reverse_masked_cosine > min_reverse_masked_cosine_evidence,
	      fragment_evidence_hits = as.integer(coverage$n_reverse_match),
	      evidence_eligible = !identical(status, "no_library_spectrum") &&
	        as.integer(coverage$n_reverse_match) > 0L,
	      evidence_eligibility_reason = if (identical(status, "no_library_spectrum")) {
	        "NO_REFERENCE_SPECTRUM"
	      } else if (as.integer(coverage$n_reverse_match) == 0L) {
	        "NO_EVIDENCE_COMPATIBLE_FRAGMENT"
	      } else {
	        "HAS_EVIDENCE_COMPATIBLE_FRAGMENT"
	      },
	      evidence_pass = as.integer(coverage$n_reverse_match) >= min_matched_peaks_evidence &
	        is.finite(reverse_masked_cosine) &
	        reverse_masked_cosine > min_reverse_masked_cosine_evidence,
	      min_matched_peaks_used = as.integer(min_matched_peaks_evidence),
	      min_reverse_masked_cosine_used = as.numeric(min_reverse_masked_cosine_evidence),
	      evidence_failure_reason = if (identical(status, "no_library_spectrum")) {
	        "missing library spectrum"
	      } else if (as.integer(coverage$n_reverse_match) < min_matched_peaks_evidence) {
	        paste0("n_matched_peaks < ", min_matched_peaks_evidence)
	      } else if (!is.finite(reverse_masked_cosine) ||
	                 reverse_masked_cosine <= min_reverse_masked_cosine_evidence) {
	        "reverse_masked_cosine <= threshold"
	      } else {
	        "pass"
	      },
	      evidence_filter_reason = if (identical(status, "no_library_spectrum")) {
	        "missing library spectrum"
	      } else if (as.integer(coverage$n_reverse_match) < min_matched_peaks_evidence) {
	        paste0("n_matched_peaks < ", min_matched_peaks_evidence)
	      } else if (!is.finite(reverse_masked_cosine) ||
	                 reverse_masked_cosine <= min_reverse_masked_cosine_evidence) {
	        "reverse_masked_cosine <= threshold"
	      } else {
        "pass"
      },
      stringsAsFactors = FALSE,
      check.names = FALSE
    )
    if (i %% progress_every == 0L || i == nrow(candidates)) {
      .rm_emit_progress(
        progress_callback = progress_callback,
        progress_path = progress_path,
        started_at = progress_started_at,
        processed_candidates = i,
        total_candidates = nrow(candidates),
        reverse_rows_so_far = i,
        pass_rows_so_far = sum(vapply(rows[seq_len(i)], function(x) {
          if (!is.data.frame(x) || !"evidence_pass" %in% names(x)) return(FALSE)
          x$evidence_pass[[1]] %in% TRUE
        }, logical(1))),
        message = "summarising reverse scores"
      )
    }
  }
  out <- do.call(rbind, rows)
  out <- .rm_add_missing_columns(out, colnames(empty))
  out <- out[, colnames(empty), drop = FALSE]
  rownames(out) <- NULL
  out
}

`compute_reverse_masked_cosine` <- 
function(
    query_peaks,
    library_peaks,
    tolerance_ppm = 10,
    tolerance_da_lowmz = 0.005,
    lowmz_threshold = 400,
    one_to_one = TRUE) {
  query_peaks <- as.data.frame(query_peaks, stringsAsFactors = FALSE, check.names = FALSE)
  library_peaks <- as.data.frame(library_peaks, stringsAsFactors = FALSE, check.names = FALSE)
  ppm <- suppressWarnings(as.numeric(tolerance_ppm[[1]]))
  if (!is.finite(ppm) || ppm < 0) ppm <- 10
  low_da <- suppressWarnings(as.numeric(tolerance_da_lowmz[[1]]))
  if (!is.finite(low_da) || low_da < 0) low_da <- 0
  lowmz_threshold <- suppressWarnings(as.numeric(lowmz_threshold[[1]]))
  if (!is.finite(lowmz_threshold) || lowmz_threshold < 0) lowmz_threshold <- 400

  lib_mz_col <- .reverse_dp_first_column(
    library_peaks, c("library_mz", "peak_mz", "fragment_mz", "mz")
  )
  lib_int_col <- .reverse_dp_first_column(
    library_peaks,
    c("library_intensity_raw", "intensity_raw", "library_intensity", "intensity_sum1",
      "fragment_intensity", "intensity")
  )
  query_mz_col <- .reverse_dp_first_column(
    query_peaks, c("query_mz", "peak_mz", "fragment_mz", "mz")
  )
  query_int_col <- .reverse_dp_first_column(
    query_peaks,
    c("query_intensity_raw", "pseudo_raw_intensity", "query_intensity", "intensity",
      "fragment_intensity")
  )

  lib <- data.frame(
    library_peak_order = seq_len(nrow(library_peaks)),
    library_peak_id = if ("library_peak_id" %in% names(library_peaks)) {
      as.character(library_peaks$library_peak_id)
    } else if (nrow(library_peaks) == 0L) {
      character()
    } else {
      paste0("library_peak_", seq_len(nrow(library_peaks)))
    },
    library_mz = suppressWarnings(as.numeric(library_peaks[[lib_mz_col]])),
    library_intensity_raw = suppressWarnings(as.numeric(library_peaks[[lib_int_col]])),
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
  query <- data.frame(
    query_peak_order = seq_len(nrow(query_peaks)),
    query_peak_id = if ("query_peak_id" %in% names(query_peaks)) {
      as.character(query_peaks$query_peak_id)
    } else if ("peak_feature_id" %in% names(query_peaks)) {
      as.character(query_peaks$peak_feature_id)
    } else if (nrow(query_peaks) == 0L) {
      character()
    } else {
      paste0("query_peak_", seq_len(nrow(query_peaks)))
    },
    query_mz = suppressWarnings(as.numeric(query_peaks[[query_mz_col]])),
    query_intensity_raw = suppressWarnings(as.numeric(query_peaks[[query_int_col]])),
    SAS = if ("SAS" %in% names(query_peaks)) {
      suppressWarnings(as.numeric(query_peaks$SAS))
    } else if ("query_sas" %in% names(query_peaks)) {
      suppressWarnings(as.numeric(query_peaks$query_sas))
    } else if ("graph_weight" %in% names(query_peaks)) {
      suppressWarnings(as.numeric(query_peaks$graph_weight))
    } else {
      rep(NA_real_, nrow(query_peaks))
    },
    stringsAsFactors = FALSE,
    check.names = FALSE
  )

  lib <- lib[is.finite(lib$library_mz) & is.finite(lib$library_intensity_raw) &
               lib$library_intensity_raw > 0, , drop = FALSE]
  query <- query[is.finite(query$query_mz) & is.finite(query$query_intensity_raw) &
                   query$query_intensity_raw > 0, , drop = FALSE]
  lib <- lib[order(lib$library_mz, lib$library_peak_order), , drop = FALSE]
  query <- query[order(query$query_mz, query$query_peak_order), , drop = FALSE]
  rownames(lib) <- NULL
  rownames(query) <- NULL
  lib$library_intensity_norm <- .reverse_dp_normalize_max100(lib$library_intensity_raw)
  query$query_intensity_norm <- .reverse_dp_normalize_max100(query$query_intensity_raw)

  empty_matched <- lib
  empty_matched$effective_tolerance_da <- numeric(nrow(empty_matched))
  if (nrow(empty_matched)) {
    ppm_da <- empty_matched$library_mz * ppm / 1e6
    empty_matched$effective_tolerance_da <- ifelse(
      empty_matched$library_mz < lowmz_threshold,
      pmax(ppm_da, low_da),
      ppm_da
    )
  }
  n_empty_matched <- nrow(empty_matched)
  empty_matched$match_pair_id <- rep(NA_character_, n_empty_matched)
  empty_matched$query_peak_id <- rep(NA_character_, n_empty_matched)
  empty_matched$query_mz <- rep(NA_real_, n_empty_matched)
  empty_matched$query_intensity_raw <- rep(0, n_empty_matched)
  empty_matched$query_intensity_norm <- rep(0, n_empty_matched)
  empty_matched$mz_error_da <- rep(NA_real_, n_empty_matched)
  empty_matched$mz_error_ppm <- rep(NA_real_, n_empty_matched)
  empty_matched$matched <- rep(FALSE, n_empty_matched)
  empty_matched$SAS <- rep(NA_real_, n_empty_matched)
  if (!nrow(lib)) {
    return(list(
      reverse_masked_cosine = 0,
      n_matched_peaks = 0L,
      matched_library_intensity_fraction = 0,
      matched_query_intensity_fraction = 0,
      matched_library_intensity = 0,
      matched_query_intensity = 0,
      total_library_intensity = 0,
      total_query_intensity = sum(query$query_intensity_raw, na.rm = TRUE),
      library_l2_coverage = 0,
      library_all_l2 = 0,
      query_masked_l2 = 0,
      query_full_l2 = sqrt(sum(query$query_intensity_raw^2, na.rm = TRUE)),
      S_SAS = NA_real_,
      matched_peaks = empty_matched,
      query_peaks = query
    ))
  }

  matched <- empty_matched
  selected <- data.frame(
    library_row = integer(), query_row = integer(), mz_error_da = numeric(),
    mz_error_ppm = numeric(), mz_error_abs = numeric(),
    query_intensity_raw = numeric(), stringsAsFactors = FALSE, check.names = FALSE
  )
  if (nrow(query)) {
    candidate_rows <- vector("list", nrow(lib))
    for (i in seq_len(nrow(lib))) {
      tolerance <- matched$effective_tolerance_da[[i]]
      if (!is.finite(tolerance)) next
      error_da <- query$query_mz - lib$library_mz[[i]]
      within <- is.finite(error_da) & abs(error_da) <= tolerance
      if (!any(within)) next
      candidate_rows[[i]] <- data.frame(
        library_row = i,
        query_row = which(within),
        mz_error_da = error_da[within],
        mz_error_ppm = error_da[within] / lib$library_mz[[i]] * 1e6,
        mz_error_abs = abs(error_da[within]),
        query_intensity_raw = query$query_intensity_raw[within],
        SAS_rank = ifelse(is.finite(query$SAS[within]), query$SAS[within], -Inf),
        stringsAsFactors = FALSE,
        check.names = FALSE
      )
    }
    candidate_table <- do.call(rbind, candidate_rows[!vapply(candidate_rows, is.null, logical(1))])
    if (is.data.frame(candidate_table) && nrow(candidate_table)) {
      candidate_table <- candidate_table[order(
        candidate_table$mz_error_abs,
        -candidate_table$query_intensity_raw,
        -candidate_table$SAS_rank,
        candidate_table$library_row,
        candidate_table$query_row
      ), , drop = FALSE]
      used_library <- rep(FALSE, nrow(lib))
      used_query <- rep(FALSE, nrow(query))
      keep <- logical(nrow(candidate_table))
      for (i in seq_len(nrow(candidate_table))) {
        li <- candidate_table$library_row[[i]]
        qi <- candidate_table$query_row[[i]]
        if (used_library[[li]]) next
        if (isTRUE(one_to_one) && used_query[[qi]]) next
        keep[[i]] <- TRUE
        used_library[[li]] <- TRUE
        if (isTRUE(one_to_one)) used_query[[qi]] <- TRUE
      }
      selected <- candidate_table[keep, , drop = FALSE]
      selected <- selected[order(selected$library_row, selected$query_row), , drop = FALSE]
      selected$match_pair_id <- sprintf("match_%03d", seq_len(nrow(selected)))
      for (i in seq_len(nrow(selected))) {
        li <- selected$library_row[[i]]
        qi <- selected$query_row[[i]]
        matched$match_pair_id[[li]] <- selected$match_pair_id[[i]]
        matched$query_peak_id[[li]] <- query$query_peak_id[[qi]]
        matched$query_mz[[li]] <- query$query_mz[[qi]]
        matched$query_intensity_raw[[li]] <- query$query_intensity_raw[[qi]]
        matched$query_intensity_norm[[li]] <- query$query_intensity_norm[[qi]]
        matched$mz_error_da[[li]] <- selected$mz_error_da[[i]]
        matched$mz_error_ppm[[li]] <- selected$mz_error_ppm[[i]]
        matched$matched[[li]] <- TRUE
        matched$SAS[[li]] <- query$SAS[[qi]]
      }
    }
  }

  L <- pmax(suppressWarnings(as.numeric(matched$library_intensity_raw)), 0)
  Qm <- pmax(suppressWarnings(as.numeric(matched$query_intensity_raw)), 0)
  denom <- sqrt(sum(L^2, na.rm = TRUE)) * sqrt(sum(Qm^2, na.rm = TRUE))
  reverse_masked_cosine <- if (is.finite(denom) && denom > 0) {
    sum(L * Qm, na.rm = TRUE) / denom
  } else {
    0
  }
  reverse_masked_cosine <- pmin(pmax(as.numeric(reverse_masked_cosine), 0), 1)
  total_library <- sum(L, na.rm = TRUE)
  total_query <- sum(query$query_intensity_raw, na.rm = TRUE)
  matched_library <- sum(L[matched$matched %in% TRUE], na.rm = TRUE)
  matched_query <- sum(Qm[matched$matched %in% TRUE], na.rm = TRUE)
  library_l2_all <- sqrt(sum(L^2, na.rm = TRUE))
  library_l2_matched <- sqrt(sum(L[matched$matched %in% TRUE]^2, na.rm = TRUE))
  sas_keep <- matched$matched %in% TRUE & is.finite(matched$SAS) & L > 0
  S_SAS <- if (any(sas_keep)) {
    stats::weighted.mean(matched$SAS[sas_keep], w = L[sas_keep])
  } else {
    NA_real_
  }
  matched$contribution_to_dot <- matched$library_intensity_raw * matched$query_intensity_raw
  matched$contribution_to_S_SAS <- ifelse(
    matched$matched & is.finite(matched$SAS),
    matched$library_intensity_raw * matched$SAS,
    0
  )
  list(
    reverse_masked_cosine = reverse_masked_cosine,
    n_matched_peaks = as.integer(sum(matched$matched %in% TRUE, na.rm = TRUE)),
    matched_library_intensity_fraction = if (total_library > 0) matched_library / total_library else 0,
    matched_query_intensity_fraction = if (total_query > 0) matched_query / total_query else 0,
    matched_library_intensity = matched_library,
    matched_query_intensity = matched_query,
    total_library_intensity = total_library,
    total_query_intensity = total_query,
    library_l2_coverage = if (library_l2_all > 0) library_l2_matched / library_l2_all else 0,
    library_all_l2 = library_l2_all,
    query_masked_l2 = sqrt(sum(Qm^2, na.rm = TRUE)),
    query_full_l2 = sqrt(sum(query$query_intensity_raw^2, na.rm = TRUE)),
    S_SAS = S_SAS,
    matched_peaks = matched,
    query_peaks = query
  )
}

