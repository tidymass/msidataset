test_that("D15 step6 through MFC reproduces the preregistered regression", {
  skip_if_not(identical(Sys.getenv("MSIDATASET_RUN_D15_MFC"), "true"),
              "Set MSIDATASET_RUN_D15_MFC=true for the large external D15 regression")
  cfg <- locked_d15()
  keys <- c("csv", "labels", "meta", "step6", "mfc_reference")
  paths <- lapply(cfg[keys], `[[`, "path")
  d15_require_files(unlist(paths))
  resource_keys <- c("adp_config", "ms1_library", "ms1_project_root", "vector_ids")
  d15_require_files(vapply(cfg[resource_keys], `[[`, character(1), "path"))
  study <- d15_fixture_root("MSIDATASET_D15_STUDY_ROOT")
  reports <- d15_fixture_root("MSIDATASET_D15_REPORTS_ROOT")
  form_path <- file.path(study, "runs/d11_kidney/report_d11/formula.tsv")
  qc_path <- file.path(reports, "D12", "formula_qc.tsv")
  d15_require_files(c(form_path, qc_path))
  meta <- data.table::fread(paths$meta, data.table = FALSE)
  x <- read_msi_data(paths$csv, mz = meta[["m/z"]] - 1.007276,
                     feature_id = sprintf("F%03d", meta$col_index))
  x <- define_roi(x, data.table::fread(paths$labels, data.table = FALSE))
  x <- detect_de_features(x)
  x <- generate_pseudo_ms2(x)
  out <- tempfile("mfc-d15-"); dir.create(out)
  x <- calculate_mfc(x, readRDS(paths$step6), out)
  z <- x@process_info$calculate_mfc$result
  expect_identical(z$status, "OK")
  observed <- z$candidates
  expected <- data.table::fread(paths$mfc_reference, data.table = FALSE)
  key <- function(y) paste(y$feature_id, y$HMDB_ID, y$best_adduct, sep = "|")
  expect_identical(key(observed), key(expected))
  expect_identical(as.integer(observed$rank_within_feature), as.integer(expected$rank_within_feature))
  expect_lt(max(abs(observed$mfc_priority_score - expected$mfc_priority_score) /
                  pmax(abs(expected$mfc_priority_score), .Machine$double.xmin)), 1e-6)
  expect_equal(length(unique(observed$feature_id)), 149)
  expect_equal(nrow(x@variable_info), 156)

  # Author formulas are read only here, after inference, for preregistered evaluation.
  formulas <- data.table::fread(form_path, data.table = FALSE)
  qc <- data.table::fread(qc_path, data.table = FALSE)
  canonical <- function(v) {
    if (is.na(v) || !nzchar(v)) return(NA_character_)
    p <- regmatches(v, gregexpr("[A-Z][a-z]?[0-9]*", v, perl = TRUE))[[1L]]
    if (!length(p) || paste0(p, collapse = "") != v) return(NA_character_)
    element <- sub("[0-9]+$", "", p); n <- sub("^[A-Za-z]+", "", p); n[n == ""] <- "1"
    counts <- tapply(as.integer(n), element, sum)
    paste0(names(counts), "=", counts, collapse = ";")
  }
  refs <- meta[match(qc$col_index[qc$formula_consistent %in% TRUE], meta$col_index), , drop = FALSE]
  targets <- vapply(refs$Formula, canonical, character(1))
  formulas <- formulas[formulas$HMDB.ID %in% observed$HMDB_ID, , drop = FALSE]
  formulas$canonical <- vapply(formulas$Formula, canonical, character(1))
  scores <- lapply(seq_len(nrow(refs)), function(i) {
    sub <- observed[observed$feature_id == sprintf("F%03d", refs$col_index[i]), , drop = FALSE]
    ids <- formulas$HMDB.ID[!is.na(formulas$canonical) & formulas$canonical == targets[i]]
    match <- sub$HMDB_ID %in% ids
    if (!any(match)) return(c(present = 0L, top1 = 0L))
    best <- max(sub$mfc_priority_score[match])
    c(present = 1L, top1 = as.integer(!any(sub$mfc_priority_score > best)))
  })
  scores <- do.call(rbind, scores)
  expect_equal(nrow(refs), 95)
  expect_equal(sum(scores[, "present"]), 89)
  expect_equal(sum(scores[, "top1"]), 88)
})
