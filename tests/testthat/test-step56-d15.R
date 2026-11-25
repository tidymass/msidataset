test_that('D15 pseudo-MS2 through step5, step6 and MFC matches the locked run', {
  skip_if_not(identical(Sys.getenv('MSIDATASET_RUN_D15_STEP56'),'true'),
    'Set MSIDATASET_RUN_D15_STEP56=true to run the large external regression')
  cfg <- locked_d15()
  keys <- c('meta', 'csv', 'labels', 'step6', 'mfc_reference',
            'ms2_library', 'adp_config', 'ms1_library', 'ms1_project_root', 'vector_ids')
  d15_require_files(vapply(cfg[keys], `[[`, character(1), 'path'))
  study <- d15_fixture_root("MSIDATASET_D15_STUDY_ROOT")
  reports <- d15_fixture_root("MSIDATASET_D15_REPORTS_ROOT")
  form_path <- file.path(study, "runs/d11_kidney/report_d11/formula.tsv")
  qc_path <- file.path(reports, "D12", "formula_qc.tsv")
  d15_require_files(c(form_path, qc_path))
  meta <- data.table::fread(cfg$meta$path, data.table=FALSE)
  x <- read_msi_data(cfg$csv$path, mz=meta[['m/z']]-1.007276, feature_id=sprintf('F%03d',meta$col_index))
  x <- define_roi(x,data.table::fread(cfg$labels$path,data.table=FALSE))
  x <- generate_pseudo_ms2(x)
  out <- tempfile('step56-'); dir.create(out)
  x <- match_spectral_library(x,list(output_dir=out))
  z <- x@process_info$match_spectral_library$result
  expect_identical(z$status,'OK')
  audit_dir <- Sys.getenv('MSIDATASET_STEP56_AUDIT_DIR')
  if (nzchar(audit_dir)) saveRDS(z$step6, file.path(audit_dir, 'reconstructed_step6.rds'))
  ref <- readRDS(cfg$step6$path)
  compare <- function(field,key,score) {
    a <- z$step6@annotation_results[[field]]; b <- ref@annotation_results[[field]]
    make_key <- function(d) do.call(paste,c(d[key],sep='|'))
    expect_equal(anyDuplicated(make_key(a)),anyDuplicated(make_key(b)),info=field)
    expect_identical(sort(make_key(a)),sort(make_key(b)),info=field)
    i <- match(make_key(b),make_key(a)); expect_false(anyNA(i),info=field)
    for (s in score) {
      av <- a[[s]][i]; bv <- b[[s]]
      expect_identical(is.na(av),is.na(bv),info=paste(field,s))
      ok <- is.finite(av) & is.finite(bv)
      err <- if(any(ok)) max(abs(av[ok]-bv[ok])/pmax(abs(bv[ok]),.Machine$double.xmin)) else 0
      expect_lt(err,1e-6)
    }
  }
  compare('flash_entropy_matched_peaks',c('peak_id','query_fragment_id','spectrum_id','library_fragment_id'),c('query_intensity','library_intensity','ppm_error'))
  compare('flash_entropy_spectrum_hits',c('query_id','spectrum_id'),c('S_FE','entropy_S'))
  compare('entropy_candidates',c('query_id','spectrum_id'),c('entropy_score','entropy_S'))
  compare('reverse_match_all',c('query_id','candidate_key'),c('S_FE','S_precursor',
    'reverse_masked_cosine','S_lib_coverage','S_SAS','S_SAS_score_for_total','S_total',
    'S_spatial_reverse','s_rev','score_ms2','score_spatial','score_total_pre_mfc'))
  compare('reverse_match',c('query_id','candidate_key'),c('reverse_masked_cosine','S_SAS'))
  compare('reverse_matched_peaks_all',c('query_id','candidate_key','library_peak_id'),c('library_intensity','query_intensity','reverse_masked_cosine'))
  x <- calculate_mfc(x,output_dir=out)
  observed <- x@process_info$calculate_mfc$result$candidates
  expected <- data.table::fread(cfg$mfc_reference$path,data.table=FALSE)
  key <- function(y) paste(y$feature_id,y$HMDB_ID,y$best_adduct,sep='|')
  expect_equal(nrow(observed),3391)
  expect_equal(length(unique(observed$feature_id)),149)
  expect_equal(nrow(x@variable_info),156)
  expect_identical(key(observed),key(expected))
  expect_identical(as.integer(observed$rank_within_feature),as.integer(expected$rank_within_feature))
  expect_lt(max(abs(observed$mfc_priority_score-expected$mfc_priority_score)/pmax(abs(expected$mfc_priority_score),.Machine$double.xmin)),1e-6)
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
