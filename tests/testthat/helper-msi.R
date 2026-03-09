make_test_msi_object <- function(n_features = 5L, n_samples = 8L, n_sections = 2L) {
  stopifnot(n_samples %% n_sections == 0)

  mz_seed <- c(130.06, 191.02, 250.10, 306.07, 401.20, 550.30, 700.40, 850.50)
  mz_vals <- mz_seed[seq_len(n_features)]
  expr <- matrix(seq_len(n_features * n_samples), nrow = n_features, ncol = n_samples)

  per_section <- n_samples / n_sections
  section_id <- rep(paste0("section_", seq_len(n_sections)), each = per_section)
  z_vals <- rep(seq_len(n_sections), each = per_section)
  x_vals <- rep(seq_len(per_section), times = n_sections)
  y_vals <- rep(seq_len(n_sections), each = per_section)

  sample_info <- data.frame(
    sample_id = paste0("spectrum_", seq_len(n_samples)),
    x = as.numeric(x_vals),
    y = as.numeric(y_vals),
    z = as.numeric(z_vals),
    subject_id = "subject_test",
    section_id = section_id,
    resolution = 50,
    matrix = "CHCA",
    stringsAsFactors = FALSE,
    check.names = FALSE
  )

  variable_info <- data.frame(
    variable_id = paste0("m", seq_len(n_features)),
    mz = mz_vals,
    ms2 = rep("no", n_features),
    has_ms2 = rep("no", n_features),
    ms2_id = rep(NA_character_, n_features),
    ms2_type = rep(NA_character_, n_features),
    stringsAsFactors = FALSE,
    check.names = FALSE
  )

  rownames(expr) <- variable_info$variable_id
  colnames(expr) <- sample_info$sample_id

  create_msi_dataset(
    expression_data = expr,
    sample_info = sample_info,
    variable_info = variable_info
  )
}

add_test_rois <- function(object) {
  ids <- object@sample_info$sample_id
  sec <- object@sample_info$section_id
  sec1 <- ids[sec == unique(sec)[1]]
  sec2 <- ids[sec == unique(sec)[2]]

  object <- add_roi(
    object = object,
    roi_id = "roi_1",
    roi_name = "ROI 1",
    sample_id = sec1[seq_len(min(2L, length(sec1)))],
    section_id = unique(sec)[1],
    overwrite = TRUE
  )

  object <- add_roi(
    object = object,
    roi_id = "roi_2",
    roi_name = "ROI 2",
    sample_id = sec2[seq_len(min(2L, length(sec2)))],
    section_id = unique(sec)[2],
    overwrite = TRUE
  )
  object
}

make_temp_mgf <- function(precursor_mz = c(191.02, 306.07),
                          scans = c("spectrum_1", "spectrum_2")) {
  stopifnot(length(precursor_mz) == length(scans))
  mgf <- tempfile(fileext = ".mgf")
  lines <- character()
  for (i in seq_along(precursor_mz)) {
    lines <- c(
      lines,
      "BEGIN IONS",
      paste0("TITLE=test_", i),
      paste0("PEPMASS=", sprintf("%.4f", precursor_mz[[i]])),
      paste0("RTINSECONDS=", 10 + i),
      paste0("SCANS=", scans[[i]]),
      "100.0 200",
      "120.0 50",
      "END IONS"
    )
  }
  writeLines(lines, mgf)
  mgf
}

last_process_fun <- function(object) {
  object@process_info[[length(object@process_info)]]@function_name
}

last_process_log <- function(object) {
  object@process_info[[length(object@process_info)]]
}

expect_valid_msi <- function(object) {
  testthat::expect_true(isTRUE(methods::validObject(object)))
  testthat::expect_s4_class(object, "msi_dataset")
}

make_targeted_pseudo_test_object <- function() {
  sample_info <- data.frame(
    sample_id = paste0("spectrum_", seq_len(8)),
    x = rep(1:2, times = 4),
    y = rep(1:2, each = 2, times = 2),
    z = rep(1:2, each = 4),
    subject_id = "subject_test",
    section_id = rep(c("section_1", "section_2"), each = 4),
    resolution = 50,
    matrix = "CHCA",
    stringsAsFactors = FALSE,
    check.names = FALSE
  )

  precursor_a <- c(10, 12, 1, 0, 9, 11, 1, 0)
  fragment_a1 <- c(9, 11, 1, 0, 8, 10, 1, 0)
  fragment_a2 <- c(8, 10, 0.5, 0, 7, 9, 0.5, 0)
  precursor_b <- c(0, 1, 9, 11, 0, 1, 10, 12)
  fragment_b1 <- c(0, 1, 8, 10, 0, 1, 9, 11)
  noise <- c(3, 1, 6, 4, 2, 7, 1, 5)

  expr <- rbind(
    precursor_a,
    fragment_a1,
    fragment_a2,
    precursor_b,
    fragment_b1,
    noise
  )
  variable_info <- data.frame(
    variable_id = paste0("m", seq_len(nrow(expr))),
    mz = c(151.02, 120.05, 90.03, 289.13, 170.08, 240.2),
    ms2 = rep("no", nrow(expr)),
    has_ms2 = rep("no", nrow(expr)),
    ms2_id = rep(NA_character_, nrow(expr)),
    ms2_type = rep(NA_character_, nrow(expr)),
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
  rownames(expr) <- variable_info$variable_id
  colnames(expr) <- sample_info$sample_id

  create_msi_dataset(
    expression_data = expr,
    sample_info = sample_info,
    variable_info = variable_info
  )
}
