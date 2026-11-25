# Resolve the split archive layout for users of the deprecated environment variables.
.msi_legacy_reference_paths <- function(resources, combined_root) {
  layouts <- list(
    MSIDATASET_D15_REPORTS_ROOT = c(
      analysis_protocol = "D15/prereg_d15.json",
      annotation_protocol = "D15r/prereg_d15r.json"
    ),
    MSIDATASET_D15_STUDY_ROOT = c(
      csv = "data_ext/samarah_liver/liver_mouse_7_ROI_2.csv",
      labels = "runs/d12_liver/region_d12/liver_mouse_7_ROI_2/labels.tsv",
      meta = "data_ext/samarah_liver/livm_meta.tsv",
      de_reference = "runs/d15_roi2/de_d15/de_features.tsv",
      de_all_reference = "runs/d15_roi2/de_d15/features_all.tsv",
      pseudo_reference = "runs/d15_resume/pseudo_d15/pseudo_ms2.tsv",
      mfc_reference = "runs/d15_resume/mfc_final_d15/mfc_candidates_filtered.tsv",
      step6 = "runs/d15_resume/mfc_retry_d15/step6.rds",
      query = "runs/d15_resume/mfc_retry_d15/pseudo_for_step5.tsv",
      features = "runs/d15_resume/mfc_retry_d15/features.tsv",
      ms2_library = "libs/ms2_neg_d15/lipidblast_hmdb_neg_merged.rds"
    ),
    MSIDATASET_D15_REFERENCE_ROOT = c(
      matching_config = "outputs/nedc_retention_corrected_20260908_153249/nedc_adp_roi1000_20260908_153256/03_resolved_config.txt",
      ms1_library = "derived_library/msi_database_ms1_flat_v2.rds",
      ms1_project_root = "."
    ),
    MSIDATASET_D15_DATA_ROOT = c(
      vector_meta = "hmdb_full_bge_v2/meta.json",
      vector_ids = "hmdb_full_bge_v2/ids.tsv",
      vector_data = "hmdb_full_bge_v2/embeddings.f32"
    )
  )
  for (variable in names(layouts)) {
    root <- .msi_deprecated_env(variable)
    for (key in intersect(names(layouts[[variable]]), names(resources))) {
      modern <- if (nzchar(combined_root)) file.path(combined_root, resources[[key]]$path) else ""
      legacy <- if (nzchar(root)) file.path(root, layouts[[variable]][[key]]) else ""
      resources[[key]]$path <- if (file.exists(legacy)) legacy else modern
    }
  }
  resources
}
