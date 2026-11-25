# Internal msi_dataset class needed to deserialize compatible step6 inputs.
methods::setClass(
  Class = "msi_dataset",
  methods::representation(
    expression_data  = "ANY",
    ms2_data         = "list",
    feature_graph    = "list",
    pseudo_ms2       = "list",
    spectral_library = "list",
    annotation_results = "list",
    archive_index    = "list",
    msi_parameters   = "list",
    annotation_table = "data.frame",
    sample_info      = "data.frame",
    variable_info    = "data.frame",
    sample_info_note     = "data.frame",
    variable_info_note   = "data.frame",
    process_info     = "list",
    other_files      = "list",
    version          = "character",
    activated        = "character"
  ),
  methods::prototype(
    expression_data      = matrix(numeric(0), nrow = 0, ncol = 0),
    ms2_data             = list(),
    feature_graph        = list(),
    pseudo_ms2           = list(),
    spectral_library     = list(),
    annotation_results   = list(),
    archive_index        = list(),
    msi_parameters       = list(),
    annotation_table     = data.frame(
      anchor_feature_id = character(), candidate_id = character(),
      candidate_name = character(), spectrum_id = character(),
      annotation_status = character(), stringsAsFactors = FALSE,
      check.names = FALSE
    ),
    sample_info          = data.frame(),
    variable_info        = data.frame(),
    sample_info_note     = data.frame(),
    variable_info_note   = data.frame(),
    process_info         = list(),
    other_files          = list(),
    version              = "0.1.0",
    activated            = "expression_data"
  )
)

methods::setValidity(
  "msi_dataset",
  function(object) {
    TRUE
  }
)
