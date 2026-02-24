# ============================================================================
# All Generic Definitions for msidataset
# ============================================================================

# ---------------------------------------------------------------------------
# Slot getters
# ---------------------------------------------------------------------------

setGeneric("expression_data",
           function(object, ...) standardGeneric("expression_data"))

setGeneric("sample_info",
           function(object, ...) standardGeneric("sample_info"))

setGeneric("variable_info",
           function(object, ...) standardGeneric("variable_info"))

setGeneric("process_info",
           function(object, ...) standardGeneric("process_info"))

setGeneric("ms2_data",
           function(object, ...) standardGeneric("ms2_data"))

setGeneric("annotation_table",
           function(object, ...) standardGeneric("annotation_table"))

setGeneric("sample_info_note",
           function(object, ...) standardGeneric("sample_info_note"))

setGeneric("variable_info_note",
           function(object, ...) standardGeneric("variable_info_note"))

# ---------------------------------------------------------------------------
# Slot setters (replacement methods)
# ---------------------------------------------------------------------------

setGeneric("expression_data<-",
           function(object, value) standardGeneric("expression_data<-"))

setGeneric("sample_info<-",
           function(object, value) standardGeneric("sample_info<-"))

setGeneric("variable_info<-",
           function(object, value) standardGeneric("variable_info<-"))

setGeneric("ms2_data<-",
           function(object, value) standardGeneric("ms2_data<-"))

setGeneric("annotation_table<-",
           function(object, value) standardGeneric("annotation_table<-"))

setGeneric("sample_info_note<-",
           function(object, value) standardGeneric("sample_info_note<-"))

setGeneric("variable_info_note<-",
           function(object, value) standardGeneric("variable_info_note<-"))

# ---------------------------------------------------------------------------
# Convenience count / summary generics
# ---------------------------------------------------------------------------

setGeneric("n_variables",
           function(object, ...) standardGeneric("n_variables"))

setGeneric("n_samples",
           function(object, ...) standardGeneric("n_samples"))

setGeneric("n_pixels",
           function(object, ...) standardGeneric("n_pixels"))

setGeneric("mz",
           function(object, ...) standardGeneric("mz"))

# ---------------------------------------------------------------------------
# Backward-compatibility generics (pixel_info alias)
# ---------------------------------------------------------------------------

setGeneric("pixel_info",
           function(object, ...) standardGeneric("pixel_info"))

setGeneric("pixel_info_note",
           function(object, ...) standardGeneric("pixel_info_note"))

# ---------------------------------------------------------------------------
# Spatial helper generics
# ---------------------------------------------------------------------------

setGeneric("get_spatial_coords",
           function(object, ...) standardGeneric("get_spatial_coords"))

setGeneric("get_section_ids",
           function(object, ...) standardGeneric("get_section_ids"))

setGeneric("get_section_pixels",
           function(object, section_id, ...) standardGeneric("get_section_pixels"))

setGeneric("get_sample_id",
           function(object, ...) standardGeneric("get_sample_id"))

setGeneric("get_variable_id",
           function(object, ...) standardGeneric("get_variable_id"))
