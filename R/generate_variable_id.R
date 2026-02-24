# ============================================================================
# Helper Function: Generate variable_id based on m/z values
# ============================================================================

#' Generate Variable IDs Based on m/z Values
#' 
#' @description
#' This function generates variable IDs based on the integer part of m/z values.
#' For example:
#' - m/z = 50.00125 -> m50_1
#' - m/z = 50.00250 -> m50_2
#' - m/z = 51.00123 -> m51_1
#' - m/z = 51.00456 -> m51_2
#' 
#' @param mz_values Numeric vector of m/z values
#' @return Character vector of variable IDs
#' 
#' @examples
#' mz_values <- c(50.00125, 50.00250, 50.00375, 51.00123, 51.00456)
#' generate_variable_id(mz_values)
#' # Returns: "m50_1" "m50_2" "m50_3" "m51_1" "m51_2"
#' 
#' @export
generate_variable_id <- function(mz_values) {
  
  # Get the integer part of m/z values
  mz_integer <- floor(mz_values)
  
  # Create a data frame for processing
  mz_df <- data.frame(
    mz = mz_values,
    mz_int = mz_integer,
    original_order = seq_along(mz_values),
    stringsAsFactors = FALSE
  )
  
  # Sort by m/z value (should already be sorted, but ensure it)
  mz_df <- mz_df[order(mz_df$mz), ]
  
  # Count occurrences within each integer mass group
  mz_df$count <- ave(
    mz_df$mz_int, 
    mz_df$mz_int, 
    FUN = seq_along
  )
  
  # Generate variable_id
  mz_df$variable_id <- paste0("m", mz_df$mz_int, "_", mz_df$count)
  
  # Restore original order
  mz_df <- mz_df[order(mz_df$original_order), ]
  
  # Return variable IDs
  return(mz_df$variable_id)
}


#' Generate Variable IDs with Optional Prefix
#' 
#' @description
#' Extended version that allows custom prefix and formatting options.
#' 
#' @param mz_values Numeric vector of m/z values
#' @param prefix Character string prefix for variable IDs (default: "m")
#' @param separator Character string separator between mass and count (default: "_")
#' @param zero_pad Logical, whether to zero-pad the count numbers (default: FALSE)
#' @param pad_width Integer, width for zero-padding if zero_pad = TRUE (default: 3)
#' 
#' @return Character vector of variable IDs
#' 
#' @examples
#' mz_values <- c(50.00125, 50.00250, 51.00123)
#' 
#' # Default format
#' generate_variable_id_extended(mz_values)
#' # Returns: "m50_1" "m50_2" "m51_1"
#' 
#' # With zero padding
#' generate_variable_id_extended(mz_values, zero_pad = TRUE)
#' # Returns: "m50_001" "m50_002" "m51_001"
#' 
#' # Custom prefix and separator
#' generate_variable_id_extended(mz_values, prefix = "mass", separator = "-")
#' # Returns: "mass50-1" "mass50-2" "mass51-1"
#' 
#' @export
generate_variable_id_extended <- function(mz_values, 
                                          prefix = "m",
                                          separator = "_",
                                          zero_pad = FALSE,
                                          pad_width = 3) {
  
  # Get the integer part of m/z values
  mz_integer <- floor(mz_values)
  
  # Create a data frame for processing
  mz_df <- data.frame(
    mz = mz_values,
    mz_int = mz_integer,
    original_order = seq_along(mz_values),
    stringsAsFactors = FALSE
  )
  
  # Sort by m/z value
  mz_df <- mz_df[order(mz_df$mz), ]
  
  # Count occurrences within each integer mass group
  mz_df$count <- ave(
    mz_df$mz_int, 
    mz_df$mz_int, 
    FUN = seq_along
  )
  
  # Format count numbers
  if (zero_pad) {
    count_str <- sprintf(paste0("%0", pad_width, "d"), mz_df$count)
  } else {
    count_str <- as.character(mz_df$count)
  }
  
  # Generate variable_id
  mz_df$variable_id <- paste0(prefix, mz_df$mz_int, separator, count_str)
  
  # Restore original order
  mz_df <- mz_df[order(mz_df$original_order), ]
  
  # Return variable IDs
  return(mz_df$variable_id)
}


#' Validate Variable ID Uniqueness
#' 
#' @description
#' Check if generated variable IDs are unique.
#' 
#' @param variable_ids Character vector of variable IDs
#' @param mz_values Numeric vector of m/z values (for informative error messages)
#' 
#' @return Logical TRUE if all unique, otherwise throws an error
#' 
#' @export
validate_variable_id <- function(variable_ids, mz_values = NULL) {
  
  if (any(duplicated(variable_ids))) {
    dup_ids <- unique(variable_ids[duplicated(variable_ids)])
    
    error_msg <- paste0(
      "Duplicate variable IDs detected: ",
      paste(head(dup_ids, 5), collapse = ", ")
    )
    
    if (!is.null(mz_values)) {
      error_msg <- paste0(
        error_msg,
        "\nThis should not happen. Please check m/z values."
      )
    }
    
    stop(error_msg)
  }
  
  return(TRUE)
}