#' Display Processing Information with Details
#'
#' @description
#' Show all processing steps applied to msi_dataset object with detailed information:
#' - Function names
#' - Timestamps
#' - Parameters used
#' - Data statistics at each step
#'
#' @param object A msi_dataset object
#' @param detailed If TRUE, show full parameter details
#' @param verbose If TRUE, print step-by-step summary
#'
#' @return Invisible data frame with processing history
#' @export
#'
#' @examples
#' \dontrun{
#' show_process_info(msi_dataset, detailed = TRUE, verbose = TRUE)
#' }

show_process_info <- function(object, detailed = TRUE, verbose = TRUE) {
  
  if (!methods::is(object, "msi_dataset"))
    stop("object must be a msi_dataset.")
  
  if (length(object@process_info) == 0) {
    cat("No processing steps recorded.\n")
    return(invisible(NULL))
  }
  
  cat("\n")
  cat("╔═══════════════════════════════════════════════════════════════════╗\n")
  cat("║                    Processing History Summary                      ║\n")
  cat("╚═══════════════════════════════════════════════════════════════════╝\n\n")
  
  # Extract process info
  process_list <- object@process_info
  process_keys <- names(process_list)
  
  # Create summary dataframe
  process_summary <- data.frame(
    Step = seq_along(process_keys),
    Timestamp = character(length(process_keys)),
    Function = character(length(process_keys)),
    Key_Parameters = character(length(process_keys)),
    stringsAsFactors = FALSE
  )
  
  for (i in seq_along(process_keys)) {
    param_obj <- process_list[[process_keys[i]]]
    
    # Extract timestamp
    timestamp <- format(param_obj@time, "%Y-%m-%d %H:%M:%S")
    process_summary$Timestamp[i] <- timestamp
    
    # Extract function name
    func_name <- param_obj@function_name
    process_summary$Function[i] <- func_name
    
    # Extract key parameters
    if (length(param_obj@parameter) > 0) {
      params <- param_obj@parameter
      # Format key parameters concisely
      key_params <- character()
      
      # Identify and format important parameters
      if (!is.null(params$zero_threshold)) {
        key_params <- c(key_params, sprintf("zero_threshold=%.1f%%", params$zero_threshold * 100))
      }
      if (!is.null(params$min_pixels)) {
        key_params <- c(key_params, sprintf("min_pixels=%d", params$min_pixels))
      }
      if (!is.null(params$min_intensity)) {
        key_params <- c(key_params, sprintf("min_intensity=%.1f", params$min_intensity))
      }
      if (!is.null(params$method)) {
        key_params <- c(key_params, sprintf("method=%s", params$method))
      }
      if (!is.null(params$num_retained)) {
        key_params <- c(key_params, sprintf("retained=%d", params$num_retained))
        if (!is.null(params$num_removed)) {
          key_params <- c(key_params, sprintf("removed=%d", params$num_removed))
        }
      }
      
      if (length(key_params) > 0) {
        process_summary$Key_Parameters[i] <- paste(key_params, collapse=" | ")
      } else {
        process_summary$Key_Parameters[i] <- "See details below"
      }
    }
  }
  
  # Print summary table
  cat("Summary Table:\n")
  cat("─────────────────────────────────────────────────────────────────────\n")
  
  for (i in 1:nrow(process_summary)) {
    cat(sprintf("%d. [%s]\n", process_summary$Step[i], process_summary$Timestamp[i]))
    cat(sprintf("   Function: %s\n", process_summary$Function[i]))
    cat(sprintf("   Parameters: %s\n", process_summary$Key_Parameters[i]))
    cat("\n")
  }
  
  # Print detailed information if requested
  if (detailed) {
    cat("\n")
    cat("╔═══════════════════════════════════════════════════════════════════╗\n")
    cat("║                    Detailed Parameter Information                 ║\n")
    cat("╚═══════════════════════════════════════════════════════════════════╝\n\n")
    
    for (i in seq_along(process_keys)) {
      param_obj <- process_list[[process_keys[i]]]
      
      cat(sprintf("Step %d: %s\n", i, param_obj@function_name))
      cat(sprintf("├─ Timestamp: %s\n", format(param_obj@time, "%Y-%m-%d %H:%M:%S")))
      cat(sprintf("├─ Package: %s\n", param_obj@pacakge_name))
      cat(sprintf("└─ Parameters:\n"))
      
      if (length(param_obj@parameter) > 0) {
        params <- param_obj@parameter
        param_names <- names(params)
        
        for (j in seq_along(param_names)) {
          param_name <- param_names[j]
          param_value <- params[[param_name]]
          
          # Format parameter value nicely
          if (is.numeric(param_value)) {
            if (param_name %in% c("zero_threshold", "frequency_threshold")) {
              formatted_value <- sprintf("%.1f%% (%.4f)", param_value * 100, param_value)
            } else {
              formatted_value <- sprintf("%d", param_value)
            }
          } else {
            formatted_value <- as.character(param_value)
          }
          
          is_last <- (j == length(param_names))
          prefix <- ifelse(is_last, "    └─ ", "    ├─ ")
          cat(sprintf("%s%s: %s\n", prefix, param_name, formatted_value))
        }
      }
      
      cat("\n")
    }
  }
  
  # Print verbose summary if requested
  if (verbose) {
    cat("\n")
    cat("╔═══════════════════════════════════════════════════════════════════╗\n")
    cat("║                    Processing Flow Analysis                       ║\n")
    cat("╚═══════════════════════════════════════════════════════════════════╝\n\n")
    
    # Calculate overall statistics
    cat("Overall Processing Summary:\n")
    cat("─────────────────────────────────────────────────────────────────────\n")
    cat(sprintf("Total steps performed: %d\n", length(process_keys)))
    
    if (nrow(process_summary) > 0) {
      start_time <- process_list[[1]]@time
      end_time <- process_list[[length(process_list)]]@time
      duration <- difftime(end_time, start_time, units = "mins")
      
      cat(sprintf("Processing duration: %.2f minutes\n", as.numeric(duration)))
      cat(sprintf("First step: %s\n", process_summary$Function[1]))
      cat(sprintf("Last step: %s\n", process_summary$Function[nrow(process_summary)]))
    }
    
    cat(sprintf("\nFinal data dimensions:\n"))
    cat(sprintf("  - m/z features: %d\n", nrow(object@variable_info)))
    cat(sprintf("  - Pixels/samples: %d\n", ncol(object@expression_data)))
    
    # Calculate retention rate
    cat("\n")
    cat("Feature Retention by Step:\n")
    cat("─────────────────────────────────────────────────────────────────────\n")
    
    for (i in seq_along(process_keys)) {
      param_obj <- process_list[[process_keys[i]]]
      params <- param_obj@parameter
      
      if (!is.null(params$num_retained) && !is.null(params$num_removed)) {
        total <- params$num_retained + params$num_removed
        retention_rate <- params$num_retained / total * 100
        
        cat(sprintf("Step %d (%s):\n", i, param_obj@function_name))
        cat(sprintf("  Retained: %d | Removed: %d | Retention: %.1f%%\n",
                    params$num_retained, params$num_removed, retention_rate))
      }
    }
  }
  
  cat("\n")
  
  # Return invisible dataframe
  invisible(process_summary)
}


# ════════════════════════════════════════════════════════════════════════════════
# 创建处理流程可视化
# ════════════════════════════════════════════════════════════════════════════════

#' Create Processing Timeline Visualization
#'
#' @description
#' Generate a visual timeline of all processing steps with timestamps
#'
#' @param object A msi_dataset object
#'
#' @return NULL (invisible). Prints a formatted timeline.
#' @export

show_process_timeline <- function(object) {
  
  if (!methods::is(object, "msi_dataset"))
    stop("object must be a msi_dataset.")
  
  if (length(object@process_info) == 0) {
    cat("No processing steps recorded.\n")
    return(invisible(NULL))
  }
  
  cat("\n")
  cat("╔═══════════════════════════════════════════════════════════════════╗\n")
  cat("║                    Processing Timeline                            ║\n")
  cat("╚═══════════════════════════════════════════════════════════════════╝\n\n")
  
  process_list <- object@process_info
  process_keys <- names(process_list)
  
  # Extract times
  times <- sapply(process_list, function(x) x@time)
  
  # Calculate elapsed time from start
  start_time <- times[1]
  
  for (i in seq_along(process_keys)) {
    param_obj <- process_list[[process_keys[i]]]
    elapsed <- difftime(param_obj@time, start_time, units = "secs")
    
    # Format output
    timestamp <- format(param_obj@time, "%H:%M:%S")
    func_short <- gsub("\\(\\)", "", param_obj@function_name)
    
    # Draw timeline
    cat(sprintf("[%02d:%02d:%02d] Step %d\n",
                as.integer(elapsed) %/% 3600,
                (as.integer(elapsed) %% 3600) %/% 60,
                as.integer(elapsed) %% 60,
                i))
    cat(sprintf("    └─ %s\n", func_short))
    cat(sprintf("       Time: %s\n", timestamp))
    
    # Add connection line
    if (i < length(process_keys)) {
      cat("    ↓\n")
    }
    
    cat("\n")
  }
  
  # Total time summary
  total_time <- difftime(times[length(times)], times[1], units = "secs")
  cat(sprintf("Total processing time: %.1f seconds\n", as.numeric(total_time)))
  cat("\n")
  
  invisible(NULL)
}


# ════════════════════════════════════════════════════════════════════════════════
# 导出处理历史为表格
# ════════════════════════════════════════════════════════════════════════════════

#' Export Processing History to Data Frame
#'
#' @description
#' Export complete processing information to a data frame for further analysis or export
#'
#' @param object A msi_dataset object
#' @param include_details If TRUE, include all parameter details
#'
#' @return Data frame with processing history
#' @export

export_process_history <- function(object, include_details = FALSE) {
  
  if (!methods::is(object, "msi_dataset"))
    stop("object must be a msi_dataset.")
  
  if (length(object@process_info) == 0) {
    cat("No processing steps recorded.\n")
    return(invisible(NULL))
  }
  
  process_list <- object@process_info
  process_keys <- names(process_list)
  
  # Create base dataframe
  history_df <- data.frame(
    Step = seq_along(process_keys),
    Timestamp = character(length(process_keys)),
    Date = character(length(process_keys)),
    Time = character(length(process_keys)),
    Function = character(length(process_keys)),
    Package = character(length(process_keys)),
    stringsAsFactors = FALSE
  )
  
  for (i in seq_along(process_keys)) {
    param_obj <- process_list[[process_keys[i]]]
    timestamp <- param_obj@time
    
    history_df$Timestamp[i] <- format(timestamp, "%Y-%m-%d %H:%M:%S")
    history_df$Date[i] <- format(timestamp, "%Y-%m-%d")
    history_df$Time[i] <- format(timestamp, "%H:%M:%S")
    history_df$Function[i] <- param_obj@function_name
    history_df$Package[i] <- param_obj@pacakge_name
    
    # Add parameter details if requested
    if (include_details) {
      params <- param_obj@parameter
      
      # Serialize parameters to JSON-like format
      if (length(params) > 0) {
        param_strings <- sapply(names(params), function(x) {
          val <- params[[x]]
          if (is.numeric(val)) {
            sprintf("%s=%.2f", x, val)
          } else {
            sprintf("%s=%s", x, as.character(val))
          }
        })
        history_df[[sprintf("Params_Step%d", i)]] <- paste(param_strings, collapse="; ")
      }
    }
  }
  
  return(history_df)
}
