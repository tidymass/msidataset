setMethod(
  f = "show",
  signature = "msi_dataset",
  definition = function(object) {
    ###check again
    check_result = check_msi_dataset(
      expression_data = object@expression_data,
      pixel_info = object@pixel_info,
      sample_info = object@sample_info,
      variable_info = object@variable_info
    )
    if (check_result != "all good.") {
      cat(crayon::red(check_result, "\n"))
      cat(crayon::red(
        "You may changed the slots, try to use update_msi_dataset().\n"
      ))
    }
    
    cat(crayon::yellow(paste(rep("-", 20), collapse = ""), "\n"))
    cat(crayon::green("msidataset version:", object@version, "\n"))
    cat(crayon::yellow(paste(rep("-", 20), collapse = ""), "\n"))
    
    # 1. expression_data
    cat(crayon::green("1.expression_data:"))
    cat(
      "[",
      nrow(object@expression_data),
      "x",
      ncol(object@expression_data),
      "data.frame]\n"
    )
    
    # 2. pixel_info
    cat(crayon::green("2.pixel_info:"))
    cat("[",
        nrow(object@pixel_info),
        "x",
        ncol(object@pixel_info),
        "data.frame]\n")
    if (nrow(object@pixel_info) <= 5) {
      pixel_id <- object@pixel_info$pixel_id
    } else{
      pixel_id <-
        c(head(object@pixel_info$pixel_id, 3), "...",
          tail(object@pixel_info$pixel_id, 2))
    }
    cat(nrow(object@pixel_info), "pixels:")
    cat(pixel_id)
    cat("\n")
    
    # 3. sample_info
    cat(crayon::green("3.sample_info:"))
    cat("[",
        nrow(object@sample_info),
        "x",
        ncol(object@sample_info),
        "data.frame]\n")
    if (nrow(object@sample_info) <= 5) {
      sample_id <- object@sample_info$sample_id
    } else{
      sample_id <-
        c(head(object@sample_info$sample_id, 3), "...",
          tail(object@sample_info$sample_id, 2))
    }
    cat(nrow(object@sample_info), "samples:")
    cat(sample_id)
    cat("\n")
    
    # 4. variable_info
    cat(crayon::green("4.variable_info:"))
    cat(
      "[",
      nrow(object@variable_info),
      "x",
      ncol(object@variable_info),
      "data.frame]\n"
    )
    if (nrow(object@variable_info) <= 5) {
      variable_id <- object@variable_info$variable_id
    } else{
      variable_id <-
        c(head(object@variable_info$variable_id, 3), "...",
          tail(object@variable_info$variable_id, 2))
    }
    cat(nrow(object@variable_info), "variables:")
    cat(variable_id)
    cat("\n")
    
    # 5. pixel_info_note
    cat(crayon::green("5.pixel_info_note:"))
    cat(
      "[",
      nrow(object@pixel_info_note),
      "x",
      ncol(object@pixel_info_note),
      "data.frame]\n"
    )
    
    # 6. sample_info_note
    cat(crayon::green("6.sample_info_note:"))
    cat(
      "[",
      nrow(object@sample_info_note),
      "x",
      ncol(object@sample_info_note),
      "data.frame]\n"
    )
    
    # 7. variable_info_note
    cat(crayon::green("7.variable_info_note:"))
    cat(
      "[",
      nrow(object@variable_info_note),
      "x",
      ncol(object@variable_info_note),
      "data.frame]\n"
    )
    
    # 8. ms2_data
    cat(crayon::green("8.ms2_data:"))
    cat("[",
        sum(unlist(lapply(object@ms2_data, function(x) {
          length(unique(x@variable_id))
        }))),
        "variables x",
        sum(unlist(lapply(object@ms2_data, function(x) {
          length(unique(x@ms2_spectra))
        }))),
        "MS2 spectra]\n")
    
    cat(crayon::yellow(paste(rep("-", 20), collapse = ""), "\n"))
    cat(crayon::green("Processing information\n"))
    if (.hasSlot(object = object, name = "process_info") &
        length(object@process_info) != 0) {
      process_info <- object@process_info
      cat(crayon::green(length(process_info), "processings in total\n"))
      if (length(process_info) > 5) {
        cat(crayon::green("Latest 3 processings show\n"))
        process_info <- tail(process_info, 3)
      }
      for (idx in seq_along(process_info)) {
        cat(crayon::green(names(process_info)[idx], paste(rep("-", 10), collapse = ""), "\n"))
        if (length(process_info[[idx]]) == 1) {
          data.frame(
            "Package" = process_info[[idx]]@pacakge_name,
            "Function used" = process_info[[idx]]@function_name,
            "Time" = process_info[[idx]]@time
          ) %>%
            print()
        } else{
          data.frame(
            "Package" = process_info[[idx]] %>% lapply(function(x)
              x@pacakge_name) %>% unlist(),
            "Function used" = process_info[[idx]] %>% lapply(function(x)
              x@function_name) %>% unlist(),
            "Time" = process_info[[idx]] %>% lapply(function(x)
              as.character(x@time)) %>% unlist()
          ) %>%
            print()
        }
      }
    } else{
      cat(crayon::red("There are no processing for your data.\n"))
    }
  }
)
