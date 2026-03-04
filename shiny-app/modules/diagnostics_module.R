diagnosticsUI <- function(id) {
  ns <- NS(id)

  fluidRow(
    column(
      width = 5,
      wellPanel(
        h4("Object Summary"),
        verbatimTextOutput(ns("summary_text"))
      )
    ),
    column(
      width = 7,
      wellPanel(
        h4("Process History (package-level)"),
        tableOutput(ns("process_table"))
      )
    )
  )
}

diagnosticsServer <- function(id, get_object, app_log) {
  moduleServer(id, function(input, output, session) {
    output$summary_text <- renderText({
      obj <- get_object()
      if (is.null(obj)) {
        return("No object loaded.")
      }

      sections <- try(msidataset::get_section_ids(obj), silent = TRUE)
      n_sections <- if (inherits(sections, "try-error") || is.null(sections)) 0 else length(unique(sections))

      paste(
        "Class: msi_dataset", "\n",
        "Dimensions:", nrow(obj@expression_data), "features x", ncol(obj@expression_data), "pixels", "\n",
        "Activated slot:", obj@activated, "\n",
        "Section count:", n_sections, "\n",
        "Process steps:", length(obj@process_info), "\n\n",
        "App activity log:\n",
        paste(utils::tail(app_log(), 15), collapse = "\n")
      )
    })

    output$process_table <- renderTable({
      obj <- get_object()
      req(obj)

      proc <- obj@process_info
      if (length(proc) == 0) {
        return(data.frame(note = "No process_info records available."))
      }

      rows <- lapply(seq_along(proc), function(i) {
        p <- proc[[i]]

        function_name <- NA_character_
        time_value <- NA_character_

        if (methods::is(p, "tidymass_parameter")) {
          function_name <- methods::slot(p, "function_name")
          time_value <- as.character(methods::slot(p, "time"))
        } else if (is.list(p)) {
          function_name <- if (!is.null(p$function_name)) as.character(p$function_name) else "unknown"
          time_value <- if (!is.null(p$time)) as.character(p$time) else NA_character_
        }

        data.frame(
          step = i,
          key = names(proc)[[i]],
          function_name = function_name,
          time = time_value,
          stringsAsFactors = FALSE
        )
      })

      do.call(rbind, rows)
    })
  })
}
