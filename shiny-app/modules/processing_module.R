processingUI <- function(id) {
  ns <- NS(id)

  fluidRow(
    column(
      width = 4,
      wellPanel(
        h4("Activation"),
        selectInput(
          ns("activate_what"),
          "Activate slot",
          choices = c("sample_info", "variable_info", "expression_data", "annotation_table"),
          selected = "expression_data"
        ),
        actionButton(ns("apply_activation"), "Apply Activation", class = "btn-primary")
      ),
      wellPanel(
        h4("Spatial Filter"),
        p("Keep pixels within x/y bounds."),
        numericInput(ns("x_min"), "x_min", value = NA),
        numericInput(ns("x_max"), "x_max", value = NA),
        numericInput(ns("y_min"), "y_min", value = NA),
        numericInput(ns("y_max"), "y_max", value = NA),
        actionButton(ns("apply_spatial_filter"), "Apply Spatial Filter", class = "btn-warning")
      )
    ),
    column(
      width = 4,
      wellPanel(
        h4("m/z Range Filter"),
        p("Keep features in selected m/z range."),
        numericInput(ns("mz_min"), "mz_min", value = NA),
        numericInput(ns("mz_max"), "mz_max", value = NA),
        actionButton(ns("apply_mz_filter"), "Apply m/z Filter", class = "btn-warning")
      )
    ),
    column(
      width = 4,
      wellPanel(
        h4("Processing Summary"),
        verbatimTextOutput(ns("processing_status"))
      )
    )
  )
}

processingServer <- function(id, get_object, set_object, append_log) {
  moduleServer(id, function(input, output, session) {
    status <- reactiveVal("No processing action yet.")

    run_processing_step <- function(step_label, fn) {
      obj <- get_object()
      if (is.null(obj)) {
        msg <- "Load data in tab 1 first."
        status(msg)
        showNotification(msg, type = "warning")
        return(invisible(FALSE))
      }

      withProgress(message = step_label, value = 0.15, {
        incProgress(0.3)
        result <- tryCatch(fn(obj), error = function(e) e)
        incProgress(0.7)

        if (inherits(result, "error")) {
          err <- conditionMessage(result)
          msg <- paste0(step_label, " failed: ", err)
          status(msg)
          append_log(msg)
          showNotification(msg, type = "error")
          return(invisible(FALSE))
        }

        if (!validate_msi_object(result)) {
          msg <- paste0(step_label, " produced an invalid object.")
          status(msg)
          append_log(msg)
          showNotification(msg, type = "error")
          return(invisible(FALSE))
        }

        set_object(result, paste0(step_label, " completed."))
        msg <- paste0(step_label, " completed successfully.")
        status(msg)
        append_log(msg)
        showNotification(msg, type = "message")
        TRUE
      })
    }

    observeEvent(get_object(), {
      obj <- get_object()
      req(obj)

      xr <- range(obj@sample_info$x, na.rm = TRUE)
      yr <- range(obj@sample_info$y, na.rm = TRUE)
      mzr <- range(obj@variable_info$mz, na.rm = TRUE)

      updateNumericInput(session, "x_min", value = xr[[1]])
      updateNumericInput(session, "x_max", value = xr[[2]])
      updateNumericInput(session, "y_min", value = yr[[1]])
      updateNumericInput(session, "y_max", value = yr[[2]])
      updateNumericInput(session, "mz_min", value = mzr[[1]])
      updateNumericInput(session, "mz_max", value = mzr[[2]])
    })

    observeEvent(input$apply_activation, {
      run_processing_step("Activation", function(obj) {
        msidataset::activate_msi_dataset(obj, what = input$activate_what)
      })
    })

    observeEvent(input$apply_spatial_filter, {
      req(is.finite(input$x_min), is.finite(input$x_max), is.finite(input$y_min), is.finite(input$y_max))
      if (input$x_min > input$x_max || input$y_min > input$y_max) {
        msg <- "Invalid spatial bounds: min cannot be greater than max."
        status(msg)
        showNotification(msg, type = "warning")
        return()
      }

      run_processing_step("Spatial filtering", function(obj) {
        msidataset::filter_msi_dataset(
          object = obj,
          target = "sample_info",
          filter_type = "spatial",
          x_min = input$x_min,
          x_max = input$x_max,
          y_min = input$y_min,
          y_max = input$y_max
        )
      })
    })

    observeEvent(input$apply_mz_filter, {
      req(is.finite(input$mz_min), is.finite(input$mz_max))
      if (input$mz_min > input$mz_max) {
        msg <- "Invalid m/z range: mz_min cannot be greater than mz_max."
        status(msg)
        showNotification(msg, type = "warning")
        return()
      }

      run_processing_step("m/z filtering", function(obj) {
        msidataset::filter_msi_dataset(
          object = obj,
          target = "variable_info",
          filter_type = "mz_range",
          mz_min = input$mz_min,
          mz_max = input$mz_max
        )
      })
    })

    output$processing_status <- renderText({
      obj <- get_object()
      if (is.null(obj)) {
        return("Load data in tab 1 first.")
      }

      paste(
        status(), "\n",
        "Current dimensions:",
        nrow(obj@expression_data), "features x",
        ncol(obj@expression_data), "pixels"
      )
    })
  })
}
