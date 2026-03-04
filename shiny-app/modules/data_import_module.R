dataImportUI <- function(id) {
  ns <- NS(id)

  fluidRow(
    column(
      width = 4,
      wellPanel(
        h4("Load RDS Data"),
        p("RDS-only mode for stability and fast startup."),
        actionButton(ns("load_demo_rds"), "Load Demo RDS", class = "btn-primary"),
        tags$hr(),
        fileInput(
          ns("rds_file"),
          "Upload .rds file",
          multiple = FALSE,
          accept = c(".rds")
        ),
        actionButton(ns("import_rds"), "Import Uploaded RDS", class = "btn-success"),
        tags$small("Tip: convert large imzML/ibd to .rds offline, then import here.")
      )
    ),
    column(
      width = 8,
      wellPanel(
        h4("Loading Monitor"),
        uiOutput(ns("loading_ui"))
      ),
      wellPanel(
        h4("Current Dataset Preview"),
        verbatimTextOutput(ns("dataset_dims")),
        h5("sample_info (head)"),
        tableOutput(ns("sample_preview")),
        h5("variable_info (head)"),
        tableOutput(ns("variable_preview"))
      ),
      wellPanel(
        h4("Status"),
        verbatimTextOutput(ns("import_status"))
      )
    )
  )
}

dataImportServer <- function(id, get_object, set_object, append_log) {
  moduleServer(id, function(input, output, session) {
    status <- reactiveVal("No dataset loaded yet.")

    loading <- reactiveValues(
      active = FALSE,
      label = "",
      stage = "Idle",
      started_at = NULL,
      last_success = "No recent load.",
      last_error = ""
    )

    update_status <- function(text, type = "default") {
      status(text)
      append_log(text)
      if (type %in% c("message", "warning", "error")) {
        showNotification(text, type = type)
      }
    }

    begin_loading <- function(label) {
      loading$active <- TRUE
      loading$label <- label
      loading$stage <- "Starting..."
      loading$started_at <- Sys.time()
      loading$last_error <- ""
    }

    end_loading <- function(success = TRUE, message = "") {
      loading$active <- FALSE
      loading$stage <- if (success) "Completed" else "Failed"
      if (success && nzchar(message)) {
        loading$last_success <- message
      }
      if (!success && nzchar(message)) {
        loading$last_error <- message
      }
    }

    output$loading_ui <- renderUI({
      if (isTRUE(loading$active)) {
        invalidateLater(500, session)
        elapsed <- max(0, as.numeric(difftime(Sys.time(), loading$started_at, units = "secs")))
        # Heuristic progress so users can see movement during large RDS unserialization.
        progress_pct <- min(95, round(8 + elapsed * 2.5))

        return(tagList(
          tags$p(tags$strong("Loading in progress"), " - ", loading$label),
          tags$div(
            class = "progress",
            tags$div(
              class = "progress-bar progress-bar-striped active",
              role = "progressbar",
              style = paste0("width: ", progress_pct, "%; min-width: 2em;"),
              paste0(progress_pct, "%")
            )
          ),
          tags$p("Stage: ", loading$stage, tags$br(), "Elapsed: ", sprintf("%.1f", elapsed), " seconds")
        ))
      }

      if (nzchar(loading$last_error)) {
        return(tagList(
          tags$p(tags$strong("Last load failed")),
          tags$p(style = "color:#b30000;", loading$last_error),
          tags$p("You can retry with another RDS file.")
        ))
      }

      tagList(
        tags$p(tags$strong("Idle")),
        tags$p(loading$last_success)
      )
    })

    load_rds_with_feedback <- function(path, success_label) {
      begin_loading(success_label)

      p <- shiny::Progress$new(session = session, min = 0, max = 1)
      on.exit(p$close(), add = TRUE)

      p$set(message = "Loading RDS dataset", value = 0.05, detail = "Preparing")
      loading$stage <- "Reading RDS from disk"

      result <- tryCatch(readRDS(path), error = function(e) e)

      p$set(value = 0.85, detail = "Validating msi_dataset object")
      loading$stage <- "Validating object"

      if (inherits(result, "error") || !validate_msi_object(result)) {
        err_msg <- if (inherits(result, "error")) conditionMessage(result) else "RDS object is not a valid msi_dataset."
        end_loading(success = FALSE, message = err_msg)
        update_status(paste0("Load failed: ", err_msg), type = "error")
        return(invisible(FALSE))
      }

      loading$stage <- "Updating app state"
      ok <- set_object(result, paste0("Loaded msi_dataset from ", success_label))
      p$set(value = 1, detail = "Done")

      if (isTRUE(ok)) {
        end_loading(success = TRUE, message = paste0("Loaded successfully: ", success_label))
        update_status(paste0("Success: ", success_label), type = "message")
      } else {
        end_loading(success = FALSE, message = "Object update was rejected.")
      }

      TRUE
    }

    observeEvent(input$load_demo_rds, {
      path <- find_demo_rds()
      if (is.null(path)) {
        update_status("Demo RDS not found in project data folder.", type = "error")
        return()
      }
      load_rds_with_feedback(path, paste0("demo RDS (", basename(path), ")"))
    })

    observeEvent(input$import_rds, {
      req(input$rds_file)
      load_rds_with_feedback(input$rds_file$datapath, paste0("uploaded RDS (", input$rds_file$name, ")"))
    })

    output$dataset_dims <- renderText({
      obj <- get_object()
      if (is.null(obj)) {
        return("No active object.")
      }

      paste(
        "Class:", class(obj), "\n",
        "Variables (rows):", nrow(obj@expression_data), "\n",
        "Pixels (columns):", ncol(obj@expression_data), "\n",
        "Activated:", obj@activated
      )
    })

    output$sample_preview <- renderTable({
      obj <- get_object()
      req(obj)
      utils::head(obj@sample_info, 8)
    })

    output$variable_preview <- renderTable({
      obj <- get_object()
      req(obj)
      utils::head(obj@variable_info, 8)
    })

    output$import_status <- renderText(status())
  })
}
