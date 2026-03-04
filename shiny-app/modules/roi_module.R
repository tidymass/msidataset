roiUI <- function(id) {
  ns <- NS(id)

  fluidRow(
    column(
      width = 5,
      wellPanel(
        h4("ROI Module"),
        p("Register a basic ROI from the current section to enable downstream ROI workflows."),
        textInput(ns("roi_id"), "ROI ID", value = "roi_1"),
        textInput(ns("roi_name"), "ROI Name", value = "Tumor Region"),
        selectInput(ns("roi_section"), "Section", choices = "All", selected = "All"),
        numericInput(ns("roi_max_pixels"), "Max pixels (0 = all)", value = 500, min = 0, step = 100),
        actionButton(ns("register_roi"), "Register ROI", class = "btn-info"),
        actionButton(ns("refresh_roi"), "Refresh ROI List", class = "btn-default")
      )
    ),
    column(
      width = 7,
      wellPanel(
        h4("Existing ROI Summary"),
        textOutput(ns("roi_count")),
        tableOutput(ns("roi_table"))
      ),
      wellPanel(
        h4("ROI Membership Preview"),
        selectInput(ns("roi_pick"), "Preview ROI", choices = character(), selected = NULL),
        tableOutput(ns("roi_membership_preview")),
        verbatimTextOutput(ns("roi_status"))
      )
    )
  )
}

roiServer <- function(id, get_object, set_object, append_log) {
  moduleServer(id, function(input, output, session) {
    status <- reactiveVal("No ROI action yet.")

    refresh_roi_views <- function(obj) {
      roi_tbl <- tryCatch(msidataset::list_roi(obj), error = function(e) data.frame())
      roi_ids <- if (nrow(roi_tbl) > 0) as.character(roi_tbl$roi_id) else character()
      updateSelectInput(
        session,
        "roi_pick",
        choices = roi_ids,
        selected = if (length(roi_ids) > 0) roi_ids[[1]] else NULL
      )
      invisible(roi_tbl)
    }

    observeEvent(get_object(), {
      obj <- get_object()
      req(obj)
      updateSelectInput(session, "roi_section", choices = safe_section_choices(obj), selected = "All")
      refresh_roi_views(obj)
    })

    observeEvent(input$refresh_roi, {
      obj <- get_object()
      req(obj)
      refresh_roi_views(obj)
      status("ROI list refreshed.")
      showNotification("ROI list refreshed.", type = "message")
    })

    observeEvent(input$register_roi, {
      obj <- get_object()
      if (is.null(obj)) {
        status("Load data first before registering ROI.")
        showNotification(status(), type = "warning")
        return()
      }

      roi_id <- trimws(input$roi_id)
      roi_name <- trimws(input$roi_name)
      if (!nzchar(roi_id)) {
        status("ROI ID is required.")
        showNotification(status(), type = "warning")
        return()
      }
      if (!nzchar(roi_name)) {
        roi_name <- roi_id
      }

      sample_tbl <- obj@sample_info
      if (!identical(input$roi_section, "All")) {
        sample_tbl <- sample_tbl[sample_tbl$section_id %in% input$roi_section, , drop = FALSE]
      }

      if (nrow(sample_tbl) == 0) {
        status("No pixels available for the selected section.")
        showNotification(status(), type = "error")
        return()
      }

      max_pixels <- as.integer(input$roi_max_pixels)
      if (is.na(max_pixels) || max_pixels < 0) {
        max_pixels <- 0L
      }

      sample_ids <- as.character(sample_tbl$sample_id)
      if (max_pixels > 0 && length(sample_ids) > max_pixels) {
        sample_ids <- sample_ids[seq_len(max_pixels)]
      }

      new_obj <- tryCatch(
        msidataset::add_roi(
          object = obj,
          roi_id = roi_id,
          sample_id = sample_ids,
          roi_name = roi_name,
          section_id = if (identical(input$roi_section, "All")) NA_character_ else as.character(input$roi_section),
          overwrite = TRUE
        ),
        error = function(e) e
      )

      if (inherits(new_obj, "error")) {
        msg <- paste0("ROI registration failed: ", conditionMessage(new_obj))
        status(msg)
        append_log(msg)
        showNotification(msg, type = "error")
        return()
      }

      set_object(new_obj, paste0("Registered ROI: ", roi_id, " (", length(sample_ids), " pixels)"))
      refresh_roi_views(new_obj)
      msg <- paste0("ROI registered: ", roi_id, " with ", length(sample_ids), " pixels.")
      status(msg)
      append_log(msg)
      showNotification(msg, type = "message")
    })

    output$roi_table <- renderTable({
      obj <- get_object()
      req(obj)

      roi_tbl <- tryCatch(msidataset::list_roi(obj), error = function(e) NULL)
      if (is.null(roi_tbl) || nrow(as.data.frame(roi_tbl)) == 0) {
        return(data.frame(note = "No ROI stored yet."))
      }
      as.data.frame(roi_tbl)
    })

    output$roi_count <- renderText({
      obj <- get_object()
      req(obj)
      roi_tbl <- tryCatch(msidataset::list_roi(obj), error = function(e) data.frame())
      paste0("ROI count: ", nrow(roi_tbl))
    })

    output$roi_membership_preview <- renderTable({
      obj <- get_object()
      req(obj)
      req(input$roi_pick)

      membership <- tryCatch(msidataset::roi_data(obj, roi_id = input$roi_pick), error = function(e) NULL)
      if (is.null(membership) || nrow(as.data.frame(membership)) == 0) {
        return(data.frame(note = "No ROI membership data available."))
      }

      utils::head(as.data.frame(membership), 20)
    })

    output$roi_status <- renderText(status())
  })
}
