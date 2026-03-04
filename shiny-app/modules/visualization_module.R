visualizationUI <- function(id) {
  ns <- NS(id)

  fluidRow(
    column(
      width = 3,
      wellPanel(
        h4("Map Settings"),
        radioButtons(
          ns("feature_mode"),
          "Feature type",
          choices = c("TIC" = "tic", "m/z" = "mz", "Variable ID" = "variable_id"),
          selected = "tic"
        ),
        numericInput(
          ns("mz_value"),
          label = "m/z value",
          value = 500
        ),
        helpText("Used when Feature type = m/z"),
        selectInput(
          ns("variable_id"),
          "Variable ID",
          choices = character(),
          selected = NULL
        ),
        selectInput(
          ns("section_id"),
          "Section",
          choices = "All",
          selected = "All"
        )
      ),
      wellPanel(
        h4("Spectrum Settings"),
        selectInput(ns("sample_id"), "Pixel sample_id", choices = character(), selected = NULL)
      )
    ),
    column(
      width = 9,
      wellPanel(
        h4("Spatial Intensity Map"),
        plotOutput(ns("spatial_plot"), height = "520px")
      ),
      fluidRow(
        column(
          width = 6,
          wellPanel(
            h4("Average Spectrum"),
            plotOutput(ns("mean_spectrum_plot"), height = "280px")
          )
        ),
        column(
          width = 6,
          wellPanel(
            h4("Selected Pixel Spectrum"),
            plotOutput(ns("pixel_spectrum_plot"), height = "280px")
          )
        )
      )
    )
  )
}

visualizationServer <- function(id, get_object, append_log) {
  moduleServer(id, function(input, output, session) {
    render_plot_error <- function(msg) {
      graphics::plot.new()
      graphics::text(0.5, 0.5, labels = msg, cex = 0.95)
    }

    observeEvent(get_object(), {
      obj <- get_object()

      updateSelectInput(session, "section_id", choices = safe_section_choices(obj), selected = "All")

      var_ids <- safe_variable_ids(obj)
      updateSelectInput(
        session,
        "variable_id",
        choices = var_ids,
        selected = if (length(var_ids) > 0) var_ids[[1]] else NULL
      )

      sample_ids <- safe_sample_ids(obj)
      updateSelectInput(
        session,
        "sample_id",
        choices = sample_ids,
        selected = if (length(sample_ids) > 0) sample_ids[[1]] else NULL
      )

      if (!is.null(obj) && nrow(obj@variable_info) > 0) {
        updateNumericInput(session, "mz_value", value = stats::median(obj@variable_info$mz, na.rm = TRUE))
      }
    })

    output$spatial_plot <- renderPlot({
      obj <- get_object()
      if (is.null(obj)) {
        render_plot_error("Load data in tab 1 first.")
        return(invisible(NULL))
      }

      if (!ensure_tidymass_parameter_class()) {
        render_plot_error("Runtime class initialization failed. Please install/load massdataset.")
        return(invisible(NULL))
      }

      selected_section <- if (identical(input$section_id, "All")) NULL else input$section_id

      p <- tryCatch(
        switch(
          input$feature_mode,
          tic = msidataset::plot_spatial(
            object = obj,
            feature_type = "tic",
            section_id = selected_section,
            facet_by = c("section")
          ),
          mz = msidataset::plot_spatial(
            object = obj,
            feature = input$mz_value,
            feature_type = "mz",
            section_id = selected_section,
            facet_by = c("section")
          ),
          variable_id = {
            req(input$variable_id)
            msidataset::plot_spatial(
              object = obj,
              feature = input$variable_id,
              feature_type = "variable_id",
              section_id = selected_section,
              facet_by = c("section")
            )
          }
        ),
        error = function(e) e
      )

      if (inherits(p, "error")) {
        msg <- paste0("Spatial map error: ", conditionMessage(p))
        append_log(msg)
        render_plot_error(msg)
        return(invisible(NULL))
      }

      print(p)
    })

    output$mean_spectrum_plot <- renderPlot({
      obj <- get_object()
      if (is.null(obj)) {
        render_plot_error("Load data in tab 1 first.")
        return(invisible(NULL))
      }

      res <- tryCatch(msidataset::mean_spectrum(obj, plot = FALSE), error = function(e) e)
      if (inherits(res, "error")) {
        msg <- paste0("Average spectrum error: ", conditionMessage(res))
        append_log(msg)
        render_plot_error(msg)
        return(invisible(NULL))
      }

      ggplot(res, aes(x = mz, y = intensity)) +
        geom_line(color = "#1f4e79", linewidth = 0.4) +
        labs(x = "m/z", y = "Mean intensity") +
        theme_minimal(base_size = 12)
    })

    output$pixel_spectrum_plot <- renderPlot({
      obj <- get_object()
      if (is.null(obj)) {
        render_plot_error("Load data in tab 1 first.")
        return(invisible(NULL))
      }
      req(input$sample_id)

      res <- tryCatch(
        msidataset::pixel_spectrum(obj, sample_id = input$sample_id, plot = FALSE),
        error = function(e) e
      )
      if (inherits(res, "error")) {
        msg <- paste0("Pixel spectrum error: ", conditionMessage(res))
        append_log(msg)
        render_plot_error(msg)
        return(invisible(NULL))
      }

      ggplot(res, aes(x = mz, y = intensity)) +
        geom_line(color = "#2b8cbe", linewidth = 0.35) +
        labs(x = "m/z", y = "Intensity", title = paste("sample_id:", input$sample_id)) +
        theme_minimal(base_size = 12)
    })

    observeEvent(input$feature_mode, {
      append_log(paste0("Visualization feature mode changed to: ", input$feature_mode))
    }, ignoreInit = TRUE)
  })
}
