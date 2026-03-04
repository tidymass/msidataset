library(shiny)
library(bslib)
library(ggplot2)
library(massdataset)
library(msidataset)
library(methods)

source(file.path("modules", "utils.R"), local = TRUE)
source(file.path("modules", "data_import_module.R"), local = TRUE)
source(file.path("modules", "processing_module.R"), local = TRUE)
source(file.path("modules", "visualization_module.R"), local = TRUE)
source(file.path("modules", "roi_module.R"), local = TRUE)
source(file.path("modules", "diagnostics_module.R"), local = TRUE)

ui <- fluidPage(
  theme = bs_theme(
    version = 5,
    bootswatch = "minty"
  ),
  titlePanel("msidataset Shiny Explorer"),
  p(
    "Beginner-friendly interface for MSI workflows.",
    tags$br(),
    tags$small("Author: Xu Xiaoyu | Email: XIAOYU021@e.ntu.edu.sg")
  ),
  tabsetPanel(
    id = "main_tabs",
    tabPanel("1) Data Import", dataImportUI("import")),
    tabPanel("2) Basic Processing", processingUI("processing")),
    tabPanel("3) Spatial Visualization", visualizationUI("viz")),
    tabPanel("4) ROI (Future-Ready)", roiUI("roi")),
    tabPanel("5) Diagnostics / Summary", diagnosticsUI("diag"))
  )
)

server <- function(input, output, session) {
  rv <- reactiveValues(
    msi_object = NULL,
    app_log = character(),
    runtime_ready = FALSE
  )

  append_log <- function(msg) {
    # Avoid creating a reactive self-dependency loop while appending logs.
    current_log <- isolate(rv$app_log)
    rv$app_log <- c(
      current_log,
      paste0(format(Sys.time(), "%Y-%m-%d %H:%M:%S"), " | ", msg)
    )
  }

  observeEvent(TRUE, {
    rt <- initialize_msidataset_runtime()
    rv$runtime_ready <- isTRUE(rt$ok)
    append_log(rt$message)

    if (!isTRUE(rt$ok)) {
      showNotification(rt$message, type = "error", duration = NULL)
    }
  }, once = TRUE)

  set_shared_object <- function(obj, message = "Object updated") {
    if (!validate_msi_object(obj)) {
      showNotification("Updated object is not a valid msi_dataset.", type = "error")
      append_log("Rejected invalid object update.")
      return(invisible(FALSE))
    }

    rv$msi_object <- obj
    append_log(message)
    showNotification(message, type = "message")
    TRUE
  }

  dataImportServer(
    id = "import",
    get_object = reactive(rv$msi_object),
    set_object = set_shared_object,
    append_log = append_log
  )

  processingServer(
    id = "processing",
    get_object = reactive(rv$msi_object),
    set_object = set_shared_object,
    append_log = append_log
  )

  visualizationServer(
    id = "viz",
    get_object = reactive(rv$msi_object),
    append_log = append_log
  )

  roiServer(
    id = "roi",
    get_object = reactive(rv$msi_object),
    set_object = set_shared_object,
    append_log = append_log
  )

  diagnosticsServer(
    id = "diag",
    get_object = reactive(rv$msi_object),
    app_log = reactive(rv$app_log)
  )
}

shinyApp(ui = ui, server = server)
