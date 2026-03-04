# msidataset Shiny App

## Folder structure

- `shiny-app/app.R`: app entrypoint
- `shiny-app/modules/utils.R`: shared helper functions
- `shiny-app/modules/data_import_module.R`: load demo/upload data + preview
- `shiny-app/modules/processing_module.R`: activation + basic filters
- `shiny-app/modules/visualization_module.R`: spatial map + spectra
- `shiny-app/modules/roi_module.R`: ROI placeholder and extension hook
- `shiny-app/modules/diagnostics_module.R`: summary and process history

## Run locally

From the package root:

```r
shiny::runApp("shiny-app")
```

## Notes

- Core computation remains in `msidataset` functions.
- App module code handles UI and reactive orchestration only.
- ROI tab is intentionally future-ready, not hard-coded to one ROI workflow.
