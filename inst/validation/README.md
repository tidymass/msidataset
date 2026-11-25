# Saved validation records

The summary and discrepancy tables retain the original Chinese descriptions and
recorded values. Source paths are relative to archive directories and paired
with SHA-256 hashes. The data and reports referenced in the tables are external.

Set `MSIDATASET_VALIDATION_ROOT` to the validation archive directory.
The script expects `study`, `package`, and `data` subdirectories. Override
their locations with `MSIDATASET_VALIDATION_STUDY_ROOT`,
`MSIDATASET_VALIDATION_PACKAGE_ROOT`, and `MSIDATASET_VALIDATION_DATA_ROOT`.
Report records default to a `reports` subdirectory under each of the study and
package roots. Override these locations with
`MSIDATASET_VALIDATION_STUDY_REPORTS_ROOT` and
`MSIDATASET_VALIDATION_PACKAGE_REPORTS_ROOT`.

From the source tree, run `Rscript inst/validation/recompute_validation.R --check`
to compare the tables with archived inputs. The package root should contain the
source copy of `inst/validation`. Use `--write` when intentionally regenerating
those tables. The script computes summary values from saved outputs.
