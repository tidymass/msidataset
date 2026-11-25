# msidataset

msidataset analyzes regional differences in mass spectrometry imaging (MSI). It
builds pseudo-MS/MS spectra from spatially correlated ion images, matches them
against a spectral library, and retains ranked metabolite-feature candidates
(MFCs) with their supporting evidence. It then uses precomputed semantic vectors
to group features into modules and retrieve related library entries. Evidence
cards support optional interpretation by a local language model. The workflow
can produce module-level clues while individual metabolite identities remain
uncertain.

## Installation

Install Cardinal from Bioconductor and massdataset from GitHub before installing msidataset.

```r
install.packages(c("BiocManager", "remotes"))
BiocManager::install("Cardinal")
remotes::install_github("tidymass/massdataset")
remotes::install_github("tidymass/msidataset")
```

## Minimal example

This example creates a small synthetic pixel CSV.
The columns contain x, y, total ion current (TIC), and two feature intensities. Both analyses record their results in the `mass_dataset` object's
`process_info`.

```r
library(msidataset)

pixels <- expand.grid(x = 1:10, y = 1:6)
signal <- seq_len(nrow(pixels))
input <- tempfile(fileext = ".csv")
write.table(cbind(pixels, tic = 100, precursor = signal * 2, fragment = signal),
            input, sep = ",", row.names = FALSE, col.names = FALSE)

msi <- read_msi_data(input, format = "csv", mz = c(100, 50))
labels <- transform(pixels, region = rep(c("periportal", "pericentral"), each = 30))
msi <- define_roi(msi, labels)
msi <- detect_de_features(msi)
msi <- generate_pseudo_ms2(msi)

msi@process_info$detect_de_features$result
msi@process_info$generate_pseudo_ms2$result
unlink(input)
```

See `vignette("msidataset")` for an example with multiple paired blocks.

## Exported functions

| Function | Purpose |
| --- | --- |
| `read_msi_data()` | Import pixel CSV or imzML into a `mass_dataset`. |
| `define_roi()` | Attach region labels to pixels. |
| `detect_de_features()` | Compare regions using paired spatial blocks. |
| `generate_pseudo_ms2()` | Reconstruct pseudo-MS/MS fragment pairs. |
| `match_spectral_library()` | Retrieve reference spectra and score fragment matches. |
| `calculate_mfc()` | Score and retain metabolite-feature candidates and exclusions. |
| `read_vector_library()` | Read rows from an external, preconverted float32 vector library. |
| `embed_features()` | Combine candidate vectors weighted by MFC priority. |
| `select_k()` | Select a cluster count using PCA, k-means, and a permutation baseline. |
| `cluster_features()` | Assign features to modules and estimate their stability. |
| `select_representatives()` | Select observed members to represent each module. |
| `retrieve_neighbors()` | Retrieve semantic neighbors and merge them by module. |
| `build_llm_evidence()` | Build evidence views with optional functional support. |
| `annotate_clusters()` | Prepare local model requests and validate model responses. |
| `run_pipeline()` | Run analysis stages through evidence construction with checkpoints. |

## External resources and local annotation

Set `MSIDATASET_REFERENCE_ROOT` to the archive containing `spectral_libraries/`,
`ms1_library/`, `vector_library/`, and `reference_runs/liver/`. Spectrum matching
and MFC scoring verify the SHA-256 hashes in `inst/config/reference_lock.yaml`.
For embedding and retrieval, pass the `vector_library/` directory containing
`meta.json`, `ids.tsv`, and `embeddings.f32`.

`match_spectral_library()` and `calculate_mfc()` also accept a named `resources`
list or a YAML file. Use `matching_config` and `ms1_library` for both functions,
`ms2_library` for spectral matching, and `vector_ids` for MFC scoring. Entries
can be paths or lists with `path` and optional `sha256`. YAML paths are relative
to the YAML file. The process record includes the resource source and file hashes.
These functions currently support only negative ion mode ([M−H]⁻, [M+Cl]⁻).
When declaring an ion mode, use `ion_mode: negative`.

`run_pipeline()` ends at evidence construction by default. Call
`annotate_clusters()` separately with an evidence view.

Local annotation is optional. Set `config$base_url` to the loopback HTTP
address and port of your llama.cpp-compatible server, and supply the model,
prompt, and response schema. `perform = FALSE` and the default
`count_tokens = FALSE` prepare a request without contacting the server.
`count_tokens = TRUE` opts in to local template and tokenizer requests.
LLM output is for review and never changes MFC ranking. Reference labels,
gold standards, author names and formulas, and evaluation truth stay outside
model input.

## Citation and license

Xiaoyu Xu and Xiaotao Shen. msidataset 0.2.1. Manuscript in preparation.

MIT. See [LICENSE](LICENSE).
