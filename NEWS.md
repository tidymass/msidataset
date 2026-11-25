# msidataset 0.2.0

This release rewrites the package since v0.1.1 for staged MSI analysis.

- Import pixel CSV and imzML into `mass_dataset` objects with recorded processing history.
- Fix imzML m/z precision by preserving measured centroids for mass matching
  and using the feature axis when no centroid is available.
- Compare spatial regions, reconstruct pseudo-MS/MS spectra, match reference
  libraries, and retain ranked metabolite-feature candidates with audit details.
- Use preconverted vector libraries for weighted feature embeddings, permutation-based
  cluster selection, stability assessment, representative selection, and neighbor retrieval.
- Construct evidence views and optionally request and validate annotations from a
  local llama.cpp-compatible server.
- Resume input- and code-hashed checkpoints and reproduce the locked validation
  workflow with external reference resources.
