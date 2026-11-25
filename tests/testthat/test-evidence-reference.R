test_that('Reconstructed base and pathway evidence matches the reference', {
  skip_if_not(identical(Sys.getenv('MSIDATASET_RUN_REFERENCE_DOWNSTREAM'),'true'),'Set MSIDATASET_RUN_REFERENCE_DOWNSTREAM=true')
  root <- reference_fixture_root()
  tsv <- function(p) data.table::fread(file.path(root,p),data.table=FALSE,showProgress=FALSE)
  reference_require_files(file.path(root, c(
    'reference_runs/liver/clustering/members.tsv',
    'reference_runs/liver/clustering/representatives.tsv',
    'reference_runs/liver/mfc/candidates.tsv',
    'reference_runs/liver/neighbors/retrieved.tsv', 'vector_library/ids.tsv',
    'reference_runs/liver/de/features.tsv', 'reference_runs/liver/mfc/features.tsv',
    'reference_runs/liver/pathways/hmdb_smpdb.tsv', 'reference_runs/liver/pathways/hmdb_reactome.tsv',
    paste0('reference_runs/liver/evidence/extended/C', 0:2, '.json'),
    paste0('reference_runs/liver/evidence/base/C', 0:2, '.json'),
    paste0('reference_runs/liver/annotation/primary/requests/C', 0:2, '.json'))))
  members <- tsv('reference_runs/liver/clustering/members.tsv')
  reps <- tsv('reference_runs/liver/clustering/representatives.tsv')
  candidates <- tsv('reference_runs/liver/mfc/candidates.tsv')
  near <- tsv('reference_runs/liver/neighbors/retrieved.tsv')
  metadata <- tsv('vector_library/ids.tsv')
  near$card_completeness <- metadata$card_completeness[match(near$HMDB_ID,metadata$HMDB_ID)]
  merged <- do.call(rbind,lapply(split(near,near$cluster_id),function(part) {
    rows <- lapply(split(part,part$HMDB_ID),function(g) {
      g <- g[order(-g$cosine,g$feature_id,g$rank),,drop=FALSE]
      b <- g[1,,drop=FALSE]; b$n_reps_hit <- length(unique(g$feature_id)); b$max_cosine <- b$cosine; b
    })
    z <- do.call(rbind,rows); head(z[order(-z$max_cosine,z$HMDB_ID),,drop=FALSE],25L)
  }))
  features <- merge(tsv('reference_runs/liver/de/features.tsv'),tsv('reference_runs/liver/mfc/features.tsv')[,c('feature_id','mz')],by='feature_id')
  d <- tempfile(); dir.create(d); on.exit(unlink(d,recursive=TRUE))
  views <- .msi_base_views(members,reps,candidates,list(neighbors=near,merged=merged),features,metadata,d)
  extended <- build_llm_evidence(candidates=candidates,representatives=reps,
    hmdb_pathways=file.path(root,'reference_runs/liver/pathways/hmdb_smpdb.tsv'),reactome_mapping=file.path(root,'reference_runs/liver/pathways/hmdb_reactome.tsv'),
    members=members,neighbors=list(neighbors=near,merged=merged),features=features,metadata=metadata,output_dir=d)$views
  for (i in 0:2) for (extended_view in c(FALSE,TRUE)) {
    path <- if(extended_view) 'reference_runs/liver/evidence/extended' else 'reference_runs/liver/evidence/base'
    expected <- jsonlite::parse_json(paste(readLines(file.path(root,path,paste0('C',i,'.json')),warn=FALSE),collapse='\n'))
    actual <- if(extended_view) extended[[as.character(i)]] else views[[as.character(i)]]
    # Actual hashes and paths must refer to fresh files, so compare scientific
    # content separately and check the new provenance independently.
    expect_identical(actual$truncation$full_evidence_sha256,digest::digest(file=actual$truncation$full_evidence_path,algo='sha256'))
    actual$truncation$full_evidence_path <- expected$truncation$full_evidence_path
    actual$truncation$full_evidence_sha256 <- expected$truncation$full_evidence_sha256
    expect_identical(actual$run_id, "msidataset_run")
    expect_identical(actual$cluster_metadata$run_id, "msidataset_run")
    actual$run_id <- expected$run_id
    actual$cluster_metadata$run_id <- expected$cluster_metadata$run_id
    expect_equal(actual,expected,tolerance=1e-6)
    if (!extended_view) expect_identical(.msi_python_json(actual), .msi_python_json(expected))
    if (extended_view) {
      request <- jsonlite::parse_json(paste(readLines(file.path(root,'reference_runs/liver/annotation/primary/requests',paste0('C',i,'.json')),warn=FALSE),collapse='\n'))
      expect_identical(.msi_python_json(actual),request$body$messages[[2]]$content)
    }
  }
})
