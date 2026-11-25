test_that('reconstructed D15 and D15r evidence agrees apart from local provenance', {
  skip_if_not(identical(Sys.getenv('MSIDATASET_RUN_D15_DOWNSTREAM'),'true'),'Set MSIDATASET_RUN_D15_DOWNSTREAM=true')
  root <- d15_fixture_root("MSIDATASET_D15_STUDY_ROOT")
  tsv <- function(p) data.table::fread(file.path(root,p),data.table=FALSE,showProgress=FALSE)
  d15_require_files(file.path(root, c(
    'runs/d15_resume/cluster_d15/cluster_members.tsv',
    'runs/d15_resume/cluster_d15/representatives.tsv',
    'runs/d15_resume/mfc_final_d15/mfc_candidates_filtered.tsv',
    'runs/d15_resume/retrieve_d15/neighbors.tsv', 'libs/hmdb_full_bge_v2/ids.tsv',
    'runs/d15_roi2/de_d15/de_features.tsv', 'runs/d15_resume/mfc_final_d15/features.tsv',
    'runs/d14c_common/pathways.tsv', 'runs/d15r_common/hmdb_reactome.tsv',
    paste0('runs/d15r_common/evidence_d15r/views/C', 0:2, '.json'),
    paste0('runs/d15_neighbors_v2/neighbors_v2_d15/views/C', 0:2, '.json'),
    paste0('runs/d15r_primary/llm_d15q/requests/C', 0:2, '.json'))))
  members <- tsv('runs/d15_resume/cluster_d15/cluster_members.tsv')
  reps <- tsv('runs/d15_resume/cluster_d15/representatives.tsv')
  candidates <- tsv('runs/d15_resume/mfc_final_d15/mfc_candidates_filtered.tsv')
  near <- tsv('runs/d15_resume/retrieve_d15/neighbors.tsv')
  metadata <- tsv('libs/hmdb_full_bge_v2/ids.tsv')
  near$card_completeness <- metadata$card_completeness[match(near$HMDB_ID,metadata$HMDB_ID)]
  merged <- do.call(rbind,lapply(split(near,near$cluster_id),function(part) {
    rows <- lapply(split(part,part$HMDB_ID),function(g) {
      g <- g[order(-g$cosine,g$feature_id,g$rank),,drop=FALSE]
      b <- g[1,,drop=FALSE]; b$n_reps_hit <- length(unique(g$feature_id)); b$max_cosine <- b$cosine; b
    })
    z <- do.call(rbind,rows); head(z[order(-z$max_cosine,z$HMDB_ID),,drop=FALSE],25L)
  }))
  features <- merge(tsv('runs/d15_roi2/de_d15/de_features.tsv'),tsv('runs/d15_resume/mfc_final_d15/features.tsv')[,c('feature_id','mz')],by='feature_id')
  d <- tempfile(); dir.create(d); on.exit(unlink(d,recursive=TRUE))
  views <- .msi_base_views(members,reps,candidates,list(neighbors=near,merged=merged),features,metadata,d)
  extended <- build_llm_evidence(candidates=candidates,representatives=reps,
    hmdb_pathways=file.path(root,'runs/d14c_common/pathways.tsv'),reactome_mapping=file.path(root,'runs/d15r_common/hmdb_reactome.tsv'),
    members=members,neighbors=list(neighbors=near,merged=merged),features=features,metadata=metadata,output_dir=d)$views
  for (i in 0:2) for (extended_view in c(FALSE,TRUE)) {
    path <- if(extended_view) 'runs/d15r_common/evidence_d15r/views' else 'runs/d15_neighbors_v2/neighbors_v2_d15/views'
    expected <- jsonlite::parse_json(paste(readLines(file.path(root,path,paste0('C',i,'.json')),warn=FALSE),collapse='\n'))
    actual <- if(extended_view) extended[[as.character(i)]] else views[[as.character(i)]]
    # Actual hashes and paths must refer to fresh files, so compare scientific
    # content separately and check the new provenance independently.
    expect_identical(actual$truncation$full_evidence_sha256,digest::digest(file=actual$truncation$full_evidence_path,algo='sha256'))
    actual$truncation$full_evidence_path <- expected$truncation$full_evidence_path
    actual$truncation$full_evidence_sha256 <- expected$truncation$full_evidence_sha256
    expect_equal(actual,expected,tolerance=1e-6)
    if (extended_view) {
      request <- jsonlite::parse_json(paste(readLines(file.path(root,'runs/d15r_primary/llm_d15q/requests',paste0('C',i,'.json')),warn=FALSE),collapse='\n'))
      expect_identical(.msi_python_json(actual),request$body$messages[[2]]$content)
    }
  }
})
