test_that('base evidence is derived from upstream rows and has local provenance', {
  d <- tempfile(); dir.create(d); on.exit(unlink(d,recursive=TRUE))
  members <- data.frame(cluster_id=0,feature_id=c('F1','F2'))
  reps <- data.frame(cluster_id=0,feature_id='F1',rank=1)
  candidates <- data.frame(feature_id='F1',HMDB_ID='H1',compound_name='A',best_adduct='[M-H]-',mfc_priority_score=0.5,rank_within_feature=1)
  features <- data.frame(feature_id='F1',mz=100,log2FC=-1)
  metadata <- data.frame(HMDB_ID='H1',taxonomy_path='Organic > Lipids')
  near <- data.frame(cluster_id=0,feature_id='F1',rank=1,HMDB_ID='H1',name='A',taxonomy_path='Organic > Lipids',card_completeness='FULL',cosine=.9,n_reps_hit=1,max_cosine=.9)
  views <- .msi_base_views(members,reps,candidates,list(neighbors=near,merged=near),features,metadata,d)
  view <- views[['0']]
  expect_equal(view$n_members,2)
  expect_equal(view$representatives[[1]]$candidates[[1]]$entity_id,'H1')
  expect_false('description' %in% names(view$semantic_neighbors[[1]]))
  expect_equal(view$semantic_neighbors[[1]]$n_reps_hit,1)
  expect_identical(view$truncation$full_evidence_sha256,digest::digest(file=view$truncation$full_evidence_path,algo='sha256'))
  expect_false(.msi_forbidden_view(view))
})

test_that('representative count follows the reference locked rule max(ceil(0.1n), 5)', {
  for (n in c(5L, 23L, 50L, 60L)) {
    vectors <- diag(n); rownames(vectors) <- sprintf('F%03d',seq_len(n))
    members <- data.frame(feature_id=rownames(vectors),cluster_id=0L)
    expect_equal(nrow(select_representatives(vectors,members)), min(n, max(5L, as.integer(ceiling(0.1*n)))))
  }
})
