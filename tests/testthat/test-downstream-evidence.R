test_that("reference pathway extension uses representative weights", {
  p <- tempfile(); dir.create(p); on.exit(unlink(p,recursive=TRUE))
  h <- data.frame(HMDB_ID=c("A","B"),pathway=c("Path A","Path B"))
  r <- data.frame(HMDB_ID=c("A","B"),pathway_id=c("R1","R2"),pathway_name=c("Reaction A","Reaction B"),pathway_type=c("metabolic","disease_or_drug"),mapping_distance=c(0,1),status="OK")
  write.table(h,file.path(p,"h.tsv"),sep="\t",row.names=FALSE,quote=FALSE)
  write.table(r,file.path(p,"r.tsv"),sep="\t",row.names=FALSE,quote=FALSE)
  cand <- data.frame(feature_id=c("f1","f1"),HMDB_ID=c("A","B"),compound_name=c("a","b"),mfc_priority_score=c(1,3))
  reps <- data.frame(cluster_id=0,feature_id="f1")
  view <- jsonlite::parse_json('{"cluster_id":"0","functional_support":{"pathways":[]}}')
  result <- build_llm_evidence(list('0'=view),cand,reps,file.path(p,"h.tsv"),file.path(p,"r.tsv"))$views[['0']]$functional_support$pathways
  expect_length(result,4)
  a <- Filter(function(x) identical(x$name,"Path A"),result)[[1]]
  expect_equal(a$weighted_support,.25)
  expect_equal(a$feature_coverage,1)
  expect_equal(a$source,"HMDB_SMPDB")
  expect_error(build_llm_evidence(list('0'=jsonlite::parse_json('{"formula":"X"}')),cand,reps,file.path(p,"h.tsv"),file.path(p,"r.tsv")),"prohibited")
})
