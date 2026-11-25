# Construct pathway-free evidence views using only upstream measurements and library
# metadata. All evidence identifiers are deterministic; full_evidence_path and
# its hash describe the freshly written artifact, not a historical artifact.
.msi_summary_evidence <- function(type, content) {
  ordered <- content[order(names(content), method='radix')]
  canonical <- .msi_python_json(list(type, ordered))
  c(list(evidence_id=paste0('E_',substr(digest::digest(charToRaw(enc2utf8(canonical)),algo='sha256',serialize=FALSE),1L,24L)),type=type),content)
}
.msi_base_views <- function(members, representatives, candidates, neighbors, features, metadata, output_dir, run_id='msidataset_run') {
  required <- list(members=c('cluster_id','feature_id'), representatives=c('cluster_id','feature_id','rank'),
                   candidates=c('feature_id','HMDB_ID','compound_name','best_adduct','mfc_priority_score','rank_within_feature'),
                   features=c('feature_id','mz','log2FC'), metadata=c('HMDB_ID','taxonomy_path'))
  inputs <- list(members=members,representatives=representatives,candidates=candidates,features=features,metadata=metadata)
  for (nm in names(required)) if (!is.data.frame(inputs[[nm]]) || !all(required[[nm]] %in% names(inputs[[nm]]))) stop('Invalid base evidence input: ',nm)
  if (!is.list(neighbors) || !is.data.frame(neighbors$neighbors) || !is.data.frame(neighbors$merged)) stop('Invalid neighbor result')
  if (anyDuplicated(features$feature_id) || anyDuplicated(metadata$HMDB_ID) || anyDuplicated(members$feature_id)) stop('Duplicated feature or metadata ID')
  dir.create(output_dir,recursive=TRUE,showWarnings=FALSE)
  if (!dir.exists(output_dir)) stop('Cannot create evidence output directory')
  output_dir <- normalizePath(output_dir)
  out <- list()
  for (cid in sort(unique(as.character(members$cluster_id)))) {
    rr <- representatives[as.character(representatives$cluster_id)==cid,,drop=FALSE]; rr <- rr[order(rr$rank),,drop=FALSE]
    if (!nrow(rr) || anyDuplicated(rr$feature_id) || anyNA(match(rr$feature_id,features$feature_id))) stop('Invalid representative mapping')
    fids <- as.character(rr$feature_id); k <- length(fids); nmembers <- sum(as.character(members$cluster_id)==cid)
    rows <- which(as.character(candidates$feature_id) %in% fids)
    cs <- lapply(rows,function(i) {
      z <- candidates[i,,drop=FALSE]; path <- metadata$taxonomy_path[match(z$HMDB_ID,metadata$HMDB_ID)]
      if (!length(path) || is.na(path)) path <- ''
      list(evidence_id=.msi_evidence_id('candidate',cid,z$feature_id,z$HMDB_ID,as.character(i-1L)),feature_id=as.character(z$feature_id),entity_id=as.character(z$HMDB_ID),name=as.character(z$compound_name),adduct=as.character(z$best_adduct),score=as.numeric(z$mfc_priority_score),rank=as.numeric(z$rank_within_feature),classes=if(nzchar(path)) list(as.character(path)) else list(),pathways=list())
    })
    seq_sum <- function(x) Reduce(`+`,x,init=0.0)
    den <- vapply(fids,function(fid) seq_sum(vapply(Filter(function(x) identical(x$feature_id,fid),cs),function(x)x$score,numeric(1))),numeric(1))
    if (any(!is.finite(den) | den <= 0)) stop('Missing or invalid representative candidate weight')
    labels <- sort(unique(unlist(lapply(cs,function(x)x$classes),use.names=FALSE)),method='radix')
    classes <- lapply(labels,function(label) {
      hits <- Filter(function(x) label %in% unlist(x$classes),cs)
      list(name=label,feature_count=length(unique(vapply(hits,`[[`,character(1),'feature_id'))),
           feature_coverage=length(unique(vapply(hits,`[[`,character(1),'feature_id')))/k,
           weighted_support=seq_sum(unlist(lapply(fids,function(fid) vapply(Filter(function(x) identical(x$feature_id,fid),hits),function(x)x$score/den[match(fid,fids)],numeric(1))),use.names=FALSE))/k,
           denominator_features=as.integer(k),evidence_id=.msi_evidence_id('mfc_class',cid,label))
    })
    unknown <- seq_sum(unlist(lapply(fids,function(fid) vapply(Filter(function(x) identical(x$feature_id,fid) && !length(x$classes),cs),function(x)x$score/den[match(fid,fids)],numeric(1))),use.names=FALSE))/k
    ordered <- if(length(classes)) order(-vapply(classes,`[[`,numeric(1),'weighted_support'),vapply(classes,`[[`,character(1),'name'),method='radix') else integer()
    class_view <- lapply(utils::head(classes[ordered],10L),function(x)x[c('name','feature_coverage','weighted_support','denominator_features','evidence_id')])
    rep_views <- lapply(fids,function(fid) {
      f <- features[match(fid,features$feature_id),,drop=FALSE]
      if (!is.finite(f$mz) || !is.finite(f$log2FC)) stop('Unresolved representative measurements')
      cands <- Filter(function(x) identical(x$feature_id,fid),cs)
      o <- order(vapply(cands,`[[`,numeric(1),'rank'),-vapply(cands,`[[`,numeric(1),'score'),vapply(cands,`[[`,character(1),'entity_id'),method='radix')
      selected <- lapply(utils::head(cands[o],5L),function(x) c(x[c('entity_id','name','adduct','rank','score','evidence_id')],list(Sub_class=x$classes)))
      content <- list(feature_id=fid,mz=as.numeric(f$mz),direction=if(f$log2FC>0) 'up' else 'down',embedding_coverage=1.0)
      c(content,list(candidates=selected,evidence_id=.msi_summary_evidence('representative',content)$evidence_id))
    })
    subset <- neighbors$neighbors[as.character(neighbors$neighbors$cluster_id)==cid,,drop=FALSE]
    terms <- trimws(unlist(strsplit(as.character(subset$taxonomy_path),' > ',fixed=TRUE),use.names=FALSE)); terms <- terms[!is.na(terms) & nzchar(terms)]
    counts <- table(terms); ord <- order(-as.integer(counts),names(counts),method='radix'); ord <- utils::head(ord,10L)
    semantics <- lapply(ord,function(i) list(evidence_id=.msi_evidence_id('semantic_class',cid,names(counts)[i]),label=names(counts)[i],count=as.integer(counts[i]),neighbor_count=as.integer(nrow(subset))))
    merged <- neighbors$merged[as.character(neighbors$merged$cluster_id)==cid,,drop=FALSE]
    nn <- .msi_evidence_neighbors(merged)
    n_omitted <- sum(vapply(fids,function(fid) max(0L,sum(vapply(cs,function(x) identical(x$feature_id,fid),logical(1)))-5L),integer(1)))
    trunc <- list(candidates_omitted=as.integer(n_omitted),classes_omitted=as.integer(max(0,length(classes)-10L)),pathways_omitted=0L,candidate_names_truncated=0L,nonrepresentative_candidates_omitted=0L)
    root <- list(evidence_id=.msi_evidence_id('cluster',cid),n_members=as.integer(nmembers),medoid=fids[1],k=as.integer(k))
    selection <- 'nearest actual member to original-embedding centroid, then nearest to that member'
    meta <- list(run_id=run_id,cluster_id=cid,module_source='PCA_KMEANS',n_members=as.integer(nmembers),medoid=fids[1],k=as.integer(k),selection_rule=selection)
    uncertainty <- as.list(c('Candidate support is not identity or probability.','Semantic neighbors are not identities.','Single-section spatial blocks are not biological replicates.'))
    other <- .msi_summary_evidence('classes_other',list(weighted_support=NULL,feature_coverage=NULL,denominator_features=as.integer(k),status='NOT_PRECOMPUTED_OVERLAPPING_LABELS'))
    unknown_summary <- .msi_summary_evidence('classes_unknown',list(weighted_support=unknown,feature_coverage=NULL,denominator_features=as.integer(k),status='SOURCE_SUMMARY'))
    support <- list(classes=class_view,classes_other=other,classes_unknown=unknown_summary,pathways=list(),pathways_other=NULL,pathways_unknown=NULL)
    # Preserve the historical cluster metadata ID for the renamed default run.
    hash_meta <- meta
    if (identical(run_id, 'msidataset_run')) hash_meta$run_id <- 'd15_resume'
    cluster_metadata <- .msi_summary_evidence('cluster_metadata', hash_meta)
    cluster_metadata$run_id <- run_id
    pack <- c(meta,list(member_feature_ids=as.list(as.character(members$feature_id[as.character(members$cluster_id)==cid])),representatives=rep_views,candidates=cs,functional_support=support,semantic_neighbors=nn,cluster_semantics=list(classes=semantics),uncertainty=uncertainty,cluster_evidence=list(root)))
    fp <- file.path(output_dir,paste0('C',cid,'.json')); writeLines(.msi_python_json(pack),fp,useBytes=TRUE)
    trunc$full_evidence_path <- fp; trunc$full_evidence_sha256 <- digest::digest(file=fp,algo='sha256')
    trunc$evidence_id <- .msi_summary_evidence('truncation',trunc[c('candidates_omitted','classes_omitted','pathways_omitted','candidate_names_truncated','nonrepresentative_candidates_omitted')])$evidence_id
    out[[cid]] <- c(meta,list(cluster_evidence=list(root),is_control=FALSE,control_variant=NULL,deterministic_followup_feature_ids=list(),representatives=rep_views,functional_support=support,independent_evidence=structure(list(),names=character()),uncertainty=uncertainty,truncation=trunc,
      core_representativeness=.msi_summary_evidence('core_representativeness',list(k=as.integer(k),n_members=as.integer(nmembers),status='CORE_VS_WHOLE_COVERAGE_NOT_PRECOMPUTED')),
      missing_resources=.msi_summary_evidence('missing_resources',list(no_card_count=NULL,no_embedding_count=NULL,status='COUNTS_NOT_AVAILABLE_IN_FULL_EVIDENCE')),
      uncertainty_evidence=lapply(uncertainty,function(s).msi_summary_evidence('uncertainty',list(text=s))),cluster_metadata=cluster_metadata,semantic_neighbors=nn,cluster_semantics=list(classes=semantics)))
  }
  out
}
