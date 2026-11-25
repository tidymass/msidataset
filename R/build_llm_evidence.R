.msi_evidence_id <- function(...) {
  parts <- c(...)
  encoded <- paste0("[", paste(vapply(parts,function(x) as.character(jsonlite::toJSON(as.character(x),auto_unbox=TRUE)),character(1)), collapse=", "), "]")
  paste0("E_", substr(digest::digest(charToRaw(enc2utf8(encoded)), algo="sha256", serialize=FALSE),1L,24L))
}

#' Build module evidence from neighbors and pathway mappings
#'
#' Derive representative-weighted HMDB/SMPDB and precomputed Reactome/ChEBI
#' pathway support, specificity adjustments, and the registered 8/2 ranking.
#' By default construct pathway-free views from upstream tables and then add
#' pathway support. Reactome/ChEBI mappings must be precomputed.
#' No reference labels, gold standard, author metadata or truth fields are used.
#' @param base_views Optional named list of pathway-free cluster views. When
#'   omitted, construct views from members, neighbors and metadata below.
#' @param candidates Data frame with feature_id, HMDB_ID, compound_name, and
#'   mfc_priority_score, in the original MFC order.
#' @param representatives Data frame with cluster_id and feature_id.
#' @param hmdb_pathways Path to the registered HMDB/SMPDB membership TSV.
#' @param reactome_mapping Path to the registered precomputed Reactome mapping TSV.
#' @param members Cluster membership data frame.
#' @param neighbors Result of [retrieve_neighbors()].
#' @param features Data frame with feature_id, mz and log2FC.
#' @param metadata Data frame or TSV path with HMDB_ID and taxonomy_path.
#' @param output_dir Writable directory for new full evidence artifacts.
#' @param run_id Run identifier included in the evidence view.
#' @return Named list of extended views and process information.
#'   If an upstream stage reports a stopping status, return it with empty views.
#' @examples
#' \dontrun{views <- build_llm_evidence(base_views, candidates, reps, hmdb, reactome)}
#' @export
build_llm_evidence <- function(base_views = NULL, candidates, representatives, hmdb_pathways, reactome_mapping,
                               members = NULL, neighbors = NULL, features = NULL, metadata = NULL,
                               output_dir = NULL, run_id = "msidataset_run") {
  status <- .msi_stage_status(neighbors)
  if (is.null(status)) status <- .msi_stage_status(members)
  if (is.null(status) && !missing(representatives)) status <- .msi_stage_status(representatives)
  if (!is.null(status)) return(list(status = status, views = list()))
  if (is.null(base_views)) {
    if (is.character(metadata) && length(metadata)==1L) metadata <- data.table::fread(metadata,data.table=FALSE,showProgress=FALSE)
    base_views <- .msi_base_views(members,representatives,candidates,neighbors,features,metadata,output_dir,run_id)
  }
  if (!all(c("feature_id","HMDB_ID","compound_name","mfc_priority_score") %in% names(candidates)) || !all(c("cluster_id","feature_id") %in% names(representatives)) || is.null(names(base_views))) stop("Invalid inputs")
  h <- data.table::fread(hmdb_pathways, data.table=FALSE, showProgress=FALSE, colClasses="character")
  r <- data.table::fread(reactome_mapping, data.table=FALSE, showProgress=FALSE, colClasses="character")
  if (!all(c("HMDB_ID","pathway") %in% names(h)) || !all(c("HMDB_ID","pathway_id","pathway_name","pathway_type","mapping_distance","status") %in% names(r))) stop("Invalid pathway mapping columns")
  h <- unique(h[toupper(h$pathway) != "[MISSING]" & nzchar(h$pathway), c("HMDB_ID","pathway"),drop=FALSE])
  r <- r[r$status == "OK",,drop=FALSE]; r$mapping_distance <- as.integer(r$mapping_distance)
  if (anyNA(r$mapping_distance)) stop("Invalid Reactome mapping distance")
  n_h <- length(unique(h$HMDB_ID)); n_r <- length(unique(r$HMDB_ID))
  hs <- table(h$pathway); rs <- table(unique(r[,c("HMDB_ID","pathway_id")])$pathway_id)
  out <- list()
  for (cid in names(base_views)) {
    view <- base_views[[cid]]
    if (.msi_forbidden_view(view)) stop("Reference or truth fields are prohibited in evidence")
    reps <- as.character(representatives$feature_id[as.character(representatives$cluster_id)==cid]); if (!length(reps) || anyDuplicated(reps)) stop("Invalid cluster representatives")
    cs <- candidates[candidates$feature_id %in% reps,,drop=FALSE]
    if (any(!is.finite(cs$mfc_priority_score) | cs$mfc_priority_score < 0)) stop("Invalid MFC priority")
    den <- vapply(reps,function(fid) Reduce(`+`,cs$mfc_priority_score[cs$feature_id==fid],init=0.0),numeric(1)); if (any(den<=0)) stop("Missing representative candidate weights")
    cs$weight <- cs$mfc_priority_score / den[match(cs$feature_id,reps)] / length(reps)
    # The pathway salts are each a frozen constant for evidence-ID reproducibility.
    make_row <- function(p, source) {
      if (source == "HMDB_SMPDB") {
        ids <- h$HMDB_ID[h$pathway==p]; sub <- cs[cs$HMDB_ID %in% ids,,drop=FALSE]; nm <- p
        typ <- if (grepl('deficiency|disease|syndrome|aciduria|uria$|emia|acidosis|cancer|oncogenic|warburg|toxicity|action pathway|drug|defect|glycogenosis|glycogen storage|intolerance|xanthomatosis|familial|congenital|hereditary|hyper[a-z]*emia|hypo[a-z]*emia|type [0-9]+[a-z]?$|type i+v?$|type [iv]+[ab]?$', p,ignore.case=TRUE,perl=TRUE)) "disease_or_drug" else "metabolic"
        size <- as.integer(hs[[p]]); dist <- NULL; eid <- .msi_evidence_id("d14c_pathway","D15r",cid,p); N <- n_h
      } else {
        x <- r[r$pathway_id==p,,drop=FALSE]; ids <- unique(x$HMDB_ID); sub <- cs[cs$HMDB_ID %in% ids,,drop=FALSE]; nm <- x$pathway_name[1]; typ <- x$pathway_type[1]; size <- as.integer(rs[[p]]); dist <- min(x$mapping_distance[x$HMDB_ID %in% sub$HMDB_ID]); eid <- .msi_evidence_id("d15r_reactome","D15r",cid,p); N <- n_r
      }
      by_id <- sub[!duplicated(sub$HMDB_ID,fromLast=TRUE),,drop=FALSE]; by_id <- by_id[order(by_id$HMDB_ID),,drop=FALSE]
      weight <- Reduce(`+`,unlist(lapply(reps,function(fid) sub$weight[sub$feature_id==fid]),use.names=FALSE),init=0.0); specificity <- log(N/size)/log(N)
      z <- list(evidence_id=eid,name=nm,pathway_type=typ,weighted_support=weight,feature_coverage=length(unique(sub$feature_id))/length(reps),n_candidates=nrow(by_id),example_candidates=as.list(utils::head(as.character(by_id$compound_name),3L)),pathway_size=size,source=source,mapping_distance_min=dist)
      if (source=="Reactome") z$pathway_id <- p
      z$specificity_weight <- specificity; z$adjusted_support <- weight*specificity; z$generic <- size>1000
      z
    }
    hp <- sort(unique(h$pathway[h$HMDB_ID %in% cs$HMDB_ID])); rp <- sort(unique(r$pathway_id[r$HMDB_ID %in% cs$HMDB_ID & r$pathway_type != "other"]))
    rows <- c(lapply(hp,make_row,source="HMDB_SMPDB"),lapply(rp,make_row,source="Reactome")); ranked <- list()
    for (typ in c("metabolic","disease_or_drug")) {
      group <- Filter(function(z) identical(z$pathway_type,typ),rows)
      if (!length(group)) next
      key <- data.frame(adjusted=vapply(group,`[[`,numeric(1),"adjusted_support"),coverage=vapply(group,`[[`,numeric(1),"feature_coverage"),count=vapply(group,`[[`,integer(1),"n_candidates"),size=vapply(group,`[[`,integer(1),"pathway_size"),name=vapply(group,`[[`,character(1),"name"),source=vapply(group,`[[`,character(1),"source"),id=vapply(group,function(z) if(is.null(z$pathway_id)) "" else z$pathway_id,character(1)))
      o <- order(-key$adjusted,-key$coverage,-key$count,key$size,key$name,key$source,key$id,method="radix"); key <- key[o,,drop=FALSE]; group <- group[o]
      for (i in seq_len(min(length(group),if(typ=="metabolic") 8L else 2L))) {
        group[[i]]$rank <- i; group[[i]]$tied_with_next <- i < length(group) && identical(unname(as.numeric(unlist(key[i,c("adjusted","coverage","count")]))),unname(as.numeric(unlist(key[i+1L,c("adjusted","coverage","count")]))))
        fields <- c("evidence_id","name","pathway_type","weighted_support","feature_coverage",
                    "n_candidates","example_candidates","rank","tied_with_next",
                    if (!is.null(group[[i]]$pathway_id)) "pathway_id", "pathway_size","source",
                    "mapping_distance_min","specificity_weight","adjusted_support","generic")
        ranked[[length(ranked)+1L]] <- group[[i]][fields]
      }
    }
    view$functional_support$pathways <- ranked; out[[cid]] <- view
  }
  info <- list(function_name = "build_llm_evidence",
               hmdb_pathways = .msi_file_record(hmdb_pathways),
               reactome_mapping = .msi_file_record(reactome_mapping),
               time = as.character(Sys.time()))
  list(views = out, process_info = info)
}
