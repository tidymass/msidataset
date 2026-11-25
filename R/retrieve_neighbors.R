#' Retrieve and merge semantic neighbors
#'
#' Search the external row-major float32 library by block matrix products.
#' Within each cluster, merge HMDB IDs by maximum cosine; ties in the chosen
#' hit use feature ID and rank, and ties across IDs use HMDB ID. Library vectors
#' are assumed to be L2 normalized as in the locked reference library.
#' @param vectors Unit matrix or result of [embed_features()].
#' @param representatives Data frame with cluster_id and feature_id.
#' @param library_path Directory understood by [read_vector_library()].
#' @param block_rows Positive integer rows per block.
#' @return List of `neighbors`, cluster `merged` neighbors, and process history.
#' @examples
#' \dontrun{near <- retrieve_neighbors(embedding, representatives, library_path)}
#' @export
retrieve_neighbors <- function(vectors, representatives, library_path, block_rows = 4096L) {
  x <- .msi_unit(vectors); ids <- if (is.list(vectors)) vectors$feature_ids else rownames(x)
  ix <- match(representatives$feature_id, ids)
  if (is.null(ids) || anyNA(ix) || anyDuplicated(ids) || !all(c("cluster_id", "feature_id") %in% names(representatives)) || length(block_rows)!=1L || is.na(block_rows) || block_rows < 1 || block_rows != floor(block_rows)) stop("Invalid representatives or block size")
  q <- x[ix,,drop=FALSE]; q <- q / sqrt(rowSums(q^2))
  meta <- jsonlite::fromJSON(file.path(library_path,"meta.json")); d <- as.integer(meta$shape[2]); nr <- as.integer(meta$shape[1])
  if (ncol(q)!=d || meta$dtype != "float32" || meta$byte_order != "little" || meta$order != "C") stop("Invalid library shape or format")
  id <- utils::read.delim(file.path(library_path,"ids.tsv"),check.names=FALSE,stringsAsFactors=FALSE,quote="",comment.char="")
  bin <- file.path(library_path,"embeddings.f32"); if(nrow(id)!=nr || anyDuplicated(id$HMDB_ID) || file.info(bin)$size != as.double(nr)*d*4) stop("Invalid library IDs or file size")
  con <- file(bin,"rb"); on.exit(close(con)); best <- lapply(seq_len(nrow(q)),function(i) data.frame(row=integer(),cosine=numeric()))
  start <- 1L
  while (start <= nr) {
    n <- min(block_rows,nr-start+1L); v <- readBin(con,"numeric",n=n*d,size=4,endian="little")
    if(length(v)!=n*d || any(!is.finite(v))) stop("Truncated or nonfinite library vector block")
    m <- matrix(v,nrow=n,byrow=TRUE); sim <- q %*% t(m); rows <- seq.int(start,length.out=n)
    for (j in seq_len(nrow(q))) {
      a <- rbind(best[[j]],data.frame(row=rows,cosine=as.numeric(sim[j,])))
      best[[j]] <- utils::head(a[order(-a$cosine,id$HMDB_ID[a$row]),,drop=FALSE],25L)
    }
    start <- start+n
  }
  near <- lapply(seq_along(best),function(j) {
    b <- best[[j]]; data.frame(cluster_id=representatives$cluster_id[j],feature_id=representatives$feature_id[j],rank=seq_len(nrow(b)),HMDB_ID=id$HMDB_ID[b$row],name=id$name[b$row],taxonomy_path=id$taxonomy_path[b$row],card_completeness=id$card_completeness[b$row],cosine=b$cosine,interpretation="semantic neighbor, not identity",stringsAsFactors=FALSE)
  }); near <- do.call(rbind,near)
  merged <- lapply(split(near,near$cluster_id),function(part) {
    groups <- split(part,part$HMDB_ID); rows <- lapply(groups,function(g) {
      g <- g[order(-g$cosine,g$feature_id,g$rank),,drop=FALSE]; best <- g[1,,drop=FALSE]; best$n_reps_hit <- length(unique(g$feature_id)); best$max_cosine <- best$cosine; best
    }); z <- do.call(rbind,rows); utils::head(z[order(-z$max_cosine,z$HMDB_ID),,drop=FALSE],25L)
  })
  list(neighbors=near,merged=do.call(rbind,merged),process_info=list(function_name="retrieve_neighbors", library_path=library_path,block_rows=block_rows,top_n=25L,time=as.character(Sys.time())))
}
