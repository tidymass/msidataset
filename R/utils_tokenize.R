# Use the server's chat template and tokenizer to measure model input. The caller
# controls whether this makes a request; no model is started by the package.
.msi_token_count <- function(base_url, messages) {
  post <- function(path, body) {
    req <- httr2::request(paste0(sub('/$', '', base_url), '/', path))
    req <- httr2::req_options(req, proxy = '')
    req <- httr2::req_timeout(req, 90)
    req <- httr2::req_body_json(req, body, auto_unbox = TRUE)
    httr2::resp_body_json(httr2::req_perform(req), simplifyVector = FALSE)
  }
  prompt <- post('apply-template', list(messages = messages))$prompt
  if (!is.character(prompt) || length(prompt) != 1L || is.na(prompt))
    stop('Invalid apply-template response')
  tokens <- post('tokenize', list(content = prompt, add_special = TRUE, parse_special = TRUE))$tokens
  if (!is.list(tokens) && !is.atomic(tokens)) stop('Invalid tokenize response')
  if (is.null(tokens)) stop('Missing tokenize response')
  length(tokens)
}

# Apply the registered evidence reduction sequence. The full evidence artifact
# remains untouched and must match the digest carried by the input view.
.msi_minimal_projection <- function(view, body, base_url) {
  path <- view$truncation$full_evidence_path; sha <- view$truncation$full_evidence_sha256
  if (!is.character(path) || length(path)!=1L || !file.exists(path) ||
      !is.character(sha) || !identical(digest::digest(file=path,algo='sha256'),sha))
    stop('Missing or mismatched full evidence artifact for token projection')
  pack <- jsonlite::parse_json(paste(readLines(path,warn=FALSE,encoding='UTF-8'),collapse='\n'))
  if (.msi_forbidden_view(pack) || !is.list(pack$candidates)) stop('Invalid full evidence artifact')
  steps <- 'minimal_projection'
  count <- function() {
    body$messages[[2]]$content <<- .msi_python_json(view)
    .msi_token_count(base_url,body$messages)
  }
  refresh <- function() {
    keys <- c('candidates_omitted','classes_omitted','pathways_omitted','candidate_names_truncated','nonrepresentative_candidates_omitted')
    view$truncation$evidence_id <<- .msi_summary_evidence('truncation',view$truncation[keys])$evidence_id
  }
  n <- count()
  if (n > 20000L) { steps <- c(steps,'N_path 10->5'); n <- count() }
  if (n > 20000L) {
    rows <- view$functional_support$classes
    view$functional_support$classes <- utils::head(rows,5L)
    view$truncation$classes_omitted <- as.integer(view$truncation$classes_omitted + max(0L,length(rows)-5L))
    refresh(); steps <- c(steps,'N_class 10->5'); n <- count()
  }
  if (n > 20000L) {
    for (i in seq_along(view$representatives)) {
      r <- view$representatives[[i]]; r$candidates <- utils::head(r$candidates,3L); view$representatives[[i]] <- r
    }
    view$truncation$candidates_omitted <- as.integer(sum(vapply(view$representatives,function(r) {
      max(0L,sum(vapply(pack$candidates,function(x) identical(x$feature_id,r$feature_id),logical(1)))-3L)
    },integer(1))))
    refresh(); steps <- c(steps,'N_cand 5->3'); n <- count()
  }
  if (n > 20000L) {
    shortened <- 0L
    for (i in seq_along(view$representatives)) {
      r <- view$representatives[[i]]
      for (j in seq_along(r$candidates)) if (nchar(r$candidates[[j]]$name,type='chars') > 60L) {
        r$candidates[[j]]$name <- substr(r$candidates[[j]]$name,1L,60L); shortened <- shortened + 1L
      }
      view$representatives[[i]] <- r
    }
    view$truncation$candidate_names_truncated <- as.integer(shortened)
    refresh(); steps <- c(steps,'candidate names truncate to 60 characters'); n <- count()
  }
  list(view=view,body=body,tokens=n,steps=steps)
}
