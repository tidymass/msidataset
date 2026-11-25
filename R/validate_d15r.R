# Semantic validator for saved model output. It reports all prose for review; it does not assert
# that prose is scientifically correct. Vocabulary is a fixed input resource.
.msi_norm_entity <- function(x) gsub('\\s+', ' ', trimws(tolower(x)), perl=TRUE)
.msi_vocabulary <- function(path) {
  con <- file(path, open='r', encoding='UTF-8'); on.exit(close(con))
  lines <- readLines(con, warn=FALSE)
  terms <- unlist(lapply(lines, function(line) {
    card <- jsonlite::parse_json(line)
    out <- card$name
    if (!is.null(card$text)) {
      rows <- strsplit(card$text, '\n', fixed=TRUE)[[1]]
      rows <- rows[grepl(': ', rows, fixed=TRUE)]
      keys <- sub(': .*$', '', rows); vals <- sub('^[^:]*: ', '', rows)
      value <- function(key) { x <- vals[keys == key]; if (length(x)) tail(x,1L) else '' }
      tax <- value('Taxonomy')
      out <- c(out, value('Name'), strsplit(tax,' > ',fixed=TRUE)[[1]], tax,
               strsplit(value('Pathways'),'; ',fixed=TRUE)[[1]])
    }
    for (tax in card$taxonomy) out <- c(out, strsplit(as.character(tax),' > ',fixed=TRUE)[[1]], tax)
    for (p in card$pathways) if (is.list(p) && !is.null(p$name)) out <- c(out,p$name)
    out
  }), use.names=FALSE)
  terms <- unique(.msi_norm_entity(terms))
  terms[nchar(terms, type='chars') >= 4L]
}
.msi_all_evidence_ids <- function(x) {
  if (!is.list(x)) return(character())
  result <- if (!is.null(names(x)) && is.character(x$evidence_id)) x$evidence_id else character()
  c(result, unlist(lapply(x, .msi_all_evidence_ids), use.names=FALSE))
}
.msi_input_numbers <- function(x, percentages=FALSE) {
  excluded <- c('evidence_id','feature_id','entity_id','cluster_id','run_id','schema_version','source','source_path','card_source','full_evidence_path','full_evidence_sha256')
  if (is.list(x)) {
    if (!is.null(names(x))) x <- x[setdiff(names(x),excluded)]
    return(unlist(lapply(x,.msi_input_numbers,percentages=percentages), use.names=FALSE))
  }
  if (is.logical(x) || is.null(x)) return(numeric())
  if (is.numeric(x)) return(if (percentages) numeric() else x)
  if (!is.character(x)) return(numeric())
  pattern <- if (percentages) '(?<![\\w.])[-+]?\\d+(?:\\.\\d+)?(?:[eE][-+]?\\d+)?(?![\\w.])\\s*%' else '(?<![\\w.])[-+]?\\d+(?:\\.\\d+)?(?:[eE][-+]?\\d+)?(?![\\w.])'
  tokens <- unlist(regmatches(x, gregexpr(pattern,x,perl=TRUE)),use.names=FALSE)
  values <- suppressWarnings(as.numeric(sub('\\s*%$', '', tokens)))
  if (percentages) values / 100 else values
}
.msi_regex_quote <- function(x) gsub('([][{}()+*^$|\\\\?.])', '\\\\\\1', x, perl=TRUE)
.msi_validate_d15r <- function(answer, view, schema, vocabulary) {
  valid <- tryCatch({.msi_validate_schema(answer,schema); TRUE}, error=function(e) FALSE)
  if (!valid) return(list(errors='SCHEMA_INVALID', review_prose=character()))
  claims <- c(list(answer$module_name_candidate_supported, answer$core_representativeness), answer$main_classes,
              answer$possible_pathway_clues, answer$conflicts_and_unknowns, answer$features_for_followup_identification)
  errors <- character()
  cited <- c(unlist(answer$supporting_evidence_ids), unlist(lapply(claims, function(x) x$evidence_ids)))
  if (any(vapply(claims,function(x) !length(x$evidence_ids),logical(1)))) errors <- c(errors,'UNCITED_CLAIM')
  if (length(setdiff(cited,.msi_all_evidence_ids(view)))) errors <- c(errors,'UNKNOWN_EVIDENCE_ID')
  nums <- .msi_input_numbers(view); pct <- .msi_input_numbers(view, TRUE)
  pattern <- '(?<![\\w.])[-+]?\\d+(?:\\.\\d+)?(?:[eE][-+]?\\d+)?(?![\\w.])'
  for (claim in claims) {
    txt <- gsub('(?m)^\\s*\\d+[.)\u3001]\\s*', '', claim$text, perl=TRUE)
    matches <- gregexpr(pattern,txt,perl=TRUE)[[1]]
    if (matches[1L] < 0L) next
    lengths <- attr(matches,'match.length')
    for (i in seq_along(matches)) {
      token <- substr(txt,matches[i],matches[i]+lengths[i]-1L); val <- as.numeric(token)
      mant <- strsplit(tolower(token),'e',fixed=TRUE)[[1]][1L]
      dec <- if (grepl('.',mant,fixed=TRUE)) nchar(sub('^[^.]*\\.','',mant)) else 0L
      candidates <- c(nums,0,1)
      after <- substring(txt,matches[i]+lengths[i])
      candidates <- c(candidates, if (grepl('^\\s*%',after,perl=TRUE)) nums*100 else pct)
      if (!any(abs(round(candidates,dec)-val) <= 1e-10,na.rm=TRUE)) errors <- c(errors,'UNSOURCED_NUMBER')
    }
  }
  text <- .msi_norm_entity(paste(vapply(claims, `[[`, character(1),'text'), collapse=' '))
  source <- .msi_norm_entity(.msi_python_json(view))
  # Fixed-string prefilter avoids compiling regular expressions for the full vocabulary.
  for (term in vocabulary) if (grepl(term,text,fixed=TRUE)) {
    pat <- paste0('(?<![A-Za-z0-9_])',.msi_regex_quote(term),'(?![A-Za-z0-9_])')
    if (grepl(pat,text,perl=TRUE) && !grepl(pat,source,perl=TRUE)) { errors <- c(errors,'OUT_OF_INPUT_ENTITY'); break }
  }
  forbidden <- c('identity','\u8eab\u4efd\u786e\u5b9a','activation','\u6fc0\u6d3b','flux','\u901a\u91cf')
  markers <- '\u672a|\u4e0d|\u65e0|\u975e|\u6ca1\u6709|\u5e76\u975e|\u4e0d\u80fd|\u4e0d\u5f97|\u4e0d\u4ee3\u8868|\u4e0d\u7b49\u4e8e|(?<![a-z])(?:not|no|without|cannot|does not)(?![a-z])'
  for (claim in claims) for (term in forbidden) {
    txt <- tolower(claim$text); pos <- gregexpr(term,txt,fixed=TRUE)[[1]]
    if (pos[1L] < 0L) next
    for (p in pos) {
      prefix <- if (p > 1L) substr(txt,1L,p-1L) else ''
      prefix <- sub('^.*[\u3002\uff1b;.!?]','',prefix,perl=TRUE)
      prefix <- substr(prefix,max(1L,nchar(prefix)-7L),nchar(prefix))
      if (!grepl(markers,prefix,perl=TRUE)) errors <- c(errors,'FORBIDDEN_TERM')
    }
  }
  support <- view$functional_support; paths <- support$pathways; classes <- support$classes
  if (answer$evidence_level=='PATHWAY_CLUE' && (!length(paths) || !length(answer$possible_pathway_clues)) ||
      answer$evidence_level=='CLASS_ONLY' && (!length(classes) || !length(answer$main_classes) || length(answer$possible_pathway_clues)) ||
      answer$evidence_level=='INSUFFICIENT' && (length(answer$main_classes) || length(answer$possible_pathway_clues))) errors <- c(errors,'EVIDENCE_LEVEL_CONFLICT')
  permitted <- unlist(view$deterministic_followup_feature_ids)
  if (any(!vapply(answer$features_for_followup_identification,function(x) x$text %in% permitted,logical(1)))) errors <- c(errors,'FOLLOWUP_NOT_PRESELECTED')
  if (isTRUE(view$is_control) && answer$evidence_level!='INSUFFICIENT') errors <- c(errors,'CONTROL_ABSTENTION')
  if (answer$evidence_level=='PATHWAY_CLUE') {
    ids <- unlist(lapply(Filter(function(x) identical(x$pathway_type,'metabolic'),paths),function(x)x$evidence_id))
    if (!length(intersect(ids,unlist(answer$module_name_candidate_supported$evidence_ids)))) errors <- c(errors,'PATHWAY_NAME_MISSING_METABOLIC_ID')
  }
  list(errors=sort(unique(errors)), review_prose=vapply(claims,`[[`,character(1),'text'))
}
