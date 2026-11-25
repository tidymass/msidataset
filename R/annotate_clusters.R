# Validate the subset of JSON Schema used by the registered model output schema.
# Unknown keywords and external references fail closed.
.msi_validate_schema <- function(value, schema, root = schema, path = "$", depth = 0L) {
  if (depth > 100L) stop("Schema nesting is too deep")
  allowed <- c("$schema", "$defs", "$ref", "type", "required", "properties", "additionalProperties", "items", "minItems", "minLength", "enum")
  unknown <- setdiff(names(schema), allowed)
  if (length(unknown)) stop("Unsupported schema keywords: ", paste(unknown, collapse=", "))
  if (!is.null(schema$`$ref`)) {
    ref <- schema$`$ref`
    if (!is.character(ref) || length(ref)!=1L || !grepl("^#/[A-Za-z0-9_$./-]+$", ref)) stop("Only local schema references are supported")
    target <- root
    for (part in strsplit(substring(ref,3L),"/",fixed=TRUE)[[1]]) { target <- target[[part]]; if (is.null(target)) stop("Unresolved schema reference") }
    return(.msi_validate_schema(value,target,root,path,depth+1L))
  }
  typ <- schema$type
  if (!is.null(typ)) {
    ok <- switch(typ, object=is.list(value) && !is.null(names(value)), array=is.list(value) && is.null(names(value)), string=is.character(value) && length(value)==1L && !is.na(value), number=is.numeric(value) && length(value)==1L && is.finite(value), integer=is.numeric(value) && length(value)==1L && is.finite(value) && value==floor(value), boolean=is.logical(value) && length(value)==1L && !is.na(value), null=is.null(value), stop("Unsupported schema type"))
    if (!ok) stop("Schema type mismatch at ",path)
  }
  if (!is.null(schema$enum) && !any(vapply(schema$enum, identical, logical(1), value))) stop("Schema enum mismatch at ",path)
  if (!is.null(schema$minLength) && nchar(value,type="chars") < schema$minLength) stop("String too short at ",path)
  if (!is.null(schema$minItems) && length(value) < schema$minItems) stop("Array too short at ",path)
  if (identical(typ,"object")) {
    required <- unlist(schema$required, use.names=FALSE); if (length(setdiff(required,names(value)))) stop("Missing required property at ",path)
    props <- schema$properties
    if (identical(schema$additionalProperties,FALSE) && length(setdiff(names(value),names(props)))) stop("Unexpected property at ",path)
    if (!is.null(schema$additionalProperties) && !identical(schema$additionalProperties,FALSE)) stop("Unsupported additionalProperties")
    for (nm in intersect(names(value),names(props))) .msi_validate_schema(value[[nm]],props[[nm]],root,paste0(path,".",nm),depth+1L)
  }
  if (identical(typ,"array") && !is.null(schema$items)) for (i in seq_along(value)) .msi_validate_schema(value[[i]],schema$items,root,paste0(path,"[",i,"]"),depth+1L)
  invisible(TRUE)
}
.msi_check_schema <- function(schema, depth = 0L) {
  if (!is.list(schema) || is.null(names(schema)) || depth > 100L) stop("Unsupported schema object")
  allowed <- c("$schema", "$defs", "$ref", "type", "required", "properties", "additionalProperties", "items", "minItems", "minLength", "enum")
  unknown <- setdiff(names(schema), allowed)
  if (length(unknown)) stop("Unsupported schema keywords: ", paste(unknown, collapse=", "))
  if (!is.null(schema$`$ref`) && length(setdiff(names(schema), "$ref"))) stop("Unsupported reference siblings")
  for (group in c("$defs", "properties")) if (!is.null(schema[[group]])) for (child in schema[[group]]) .msi_check_schema(child,depth+1L)
  if (!is.null(schema$items)) .msi_check_schema(schema$items,depth+1L)
  invisible(TRUE)
}
.msi_forbidden_view <- function(x) {
  forbidden <- c("formula", "author_name", "author_formula", "reference_label", "reference_labels", "gold_standard", "pathway_truth", "truth_mapped", "majority_pathway", "expected_higher_in", "basis_quote", "d14_pathway_truth", "hit_vocabulary")
  if (!is.list(x)) return(FALSE)
  if (!is.null(names(x)) && any(tolower(names(x)) %in% forbidden)) return(TRUE)
  any(vapply(x,.msi_forbidden_view,logical(1)))
}

#' Construct or send a local cluster annotation request
#'
#' Construct the OpenAI-compatible request using a configured prompt and output
#' schema. Only loopback HTTP endpoints are accepted. Response validation uses
#' the JSON Schema subset used by the registered output schema. With `vocabulary_cards` in config,
#' responses also receive semantic evidence checks. All prose remains for review.
#' Optional token counting uses the local server's actual tokenizer.
#' @param view An evidence view parsed with `jsonlite::parse_json()`.
#' @param config List with `base_url`, `model`, `system_prompt`, `output_schema`,
#'   `temperature`, `top_k`, `seed`, and `max_tokens`. Set `base_url` to the
#'   loopback HTTP address and port of your llama.cpp-compatible server.
#' @param response Optional parsed model answer to validate offline without sending
#'   a request. Requires `config$vocabulary_cards`.
#' @param count_tokens If `TRUE`, call the local server's chat template and
#'   tokenizer endpoints to measure the input. This never starts a server.
#' @param perform If `FALSE` (default), prepare the request body. With
#'   `count_tokens = FALSE` this makes no HTTP requests. If `TRUE`, send a
#'   completion request to the configured local server.
#' @return The request body, and when sent, the schema-validated response and
#'   process information. No candidate ranking is modified.
#' @examples
#' \dontrun{request <- annotate_clusters(view, config, perform = FALSE)}
#' @export
annotate_clusters <- function(view, config, perform = FALSE, count_tokens = FALSE, response = NULL) {
  fields <- c("base_url","model","system_prompt","output_schema","temperature","top_k","seed","max_tokens")
  if (!is.list(config) || !all(fields %in% names(config)) || any(vapply(config[fields],is.null,logical(1)))) stop("Incomplete annotation configuration")
  url <- config$base_url
  if (!is.character(url) || length(url)!=1L || !grepl("^http://(localhost|127\\.0\\.0\\.1|\\[::1\\])(:[0-9]+)?/?$",url,ignore.case=TRUE)) stop("Only a loopback HTTP base URL is permitted")
  if (.msi_forbidden_view(view)) stop("Reference or truth fields are prohibited in evidence")
  if (!is.list(view) || is.null(names(view))) stop("Evidence view must be a JSON object")
  schema <- jsonlite::parse_json(paste(readLines(config$output_schema,warn=FALSE,encoding="UTF-8"),collapse="\n"))
  .msi_check_schema(schema)
  prompt <- readChar(config$system_prompt,file.info(config$system_prompt)$size,useBytes=TRUE)
  body <- list(model=config$model,temperature=config$temperature,max_tokens=config$max_tokens,
    response_format=list(type="json_schema",json_schema=list(name="module_annotation",strict=TRUE,schema=schema)),
    messages=list(list(role="system",content=prompt),list(role="user",content=.msi_python_json(view))),top_k=config$top_k,seed=config$seed)
  info <- list(function_name="annotate_clusters",model=config$model,base_url=url,schema=.msi_file_record(config$output_schema),performed=isTRUE(perform),validation="SCHEMA_ONLY",accepted=FALSE,time=as.character(Sys.time()))
  if (isTRUE(count_tokens)) {
    count <- function() .msi_token_count(url, body$messages)
    update <- function() { body$messages[[2]]$content <<- .msi_python_json(view) }
    info$projection_steps <- character()
    info$input_tokens <- count()
    if (info$input_tokens > 20000L && length(view$functional_support$pathways)) {
      for (i in seq_along(view$functional_support$pathways)) {
        p <- view$functional_support$pathways[[i]]
        p$example_candidates <- utils::head(p$example_candidates, 1L)
        view$functional_support$pathways[[i]] <- p
      }
      info$projection_steps <- c(info$projection_steps, "examples_3_to_1")
      update(); info$input_tokens <- count()
    }
    if (info$input_tokens > 20000L && length(view$functional_support$pathways)) {
      view$functional_support$pathways <- Filter(function(p) {
        !is.null(p$rank) && !is.null(p$pathway_type) &&
          p$rank <= if (identical(p$pathway_type, "metabolic")) 5L else 1L
      }, view$functional_support$pathways)
      info$projection_steps <- c(info$projection_steps, "pathways_8_2_to_5_1")
      update(); info$input_tokens <- count()
    }
    if (info$input_tokens > 20000L && !isTRUE(view$is_control)) {
      projected <- .msi_minimal_projection(view, body, url)
      view <- projected$view; body <- projected$body; info$input_tokens <- projected$tokens
      info$projection_steps <- c(info$projection_steps, projected$steps)
    }
    info$token_projection_status <- if (info$input_tokens <= 20000L) "OK" else "OVER_20000_TOKEN_CAP_AFTER_MINIMAL_PROJECTION"
    if (info$input_tokens > 20000L && isTRUE(perform)) stop(info$token_projection_status)
  }
  if (!is.null(response)) {
    if (isTRUE(perform)) stop("Cannot provide response and perform a request")
    if (is.null(config$vocabulary_cards)) stop("vocabulary_cards is required for semantic validation")
    validation <- .msi_validate_annotation(response, view, schema, .msi_vocabulary(config$vocabulary_cards))
    info$validation <- "SEMANTIC_V2"
    info$semantic_valid <- !length(validation$errors)
    info$accepted <- info$semantic_valid && !is.null(info$input_tokens) && info$input_tokens <= 20000L
    return(list(request=body,answer=response,validation=validation,process_info=info))
  }
  if (!isTRUE(perform)) return(list(request=body,process_info=info))
  request <- httr2::request(paste0(sub("/$","",url),"/v1/chat/completions"))
  request <- httr2::req_options(request,proxy="")
  request <- httr2::req_body_json(request,body,auto_unbox=TRUE,null="null",digits=NA)
  response <- httr2::resp_body_json(httr2::req_perform(request),simplifyVector=FALSE)
  choice <- response$choices[[1]]
  if (identical(choice$finish_reason,"length") && is.null(config$vocabulary_cards)) stop("Truncated model response")
  content <- choice$message$content
  if (is.null(config$vocabulary_cards)) {
    if (!is.character(content) || length(content)!=1L) stop("Missing response content")
    answer <- jsonlite::parse_json(content)
    .msi_validate_schema(answer,schema)
    return(list(request=body,answer=answer,validation=NULL,process_info=info))
  }
  parse_failed <- !is.character(content) || length(content)!=1L || is.na(content)
  answer <- if (!parse_failed) tryCatch(jsonlite::parse_json(content),error=function(e) { parse_failed <<- TRUE; NULL }) else NULL
  validation <- if (parse_failed) list(errors="INVALID_RESPONSE",review_prose=character()) else
    .msi_validate_annotation(answer, view, schema, .msi_vocabulary(config$vocabulary_cards))
  if (identical(choice$finish_reason,"length")) validation$errors <- c(validation$errors,"TRUNCATED_LENGTH")
  info$validation <- "SEMANTIC_V2"
  info$semantic_valid <- !length(validation$errors)
  info$accepted <- info$semantic_valid && !is.null(info$input_tokens) && info$input_tokens <= 20000L
  list(request=body,answer=answer,validation=validation,process_info=info)
}
