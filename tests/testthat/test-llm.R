test_that("local annotation builds requests and validates mocked responses", {
  dir <- tempfile(); dir.create(dir); on.exit(unlink(dir,recursive=TRUE))
  prompt <- file.path(dir,"prompt.txt"); writeLines("system",prompt)
  schema <- file.path(dir,"schema.json")
  writeLines('{"type":"object","required":["answer"],"additionalProperties":false,"properties":{"answer":{"type":"string","minLength":1}}}',schema)
  cfg <- list(base_url="http://localhost",model="local",system_prompt=prompt,output_schema=schema,temperature=0,top_k=1,seed=20261105,max_tokens=1500)
  view <- jsonlite::parse_json('{"cluster_id":"0","representatives":[]}')
  body <- annotate_clusters(view,cfg)$request
  expect_equal(body$model,"local")
  expect_equal(jsonlite::parse_json(body$messages[[2]]$content),view)
  mock <- function(req) httr2::response(status_code=200,headers=list('Content-Type'='application/json'),body=charToRaw('{"choices":[{"finish_reason":"stop","message":{"content":"{\\"answer\\":\\"OK\\"}"}}]}'))
  ans <- httr2::with_mocked_responses(mock, annotate_clusters(view,cfg,perform=TRUE))
  expect_equal(ans$answer$answer,"OK")
  bad <- function(req) httr2::response(status_code=200,headers=list('Content-Type'='application/json'),body=charToRaw('{"choices":[{"finish_reason":"stop","message":{"content":"{\\"answer\\":\\"\\"}"}}]}'))
  expect_error(httr2::with_mocked_responses(bad, annotate_clusters(view,cfg,perform=TRUE)),"too short")
  cfg$base_url <- "http://example.com"
  expect_error(annotate_clusters(view,cfg),"loopback")
  cfg$base_url <- "http://[::1]"
  expect_error(annotate_clusters(jsonlite::parse_json('{"formula":"C1"}'),cfg),"prohibited")
})

test_that("reference user JSON preserves Python float and Unicode conventions", {
  x <- jsonlite::parse_json('{"integral_float":1.0,"integer":1,"small":1e-7,"large":1e16,"text":"\u00e9","object":{},"array":[]}')
  expect_identical(.msi_python_json(x), '{"integral_float":1.0,"integer":1,"small":1e-07,"large":1e+16,"text":"\u00e9","object":{},"array":[]}')
})

test_that("reference recorded request user messages match byte for byte", {
  root <- reference_fixture_root()
  prompt <- file.path(root, "reference_runs/liver/annotation/system_prompt.txt")
  schema <- file.path(root, "reference_runs/liver/annotation/response.schema.json")
  paths <- file.path(root, "reference_runs/liver/evidence/extended", paste0("C", 0:2, ".json"))
  requests <- file.path(root, "reference_runs/liver/annotation/primary/requests", paste0("C", 0:2, ".json"))
  skip_if_not(all(file.exists(c(prompt, schema, paths, requests))), "reference reference data unavailable")
  cfg <- list(base_url = "http://localhost", model = "gpt-oss-120b", system_prompt = prompt,
              output_schema = schema, temperature = 0, top_k = 1, seed = 20261105, max_tokens = 1500)
  for (i in seq_along(paths)) {
    view <- jsonlite::parse_json(paste(readLines(paths[i], warn = FALSE, encoding = "UTF-8"), collapse = "\n"))
    expected <- jsonlite::fromJSON(requests[i], simplifyVector = FALSE)$body$messages[[2]]$content
    expect_identical(annotate_clusters(view, cfg)$request$messages[[2]]$content, expected)
  }
})

test_that("token counting uses chat template then llama tokenizer", {
  d <- tempfile(); dir.create(d); on.exit(unlink(d, recursive = TRUE))
  writeLines('system', file.path(d, 'prompt'))
  writeLines('{"type":"object"}', file.path(d, 'schema'))
  cfg <- list(base_url='http://localhost', model='local',
              system_prompt=file.path(d, 'prompt'), output_schema=file.path(d, 'schema'),
              temperature=0, top_k=1, seed=1, max_tokens=1500)
  calls <- character()
  mock <- function(req) {
    calls <<- c(calls, req$url)
    body <- if (endsWith(req$url, '/apply-template')) '{"prompt":"formatted"}' else '{"tokens":[1,2,3]}'
    httr2::response(status_code=200, headers=list('Content-Type'='application/json'), body=charToRaw(body))
  }
  result <- httr2::with_mocked_responses(mock, annotate_clusters(list(cluster_id='0'), cfg, count_tokens=TRUE))
  expect_identical(result$process_info$input_tokens, 3L)
  expect_match(calls[1], '/apply-template$')
  expect_match(calls[2], '/tokenize$')
})

test_that('token projection preserves full evidence and applies registered reductions', {
  d <- tempfile(); dir.create(d); on.exit(unlink(d,recursive=TRUE))
  writeLines('system',file.path(d,'prompt'))
  writeLines('{"type":"object"}',file.path(d,'schema'))
  full <- file.path(d,'full.json')
  full_candidates <- lapply(1:5,function(i) list(feature_id='F1',name=paste(rep('a',61),collapse='')))
  writeLines(.msi_python_json(list(candidates=full_candidates)),full)
  cfg <- list(base_url='http://localhost',model='local',system_prompt=file.path(d,'prompt'),
    output_schema=file.path(d,'schema'),temperature=0,top_k=1,seed=1,max_tokens=1500)
  pathways <- list(list(rank=1L,pathway_type='metabolic',example_candidates=list('a','b','c')))
  view <- list(cluster_id='0',is_control=FALSE,representatives=list(list(feature_id='F1',candidates=full_candidates)),
    functional_support=list(classes=lapply(1:6,function(i) list(name=as.character(i))),pathways=pathways),
    truncation=list(candidates_omitted=0L,classes_omitted=0L,pathways_omitted=0L,candidate_names_truncated=0L,
      nonrepresentative_candidates_omitted=0L,full_evidence_path=full,
      full_evidence_sha256=digest::digest(file=full,algo='sha256')))
  calls <- 0L; long <- paste0('{"tokens":[',paste(rep('1',20001),collapse=','),']}')
  mock <- function(req) {
    if (endsWith(req$url,'/apply-template')) body <- '{"prompt":"formatted"}' else {
      calls <<- calls+1L; body <- if (calls<8L) long else '{"tokens":[1,2]}'
    }
    httr2::response(status_code=200,headers=list('Content-Type'='application/json'),body=charToRaw(body))
  }
  result <- httr2::with_mocked_responses(mock,annotate_clusters(view,cfg,count_tokens=TRUE))
  expect_equal(result$process_info$input_tokens,2L)
  expect_identical(result$process_info$projection_steps,c('examples_3_to_1','pathways_8_2_to_5_1','minimal_projection','N_path 10->5','N_class 10->5','N_cand 5->3','candidate names truncate to 60 characters'))
  projected <- jsonlite::parse_json(result$request$messages[[2]]$content)
  expect_length(projected$functional_support$classes,5)
  expect_equal(projected$truncation$classes_omitted,1)
  expect_equal(projected$truncation$candidates_omitted,2)
  expect_equal(projected$truncation$candidate_names_truncated,3)
  expect_length(projected$representatives[[1]]$candidates,3)
  expect_identical(view$truncation$full_evidence_sha256,digest::digest(file=full,algo='sha256'))
})

test_that('semantic validation reports invalid and truncated responses without accepting them', {
  d <- tempfile(); dir.create(d); on.exit(unlink(d,recursive=TRUE))
  writeLines('system',file.path(d,'prompt')); writeLines('{"type":"object"}',file.path(d,'schema'))
  writeLines('{"name":"example"}',file.path(d,'vocab'))
  cfg <- list(base_url='http://localhost',model='local',system_prompt=file.path(d,'prompt'),
    output_schema=file.path(d,'schema'),vocabulary_cards=file.path(d,'vocab'),temperature=0,top_k=1,seed=1,max_tokens=1500)
  mock <- function(req) httr2::response(status_code=200,headers=list('Content-Type'='application/json'),
    body=charToRaw('{"choices":[{"finish_reason":"length","message":{"content":"{"}}]}'))
  result <- httr2::with_mocked_responses(mock,annotate_clusters(list(cluster_id='0'),cfg,perform=TRUE))
  expect_identical(result$validation$errors,c('INVALID_RESPONSE','TRUNCATED_LENGTH'))
  expect_false(result$process_info$accepted)
  expect_identical(result$process_info$validation, "SEMANTIC_V2")
})

test_that('JSON null response is a schema failure, not malformed JSON', {
  d <- tempfile(); dir.create(d); on.exit(unlink(d,recursive=TRUE))
  writeLines('system',file.path(d,'prompt')); writeLines('{"type":"object"}',file.path(d,'schema'))
  writeLines('{"name":"example"}',file.path(d,'vocab'))
  cfg <- list(base_url='http://localhost',model='local',system_prompt=file.path(d,'prompt'),
    output_schema=file.path(d,'schema'),vocabulary_cards=file.path(d,'vocab'),temperature=0,top_k=1,seed=1,max_tokens=1500)
  mock <- function(req) httr2::response(status_code=200,headers=list('Content-Type'='application/json'),
    body=charToRaw('{"choices":[{"finish_reason":"stop","message":{"content":"null"}}]}'))
  result <- httr2::with_mocked_responses(mock,annotate_clusters(list(cluster_id='0'),cfg,perform=TRUE))
  expect_identical(result$validation$errors,'SCHEMA_INVALID')
})
