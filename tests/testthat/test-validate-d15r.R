test_that('D15r v2p2 checks evidence citations and metabolic pathway names', {
  schema <- jsonlite::parse_json('{"type":"object"}')
  claim <- function(text, ids) list(text=text,evidence_ids=as.list(ids))
  view <- list(is_control=FALSE, functional_support=list(classes=list(), pathways=list(list(evidence_id='E_path',pathway_type='metabolic'))), deterministic_followup_feature_ids=list())
  answer <- list(module_name_candidate_supported=claim('clue','E_missing'), core_representativeness=claim('core','E_path'), main_classes=list(), possible_pathway_clues=list(claim('clue','E_path')), conflicts_and_unknowns=list(), features_for_followup_identification=list(), supporting_evidence_ids=list(), evidence_level='PATHWAY_CLUE')
  result <- .msi_validate_d15r(answer,view,schema,character())
  expect_setequal(result$errors,c('UNKNOWN_EVIDENCE_ID','PATHWAY_NAME_MISSING_METABOLIC_ID'))
  answer$module_name_candidate_supported$evidence_ids <- list('E_path')
  expect_length(.msi_validate_d15r(answer,view,schema,character())$errors,0)
  answer$core_representativeness$text <- 'activation'
  expect_true('FORBIDDEN_TERM' %in% .msi_validate_d15r(answer,view,schema,character())$errors)
  answer$core_representativeness$text <- 'not activation'
  expect_false('FORBIDDEN_TERM' %in% .msi_validate_d15r(answer,view,schema,character())$errors)
})

test_that('saved D15r responses match v2p2 flags and review text', {
  root <- d15_fixture_root("MSIDATASET_D15_STUDY_ROOT")
  cards <- Sys.getenv("MSIDATASET_D15_VOCABULARY_CARDS")
  d15_require_files(c(cards, file.path(root, c('schemas/llm_output.schema.json',
    paste0('runs/d15r_common/evidence_d15r/views/C', 0:2, '.json'),
    paste0('runs/d15r_primary/llm_d15q/responses/C', 0:2, '.attempt_1.json'),
    paste0('runs/d15r_cross/llm_d15q/responses/C', 0:2, '.attempt_1.json')))))
  parse <- function(p) jsonlite::parse_json(paste(readLines(p,warn=FALSE,encoding='UTF-8'),collapse='\n'))
  schema <- parse(file.path(root,'schemas/llm_output.schema.json'))
  vocab <- .msi_vocabulary(cards)
  for (model in c('primary','cross')) for (i in 0:2) {
    view <- parse(file.path(root,'runs/d15r_common/evidence_d15r/views',paste0('C',i,'.json')))
    saved <- parse(file.path(root,paste0('runs/d15r_',model,'/llm_d15q/responses/C',i,'.attempt_1.json')))
    answer <- jsonlite::parse_json(saved$raw_response$choices[[1]]$message$content)
    actual <- .msi_validate_d15r(answer,view,schema,vocab)
    expect_identical(as.character(actual$errors),as.character(unlist(saved$validation_errors)))
    expect_identical(as.character(actual$review_prose),as.character(unlist(saved$review_prose)))
  }
})
