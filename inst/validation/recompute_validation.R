# Read saved analysis outputs. This script does not rerun discovery or inference.
validation_roots <- function() {
  base <- Sys.getenv("MSIDATASET_VALIDATION_ROOT")
  if (!nzchar(base) || !dir.exists(base)) {
    stop("Set MSIDATASET_VALIDATION_ROOT to an existing validation archive directory.",
         call. = FALSE)
  }
  study <- Sys.getenv("MSIDATASET_VALIDATION_STUDY_ROOT", file.path(base, "study"))
  package <- Sys.getenv("MSIDATASET_VALIDATION_PACKAGE_ROOT", file.path(base, "package"))
  c(WEEK = study, PKG = package,
    DATA = Sys.getenv("MSIDATASET_VALIDATION_DATA_ROOT", file.path(base, "data")),
    WEEK_REPORTS = Sys.getenv("MSIDATASET_VALIDATION_STUDY_REPORTS_ROOT", file.path(study, "reports")),
    PKG_REPORTS = Sys.getenv("MSIDATASET_VALIDATION_PACKAGE_REPORTS_ROOT", file.path(package, "reports")))
}

validation_paths <- function(sources, roots) {
  sources <- sub("^(WEEK|PKG)/reports/", "\\1_REPORTS/", sources)
  root <- sub("/.*", "", sources)
  stopifnot(all(root %in% names(roots)))
  file.path(roots[root], sub("^[^/]*/", "", sources))
}

validation_tsv <- function(path) {
  if (file.info(path)$size == 0L) return(data.frame())
  data.table::fread(path, sep = "\t", header = TRUE, na.strings = NULL,
                    encoding = "UTF-8", data.table = FALSE)
}

validation_read_summary <- function(path) {
  read.delim(path, colClasses = "character", na.strings = NULL, quote = "",
             comment.char = "", check.names = FALSE, fileEncoding = "UTF-8")
}

validation_optional_sources <- function() {
  c("WEEK/runs/d15r_common/hmdb_reactome.tsv",
    "WEEK/runs/d12_liver/region_d12/liver_mouse_7_ROI_2/labels.tsv",
    "DATA/pk2c/pipeline-release.rds")
}

validation_text <- function(value) {
  if (!length(value)) return("")
  if (is.numeric(value)) {
    return(paste(formatC(value, digits = 12, format = "g", flag = "#") |>
                   trimws() |> sub("(\\.[0-9]*?)0+(e|$)", "\\1\\2", x = _, perl = TRUE) |>
                   sub("\\.(e|$)", "\\1", x = _), collapse = " / "))
  }
  paste(value, collapse = " / ")
}

validation_write <- function(data, path) {
  stopifnot(!any(grepl("[\t\r\n]", unlist(data))))
  write.table(data, path, sep = "\t", row.names = FALSE, quote = FALSE,
              na = "", fileEncoding = "UTF-8")
}

validation_ari <- function(counts) {
  pairs <- function(x) sum(choose(x, 2))
  rows <- pairs(rowSums(counts))
  cols <- pairs(colSums(counts))
  expected <- rows * cols / pairs(sum(counts))
  denominator <- (rows + cols) / 2 - expected
  if (denominator == 0) return(1)
  (pairs(counts) - expected) / denominator
}

recompute_validation <- function(roots = validation_roots(),
                                 omit = character()) {
  cache <- new.env(parent = emptyenv())
  hashes <- new.env(parent = emptyenv())
  read <- function(source) {
    if (!exists(source, cache, inherits = FALSE)) {
      path <- validation_paths(source, roots)
      value <- switch(tools::file_ext(path),
                      tsv = validation_tsv(path),
                      json = jsonlite::fromJSON(path, simplifyVector = FALSE),
                      rds = readRDS(path),
                      stop("Unsupported source: ", source))
      assign(source, value, cache)
    }
    cache[[source]]
  }
  tsv <- function(path) read(paste0("WEEK/", path))
  j <- function(path) read(paste0("WEEK/", path))
  yes <- function(x) toupper(as.character(x)) == "TRUE"
  num <- function(x) as.numeric(x)
  frac <- function(n, d) paste0(n, "/", d)
  size <- function(x) as.integer(table(x))
  distinct <- function(x) length(unique(x[!is.na(x) & nzchar(x)]))
  field <- function(obj, path) Reduce(function(x, key) x[[key]], path, init = obj)
  compare <- function(x, reference, limit = 1e-10) {
    a <- suppressWarnings(as.numeric(x))
    b <- suppressWarnings(as.numeric(reference))
    if (is.finite(a) && is.finite(b)) return(abs(a - b) <= limit)
    identical(x, reference)
  }
  rows <- list()
  references <- list()
  add <- function(id, stage, purpose, metric, source, compute,
                  criterion = "无验收阈值", status = "仅描述",
                  ref = "", historical = "", caveat = "", method = "重算",
                  reference_source = "") {
    if (id %in% omit) return(invisible(NULL))
    stopifnot(!id %in% names(rows), length(source) > 0L)
    source <- unique(c(source, reference_source[nzchar(reference_source)]))
    paths <- validation_paths(source, roots)
    stopifnot(all(file.exists(paths)))
    for (i in seq_along(source)) {
      if (!exists(source[i], hashes, inherits = FALSE)) {
        hashes[[source[i]]] <- digest::digest(file = paths[i], algo = "sha256")
      }
    }
    value <- validation_text(compute)
    stopifnot(nzchar(value), !anyNA(compute))
    result <- if (is.function(status)) status(compute) else status
    rows[[id]] <<- data.frame(
      validation_id = id, section = section, dataset = dataset, stage = stage,
      purpose = purpose, metric = metric, value = value,
      reference_value = validation_text(ref), pass_criterion = criterion,
      status = result, independent = independent,
      source_file = paste(source, collapse = "; "),
      source_sha256 = paste(vapply(source, function(x) hashes[[x]], ""), collapse = "; "),
      commit = "未核实", recompute = method, caveat = caveat,
      stringsAsFactors = FALSE
    )
    if (nzchar(historical)) {
      references[[id]] <<- data.frame(validation_id = id,
        historical_value = validation_text(ref), historical_source = historical,
        recomputed_value = value, source_file = rows[[id]]$source_file,
        interpretation = caveat)
    }
    invisible(NULL)
  }
  eq <- function(expected) function(x) if (isTRUE(all.equal(x, expected))) "PASS" else "FAIL"
  confirm <- function(x) if (isTRUE(x)) "PASS" else "FAIL"
  wk <- function(x) paste0("WEEK/", x)
  pk <- function(x) paste0("PKG/", x)

  section <- "外部验证"
  dataset <- "Samarah 肝脏 鼠7 ROI2（D15）"
  independent <- "外部数据来源；鼠7已用于D12/D13，ROI与参数沿用开发阶段"
  cluster <- "runs/d15_resume/cluster_d15/"
  evaluate <- "runs/d15_resume/evaluate_d15/"
  d15r <- "runs/d15r_evaluation/evaluate_d15r/"
  de <- "runs/d15_roi2/de_d15/de_features.tsv"
  all_de <- "runs/d15_roi2/de_d15/features_all.tsv"
  locked <- "WEEK/deliverables/LOCKED_D15/key_numbers.tsv"
  old15 <- "WEEK/deliverables/D15/REPORT_D15_ZH.md"
  old15r <- "WEEK/deliverables/D15r/REPORT_D15r_ZH.md"
  d15 <- function(id, stage, metric, source, value, ref = "",
                  caveat = "", method = "重算", status = "仅描述",
                  criterion = "无验收阈值", historical = locked) {
    add(id, stage, stage, metric, wk(source), value, criterion, status, ref,
        if (length(ref) && nzchar(validation_text(ref))) historical else "",
        caveat, method)
  }
  d15("de_total", "差异分析", "入选 feature", de, nrow(tsv(de)), 58,
      "单切片方块是配对筛选单元")
  d15("de_pc", "差异分析", "PC", de, sum(tsv(de)$direction == "pericentral"), 38)
  d15("d15_roi_ion_feature", "差异分析", "图示 PC 最大效应 feature / log2FC", de, {
    x <- subset(tsv(de), direction == "pericentral")
    i <- which.max(x$log2FC)
    c(x$feature_id[i], validation_text(x$log2FC[i]))
  }, caveat = "只用于选图，无阈值调整")
  d15("de_pp", "差异分析", "PP", de, sum(tsv(de)$direction == "periportal"), 20)
  tiles <- "runs/d15_roi2/de_d15/tiles.tsv"
  d15("d15_tiles", "空间划分", "合格方块", tiles, nrow(tsv(tiles)), 393, historical = old15)
  roi <- "runs/d12_liver/region_d12/liver_mouse_7_ROI_2/labels.tsv"
  d15("d15_roi_regions", "空间划分", "PP / PC / 血管像素", roi,
      as.integer(table(tsv(roi)$region)[c("periportal", "pericentral", "VESSEL")]),
      c(78425, 85740, 8641), historical = old15)
  overlap <- paste0(evaluate, "de_overlap.tsv")
  independent <- "否；D13八鼠分析已包含鼠7，参考来自同一外部研究"
  d15("de_overlap", "方向复现", "D15 与 D13 差异交集", overlap, nrow(tsv(overlap)), 47)
  d15("de_overlap_same_direction", "方向复现", "交集方向一致", overlap,
      with(tsv(overlap), sum(direction_d15 == direction_d13)), 47)
  d15("d15_overlap_rate", "方向复现", "交集方向一致 / 交集",
      overlap, with(tsv(overlap), frac(sum(direction_d15 == direction_d13), nrow(tsv(overlap)))),
      "47/47", historical = old15)
  d12de <- "runs/d12_liver/de_d12/de_features.tsv"
  d15("d15_overlap_recount", "方向复现", "从两套 DE 结果重做交集和方向比对", c(de, d12de), {
    x <- merge(tsv(de)[c("feature_id", "direction")], tsv(d12de)[c("feature_id", "direction")],
               by = "feature_id")
    frac(sum(x$direction.x == x$direction.y), nrow(x))
  }, "47/47", historical = old15,
  caveat = "D13 沿用 D12 的八鼠差异筛选，包含鼠7")
  independent <- "外部数据来源；同切片产生离子和伪谱，方法沿用D12/D13"
  pseudo <- "runs/d15_resume/pseudo_d15/pseudo_ms2.tsv"
  d15("pseudo_pairs", "伪二级", "碎片关系", pseudo, nrow(tsv(pseudo)), 144)
  d15("pseudo_precursors", "伪二级", "有碎片前体", pseudo, distinct(tsv(pseudo)$precursor_feature_id), 28)
  reverse <- "runs/d15_resume/coverage_d15/reverse_match_all.tsv"
  qualified <- function() subset(tsv(reverse), yes(precursor_gate_pass) & n_matched_peaks >= 2)
  d15("pseudo_spectrum_matches", "谱库匹配", "合格谱级匹配", reverse, nrow(qualified()), 9)
  d15("pseudo_spectrum_features", "谱库匹配", "有合格匹配的 feature", reverse, distinct(qualified()$feature_id), 7)
  library <- "libs/ms2_neg_d15/conversion_counts.json"
  d15("d15_lipidblast_conversion", "谱库匹配", "负离子源记录 / 展开库谱", library,
      with(j(library)$counts, c(negative_source_records, converted_spectra)),
      c(792757, 794681), method = "读转换输出计数", historical = old15)
  d15("d15_spectrum_unmapped", "谱库匹配", "合格匹配中无 HMDB 映射", reverse,
      sum(is.na(qualified()$target_hmdb_id) | qualified()$target_hmdb_id == ""), 2, historical = old15)
  mfc <- "runs/d15_resume/mfc_final_d15/mfc_candidates_filtered.tsv"
  d15("mfc_candidates", "MFC", "过滤后候选行", mfc, nrow(tsv(mfc)), 3391)
  d15("mfc_features_with_candidates", "MFC", "有候选 feature", mfc, distinct(tsv(mfc)$feature_id), 149)
  d15("mfc_features_total", "MFC", "输入 feature", all_de, nrow(tsv(all_de)), 156)
  d15("d15_mfc_coverage", "MFC", "有候选 / 全部 feature", c(mfc, all_de),
      frac(distinct(tsv(mfc)$feature_id), nrow(tsv(all_de))), "149/156", historical = old15)
  d15("mfc_tier1_pseudo_ms1", "MFC", "同时有伪二级与 MS1 的候选行", mfc,
      with(tsv(mfc), sum(yes(has_pseudo_ms2) & yes(has_ms1))), 7)
  coverage <- "runs/d15_resume/coverage_d15/candidate_coverage.tsv"
  d15("d15_candidates_with_spectra", "谱库匹配", "有谱 MFC 候选行 / MFC 候选行", c(mfc, coverage), {
    joined <- merge(tsv(mfc)[c("feature_id", "HMDB_ID")], tsv(coverage),
                    by = c("feature_id", "HMDB_ID"), all.x = TRUE)
    frac(sum(yes(joined$has_library_spectrum), na.rm = TRUE), nrow(tsv(mfc)))
  }, "435/3391", historical = old15, caveat = "谱库中有记录与伪二级匹配是两种口径")
  independent <- "外部作者分子式仅用于事后评价；同一研究的预选峰，鼠7曾参与流程开发"
  formula <- paste0(evaluate, "formula_evaluation.tsv")
  ranks <- function() num(tsv(formula)$correct_formula_rank)
  d15("formula_evaluable", "分子式评价", "可评估", formula, nrow(tsv(formula)), 95)
  d15("formula_candidate_present", "分子式评价", "正确式出现", formula, sum(yes(tsv(formula)$correct_formula_present)), 89)
  d15("formula_top1", "分子式评价", "正确式首位（含并列）", formula, sum(ranks() <= 1, na.rm = TRUE), 88)
  d15("formula_top5", "分子式评价", "正确式前五（含并列）", formula, sum(ranks() <= 5, na.rm = TRUE), 89)
  d15("d15_formula_rank_bins", "分子式评价", "首位 / 第2–5位 / 未出现", formula,
      c(sum(ranks() == 1, na.rm = TRUE), sum(ranks() > 1 & ranks() <= 5, na.rm = TRUE),
        sum(!yes(tsv(formula)$correct_formula_present))), caveat = "同分排名计入首位")
  formula_comparison <- paste0(evaluate, "mfc_reference_comparison.tsv")
  for (version in c("d13_mouse1", "d13_filtered_mouse8")) {
    d15(paste0("d15_compare_", version), "分子式对照", paste(version, "正确式出现 / 首位 / 前五 / 可评估"),
        formula_comparison, {
          val <- as.numeric(tsv(formula_comparison)[[paste0(version, "_rank")]])
          c(sum(!is.na(val)), sum(val <= 1, na.rm = TRUE), sum(val <= 5, na.rm = TRUE),
            nrow(tsv(formula_comparison)))
        }, c(89, 88, 89, 95), historical = old15,
        caveat = "同一研究中可复用的两张切片；其他切片缺可比 MFC 参考")
  }
  independent <- "否；同一ROI的筛选结果及候选向量用于拟合与评估"
  members <- paste0(cluster, "cluster_members.tsv")
  stability <- paste0(cluster, "stability.tsv")
  reps <- paste0(cluster, "representatives.tsv")
  ktable <- paste0(cluster, "k_selection.tsv")
  retrieve <- "runs/d15_resume/retrieve_d15/neighbors.tsv"
  neighbors <- "runs/d15_neighbors_v2/neighbors_v2_d15/"
  d15("cluster_vectors", "聚类", "进入聚类的向量", members, nrow(tsv(members)), 55)
  d15("selected_k", "聚类", "合格 k 中差值最大的 k", ktable, {
    k <- subset(tsv(ktable), min_cluster_size >= 5)
    k$k[which.max(k$silhouette - k$null_mean)]
  }, 3, criterion = "k=2..8、每簇≥5、观测减置换均值最大；与保存选择一致",
  status = function(x) confirm(identical(x, tsv(ktable)$k[yes(tsv(ktable)$selected)])))
  qc <- paste0(cluster, "cluster_qc.json")
  d15("d15_cluster_qc", "聚类", "PCA20 累计解释率（%） / Davies–Bouldin", qc,
      with(j(qc), c(round(pca_explained_variance_cumulative * 100, 2), round(davies_bouldin, 3))),
      c(86.96, 2.074), method = "读保存统计量并按历史显示精度取值", historical = old15)
  d15("d15_k_delta", "聚类", "所选 k 的观测减置换均值", ktable,
      with(subset(tsv(ktable), yes(selected)), silhouette - null_mean))
  d15("d15_unembedded", "嵌入", "入选但未入聚类 feature", c(de, members),
      length(setdiff(tsv(de)$feature_id, tsv(members)$feature_id)), 3, historical = old15)
  for (i in 0:2) {
    prefix <- paste0("cluster_c", i)
    d15(paste0(prefix, "_size"), "聚类", paste0("C", i, " 成员"), members,
        sum(tsv(members)$cluster_id == i), c(14, 18, 23)[i + 1])
    d15(paste0(prefix, "_jaccard"), "稳定性", paste0("C", i, " Jaccard 中位数"), stability,
        subset(tsv(stability), cluster_id == i)$median_jaccard,
        c(0.666666666667, 1, 0.843434343434)[i + 1],
        "沿用历史 50 次 90% 子抽样", "读保存统计量")
    d15(paste0(prefix, "_representatives"), "代表", paste0("C", i, " 代表数"), reps,
        sum(tsv(reps)$cluster_id == i), 5, criterion = "按冻结规则应取5个",
        status = eq(5))
    view <- paste0(neighbors, "full_evidence/C", i, ".json")
    d15(paste0(prefix, "_neighbors"), "neighbors_v2", paste0("C", i, " 原始近邻 / 去重 / 合并保留"),
        c(retrieve, view), {
      n <- subset(tsv(retrieve), cluster_id == i)
      c(nrow(n), distinct(n$HMDB_ID), length(j(view)$semantic_neighbors))
    }, c(125, c(117, 78, 122)[i + 1], 25), historical = old15)
    d15(paste0(prefix, "_neighbors_order"), "neighbors_v2", paste0("C", i, " HMDB ID 排序一致"),
        c(retrieve, view), {
      n <- subset(tsv(retrieve), cluster_id == i)
      top <- aggregate(cosine ~ HMDB_ID, n, max)
      top <- top[order(-top$cosine, top$HMDB_ID, method = "radix"), ]
      identical(head(top$HMDB_ID, 25), vapply(j(view)$semantic_neighbors, `[[`, "", "HMDB_ID"))
    }, criterion = "最大余弦降序，平分按 HMDB ID，前25逐位一致", status = confirm)
  }

  independent <- "外部预登记参考；只作末端评价，仍是同一组织和同一研究"
  marker <- "runs/d15_neighbors_resume_markers/markers_d15/marker_annotations.tsv"
  clustered <- function() subset(tsv(marker), status == "CLUSTERED")
  d15("gold_markers_total", "标志物", "预登记条目", marker, nrow(tsv(marker)), 22)
  d15("gold_markers_clustered", "标志物", "进入聚类", marker, nrow(clustered()), 11)
  d15("gold_clustered_in_c2", "标志物", "归入 C2", marker, sum(clustered()$cluster_id == 2), 11)
  d15("d15_marker_formula", "标志物", "已入簇条目可评估 / 正确式出现 / 伪二级支持",
      marker, with(clustered(), c(sum(yes(formula_evaluable)), sum(yes(correct_formula_present), na.rm = TRUE),
                                   sum(yes(correct_formula_pseudo_ms2), na.rm = TRUE))),
      c(10, 9, 0), historical = old15)
  independent <- "否；D15p/q/r 围绕同一ROI与簇开展，事后补充证据"
  annotations <- paste0(d15r, "annotations.tsv")
  marker_compare <- paste0(d15r, "marker_comparison.tsv")
  marker_r <- paste0(d15r, "marker_concordance.tsv")
  for (model_name in c("primary", "cross")) {
    for (i in 0:2) {
      a <- subset(tsv(annotations), cluster_id == i & model == model_name)
      id <- paste0("name_c", i, "_", model_name)
      d15(id, "D15r 注释", paste0("C", i, " ", model_name, " 名称"), annotations,
          a$name, caveat = paste(a$evidence_level, a$status, "；模型名保留原文"), method = "读保存输出")
    }
    for (version in c("D15", "D15p", "D15q", "D15r")) {
      value <- subset(tsv(marker_compare), model == model_name)[[version]]
      d15(paste0(tolower(version), "_marker_", model_name), "标志物严格词匹配",
          paste(version, model_name, "匹配 / 评价"), marker_compare,
          frac(sum(value == "MATCH"), length(value)),
          paste0(c(D15 = 0, D15p = 3, D15q = 3, D15r = 0)[[version]], "/11"),
          historical = old15r, caveat = "D15 用类别词，p/q/r 用通路词；口径不同")
    }
    r <- subset(tsv(marker_r), model == model_name)
    d15(paste0("d15r_marker_compared_", model_name), "D15r 标志物", paste(model_name, "评价条目"), marker_r, nrow(r), 11)
    d15(paste0("d15r_marker_match_", model_name), "D15r 标志物", paste(model_name, "严格通路词匹配"), marker_r,
        sum(r$with_pathway == "MATCH"), 0)
  }
  pcov <- paste0(d15r, "candidate_coverage.tsv")
  for (i in 0:2) {
    for (source_name in c("HMDB_SMPDB", "Reactome")) {
      d15(paste0("pathway_usable_c", i, "_", tolower(source_name)), "D15r 通路证据",
          paste0("C", i, " ", source_name, " 可用候选行 / 代表候选行"), pcov,
          with(subset(tsv(pcov), cluster_id == i & source == source_name),
               frac(usable_evidence, candidate_rows)), method = "读保存计数并组合分子分母")
    }
  }
  req <- paste0("runs/d15_neighbors_v2_", c("primary", "cross"), "/llm_d15/request_index.json")
  requests <- function() c(j(req[1]), j(req[2]))
  d15("d15_model_controls_valid", "D15 双模型", "通过校验 / 原始簇及缺证据对照响应", req, {
    r <- requests()
    frac(sum(vapply(r, function(x) x$validation == "OK", TRUE)), length(r))
  }, "8/8", historical = old15)
  d15("d15_model_tokens", "D15 双模型", "实际最大输入 token", req,
      max(vapply(requests(), function(x) as.numeric(x$input_tokens), 0)),
      criterion = "严格低于 20000", status = function(x) confirm(x < 20000))
  rs <- paste0(d15r, "summary.json")
  for (m in c("primary", "cross")) {
    d15(paste0("d15r_control_", m), "D15r 证据对照",
        paste(m, "缺证据不足 / 采用换入通路 / 冲突引用换入证据"), rs, {
      r <- j(rs)$controls[[m]]
      c(r$missing_insufficient, r$adopted_injected_pathway, r$conflict_cites_injected)
    }, ref = if (m == "primary") c(TRUE, TRUE, FALSE) else c(TRUE, FALSE, TRUE),
    method = "读保存对照判定", historical = old15r,
    caveat = "对照探查输入依赖和校验器盲点")
  }
  d15("d15r_pathway_agreement", "D15r 双模型", "通路名一致 / 簇", annotations, {
    a <- subset(tsv(annotations), status == "OK")
    p <- merge(subset(a, model == "primary"), subset(a, model == "cross"), by = "cluster_id")
    frac(sum(nzchar(p$named_pathway.x) & p$named_pathway.x == p$named_pathway.y), length(unique(a$cluster_id)))
  }, "1/3", historical = old15r, caveat = "C2；交叉模型原证据等级为 CLASS_ONLY")
  d15("d15r_category_agreement", "D15r 双模型", "类别一致 / 簇", rs,
      {
        saved <- j(rs)$clusters
        frac(sum(vapply(saved, function(x) isTRUE(x$agreement$category), TRUE)), length(saved))
      },
      "0/3", method = "读保存统计量", historical = old15r)
  d15("d15r_annotation_valid", "D15r 双模型", "输出校验通过 / 簇模型输出", annotations,
      with(tsv(annotations), frac(sum(status == "OK"), nrow(tsv(annotations)))),
      caveat = "校验证据引用与词汇范围")
  reactome <- "runs/d15r_common/hmdb_reactome.tsv"
  d15("reactome_inchikey_matched", "D15r 通路映射", "全键匹配 HMDB ID",
      reactome, with(tsv(reactome), distinct(HMDB_ID[nzchar(matched_CHEBI)])), 25279)
  d15("reactome_mapped_hmdb", "D15r 通路映射", "有小鼠通路的 HMDB ID",
      reactome, with(tsv(reactome), distinct(HMDB_ID[nzchar(pathway_id)])), 10345)
  offline <- "runs/d15r_offline/offline_d15r/d14e_reranking.tsv"
  d15("d15r_offline_reranking", "D15r 离线对照", "挤出前八 / 首位改变 / Reactome 首位 / 簇", offline,
      with(tsv(offline), c(sum(yes(displaced_from_top8)), sum(yes(top1_changed)),
                            sum(yes(top1_reactome)), nrow(tsv(offline)))),
      c(0, 3, 3, 13), historical = old15r, caveat = "用 D14e 模拟数据，非独立生物验证")

  manifest <- "deliverables/LOCKED_D15/manifest.tsv"
  manifest_rows <- tsv(manifest)
  manifest_rows$ref <- wk(sub("^.*/msidataset_week_20261105/", "", manifest_rows$source))
  registered <- manifest_rows[manifest_rows$ref %in% ls(hashes, all.names = TRUE), ]
  add("d15_manifest_hashes", "来源完整性", "核对锁定清单中本次使用的原始来源",
      "SHA-256 一致 / 共有来源", c(wk(manifest), registered$ref), {
        observed <- vapply(registered$ref, function(x) hashes[[x]], "")
        frac(sum(observed == registered$sha256), nrow(registered))
      }, "本次可比原始来源的 SHA-256 均等于 LOCKED_D15 清单",
      function(x) confirm(x == frac(nrow(registered), nrow(registered))),
      caveat = "仅检查锁定清单与本表直接使用的 D15 来源交集")

  section <- "内部复现"
  dataset <- "ADP（Task2 / T1 / T6）"
  independent <- "否；同一 ADP 标准品切片和黄金划分用于开发回归"
  golden <- "../msidataset_task2_runs/20261016_01/MSIDataset_Task2_20261015/outputs/adp_core_closeout_20261016_02/all_cluster_members_distance.tsv"
  previous <- "WEEK/deliverables/REPORT_ZH.md"
  labels <- function(run) paste0("runs/", run, "/cluster/", c("cluster_members.tsv", "noise.tsv"))
  partition <- function(files) {
    pieces <- lapply(files, tsv)
    typed <- lapply(pieces, function(x) x[intersect(names(x), c("feature_id", "cluster_id"))])
    do.call(rbind, typed[vapply(typed, nrow, 0L) > 0L])
  }
  for (run in c("t1_adp", "t6_repro_adp")) {
    tag <- if (run == "t1_adp") "t1" else "t6"
    files <- labels(run)
    if (run == "t6_repro_adp") {
      add("adp_t6_partition", "HDBSCAN", "黄金划分规模", "向量 / C0 / C1 / C2 / noise",
          wk(files), c(nrow(partition(files)), size(tsv(files[1])$cluster_id), nrow(tsv(files[2]))),
          "HDBSCAN 0.8.44、15/5/EOM", ref = c(513, 69, 64, 105, 275), historical = previous)
    }
    add(paste0("adp_", tag, "_golden_ari"), "HDBSCAN", "黄金划分逐 feature 比对", "ARI（含噪声）",
        wk(c(golden, files)), {
      p <- partition(files)
      g <- tsv(golden)
      stopifnot(!anyDuplicated(g$feature_id), !anyDuplicated(p$feature_id),
                all(g$feature_id %in% p$feature_id), all(p$cluster_id[p$feature_id %in% g$feature_id] >= 0))
      truth <- g$cluster_id[match(p$feature_id, g$feature_id)]
      truth[is.na(truth)] <- -1
      validation_ari(table(p$cluster_id, truth))
    }, "ARI=1", eq(1), 1, "WEEK/reports/T6/reproduction_check.tsv",
    "黄金文件记录已入簇成员；同一 513 向量集合中其他 feature 作为噪声")
  }
  adp_reps <- paste0("runs/", c("t1_adp", "t6_repro_adp"), "/represent/cluster_representatives.tsv")
  add("adp_reps", "代表", "代表表逐字段复现", "C0 / C1 / C2 代表数",
      wk(adp_reps), size(tsv(adp_reps[2])$cluster_id),
      "与 T1 表逐字段一致", function(x) confirm(identical(tsv(adp_reps[1]), tsv(adp_reps[2]))),
      c(7, 7, 11), previous)
  standards <- "reports/T5a/five_standards_recall.tsv"
  add("adp_standards", "MFC", "五标准品召回回归", "Asp / Glu / Citrate / AMP / ADP 参考排名",
      wk(standards), tsv(standards)$reference_ranks,
      "五标准品候选及分数并列按保存基准一致",
      function(x) confirm(all(yes(tsv(standards)$reference_candidate_recalled) &
                                tsv(standards)$comparison_status == "PASS")),
      c(6, 4, 1, 10, 1), previous, "Glu、AMP 参考候选存在并列", "读保存的逐项秩")
  dataset <- "E1 膀胱（T6）"
  independent <- "否；E1 单切片开发集，与 T4 及 T3 的产物交叉比较"
  for (file in c("ingest/peaks.tsv", "de/de_features.tsv", "mfc/mfc_candidates.tsv",
                 "spatial_modules/spatial_modules.tsv", "spatial_modules/module_members.tsv")) {
    original <- if (startsWith(file, "spatial_modules")) "t4_bladder_spatial_retry" else "t4_bladder"
    old <- paste0("runs/", original, "/", file)
    current <- paste0("runs/t6_repro_bladder/", file)
    id <- paste0("t6_", tools::file_path_sans_ext(basename(file)))
    add(id, sub("/.*", "", file), "同一数据既有结果逐字段复现", "重跑产物行数", wk(c(old, current)),
        nrow(tsv(current)), "与原表逐字段一致",
        function(x) confirm(identical(tsv(old), tsv(current))),
        nrow(tsv(old)), reference_source = wk(old), caveat = if (grepl("mfc", file)) "复用候选缓存" else "")
  }
  t6 <- "reports/T6/reproduction_check.tsv"
  for (spec in list(c("formal", "^acceptance query", "正式卡及完整卡逐字段相等", 8, 8),
                    c("release", "^query ", "早期卡及完整卡逐字段相等", 0, 8),
                    c("view_sha", "^view SHA ", "视图逐字节 SHA 一致", 0, 15),
                    c("full_sha", "^full evidence SHA ", "全量证据 SHA 一致", 15, 15))) {
    checks <- tsv(t6)
    checks <- checks[grepl(spec[2], checks$item), ]
    files <- strsplit(checks$source_paths, "; ", fixed = TRUE)
    files <- lapply(files, function(pair) {
      stopifnot(all(grepl("(^|/)runs/", pair)))
      wk(sub("^.*/msidataset_week_20261105/", "", pair))
    })
    results <- function() {
      vapply(files, function(pair) {
        stopifnot(length(pair) == 2L)
        if (endsWith(spec[1], "sha")) return(identical(hashes[[pair[1]]], hashes[[pair[2]]]))
        identical(read(pair[1]), read(pair[2]))
      }, TRUE)
    }
    add(paste0("t6_", spec[1]), "T6 复现审计", spec[3], "一致 / 比较文件",
        c(wk(t6), unlist(files)), frac(sum(results()), length(files)),
        "每个文件通过原比对", function(x) confirm(all(results())), frac(spec[4], spec[5]),
        historical = wk(t6), caveat = switch(spec[1], release = "新增空 annotation_status 字段",
                                       view_sha = "完全证据路径字段含新的绝对路径", ""),
        method = if (endsWith(spec[1], "sha")) "重算文件 SHA-256" else "重比 JSON 字段")
  }
  dataset <- "ADP 与 E1（T6）"
  independent <- "否；既有原始模型响应及新旧同源证据视图"
  revalidated <- "runs/t6_repro_audit/report/revalidated_responses.tsv"
  add("t6_validator", "模型校验复现", "逐条复核缓存", "错误判定一致 / 响应",
      wk(revalidated), with(tsv(revalidated), frac(sum(same == "PASS"), length(same))),
      "错误集合逐项相等", function(x) confirm(all(tsv(revalidated)$same == "PASS")),
      "30/30", previous, "其中有失败响应，另行统计")
  add("t6_validator_failed", "模型校验复现", "区分校验失败与回归一致", "校验失败响应",
      wk(revalidated), sum(tsv(revalidated)$status != "OK"), ref = 4, historical = previous)
  consensus <- "deliverables/two_model_consensus.tsv"
  for (module_name in c("M036", "M020")) {
    add(paste0("t6_", tolower(module_name), "_levels"), "双模型复核",
        "核对开发阶段证据层级", paste(module_name, "本地 / API 证据层级"), wk(consensus), {
          rec <- tsv(consensus)
          rec <- rec[rec$module == module_name, ]
          c(rec$local_level, rec$api_level)
        }, ref = if (module_name == "M036") rep("CLASS_ONLY", 2) else c("CLASS_ONLY", "PATHWAY_CLUE"),
        historical = previous, caveat = "使用已保存的响应审计表", method = "读机器逐模块记录")
  }
  add("t6_model_consensus", "双模型复核", "保存的共同类别裁定", "逐模块审计条目",
      wk(consensus), nrow(tsv(consensus)), caveat = "M036 两条均 CLASS_ONLY；M020 层级有别")

  section <- "其他数据集"
  dataset <- "E1 膀胱（D2–D8a）"
  independent <- "否；E1 单切片开发集，D5/D7 聚类尺度有事后调整"
  labels_compare <- function(a, b) {
    p <- merge(partition(a), partition(b), by = "feature_id")
    t <- table(p$cluster_id.x, p$cluster_id.y)
    c(nrow(p), validation_ari(t))
  }
  d5 <- labels("d5_bladder")
  d7 <- labels("d7_bladder")
  d7report <- "WEEK/deliverables/D7/REPORT_D7_ZH.md"
  add("d4_cluster_robustness", "嵌入库版本对照", "固定 HDBSCAN 对卡片描述的敏感性",
      "v1 簇数 / v2 簇数", wk("reports/D4b/compare_v1_v2.tsv"), {
    x <- subset(tsv("reports/D4b/compare_v1_v2.tsv"), section == "cluster")
    c(sum(x$version == "v1"), sum(x$version == "v2"))
  }, ref = c(2, 0), historical = "WEEK/deliverables/D4/REPORT_D4_ZH.md")
  for (run in c("d2_bladder", "d4_bladder", "d5_bladder")) {
    f <- labels(run)
    stage_id <- sub("_bladder$", "", run)
    ref <- switch(stage_id, d2 = c(2, 64, 118), d4 = c(0, 118, 118), d5 = c(2, 68, 118))
    history <- switch(stage_id, d2 = "D1/REPORT_D1_ZH.md", d4 = "D4/REPORT_D4_ZH.md",
                      d5 = "D5/REPORT_D5_ZH.md")
    add(paste0(stage_id, "_bladder_labels"), "膀胱 HDBSCAN", "同切片聚类敏感性",
        "簇数 / 噪声 / 输入", wk(f),
        c(distinct(tsv(f[1])$cluster_id), nrow(tsv(f[2])), nrow(partition(f))),
        ref = ref, historical = paste0("WEEK/deliverables/", history),
        caveat = if (stage_id == "d5") "D5 事后调整最小规模至 10" else "使用冻结 15/5/EOM")
  }
  f <- c(labels("d2_bladder"), labels("d4_bladder"))
  add("d4_bladder_ari", "膀胱 HDBSCAN", "固定参数使用两个嵌入库",
      "共同 feature / ARI", wk(f), labels_compare(f[1:2], f[3:4]),
      ref = c(118, 0), historical = "WEEK/deliverables/D4/REPORT_D4_ZH.md")
  f <- c(labels("d2_bladder"), d5)
  add("d5_bladder_ari", "膀胱 HDBSCAN", "同时变化的库和最小簇规模",
      "共同 feature / ARI", wk(f), {
        a <- labels_compare(f[1:2], f[3:4])
        c(a[1], round(a[2], 3))
      }, ref = c(118, 0.808), historical = "WEEK/deliverables/D5/REPORT_D5_ZH.md")
  for (task in c("D2", "D5", "D7")) {
    control <- paste0("reports/", task, "/negative_control_comparison.tsv")
    if (task %in% c("D2", "D7")) {
      add(paste0(tolower(task), "_control"), "缺证据对照", "记录模型面对证据删除的结果",
          "保存的对照行", wk(control), nrow(tsv(control)),
          caveat = "删除对照与完整证据使用同一模块，逐模型原判定见源表")
    }
  }
  add("d7_pseudo", "伪二级", "扩充碎片池后检查查询", "查询谱 / 碎片行",
      wk("reports/D7/pseudo_fragments.tsv"),
      with(tsv("reports/D7/pseudo_fragments.tsv"), c(distinct(anchor_feature_id), length(anchor_feature_id))),
      ref = c(40, 276), historical = d7report)
  add("d7_precursor", "谱库匹配", "按冻结前体门统计", "前体查询 / 合格查询–库谱对",
      wk("reports/D7/precursor_coverage.tsv"),
      with(tsv("reports/D7/precursor_coverage.tsv"),
           c(sum(eligible_library_spectra > 0), sum(eligible_library_spectra))),
      ref = c(7, 316), historical = d7report)
  add("d7_spectral_support", "谱库匹配", "评估扩库是否补回谱支持", "有谱图支持 feature",
      wk("reports/D7/ms2_supported_features.tsv"), nrow(tsv("reports/D7/ms2_supported_features.tsv")),
      ref = 0, historical = d7report, caveat = "前体通过不等于碎片匹配")
  bladder_mfc <- "runs/d7_bladder/mfc/mfc_candidates.tsv"
  add("d7_mfc", "MFC", "扩库后候选覆盖", "候选行 / 有候选 feature", wk(bladder_mfc),
      with(tsv(bladder_mfc), c(nrow(tsv(bladder_mfc)), distinct(feature_id))),
      ref = c(4226, 138), historical = d7report)
  add("d7_d5_ari", "聚类", "扩库后的划分回归", "共同 feature / ARI（噪声作为一类）",
      wk(c(d5, d7)), labels_compare(d5, d7),
      "ARI=1", function(x) confirm(x[2] == 1), c(118, 1), d7report)
  add("d7_stability", "稳定性", "稳定性描述", "C0 / C1 达标次数与分母",
      wk("reports/D7/stability.tsv"), with(tsv("reports/D7/stability.tsv"), frac(reproduced, subsamples)),
      ref = c("40/50", "39/50"), historical = d7report, method = "读保存的子抽样计数")
  add("d7_models", "双模型", "类别结论与两模型校验", "一致簇数 / 簇数",
      wk("reports/D7/model_summary.json"),
      with(j("reports/D7/model_summary.json"), frac(counts$agreement, clusters)),
      method = "读保存统计量", caveat = "一致为候选类别级")
  dataset <- "ADP（D8b）"
  independent <- "否；复用 ADP 开发样本、谱库及原有验证结果"
  add("d8b_tiers", "谱证据展示复核", "三层候选证据计数", "Tier 1 / 2 / 3 候选行",
      wk("reports/D8b/tier_counts.tsv"), tsv("reports/D8b/tier_counts.tsv")$candidate_rows,
      ref = c(40, 31, 15543), historical = "WEEK/reports/D8b/gate.json", method = "读保存计数")

  dataset <- "Wang 肾脏（D9–D13）"
  independent <- "外部数据来源；同切片在D9–D13反复使用，区域及筛选流程经过开发"
  for (task in c("D9", "D10")) {
    f <- paste0("reports/", task, "/gold_standard_eval.tsv")
    add(paste0(tolower(task), "_gold_detected"), "肾脏参考评价", "金标准检出",
        "检出 / 参考", wk(f), with(tsv(f), frac(sum(yes(detected)), length(detected))),
        ref = "13/14", historical = paste0("WEEK/deliverables/", task, "/REPORT_", task, "_ZH.md"),
        caveat = "检出含串链 bin；后续方向评价剔除")
    add(paste0(tolower(task), "_gold_de_top1"), "肾脏参考评价", "通过差异筛选的参考分子式回归",
        "正确式首位（含并列） / 已筛选参考", wk(f), {
          x <- subset(tsv(f), status == "OK" & de_status == "DE")
          frac(sum(num(x$correct_formula_rank) <= 1, na.rm = TRUE), nrow(x))
        }, caveat = "分母只包括通过差异筛选的参考物")
    add(paste0(tolower(task), "_gold_direction"), "肾脏参考评价", "可评价参考方向一致",
        "一致 / 可评价（含非差异）", wk(f), {
      x <- subset(tsv(f), status == "OK")
      frac(sum(yes(x$direction_agrees)), nrow(x))
    }, caveat = "探索性单切片")
  }
  control <- "reports/D10b/negative_control_comparison.tsv"
  add("d10b_control", "注释对照", "缺证据输入检查", "证据不足且通过 / 模型",
      wk(control), with(tsv(control),
       frac(sum(control_validation == "OK" & control_evidence_level == "INSUFFICIENT"), nrow(tsv(control)))),
      caveat = "同一已有簇的证据移除对照")
  inv <- "reports/D11/mfc_invariance.tsv"
  add("d11_mfc_invariance", "MFC", "扩大特征范围后的分数和排名变化", "共享 / 共享分数相等 / 共享排名相等",
      wk(inv), with(subset(tsv(inv), `_merge` == "both"),
                         c(length(priority_equal), sum(yes(priority_equal)), sum(yes(rank_equal)))),
      caveat = "分数变化审计，未覆盖共享行以外的候选")
  for (organ in c("kidney", "liver")) {
    liver <- organ == "liver"
    dataset <- if (liver) "Samarah 八鼠肝脏（D12/D13）" else "Wang 肾脏（D13）"
    independent <- if (liver) "否；同源八鼠开发批次，包含D15鼠7" else
      "外部数据来源；D10与D13复用同一肾脏切片"
    d13 <- paste0("runs/d13_", organ, "/")
    report <- "WEEK/deliverables/D13/REPORT_D13_ZH.md"
    ident <- paste0(d13, "identity_d13/identity.tsv")
    tag <- paste0("d13_", organ)
    add(paste0(tag, "_formula"), "候选过滤", "评估过滤导致的正确式损失",
        "评价 / 过滤前出现 / 过滤后出现 / 损失 / 过滤后首位 / 前五", wk(ident),
        with(tsv(ident), c(length(feature_id), sum(yes(correct_formula_before)),
                          sum(yes(correct_formula_present)), sum(yes(lost_correct_formula)),
                          sum(yes(top1_including_ties)), sum(yes(top5_including_ties)))),
        ref = if (liver) c(95, 89, 89, 0, 88, 89) else c(11, 11, 11, 0, 11, 11),
        historical = report, caveat = "第一项是可评价分母，非覆盖分母；排名含并列")
    removed <- paste0(d13, "mfc/removed_candidates.tsv")
    final <- paste0(d13, "mfc/mfc_candidates_filtered.tsv")
    add(paste0(tag, "_removed"), "候选过滤", "无嵌入候选剔除", "删除 / 原候选行",
        wk(c(removed, final)), frac(nrow(tsv(removed)), nrow(tsv(removed)) + nrow(tsv(final))),
        ref = if (liver) "186/4687" else "786/7329", historical = report)
    fs <- paste0(d13, "mfc/feature_status.tsv")
    add(paste0(tag, "_de_removed"), "候选过滤", "差异 feature 的候选删除与缺失变化",
        "删除 / 原候选行；受影响 / 新无候选 / 总无候选", wk(fs), {
          x <- subset(tsv(fs), yes(is_de))
          with(x, c(frac(sum(removed), sum(candidates_before)), sum(removed > 0),
                    sum(candidates_before > 0 & candidates_after == 0), sum(candidates_after == 0)))
        }, ref = if (liver) c("71/2074", 25, 0, 0) else c("317/3003", 175, 14, 142),
        historical = report)
    cv <- paste0(d13, "embed_d13/coverage.tsv")
    mem <- paste0(d13, "cluster_d13/cluster_members.tsv")
    add(paste0(tag, "_cluster"), "聚类", "过滤后的聚类覆盖", "差异 / 聚类 / 簇数",
        wk(c(cv, mem)),
        c(nrow(tsv(cv)), nrow(tsv(mem)), length(unique(tsv(mem)$cluster_id))),
        ref = if (liver) c(57, 57, 4) else c(477, 335, 6), historical = report)
    cont <- paste0(d13, "cluster_d13/cluster_correspondence.tsv")
    add(paste0(tag, "_ari"), "聚类比较", "与前序划分在共同 feature 上比对", "共同 feature / ARI",
        wk(cont), {
      tab <- as.matrix(tsv(cont)[-1])
      c(sum(tab), round(validation_ari(tab), 3))
    }, ref = if (liver) c(40, 0.758) else c(217, 0.607), historical = report,
    caveat = "肝脏对 D12b，肾脏对 D10，分母限制在共同 feature")
    anno <- paste0(d13, "report_d13/cluster_annotations.tsv")
    add(paste0(tag, "_models"), "双模型", "同首类证据 ID 且无首位并列", "一致 / 全簇", wk(anno),
        with(tsv(anno), frac(sum(yes(models_agree)), length(models_agree))),
        ref = if (liver) "3/4" else "2/6", historical = report,
        caveat = "分母保留交叉模型失败；类别匹配与身份不同")
  }
  d12de <- "runs/d12_liver/de_d12/de_features.tsv"
  add("d12_de", "八鼠差异筛选", "四舍五入前的保存筛选计数", "差异 / PC / PP",
      wk(d12de), with(tsv(d12de), c(length(direction), sum(direction == "pericentral"),
                                    sum(direction == "periportal"))),
      ref = c(57, 39, 18), historical = "WEEK/deliverables/D12/REPORT_D12_ZH.md")
  d12id <- "reports/D12/mfc_identity_all.tsv"
  add("d12_formula", "同研究参考复核", "参考式可评价时的 MFC 排名",
      "正确式出现 / 首位 / 前五 / 可评估", wk(d12id), {
    x <- subset(tsv(d12id), reference_status == "reliable")
    with(x, c(sum(yes(correct_formula_present)), sum(yes(top1_including_ties)),
              sum(yes(top5_including_ties)), nrow(x)))
  }, ref = c(89, 89, 89, 95), historical = "WEEK/deliverables/D12/REPORT_D12_ZH.md",
      caveat = "参考仅用于末端；依赖同一研究预选离子")
  gold12 <- "reports/D12/gold_standard_eval.tsv"
  add("d12_marker_formula", "标志物", "原研究参考式", "首位 / 前五 / 可评估 / 预登记",
      wk(gold12), {
        x <- tsv(gold12)
        rank <- num(x$correct_formula_rank[x$reference_status == "reliable"])
        c(sum(rank <= 1, na.rm = TRUE), sum(rank <= 5, na.rm = TRUE), length(rank), nrow(x))
      }, ref = c(20, 20, 21, 22), historical = "WEEK/deliverables/D12/REPORT_D12_ZH.md")
  for (branch in c("liver", "kidney")) {
    f <- if (branch == "liver") "runs/d12_liver/cluster_d12b/k_selection.tsv" else
      "runs/d12_liver/cluster_d12b/kidney_k_selection.tsv"
    dataset <- if (branch == "liver") "Samarah 八鼠肝脏（D12b）" else "Wang 肾脏（D12b）"
    add(paste0("d12b_", branch, "_k"), "事后敏感性", "最小簇约束下重选 k", "所选 k",
        wk(f), tsv(f)$k[yes(tsv(f)$selected)], ref = if (branch == "liver") 4 else 5,
        historical = "WEEK/deliverables/D12/REPORT_D12_ZH.md")
  }

  dataset <- "五通路模拟（D14–D14e）"
  independent <- "否；人工构造的已知标签，后续证据规则在同一模拟集上迭代"
  for (arm in c("d14_A", "d14_B", "d14b_A")) {
    folder <- if (arm == "d14b_A") "evaluate_d14b" else "evaluate_d14"
    f <- paste0("runs/", arm, "/", folder, "/assignments_evaluation.tsv")
    add(paste0(tolower(arm), "_partition"), "模拟划分评价", "与五通路真值比较",
        "feature / 簇 / ARI / NMI / 加权纯度", wk(f), {
      x <- tsv(f)
      tab <- table(x$cluster_id, x$pathway)
      prob <- tab / sum(tab)
      expected <- outer(rowSums(prob), colSums(prob))
      entropy <- function(x) -sum(x * log(x))
      nmi <- sum(prob[prob > 0] * log((prob / expected)[prob > 0])) /
        mean(c(entropy(rowSums(prob)), entropy(colSums(prob))))
      c(nrow(x), nrow(tab), round(validation_ari(tab), 3), round(nmi, 3),
        round(sum(apply(tab, 1L, max)) / sum(tab), 3))
    }, ref = switch(arm, d14_A = c(105, 8, 0.443, 0.585, 0.752),
                    d14_B = c(108, 6, 0.814, 0.839, 0.917), d14b_A = c(108, 7, 0.498, 0.686, 0.815)),
    historical = "WEEK/deliverables/D14b/REPORT_D14b_ZH.md",
    caveat = "B 为提供理想候选的上限参照；D14b 事后增加失水加合物")
    mf <- sub("assignments_evaluation", "mfc_feature_evaluation", f)
    add(paste0(tolower(arm), "_recall"), "模拟分子式回收", "按加合物统计正确式",
        "[M-H]- / [M+Cl]- / [M-H2O-H]- 正确式出现", wk(mf), {
      x <- tsv(mf)
      vapply(c("[M-H]-", "[M+Cl]-", "[M-H2O-H]-"), function(adduct) {
        y <- x[x$adduct == adduct, ]
        frac(sum(yes(y$correct_formula_present)), nrow(y))
      }, "")
    }, ref = switch(arm, d14_A = c("36/36", "36/36", "0/36"),
                    d14_B = rep("36/36", 3), d14b_A = c("36/36", "36/36", "36/36")),
    historical = "WEEK/deliverables/D14b/REPORT_D14b_ZH.md")
  }
  for (arm in c("A", "B", "b_A")) {
    task <- if (arm == "b_A") "D14b" else "D14"
    f <- paste0("reports/", task, "/annotation_summary.json")
    add(paste0("d14_models_", arm), "模拟注释评价", "已有簇多数通路的词匹配",
        paste(arm, "主模型 / 交叉模型命中数与有效分母"), wk(f), {
          x <- j(f)[[arm]]
          c(frac(x$primary$majority_hits, x$primary$valid), frac(x$cross$majority_hits, x$cross$valid))
        }, ref = switch(arm, A = c("3/8", "6/7"), B = c("2/4", "4/4"), b_A = c("6/7", "6/7")),
        historical = "WEEK/deliverables/D14b/REPORT_D14b_ZH.md", method = "读保存的独立模型计数")
  }
  for (task in c("d14a", "d14b")) {
    f <- if (task == "d14a") "runs/d14_A/evaluate_d14/mfc_feature_evaluation.tsv" else
      "runs/d14b_A/evaluate_d14b/mfc_feature_evaluation.tsv"
    add(paste0(task, "_adduct_top1"), "模拟分子式回收", "按加合物计算首位（含并列）", "三个加合物首位比",
        wk(f), {
          x <- tsv(f)
          vapply(c("[M-H]-", "[M+Cl]-", "[M-H2O-H]-"), function(adduct) {
            y <- x[x$adduct == adduct, ]
            frac(sum(yes(y$top1_including_ties)), nrow(y))
          }, "")
        }, ref = if (task == "d14a") c("31/36", "23/36", "0/36") else c("24/36", "17/36", "20/36"),
        historical = "WEEK/deliverables/D14b/REPORT_D14b_ZH.md")
  }
  for (task in c("d14c", "d14d", "d14e")) {
    f <- paste0("runs/", task, "_evaluation/evaluate_", task, "/annotations.tsv")
    for (arm in c("A", "B")) for (m in c("primary", "cross")) {
      add(paste(task, arm, m, sep = "_"), "模拟注释评价", "保存响应的通路词命中",
          paste(arm, m, "宽松 / 名称严格 / 有效 / 全簇"), wk(f), {
        x <- tsv(f)
        ver <- if (task == "d14c") x$setting == "pathway" else
          tolower(x$version) == task
        whole <- x[ver & x$arm == arm & x$model == m, ]
        valid <- subset(whole, status == "OK")
        strict <- if (task == "d14c") valid$strict_hit else valid$name_strict
        c(sum(yes(valid$loose_hit)), sum(yes(strict)), nrow(valid), nrow(whole))
      }, caveat = "同一模拟集上的事后证据/提示变化；分母保留失败记录")
    }
  }

  section <- "R 包复现"
  independent <- "否；用 D15 既有产物作 Python/R 或包内技术回归"
  pkg_metric <- function(task, id, stage, metric, key, criterion, status,
                         ref = "", caveat = "", source = paste0("reports/", task, "/gate.json")) {
    add(id, stage, paste(task, stage, "验证"), metric, pk(source),
        field(read(pk(source)), key), criterion, status, ref, caveat = caveat,
        method = "读保存的机器回归记录")
  }
  for (task in c("PK1", "PK1b", "PK1c", "PK2", "PK2b", "PK2c")) {
    dataset <- paste("D15 R 包", task)
    f <- paste0("reports/", task, "/gate.json")
    add(paste0(tolower(task), "_gate"), "阶段验收", "保留阶段原验收状态", "原 gate 状态",
        pk(f), read(pk(f))$status, "阶段各自的验收标准", read(pk(f))$status,
        method = "读机器 gate", caveat = "阶段状态沿用原记录，细项见同阶段行")
    if (task == "PK1") {
      for (m in c("selected", "pericentral", "periportal", "n_tiles",
                  "log2FC_max_relative_error", "p_max_relative_error", "fdr_max_relative_error")) {
        quality <- endsWith(m, "_error")
        pkg_metric(task, paste0("pk1_de_", m), "差异回归", m, c("regression", "differential", m),
                   if (quality) "相对误差 <1e-6" else "计数与锁定 D15 一致",
                   if (quality) function(x) confirm(x < 1e-6) else
                     eq(switch(m, selected = 58L, pericentral = 38L, periportal = 20L, n_tiles = 393L)))
      }
      pkg_metric(task, "pk1_pseudo_pairs", "伪二级", "成对键集合一致",
                 c("regression", "pseudo_ms2", "pair_set_identical"), "成对键集合完全相同", confirm)
      pkg_metric(task, "pk1_vectors", "向量转换", "随机1000行 Float32逐值一致",
                 c("regression", "vector_conversion", "random_1000_exact"), "完全相同", confirm)
    }
    if (task == "PK1b") {
      for (m in c("keys_identical", "key_multiset_equal", "rank_identical", "max_score_relative_error")) {
        error <- m == "max_score_relative_error"
        pkg_metric(task, paste0("pk1b_mfc_", m), "MFC", m, m,
          if (error) "相对误差 <1e-6" else "完全一致",
          if (error) function(x) confirm(x < 1e-6) else confirm,
          caveat = "包副本连接了锁定的预计算 step6",
          source = "reports/PK1b/connected_regression.json")
      }
      for (m in c("evaluated", "present", "top1_including_ties")) {
        expected <- c(evaluated = 95, present = 89, top1_including_ties = 88)[[m]]
        pkg_metric(task, paste0("pk1b_formula_", m), "分子式评估", m, m,
                   "与 D15 保存评价一致", eq(expected), expected,
                   "作者式只用于末端评估", "reports/PK1b/formula_regression.json")
      }
    }
    if (task == "PK1c") {
      for (m in c("matching_pair_multisets_equal", "maximum_step56_relative_error",
                  "mfc_keys_order_ranks_identical")) {
        error <- startsWith(m, "maximum")
        pkg_metric(task, paste0("pk1c_", m), "step5/6 与 MFC", m, c("regression", m),
                   if (error) "相对误差 <1e-6" else "冻结回归检查为真",
                   if (error) function(x) confirm(x < 1e-6) else confirm,
                   caveat = "隔离副本内由原始输入重建谱匹配，不依赖 step6 缓存")
      }
    }
    if (task == "PK2") {
      kcomp <- "reports/PK2/k_selection_comparison.tsv"
      kt <- function() read(pk(kcomp))
      add("pk2_selected_k", "select_k", "R 对 Python 所选 k", "R 选 k",
          pk(kcomp), kt()$k[yes(kt()$selected_R)], "所选 k 相同",
          function(x) confirm(identical(x, kt()$k[yes(kt()$selected_D15)])), 3)
      add("pk2_k_observed_fail", "select_k", "逐 k 观测 silhouette 原验收", "超差 k",
          pk(kcomp), with(kt(), k[abs(silhouette_R - silhouette_D15) >= 1e-6]),
          "所有 k 误差 <1e-6", function(x) if (length(x)) "FAIL" else "PASS",
          caveat = "k=4..8 在两套实现中均不合格")
      for (k in 2:8) {
        a <- kt()[kt()$k == k, ]
        old <- tsv(ktable)[tsv(ktable)$k == k, ]
        add(paste0("pk2_k_delta_", k), "select_k", "绘图用 R 版选择曲线", paste("k", k, "观测减基线"),
            c(pk(kcomp), wk(ktable)), a$silhouette_R - a$null_mean_R,
            ref = old$silhouette - old$null_mean, historical = wk(ktable),
            caveat = "历史基准差异；PKG-D1采用R基准，最小簇规则仍生效", reference_source = wk(ktable))
        add(paste0("pk2_k_observed_", k), "select_k", "逐 k 与锁定观测 silhouette 对比",
            paste("k", k, "观测 silhouette"), c(pk(kcomp), wk(ktable)), a$silhouette_R,
            "原登记误差<1e-6", function(x) confirm(abs(x - old$silhouette) < 1e-6),
            old$silhouette, wk(ktable), "历史基准差异；k=4–8 在两种实现中均不合格",
            reference_source = wk(ktable))
      }
      pkg_metric(task, "pk2_embedding", "嵌入", "向量最大绝对误差",
                 c("regression", "embedding", "max_abs_error"), "<1e-6",
                 function(x) confirm(x < 1e-6))
      pkg_metric(task, "pk2_ari", "聚类", "重新标簇后的 ARI", c("regression", "clusters", "ari"), "ARI=1", eq(1))
      add("pk2_stability_original_fail", "稳定性", "原验收的 C0 中位 Jaccard 差",
          "R 与 Python 的 C0 差值（gate 四位精度）", pk(f), {
        x <- read(pk(f))$regression$clusters
        abs(x$stability_original_cluster_0_new - x$stability_original_cluster_0_reference)
      }, "差 ≤0.1", function(x) confirm(x <= 0.1), caveat = "原验收未通过；PKG-D1 改为各簇 >0.5",
      method = "以保存的回归统计量求差")
      pkg_metric(task, "pk2_neighbors_original_fail", "近邻", "逐代表 top25 错位条数（共375位）",
                 c("regression", "neighbors", "positions_mismatched"), "ID 与顺序逐位一致",
                 eq(0), 0, "F139 第4、5位互换；改用合并后检查见 PKG-D1")
      pkg_metric(task, "pk2_representatives", "代表", "对齐簇后 ID、角色及顺序一致",
                 c("regression", "representatives", "feature_role_order_identical_after_relabel"),
                 "对齐后逐位一致", confirm)
      pkg_metric(task, "pk2_request", "请求", "原始请求用户消息逐字节一致",
                 c("regression", "request", "user_message_string_exact"), "逐字节一致", confirm)
    }
    if (task == "PK2b") {
      for (m in c("mfc_rows_order_and_ranks_identical", "all_stabilities_gt_0_5",
                  "representatives_and_roles_equal_after_cluster_alignment",
                  "merged_neighbor_ids_and_order_equal", "request_user_message_exact_C0_C1_C2")) {
        pkg_metric(task, paste0("pk2b_", m), "连接及回归", m, c("regression", m),
                   "按 PKG-D1 或同输入历史请求口径", confirm, caveat = "使用锁定的基础视图和已有映射")
      }
      pkg_metric(task, "pk2b_base_evidence", "证据", "新基础视图语义比较",
                 c("regression", "new_base_evidence_semantics"), "基础视图从头重建", "PARTIAL",
                 caveat = "到 PK2c 才实现")
    }
    if (task == "PK2c") {
      pkg_metric(task, "pk2c_mfc_error", "MFC", "最大相对误差", c("regression", "mfc_max_relative_error"),
                 "<1e-6", function(x) confirm(x < 1e-6))
      pkg_metric(task, "pk2c_ari", "聚类", "ARI", c("regression", "adjusted_rand_index"), "ARI=1", eq(1), 1)
      pkg_metric(task, "pk2c_k", "聚类", "所选 k", c("regression", "selected_k"), "所选 k 一致", eq(3), 3)
      pkg_metric(task, "pk2c_stability", "聚类", "所有簇稳定性中位数 >0.5",
                 c("regression", "stability_all_greater_than_0_5"), "每簇 >0.5", confirm)
      ev <- "reports/PK2c/evidence_comparison.tsv"
      for (m in c("base_semantic_except_provenance", "d15r_semantic_except_provenance",
                  "request_bytes_after_provenance_substitution", "actual_new_request_identical")) {
        success <- function(x) x == "PASS" | yes(x)
        add(paste0("pk2c_", m), "证据与请求", "从原始 CSV 到证据视图", paste(m, "通过 / 簇"),
            pk(ev), {
          x <- read(pk(ev))[[m]]
          frac(sum(success(x)), length(x))
        }, "各簇一致（按本行口径）",
        function(x) confirm(all(success(read(pk(ev))[[m]]))),
        caveat = "新请求真实路径与哈希不同；完整链使用隔离副本和共享补丁")
      }
      pkg_metric(task, "pk2c_saved_responses", "缓存响应", "六条 D15r 响应校验/审阅",
                 c("regression", "saved_d15r_v2p2_errors_and_review"), "匹配两模型 C0–C2", "PASS",
                 caveat = "只复用已保存响应")
      pkg_metric(task, "pk2c_token", "分词", "真实 tokenizer 回归",
                 c("regression", "tokenization"), "实际 tokenizer 数量可比", "PARTIAL",
                 caveat = "只有 mock 与投影验证")
    }
  }
  dataset <- "D15 R 包 PKG-D1"
  add("pkg_d1_posthoc", "验收变更", "记录验收标准是在看过结果后修改", "post_hoc",
      pk("DECISIONS_PKG.json"), read(pk("DECISIONS_PKG.json"))$post_hoc,
      caveat = "原逐 k 观测差 <1e-6、稳定性差 ≤0.1、逐代表 top25 逐位一致；新标准保留同 k/ARI=1，要求每簇稳定性 >0.5 且按簇合并后的25条顺序一致",
      method = "读决策记录")

  saved <- "DATA/pk2c/pipeline-release.rds"
  release <- function() read(saved)
  cluster_map <- function() {
    p <- merge(release()$clusters$members, tsv(members)[c("feature_id", "cluster_id")],
               by = "feature_id")
    stopifnot(nrow(p) == nrow(tsv(members)), nrow(p) == nrow(release()$clusters$members))
    p
  }
  align_clusters <- function(x) {
    p <- cluster_map()
    tab <- table(p$cluster_id.x, p$cluster_id.y)
    stopifnot(all(rowSums(tab > 0) == 1L), all(colSums(tab > 0) == 1L))
    as.integer(colnames(tab)[max.col(tab)])[match(x, as.integer(rownames(tab)))]
  }
  dataset <- "D15 R 包 PK2c 保存对象"
  add("pk2c_ari_recount", "聚类", "按 feature ID 对齐 R 与锁定簇", "ARI",
      c(saved, wk(members)), {
        p <- cluster_map()
        validation_ari(table(p$cluster_id.x, p$cluster_id.y))
      }, "ARI=1", eq(1), 1)
  for (i in 0:2) {
    add(paste0("pk2c_stability_c", i), "稳定性", "按历史簇号对齐子抽样中位数",
        paste0("C", i, " Jaccard 中位数"), c(saved, wk(c(members, stability))), {
          x <- release()$clusters$stability
          x$median_jaccard[align_clusters(x$cluster_id) == i]
        }, "PKG-D1 每簇>0.5", function(x) confirm(x > 0.5),
        tsv(stability)$median_jaccard[tsv(stability)$cluster_id == i], wk(stability),
        "历史基准差异；PKG-D1 是事后验收标准", "读保存统计量并对齐簇",
        reference_source = wk(stability))
  }
  add("pk2c_mfc_keys_recount", "MFC", "逐行比对保存的 R 与 Python MFC",
      "候选键及排名一致 / 全部候选", c(saved, wk(mfc)), {
        x <- release()$candidates
        y <- tsv(mfc)
        stopifnot(nrow(x) == nrow(y))
        fields <- c("feature_id", "HMDB_ID", "best_adduct", "rank_within_feature")
        good <- Reduce(`&`, lapply(fields, function(k) as.character(x[[k]]) == as.character(y[[k]])))
        frac(sum(good), nrow(y))
      }, "逐行候选键与排名完全相同", eq("3391/3391"), "3391/3391")
  add("pk2c_mfc_error_recount", "MFC", "比对保存的 R 与 Python MFC 分数",
      "最大相对误差", c(saved, wk(mfc)), {
        x <- release()$candidates$mfc_priority_score
        y <- tsv(mfc)$mfc_priority_score
        stopifnot(length(x) == length(y))
        max(abs(x - y) / pmax(abs(y), .Machine$double.xmin))
      }, "<1e-6", function(x) confirm(x < 1e-6))
  add("pk2c_representatives_recount", "代表", "簇对齐后比对代表 ID 与角色",
      "一致 / 全部代表", c(saved, wk(c(members, reps))), {
        x <- release()$representatives
        x$cluster_id <- align_clusters(x$cluster_id)
        y <- tsv(reps)
        p <- merge(x, y, by = c("cluster_id", "rank"))
        stopifnot(nrow(x) == nrow(y))
        frac(sum(p$feature_id.x == p$feature_id.y & p$role.x == p$role.y), nrow(y))
      }, "选点、排序和角色逐项一致", eq("15/15"), "15/15")
  views <- paste0(neighbors, "full_evidence/C", 0:2, ".json")
  add("pk2c_merged_recount", "近邻", "簇对齐后重检合并的历史近邻",
      "逐位一致 / 全部合并近邻", c(saved, wk(c(members, views))), {
        x <- release()$neighbors$merged
        x$cluster_id <- align_clusters(x$cluster_id)
        matches <- vapply(0:2, function(i) {
          y <- subset(x, cluster_id == i)
          ids <- y$HMDB_ID  # rank retains the winning representative rank; rows carry merged order.
          old <- vapply(j(views[i + 1])$semantic_neighbors, `[[`, "", "HMDB_ID")
          stopifnot(length(ids) == length(old))
          sum(ids == old)
        }, 0L)
        frac(sum(matches), nrow(x))
      }, "每簇25条 HMDB ID 按位相同", eq("75/75"), "75/75")
  add("pk2c_per_rep_neighbors_recount", "近邻", "对比逐代表的旧验收口径",
      "逐位一致 / 全部近邻", c(saved, wk(retrieve)), {
        x <- release()$neighbors$neighbors
        y <- tsv(retrieve)
        p <- merge(x, y, by = c("feature_id", "rank"))
        stopifnot(nrow(p) == nrow(x), nrow(x) == nrow(y))
        frac(sum(p$HMDB_ID.x == p$HMDB_ID.y), nrow(y))
      }, "旧标准逐代表 top25 按位相同", eq("375/375"), "375/375", wk(retrieve),
      "历史基准差异；PKG-D1 采用按簇合并的近邻")

  summary <- do.call(rbind, unname(rows))
  rownames(summary) <- NULL
  locked_table <- read(locked)
  old <- locked_table[locked_table$metric %in% summary$validation_id, ]
  matched <- match(old$metric, summary$validation_id)
  has_previous_number <- nzchar(old$prior_report)
  summary$reference_value[matched[has_previous_number]] <- old$value[has_previous_number]
  for (i in seq_len(nrow(old))) {
    id <- old$metric[i]
    references[[paste0("locked:", id)]] <- data.frame(validation_id = id,
      historical_value = old$value[i], historical_source = locked,
      recomputed_value = summary$value[matched[i]], source_file = summary$source_file[matched[i]],
      interpretation = summary$caveat[matched[i]])
  }
  refs <- unique(do.call(rbind, unname(references)))
  discrepancies <- refs[!mapply(compare, refs$recomputed_value, refs$historical_value), ]
  list(summary = summary, discrepancies = discrepancies)
}

if (sys.nframe() == 0L) {
  args <- commandArgs(trailingOnly = TRUE)
  stopifnot(length(args) == 1L, args[1] %in% c("--check", "--write"))
  roots <- validation_roots()
  output <- file.path(roots[["PKG"]], "inst/validation")
  result <- recompute_validation(roots)
  filenames <- c(summary = "validation_summary.tsv", discrepancies = "discrepancies.tsv")
  for (item in names(filenames)) {
    path <- file.path(output, filenames[[item]])
    if (args[1] == "--write") validation_write(result[[item]], path) else {
      expected <- validation_read_summary(path)
      stopifnot(isTRUE(all.equal(expected, result[[item]], check.attributes = FALSE)))
    }
  }
  cat("Validation rows:", nrow(result$summary), "; discrepancies:", nrow(result$discrepancies), "\n")
}
