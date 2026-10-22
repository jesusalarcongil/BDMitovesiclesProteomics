createMetaData = function(data, col_range, row_range, sampleNameColumns, metaColNames) {

  meta_data = data.frame(data[row_range, col_range, drop = FALSE])
  
  if(length(metaColNames) != length(col_range)) {
    stop("Number of provided column names (metaColNames) must match the selected column range.")
  }
  
  if (is.numeric(sampleNameColumns)) {
    rownames(meta_data) = meta_data[[sampleNameColumns]]
  } else {
    rownames(meta_data) = meta_data[, sampleNameColumns]
  }
  
  colnames(meta_data) = metaColNames
  
  return(meta_data)
}


relevelMetaData = function(meta_data, ref_levels) {
  for (col in names(ref_levels)) {
    if (col %in% colnames(meta_data)) {
      meta_data[[col]] = as.factor(meta_data[[col]])
      meta_data[[col]] = relevel(meta_data[[col]], ref = ref_levels[[col]])
    } else {
      warning(paste("Column", col, "not found in meta_data. Skipping..."))
    }
  }
  
  return(meta_data)
}


createExpDesign = function(meta_data, ref_levels) {

  selected_cols = names(ref_levels)
  
  missing_cols = setdiff(selected_cols, colnames(meta_data))
  if (length(missing_cols) > 0)
    warning("Missing columns in meta_data: ", paste(missing_cols, collapse = ", "))
  
  expDesign = data.frame(label = rownames(meta_data), stringsAsFactors = FALSE)
  
  for (col in selected_cols) {
    if (col %in% colnames(meta_data)) {
      vals = as.character(meta_data[[col]])
      expDesign[[col]] = relevel(factor(vals), ref = ref_levels[[col]])
    }
  }
  return(expDesign)
}


createProteinData = function(data, proteinNameColumns, col_range) {
  
  data_t = t(data)
  
  colnames(data_t) = data_t[proteinNameColumns, ]
  
  clean_data = data.frame(data_t[-col_range, ], stringsAsFactors = FALSE)
  
  return(clean_data)
}

# hybridProteomicsAnalysis: Executes the hybrid Proteomics pipeline
# (Continuous DEP + Fisher Presence/Absence MNAR)
#
# proteins_unique: Dataframe containing formatted and unique protein data.
# dataColumns: Indices or names of the columns containing intensity data.
# expDesign: Dataframe containing the experimental design.
# relevelVariables: List or dataframe with the variables used to relevel
#                   the SummarizedExperiment (SE).
# formulades: Formula for the linear model.
# contrastFormula: List containing the contrast formulas.
# contrastList: Character vector with the names of the contrasts.
# condition_cols: Character vector with the column names in expDesign that
#                 define the biological condition.
# alphaSig: FDR significance threshold.
# fcThreshold: Linear Fold Change threshold.
# min_valid_samples: Minimum number of presences in one group vs 0 in
#                    the other to consider a PA contrast significant.
#
# A list containing 4 elements: SE_Combined, SE_Pure, SE_PA_Only, and the
# final results table.

hybridProteomicsAnalysis =  function(proteins_unique,
                                     dataColumns,
                                     expDesign,
                                     relevelVariables,
                                     formulades, 
                                     contrastFormula, 
                                     contrastList,
                                     condition_cols,
                                     alphaSig = 0.05, 
                                     fcThreshold = 1.5,
                                     min_valid_samples = 4) {
  
  # STANDARD DEP ANALYSIS
  SEData = make_se(proteins_unique, dataColumns, expDesign, colnames(relevelVariables))
  noMissVal = rowSums(!is.na(assay(SEData)))
  SEData = SEData[noMissVal >= 4, ]
  SEDataImp = impute(SEData, fun = "MinProb", q = 0.01, MARGIN = 2L)
  
  vsn_fit_modelo_DEP = vsn::vsnMatrix(2 ^ assay(SEDataImp))
  
  SEDataImpNorm = SEDataImp
  assay(SEDataImpNorm) = vsn::predict(vsn_fit_modelo_DEP, 2 ^ assay(SEDataImp))
  
  SEDataImpNormAnalyzed = test_diff(SEDataImpNorm,
                                     type = "multiple",
                                     design_formula = formulades,
                                     contrastFormula = contrastFormula,
                                     contrastList = contrastList)
  
  sigProt_DEP = add_rejections(SEDataImpNormAnalyzed,
                                alpha = alphaSig,
                                lfc = log2(fcThreshold))
  
  # IDENTIFY MNAR PROTEINS (Presence/Absence)
  # Create a single biological condition column based on the chosen factors
  expDesign$condition = apply(expDesign[, condition_cols, drop = FALSE], 1,
                              paste, collapse = "_")
  unique_conditions = unique(expDesign$condition)
  
  # Logical and count matrices
  mnar_matrix = matrix(FALSE, nrow = nrow(proteins_unique),
                        ncol = length(unique_conditions))
  rownames(mnar_matrix) = proteins_unique$name
  colnames(mnar_matrix) = unique_conditions
  
  presence_matrix = matrix(0, nrow = nrow(proteins_unique),
                           ncol = length(unique_conditions))
  rownames(presence_matrix) = proteins_unique$name
  colnames(presence_matrix) = unique_conditions
  
  for(cond in unique_conditions) {
    samples_in_cond = expDesign$label[expDesign$condition == cond]
    col_indices = match(samples_in_cond, colnames(proteins_unique))
    
    mnar_matrix[, cond] =
      rowSums(is.na(proteins_unique[, col_indices])) == length(samples_in_cond)
    presence_matrix[, cond] = rowSums(!is.na(proteins_unique[, col_indices]))
  }
  
  # Keep only proteins missing entirely in at least one group
  mnar_any_group = rowSums(mnar_matrix) > 0
  PAProteins = proteins_unique[mnar_any_group, ]
  
  # MATRIX PREPARATION FOR FOLD CHANGE (LOD Imputation + VSN)
  PA_binary = !is.na(PAProteins[, dataColumns]) & PAProteins[, dataColumns] > 0
  PA_binary = as.matrix(PA_binary) * 1
  
  # Calculate global minimum (Limit of Detection)
  global_min = min(proteins_unique[proteins_unique > 0], na.rm = TRUE)
  raw_PA = PAProteins
  
  # Impute NAs and zeros with the global minimum for qualitative proteins
  raw_PA[is.na(PAProteins) | PAProteins <= 0] = as.numeric(global_min)
  
  SEData_PA = make_se(raw_PA, dataColumns, expDesign, colnames(relevelVariables))
  
  # Apply VSN: using the previously generated 'vsn_fit_modelo_DEP' model
  vsn_PA_matrix = vsn::predict(vsn_fit_modelo_DEP, 2 ^ assay(SEData_PA))
  vsn_PA_data = assay(vsn_PA_matrix)
  
  # PAIRWISE FISHER'S EXACT TEST FOR THE NA SUBSET
  # Save results in a list indexed by the real protein name
  PA_results = list()
  PA_names = rownames(assay(vsn_PA_matrix))
  
  for (i in seq_len(nrow(PA_binary))) {
    prot_name = PA_names[i]
    PA_results[[prot_name]] = list()
    
    # Iterate contrast by contrast
    for(contrast_name in contrastList) {
      
      # Extract the two groups involved
      groups = strsplit(contrast_name, "_vs_")[[1]]
      g1 = groups[1] 
      g2 = groups[2] 
      
      # MNAR check
      is_mnar_g1 = mnar_matrix[prot_name, g1]
      is_mnar_g2 = mnar_matrix[prot_name, g2]
      
      # If it's not completely missing in either, it's quantitative (DEP)
      if(!(is_mnar_g1 | is_mnar_g2)) {
        PA_results[[prot_name]][[paste0(contrast_name, "_p.val")]] = NA
        PA_results[[prot_name]][[paste0(contrast_name, "_diff")]]  = NA
        next
      }
      # If completely missing in both (0 vs 0), there are no differences
      if(is_mnar_g1 & is_mnar_g2) {
        PA_results[[prot_name]][[paste0(contrast_name, "_p.val")]] = 1
        PA_results[[prot_name]][[paste0(contrast_name, "_diff")]]  = 0
        next
      }
      
      # Pairwise Fisher's Exact Test
      idx_g1 = which(expDesign$condition == g1)
      idx_g2 = which(expDesign$condition == g2)
      
      pres_g1 = sum(PA_binary[i, idx_g1])
      pres_g2 = sum(PA_binary[i, idx_g2])
      
      n_g1 = length(idx_g1)
      n_g2 = length(idx_g2)
      
      # Build 2x2 contingency matrix
      m = matrix(c(pres_g1, n_g1 - pres_g1, pres_g2, n_g2 - pres_g2), nrow = 2)
      pval = fisher.test(m)$p.value
      
      # Fold Change Calculation (Using imputed VSN matrix)
      vals_g1 = vsn_PA_data[i, idx_g1]
      vals_g2 = vsn_PA_data[i, idx_g2]
      diff_VSN = mean(data.matrix(vals_g1)) - mean(data.matrix(vals_g2))
      
      PA_results[[prot_name]][[paste0(contrast_name, "_p.val")]] = pval
      PA_results[[prot_name]][[paste0(contrast_name, "_diff")]]  = diff_VSN
    }
  }
  
  # TABLE RECONSTRUCTION AND FINAL FILTERING
  firth_df = purrr::imap_dfr(PA_results, function(prot_res, prot_name) {
    tibble(
      protein = prot_name,
      metric = names(prot_res),
      value = unlist(prot_res, use.names = FALSE) 
    )
  })
  
  firth_df = firth_df %>%
    mutate(
      metric_type = case_when(
        grepl("_p\\.val$", metric) ~ "p.val",
        grepl("_diff$", metric) ~ "diff"
      ),
      contrast = gsub("_p\\.val$|_diff$", "", metric)
    ) %>%
    dplyr::select(-metric) %>%
    pivot_wider(names_from = metric_type, values_from = value) %>%
    separate(contrast, into = c("group1", "group2"), sep = "_vs_",
             remove = FALSE) %>%
    filter(!is.na(p.val))
  
  mnar_long = as.data.frame(mnar_matrix) %>%
    tibble::rownames_to_column("protein") %>%
    pivot_longer(-protein, names_to = "group", values_to = "is_mnar")
  
  presence_long = as.data.frame(presence_matrix) %>%
    tibble::rownames_to_column("protein") %>%
    pivot_longer(-protein, names_to = "group", values_to = "presence_count")
  
  firth_results_final = firth_df %>%
    left_join(mnar_long, by = c("protein" = "protein",
                                "group1" = "group")) %>%
    rename(is_mnar_g1 = is_mnar) %>%
    left_join(mnar_long, by = c("protein" = "protein",
                                "group2" = "group")) %>%
    rename(is_mnar_g2 = is_mnar) %>%
    left_join(presence_long, by = c("protein" = "protein",
                                    "group1" = "group")) %>%
    rename(count_g1 = presence_count) %>%
    left_join(presence_long, by = c("protein" = "protein",
                                    "group2" = "group")) %>%
    rename(count_g2 = presence_count) %>%
    mutate(
      valid_for_PA = (is_mnar_g1 | is_mnar_g2) & !(is_mnar_g1 & is_mnar_g2),
    ) %>%
    filter(valid_for_PA == TRUE) %>%
    
    # Calculate N per group and establish significance according to presence
    rowwise() %>%
    mutate(
      n_g1 = sum(expDesign$condition == group1),
      n_g2 = sum(expDesign$condition == group2),
      significant = ((count_g1 >= min_valid_samples & count_g2 == 0) | 
                       (count_g2 >= min_valid_samples & count_g1 == 0)) &
        abs(diff) > log2(fcThreshold)
    ) %>%
    ungroup() %>%
    
    # Maintain classic FDR adjustment
    group_by(contrast) %>%
    mutate(p.adj = p.adjust(p.val, method = "BH")) %>%
    ungroup() %>%
    
    dplyr::select(protein, contrast, group1, group2, diff, p.val, p.adj, significant)
  
  # CREATE PURE DEP OBJECT
  firth_wide = firth_results_final %>%
    pivot_wider(
      id_cols = protein,
      names_from = contrast,
      values_from = c(diff, p.val, p.adj, significant),
      names_glue = "{contrast}_{.value}"
    ) %>%
    rename(name = protein) %>%
    as.data.frame()
  
  rd_dep_pure = as.data.frame(rowData(sigProt_DEP), row.names = NULL)
  sigProt_DEP_Pure = sigProt_DEP
  
  for (i in 1:nrow(firth_results_final)) {
    prot = firth_results_final$protein[i]
    contrast_name = firth_results_final$contrast[i]
    
    pval_col = paste0(contrast_name, "_p.val")
    diff_col = paste0(contrast_name, "_diff")
    padj_col = paste0(contrast_name, "_p.adj")
    sig_col  = paste0(contrast_name, "_significant")
    ci_l_col = paste0(contrast_name, "_CI.L")
    ci_r_col = paste0(contrast_name, "_CI.R")
    
    if (prot %in% rd_dep_pure$name) {
      if (pval_col %in% colnames(rd_dep_pure)) rd_dep_pure[rd_dep_pure$name == prot, pval_col] = NA
      if (diff_col %in% colnames(rd_dep_pure)) rd_dep_pure[rd_dep_pure$name == prot, diff_col] = NA
      if (padj_col %in% colnames(rd_dep_pure)) rd_dep_pure[rd_dep_pure$name == prot, padj_col] = NA
      if (sig_col %in% colnames(rd_dep_pure))  rd_dep_pure[rd_dep_pure$name == prot, sig_col]  = FALSE
      if (ci_l_col %in% colnames(rd_dep_pure)) rd_dep_pure[rd_dep_pure$name == prot, ci_l_col] = NA
      if (ci_r_col %in% colnames(rd_dep_pure)) rd_dep_pure[rd_dep_pure$name == prot, ci_r_col] = NA
    }
  }
  
  # FDR readjustment
  for (contrast_name in contrastList) {
    pval_col = paste0(contrast_name, "_p.val")
    padj_col = paste0(contrast_name, "_p.adj")
    diff_col = paste0(contrast_name, "_diff")
    sig_col  = paste0(contrast_name, "_significant")
    
    sig_type_col = paste0(contrast_name, "_sig_type")
    rd_dep_pure[[sig_type_col]] = "Not Significant"
    
    if (pval_col %in% colnames(rd_dep_pure)) {
      valid_idx = which(!is.na(rd_dep_pure[[pval_col]]))
      pvals_valid = rd_dep_pure[[pval_col]][valid_idx]
      
      if (length(pvals_valid) > 0) {
        new_padj = tryCatch({
          fdr_res = fdrtool::fdrtool(pvals_valid, statistic = "pvalue",
                                     plot = FALSE, verbose = FALSE)
          fdr_res$qval
        }, error = function(e) {
          p.adjust(pvals_valid, method = "BH")
        })
        
        rd_dep_pure[[padj_col]][valid_idx] = new_padj
        
        rd_dep_pure[[sig_col]][valid_idx] = rd_dep_pure[[padj_col]][valid_idx] < alphaSig & 
          abs(rd_dep_pure[[diff_col]][valid_idx]) > log2(fcThreshold)
        
        rd_dep_pure[[sig_type_col]][valid_idx] = ifelse(rd_dep_pure[[sig_col]][valid_idx], "DEP", "Not Significant")
      }
    }
  }
  
  sig_cols = grep("_significant$", colnames(rd_dep_pure), value = TRUE)
  rd_dep_pure$significant = rowSums(rd_dep_pure[, sig_cols, drop = FALSE], na.rm = TRUE) > 0
  rowData(sigProt_DEP_Pure) = DataFrame(rd_dep_pure)
  resultsProt_DEP = get_results(sigProt_DEP_Pure)
  
  # CREATE EXCLUSIVE PRESENCE/ABSENCE DEP OBJECT
  idx_pa = match(firth_wide$name, PAProteins$name)
  rd_pa_base = data.frame(
    name              = firth_wide$name,
    ID                = PAProteins$ID[idx_pa],
    Primary.Gene.name = PAProteins$Primary.Gene.name[idx_pa],
    imputed           = TRUE 
  )
  rd_pa_full = merge(rd_pa_base, firth_wide, by = "name")
  rownames(rd_pa_full) = rd_pa_full$name
  sig_cols_pa = grep("_significant$", colnames(rd_pa_full), value = TRUE)
  rd_pa_full$significant = rowSums(rd_pa_full[, sig_cols_pa, drop = FALSE], na.rm = TRUE) > 0
  assay_pa_only = vsn_PA_data[firth_wide$name, , drop = FALSE]
  SE_PA_Only = SummarizedExperiment(
    assays  = list(counts = assay_pa_only),
    rowData = DataFrame(rd_pa_full),
    colData = colData(sigProt_DEP)
  )
  resultsProt_PA = get_results(SE_PA_Only)
  
  # CREATE COMBINED DEP OBJECT
  rd_dep_combined = rd_dep_pure
  
  for (prot in firth_wide$name) {
    for (contrast_name in contrastList) {
      pval_col = paste0(contrast_name, "_p.val")
      diff_col = paste0(contrast_name, "_diff")
      padj_col = paste0(contrast_name, "_p.adj")
      sig_col  = paste0(contrast_name, "_significant")
      sig_type_col = paste0(contrast_name, "_sig_type")
      
      if (pval_col %in% colnames(firth_wide) &&
          !is.na(firth_wide[firth_wide$name == prot, pval_col])) {
        
        rd_dep_combined[rd_dep_combined$name == prot, pval_col] = firth_wide[firth_wide$name == prot, pval_col]
        rd_dep_combined[rd_dep_combined$name == prot, diff_col] = firth_wide[firth_wide$name == prot, diff_col]
        rd_dep_combined[rd_dep_combined$name == prot, padj_col] = firth_wide[firth_wide$name == prot, padj_col] 
        rd_dep_combined[rd_dep_combined$name == prot, sig_col]  = firth_wide[firth_wide$name == prot, sig_col]
        
        is_sig_PA = firth_wide[firth_wide$name == prot, sig_col]
        rd_dep_combined[rd_dep_combined$name == prot, sig_type_col] = ifelse(is_sig_PA, "PA", "Not Significant")
      }
    }
  }
  
  sig_cols_comb = grep("_significant$", colnames(rd_dep_combined), value = TRUE)
  rd_dep_combined$significant = rowSums(rd_dep_combined[, sig_cols_comb, drop = FALSE], na.rm = TRUE) > 0
  
  sigProt_DEP_Combined = sigProt_DEP
  rowData(sigProt_DEP_Combined) = DataFrame(rd_dep_combined)
  resultsProt_Combined = get_results(sigProt_DEP_Combined)
  
  return(list(
    SE_Combined = sigProt_DEP_Combined,
    SE_DEP = sigProt_DEP_Pure,
    SE_PA = SE_PA_Only,
    results_DEP = resultsProt_DEP,
    Results_PA = resultsProt_PA,
    Results_Combined = resultsProt_Combined
  ))
}


# hybridMixedProteomicsAnalysis: Executes the hybrid Proteomics pipeline
# in mixed data (Continuous DEP + Fisher Presence/Absence MNAR)
#
# proteins_unique_1: Dataframe containing formatted and unique protein data.
# proteins_unique_2: Dataframe containing formatted and unique protein data to combine
# dataColumns_1: Indices/names of columns containing intensity data for dataset 1.
# dataColumns_2: Indices/names of columns containing intensity data for dataset 2.
# expDesign_1: Dataframe containing the experimental design 1.
# expDesign_2: Dataframe containing the experimental design 2.
# relevelVariables_1: Variables used to relevel dataset 1.
# relevelVariables_2: Variables used to relevel dataset 2.
# formulades: Formula for the linear model.
# contrastFormula: List containing the contrast formulas.
# contrastList: Character vector with the names of the contrasts.
# condition_cols: Character vector with the column names in expDesign that define the biological condition.
# alphaSig: FDR significance threshold.
# fcThreshold: Linear Fold Change threshold.
# min_valid_samples: Minimum number of presences in one group vs 0 in the other.
#
# Returns: A list containing SE_Combined, SE_DEP, SE_PA, and their tables.

# hybridMixedProteomicsAnalysis: Executes the hybrid Proteomics pipeline
# in mixed data (Continuous DEP + Fisher Presence/Absence MNAR)

hybridMixedProteomicsAnalysis =  function(proteins_unique_1,
                                          proteins_unique_2,
                                          dataColumns_1,
                                          dataColumns_2,
                                          expDesign_1,
                                          expDesign_2,
                                          relevelVariables_1,
                                          relevelVariables_2,
                                          formulades, 
                                          contrastFormula, 
                                          contrastList,
                                          condition_cols,
                                          alphaSig = 0.05, 
                                          fcThreshold = 1.5,
                                          min_valid_samples = 6) {
  
  # STANDARD DEP ANALYSIS 
  SEData_1 = make_se(proteins_unique_1, dataColumns_1, expDesign_1, colnames(relevelVariables_1))
  SEData_2 = make_se(proteins_unique_2, dataColumns_2, expDesign_2, colnames(relevelVariables_2))
  
  noMissVal_1 = rowSums(!is.na(assay(SEData_1)))
  SEData_1 = SEData_1[noMissVal_1 >= 6, ]
  noMissVal_2 = rowSums(!is.na(assay(SEData_2)))
  SEData_2 = SEData_2[noMissVal_2 >= 4, ]
  
  SEDataImp_1 = impute(SEData_1, fun = "MinProb", q = 0.01, MARGIN = 2L)
  SEDataImp_2 = impute(SEData_2, fun = "MinProb", q = 0.01, MARGIN = 2L)
  
  SEDataImp = combineSEData(SEDataImp_2, SEDataImp_1, relevelVariables_1, imputeOpt = TRUE)
  SEDataRaw = combineSEData(SEData_2, SEData_1, relevelVariables_1, imputeOpt = FALSE)
  
  vsn_fit_modelo_DEP = vsn::vsnMatrix(2 ^ assay(SEDataImp))
  
  SEDataImpNorm = SEDataImp
  assay(SEDataImpNorm) = vsn::predict(vsn_fit_modelo_DEP, 2 ^ assay(SEDataImp))
  
  SEDataImpNormAnalyzed = test_diff(SEDataImpNorm,
                                    type = "multiple",
                                    design_formula = formulades,
                                    contrastFormula = contrastFormula,
                                    contrastList = contrastList)
  
  sigProt_DEP = add_rejections(SEDataImpNormAnalyzed,
                               alpha = alphaSig,
                               lfc = log2(fcThreshold))
  
  # IDENTIFY MNAR PROTEINS
  expDesign = rbind(expDesign_1, expDesign_2)
  expDesign$condition = apply(expDesign[, condition_cols, drop = FALSE], 1, paste, collapse = "_")
  
  all_groups = unique(unlist(strsplit(contrastList, "_vs_")))
  design_matrix = model.matrix(formulades, data = expDesign)
  
  assay_raw = assay(SEDataRaw)
  prot_names = rownames(assay_raw)
  
  is_missing_mat = is.na(assay_raw) | assay_raw <= 0
  
  mnar_matrix = matrix(FALSE, nrow = nrow(assay_raw), ncol = length(all_groups))
  rownames(mnar_matrix) = prot_names
  colnames(mnar_matrix) = all_groups
  
  presence_matrix = matrix(0, nrow = nrow(assay_raw), ncol = length(all_groups))
  rownames(presence_matrix) = prot_names
  colnames(presence_matrix) = all_groups
  
  group_indices = list()
  
  for(g in all_groups) {
    if (g %in% expDesign$condition) {
      idx = which(expDesign$condition == g)
    } else if (g %in% colnames(design_matrix)) {
      idx = which(design_matrix[, g] == 1)
    } else {
      stop(paste0("Error: Group '", g, "' not found in conditions or Limma design matrix."))
    }
    
    group_indices[[g]] = idx
    mnar_matrix[, g] = rowSums(is_missing_mat[, idx, drop = FALSE]) == length(idx)
    presence_matrix[, g] = rowSums(!is_missing_mat[, idx, drop = FALSE])
  }
  
  mnar_any_group = rowSums(mnar_matrix) > 0
  PA_prot_names = prot_names[mnar_any_group]
  
  # MATRIX PREPARATION FOR FOLD CHANGE
  assay_raw_PA = assay_raw[mnar_any_group, , drop = FALSE]
  
  PA_binary = !is_missing_mat[mnar_any_group, , drop = FALSE]
  PA_binary = as.matrix(PA_binary) * 1
  
  global_min = min(assay_raw[assay_raw > 0], na.rm = TRUE)
  assay_raw_PA[is.na(assay_raw_PA) | assay_raw_PA <= 0] = as.numeric(global_min)
  
  SEData_PA = SEDataRaw[mnar_any_group, ]
  assay(SEData_PA) = assay_raw_PA
  
  vsn_PA_matrix = vsn::predict(vsn_fit_modelo_DEP, 2 ^ assay(SEData_PA))
  vsn_PA_data = assay(vsn_PA_matrix)
  
  # PAIRWISE FISHER'S EXACT TEST FOR THE NA SUBSET
  PA_results = list()
  
  for (i in seq_len(nrow(PA_binary))) {
    prot_name = PA_prot_names[i]
    PA_results[[prot_name]] = list()
    
    for(contrast_name in contrastList) {
      groups = strsplit(contrast_name, "_vs_")[[1]]
      g1 = groups[1] 
      g2 = groups[2] 
      
      is_mnar_g1 = mnar_matrix[prot_name, g1]
      is_mnar_g2 = mnar_matrix[prot_name, g2]
      
      if(!(is_mnar_g1 | is_mnar_g2)) {
        PA_results[[prot_name]][[paste0(contrast_name, "_p.val")]] = NA
        PA_results[[prot_name]][[paste0(contrast_name, "_diff")]]  = NA
        next
      }
      if(is_mnar_g1 & is_mnar_g2) {
        PA_results[[prot_name]][[paste0(contrast_name, "_p.val")]] = 1
        PA_results[[prot_name]][[paste0(contrast_name, "_diff")]]  = 0
        next
      }
      
      idx_g1 = group_indices[[g1]]
      idx_g2 = group_indices[[g2]]
      
      pres_g1 = sum(PA_binary[i, idx_g1])
      pres_g2 = sum(PA_binary[i, idx_g2])
      
      n_g1 = length(idx_g1)
      n_g2 = length(idx_g2)
      
      m = matrix(c(pres_g1, n_g1 - pres_g1, pres_g2, n_g2 - pres_g2), nrow = 2)
      pval = fisher.test(m)$p.value
      
      vals_g1 = vsn_PA_data[i, idx_g1]
      vals_g2 = vsn_PA_data[i, idx_g2]
      diff_VSN = mean(data.matrix(vals_g1)) - mean(data.matrix(vals_g2))
      
      PA_results[[prot_name]][[paste0(contrast_name, "_p.val")]] = pval
      PA_results[[prot_name]][[paste0(contrast_name, "_diff")]]  = diff_VSN
    }
  }
  
  # TABLE RECONSTRUCTION AND FINAL FILTERING
  firth_df = purrr::imap_dfr(PA_results, function(prot_res, prot_name) {
    tibble(
      protein = prot_name,
      metric = names(prot_res),
      value = unlist(prot_res, use.names = FALSE) 
    )
  })
  
  firth_df = firth_df %>%
    mutate(
      metric_type = case_when(
        grepl("_p\\.val$", metric) ~ "p.val",
        grepl("_diff$", metric) ~ "diff"
      ),
      contrast = gsub("_p\\.val$|_diff$", "", metric)
    ) %>%
    dplyr::select(-metric) %>%
    pivot_wider(names_from = metric_type, values_from = value) %>%
    separate(contrast, into = c("group1", "group2"), sep = "_vs_", remove = FALSE) %>%
    filter(!is.na(p.val))
  
  mnar_long = as.data.frame(mnar_matrix) %>%
    tibble::rownames_to_column("protein") %>%
    pivot_longer(-protein, names_to = "group", values_to = "is_mnar")
  
  presence_long = as.data.frame(presence_matrix) %>%
    tibble::rownames_to_column("protein") %>%
    pivot_longer(-protein, names_to = "group", values_to = "presence_count")
  
  group_sizes = data.frame(
    group = names(group_indices),
    n_group = sapply(group_indices, length),
    stringsAsFactors = FALSE
  )
  
  firth_results_final = firth_df %>%
    left_join(mnar_long, by = c("protein" = "protein", "group1" = "group")) %>%
    rename(is_mnar_g1 = is_mnar) %>%
    left_join(mnar_long, by = c("protein" = "protein", "group2" = "group")) %>%
    rename(is_mnar_g2 = is_mnar) %>%
    left_join(presence_long, by = c("protein" = "protein", "group1" = "group")) %>%
    rename(count_g1 = presence_count) %>%
    left_join(presence_long, by = c("protein" = "protein", "group2" = "group")) %>%
    rename(count_g2 = presence_count) %>%
    left_join(group_sizes, by = c("group1" = "group")) %>%
    rename(n_g1 = n_group) %>%
    left_join(group_sizes, by = c("group2" = "group")) %>%
    rename(n_g2 = n_group) %>%
    mutate(
      valid_for_PA = (is_mnar_g1 | is_mnar_g2) & !(is_mnar_g1 & is_mnar_g2)
    ) %>%
    filter(valid_for_PA == TRUE) %>%
    rowwise() %>%
    mutate(
      req_g1 = round(n_g1 * 0.8),
      req_g2 = round(n_g2 * 0.8),
      
      significant = ((count_g1 >= req_g1 & count_g2 == 0) | 
                       (count_g2 >= req_g2 & count_g1 == 0)) &
        abs(diff) > log2(fcThreshold)
    ) %>%
    ungroup() %>%
    group_by(contrast) %>%
    mutate(p.adj = p.adjust(p.val, method = "BH")) %>%
    ungroup() %>%
    dplyr::select(protein, contrast, group1, group2, diff, p.val, p.adj, significant)
  
  # CREATE PURE DEP OBJECT
  firth_wide = firth_results_final %>%
    pivot_wider(
      id_cols = protein,
      names_from = contrast,
      values_from = c(diff, p.val, p.adj, significant),
      names_glue = "{contrast}_{.value}"
    ) %>%
    rename(name = protein) %>%
    as.data.frame()
  
  rd_dep_pure = as.data.frame(rowData(sigProt_DEP), row.names = NULL)
  sigProt_DEP_Pure = sigProt_DEP
  
  idx_match_6 = match(firth_results_final$protein, rd_dep_pure$name)
  for (i in 1:nrow(firth_results_final)) {
    if (!is.na(idx_match_6[i])) {
      c_name = firth_results_final$contrast[i]
      rd_dep_pure[idx_match_6[i], paste0(c_name, "_p.val")] = NA
      rd_dep_pure[idx_match_6[i], paste0(c_name, "_diff")] = NA
      rd_dep_pure[idx_match_6[i], paste0(c_name, "_p.adj")] = NA
      rd_dep_pure[idx_match_6[i], paste0(c_name, "_significant")] = FALSE
    }
  }
  
  for (contrast_name in contrastList) {
    pval_col = paste0(contrast_name, "_p.val")
    padj_col = paste0(contrast_name, "_p.adj")
    diff_col = paste0(contrast_name, "_diff")
    sig_col  = paste0(contrast_name, "_significant")
    
    sig_type_col = paste0(contrast_name, "_sig_type")
    rd_dep_pure[[sig_type_col]] = "Not Significant"
    
    if (pval_col %in% colnames(rd_dep_pure)) {
      valid_idx = which(!is.na(rd_dep_pure[[pval_col]]))
      pvals_valid = rd_dep_pure[[pval_col]][valid_idx]
      
      if (length(pvals_valid) > 0) {
        new_padj = tryCatch({
          fdr_res = fdrtool::fdrtool(pvals_valid, statistic = "pvalue",
                                     plot = FALSE, verbose = FALSE)
          fdr_res$qval
        }, error = function(e) {
          p.adjust(pvals_valid, method = "BH")
        })
        
        rd_dep_pure[[padj_col]][valid_idx] = new_padj
        rd_dep_pure[[sig_col]][valid_idx] = rd_dep_pure[[padj_col]][valid_idx] < alphaSig & 
          abs(rd_dep_pure[[diff_col]][valid_idx]) > log2(fcThreshold)
        
        rd_dep_pure[[sig_type_col]][valid_idx] = ifelse(rd_dep_pure[[sig_col]][valid_idx], "DEP", "Not Significant")
      }
    }
  }
  
  sig_cols = grep("_significant$", colnames(rd_dep_pure), value = TRUE)
  rd_dep_pure$significant = rowSums(rd_dep_pure[, sig_cols, drop = FALSE], na.rm = TRUE) > 0
  rowData(sigProt_DEP_Pure) = DataFrame(rd_dep_pure)
  resultsProt_DEP = get_results(sigProt_DEP_Pure)
  
  # CREATE EXCLUSIVE PRESENCE/ABSENCE DEP OBJECT
  idx_pa = match(firth_wide$name, rd_dep_pure$name)
  rd_pa_base = data.frame(
    name    = firth_wide$name,
    ID      = rd_dep_pure$ID[idx_pa],
    imputed = TRUE 
  )
  if ("Primary.Gene.name" %in% colnames(rd_dep_pure)) {
    rd_pa_base$Primary.Gene.name = rd_dep_pure$Primary.Gene.name[idx_pa]
  }
  
  rd_pa_full = merge(rd_pa_base, firth_wide, by = "name")
  rownames(rd_pa_full) = rd_pa_full$name
  sig_cols_pa = grep("_significant$", colnames(rd_pa_full), value = TRUE)
  rd_pa_full$significant = rowSums(rd_pa_full[, sig_cols_pa, drop = FALSE], na.rm = TRUE) > 0
  assay_pa_only = vsn_PA_data[firth_wide$name, , drop = FALSE]
  SE_PA_Only = SummarizedExperiment(
    assays  = list(counts = assay_pa_only),
    rowData = DataFrame(rd_pa_full),
    colData = colData(sigProt_DEP)
  )
  resultsProt_PA = get_results(SE_PA_Only)
  
  # CREATE COMBINED DEP OBJECT
  rd_dep_combined = rd_dep_pure
  
  idx_match_8 = match(firth_wide$name, rd_dep_combined$name)
  valid_rows = !is.na(idx_match_8)
  
  for (contrast_name in contrastList) {
    pval_col = paste0(contrast_name, "_p.val")
    diff_col = paste0(contrast_name, "_diff")
    padj_col = paste0(contrast_name, "_p.adj")
    sig_col  = paste0(contrast_name, "_significant")
    sig_type_col = paste0(contrast_name, "_sig_type")
    
    if (pval_col %in% colnames(firth_wide)) {
      has_val = !is.na(firth_wide[[pval_col]])
      update_idx = valid_rows & has_val
      
      target_idx = idx_match_8[update_idx]
      source_idx = update_idx
      
      rd_dep_combined[target_idx, pval_col] = firth_wide[source_idx, pval_col]
      rd_dep_combined[target_idx, diff_col] = firth_wide[source_idx, diff_col]
      rd_dep_combined[target_idx, padj_col] = firth_wide[source_idx, padj_col] 
      rd_dep_combined[target_idx, sig_col]  = firth_wide[source_idx, sig_col]
      
      is_sig_PA = firth_wide[source_idx, sig_col]
      rd_dep_combined[target_idx, sig_type_col] = ifelse(is_sig_PA, "PA", "Not Significant")
    }
  }
  
  sig_cols_comb = grep("_significant$", colnames(rd_dep_combined), value = TRUE)
  rd_dep_combined$significant = rowSums(rd_dep_combined[, sig_cols_comb, drop = FALSE], na.rm = TRUE) > 0
  
  sigProt_DEP_Combined = sigProt_DEP
  rowData(sigProt_DEP_Combined) = DataFrame(rd_dep_combined)
  resultsProt_Combined = get_results(sigProt_DEP_Combined)
  
  return(list(
    SE_Combined = sigProt_DEP_Combined,
    SE_DEP = sigProt_DEP_Pure,
    SE_PA = SE_PA_Only,
    Results_DEP = resultsProt_DEP,
    Results_PA = resultsProt_PA,
    Results_Combined = resultsProt_Combined
  ))
}

# Combine Two SummarizedExperiment Objects by Matching Proteins and Relevel Metadata
#
# combineSEData() combines two SummarizedExperiment objects by doing the following:
# 1. Extracts protein IDs (assumed to be the row names) from the main data object.
# 2. Extracts the assay data from the secondary object.
# 3. Creates an empty matrix with the protein IDs as rows and copies over the assay data
#    for proteins found in common (ignoring case).
# 4. Creates a new SummarizedExperiment from the subset of the secondary data using its colData.
# 5. Combines the main and the new secondary SummarizedExperiment via \code{cbind}.
# 6. Relevels specified metadata fields within the combined \code{colData} using the provided
#    named list \code{relevelVar}.
#
# -mainSE: A SummarizedExperiment object containing the main proteomics data.
#   It is assumed that the protein identifiers are provided as row names.
# -otherSE: A SummarizedExperiment object containing additional data.
#   Its assay is used for extracting the corresponding protein measurements.
# -relevelVar: A named list where the names correspond to metadata fields in the colData of the
#   combined experiment, and the values are the desired reference levels.
#
# -Output: A combined SummarizedExperiment object after aligning the proteins, merging the assays,
#   and releveling the specified metadata factors.

combineSEDataOld <- function(mainSE, otherSE, relevelVar) {
  # Extract protein IDs from the main dataset (assumed to be row names)
  proteinIDs <- rownames(mainSE)
  
  # Extract the assay data from the other dataset
  otherAssayData <- assay(otherSE)
  
  # Create an empty matrix for the subset with rows equal to proteinIDs and matching columns to otherAssayData
  matSubset <- matrix(NA, nrow = length(proteinIDs), ncol = ncol(otherAssayData))
  rownames(matSubset) <- proteinIDs
  colnames(matSubset) <- colnames(otherAssayData)
  
  # Find common proteins (ignoring case) between the main data and the other assay data
  commonIndices <- match(tolower(proteinIDs), tolower(rownames(otherAssayData)))
  keepProteins <- !is.na(commonIndices)
  
  # Fill the subset matrix with the corresponding data from otherAssayData for proteins found in common
  matSubset[keepProteins, ] <- otherAssayData[commonIndices[keepProteins], ]
  
  # Create a SummarizedExperiment for the other subset using the colData from otherSE
  otherPartSE <- SummarizedExperiment(assays = list(matSubset), colData = colData(otherSE))
  
  # Combine the main dataset with the otherPartSE via cbind
  combinedSE <- cbind(mainSE, otherPartSE)
  
  # Relevel the specified metadata in the combined colData according to relevelVar
  for (col in names(relevelVar)) {
    if (col %in% names(combinedSE@colData@listData)) {
      combinedSE@colData@listData[[col]] <- relevel(combinedSE@colData@listData[[col]], ref = relevelVar[[col]])
    } else {
      warning(paste("Column", col, "not found in colData. Skipping..."))
    }
  }
  
  return(combinedSE)
}

combineSEData <- function(mainSE, otherSE, relevelVar, imputeOpt) {
  
  allProteins <- union(rownames(mainSE), rownames(otherSE))
  
  mat1 <- as.matrix(assay(mainSE))
  mat2 <- as.matrix(assay(otherSE))
  
  n_total_samples <- ncol(mat1) + ncol(mat2)
  combined_mat <- matrix(NA, nrow = length(allProteins), ncol = n_total_samples)
  rownames(combined_mat) <- allProteins
  colnames(combined_mat) <- c(colnames(mat1), colnames(mat2))
  
  combined_mat[rownames(mat1), colnames(mat1)] <- mat1
  combined_mat[rownames(mat2), colnames(mat2)] <- mat2
  
  cd1 <- as.data.frame(colData(mainSE))
  cd2 <- as.data.frame(colData(otherSE))
  
  combined_df <- dplyr::bind_rows(cd1, cd2)
  combined_colData <- S4Vectors::DataFrame(combined_df)
  
  rd1 <- as.data.frame(rowData(mainSE))
  rd1$row_id <- rownames(mainSE)
  
  rd2 <- as.data.frame(rowData(otherSE))
  rd2$row_id <- rownames(otherSE)
  
  combined_rd_df <- dplyr::bind_rows(rd1, rd2) %>%
    dplyr::group_by(row_id) %>%
    dplyr::summarise(
      
      dplyr::across(dplyr::any_of("imputed"), ~ any(.x, na.rm = TRUE)),
      dplyr::across(dplyr::any_of("num_NAs"), ~ sum(.x, na.rm = TRUE)),
      
      dplyr::across(!dplyr::any_of(c("imputed", "num_NAs")), ~ {
        valores_unicos <- unique(.x[!is.na(.x)])
        
        if (length(valores_unicos) > 1) {
          stop(sprintf(
            "Error! Protein '%s'has contradictory data '%s'. Values: %s",
            dplyr::cur_group()$row_id,
            dplyr::cur_column(),
            paste(valores_unicos, collapse = " vs ")
          ))
        }
        
        if (length(valores_unicos) == 1) {
          valores_unicos
        } else {
          .x[1] 
        }
      }),
      .groups = "drop"
    ) %>%
    as.data.frame()
  
  rownames(combined_rd_df) <- combined_rd_df$row_id
  combined_rd_df <- combined_rd_df[allProteins, ]
  combined_rd_df$row_id <- NULL 
  
  combined_rowData <- S4Vectors::DataFrame(combined_rd_df)

  combinedSE <- SummarizedExperiment(
    assays  = list(counts = combined_mat), 
    colData = combined_colData,
    rowData = combined_rowData
  )
  
  if (imputeOpt == TRUE){
    combinedSE = impute(combinedSE, fun = "MinProb", q = 0.01, MARGIN = 2L)
  }
  
  for (col in names(relevelVar)) {
    if (col %in% colnames(colData(combinedSE))) {
      combinedSE[[col]] <- relevel(factor(combinedSE[[col]]), ref = relevelVar[[col]])
    } else {
      warning(paste("Column", col, "not found in colData. Skipping..."))
    }
  }
  
  return(combinedSE)
}

generate_donut_plots_per_contrast <- function(results_DEP, results_PA, contrastList, outputDir) {
  
  for (contrast_name in contrastList) {
    
    col_pval <- paste0(contrast_name, "_p.val")
    
    if (!col_pval %in% colnames(results_DEP) || !col_pval %in% colnames(results_PA)) {
      warning(paste("Contrast", contrast_name, "has no p-values. Skipping..."))
      next
    }
    
    n_both <- sum(!is.na(results_DEP[[col_pval]]))
    n_absent <- sum(!is.na(results_PA[[col_pval]]))
    
    if (n_both == 0 && n_absent == 0) next
    
    proteinsPercentage <- data.frame(
      type = factor(c("Present in both groups", "Absent from one of the groups"),
                    levels = c("Present in both groups", "Absent from one of the groups")),
      number = c(n_both, n_absent)
    )
    
    donut_plot <- ggplot(proteinsPercentage, aes(x = 2, y = number, fill = type)) +
      geom_bar(stat = "identity", width = 1, color = "black") +
      coord_polar(theta = "y", start = 0) +
      xlim(0.1, 3) +
      theme_void() +
      theme(legend.position = "top",
            legend.text = element_text(size = 60),
            legend.title = element_text(size = 70, face = "bold",
                                        margin = margin(b = 40),
                                        hjust = 0.5), 
            legend.box.margin = margin(20, 0, 0, 0),
            legend.box.spacing = unit(-0.5, "cm"),
            plot.title = element_text(size = 65, face = "bold", hjust = 0.5, margin = margin(b = -80))) +
      labs(fill = "Total identified proteins") +
      geom_text(aes(label = number, x = 2.9),
                position = position_stack(vjust = 0.5),
                color = "black", size = 30) +
      scale_fill_manual(values = c("Absent from one of the groups" = "#870052",
                                   "Present in both groups" = "dimgray")) +
      guides(fill = guide_legend(direction = "vertical"))
    
    dir2save <- file.path(outputDir, paste0("DonutPlot_", contrast_name, ".tiff"))
    tiff(filename = dir2save, width = 4500, height = 4250, res = 300)
    print(donut_plot)
    dev.off()
  }
}