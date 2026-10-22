calculateES = function(df, ident, contrastList) {

  ident = deparse(substitute(ident))
  
  enrichmentScores = list()
  
  for (x in contrastList) {
    contrastES = paste0(x, "_ES")
    contrastRatio = paste0(x, "_ratio")
    contrastPval = paste0(x, "_p.val")
    
    foldchange = df[[contrastRatio]]
    pval = df[[contrastPval]]
    pval = ifelse(pval == 0 | pval < 1e-300, 1e-300, pval)
    enrichmentMetric = foldchange * -log10(pval)
    
    names(enrichmentMetric) = df[[ident]]
    
    enrichmentScores[[contrastES]] = sort(enrichmentMetric, decreasing = TRUE)
  }
  
  return(enrichmentScores)
}


translate_ID = function(input_ids, source_type, target_type, species_db = NULL, contrastList = NULL,
                         contrastESData = NULL, oriESData = NULL, source_taxon = NULL,
                         target_taxon = NULL, manual_map = NULL, manual_map_inter = NULL) {
  
  input_ids = as.character(input_ids)
  mapped = NULL
  
  if (!is.null(source_taxon) && !is.null(target_taxon) && (source_taxon != target_taxon)) {
    message("Using homologene() for interspecies mapping...")
    numeric_ids = suppressWarnings(as.numeric(input_ids))
    homolog_data = homologene(numeric_ids, inTax = source_taxon, outTax = target_taxon)
    source_key_col = paste0(source_taxon, "_ID")
    if (target_type == "ENTREZID") {
      target_key_col = paste0(target_taxon, "_ID")
    } else if (target_type == "SYMBOL") {
      target_key_col = as.character(target_taxon)
    } else {
      stop("For interspecies mapping, target_type must be either 'ENTREZID' or 'SYMBOL'")
    }
    mapped = homolog_data[[ target_key_col ]][match(numeric_ids, homolog_data[[ source_key_col ]])]
    if (!is.null(manual_map_inter)) {
      missing = is.na(mapped)
      mapped[missing] = manual_map_inter[[target_key_col]][match(input_ids[missing], manual_map_inter[[source_key_col]])]
    }
    
  } else {
    if (is.null(species_db)) {
      stop("For intraspecies mapping you must provide 'species_db' (e.g., org.Mm.eg.db).", call. = FALSE)
    }
    
    message("Using mapIds() for intraspecies mapping...")
    mapped = mapIds(species_db, keys = input_ids, column = target_type, keytype = source_type, multiVals = "first")
    
    mapped = as.character(mapped)
    
    if (!is.null(manual_map)) {
      missing = is.na(mapped)
      mapped[missing] = manual_map[[target_type]][match(input_ids[missing], manual_map[[source_type]])]
    }
  }
  
  if (!is.null(contrastList)){
    translated_results = list()
    for (contrast in contrastList) {
      contrastES = paste0(contrast, oriESData)
      contrastTarget = paste0(contrast, oriESData, target_type)
      ESTemp = getElement(contrastESData,contrastES)
      entrezIDNamesTemp = mapped[match(names(ESTemp), input_ids)]
      names(ESTemp) = entrezIDNamesTemp
      ESTemp = ESTemp[!is.na(names(ESTemp))]
        
      translated_results[[contrastTarget]] = ESTemp
    }
  }
  
  if (!is.null(contrastList)) {
    return(translated_results)
  } else {
    names(mapped) = input_ids
    out_df = data.frame(original = input_ids,
                         translated = mapped,
                         stringsAsFactors = FALSE)
    if (!is.null(source_taxon) && !is.null(target_taxon) && (source_taxon != target_taxon)){
      colnames(out_df)[2] = paste0(as.character(target_taxon), "_", target_type)
      return(out_df)
    } else {
      organismTAXIDData = metadata(species_db)
      organismTAXID = organismTAXIDData$value[11]
      colnames(out_df)[2] = paste0(organismTAXID, "_", target_type)
      return(out_df)
    }
  }
}


runEnrichmentAnalysis = function(geneListData, databases = c("GO_BP", "GO_MF", "GO_CC", "KEGG", "Reactome", "WikiPathways"),
                                  out_dir, OrgDb, organism, keyType, minGSSize, maxGSSize, pAdjustMethod, pvalueCutoff,
                                  verbose = FALSE, nPermSimple = 100000, eps = 0, adjpval = TRUE,
                                  direction = c("both", "enriched", "depleted"), widthValue = widthValue, 
                                  heightValue = heightValue, base_size = base_size) {

  if(!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)
  
  direction = match.arg(direction)

  if(!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)
  
  resultList = list()
  
  for(contrast in names(geneListData)) {
    geneList = geneListData[[contrast]]
    resultList[[contrast]] = list()
    
    for(db in databases) {
      
      enrichmentResult = NULL
      plotTitle = NULL
      
      if(db == "KEGG") {
        enrichmentResult = gseKEGG(geneList = geneList, organism = organism, 
                                    minGSSize = minGSSize, maxGSSize = maxGSSize, 
                                    pvalueCutoff = pvalueCutoff, pAdjustMethod = pAdjustMethod,
                                    keyType = keyType, nPermSimple = nPermSimple, eps = eps)
        plotTitle = "KEGG Terms"
      } else if(db %in% c("GO_BP", "GO_MF", "GO_CC")) {
        ont_value = switch(db,
                            "GO_BP" = "BP",
                            "GO_MF" = "MF",
                            "GO_CC" = "CC")
        enrichmentResult = gseGO(geneList = geneList, ont = ont_value, OrgDb = OrgDb, 
                                  keyType = keyType, minGSSize = minGSSize, maxGSSize = maxGSSize, 
                                  pAdjustMethod = pAdjustMethod, pvalueCutoff = pvalueCutoff, 
                                  verbose = verbose, nPermSimple = nPermSimple, eps = eps)
        plotTitle = paste0("GO-", ont_value, " Terms")
      } else if(db == "Reactome") {
        enrichmentResult = gsePathway(geneList = geneList, organism = organism, 
                                       minGSSize = minGSSize, maxGSSize = maxGSSize, 
                                       pvalueCutoff = pvalueCutoff, pAdjustMethod = pAdjustMethod,
                                       verbose = verbose, nPermSimple = nPermSimple, eps = eps)
        plotTitle = "Reactome Terms"
      } else if(db == "WikiPathways") {
        enrichmentResult = gseWP(geneList = geneList, organism = organism, 
                                  minGSSize = minGSSize, maxGSSize = maxGSSize, 
                                  pvalueCutoff = pvalueCutoff, pAdjustMethod = pAdjustMethod,
                                  nPermSimple = nPermSimple, eps = eps)
        plotTitle = "WikiPathways Terms"
      } else {
        warning("Database '", db, "' is not supported. Skipping...")
        next
      }
      
      df_enrichment = as.data.frame(enrichmentResult@result)
      if(nrow(df_enrichment) == 0) next
      
      csv_file = file.path(out_dir, paste0("Enrichment_", db, "_", contrast, ".csv"))
      write.csv(df_enrichment, csv_file, row.names = FALSE)
      
      if(direction == "enriched") {
        df_enrichment = df_enrichment[df_enrichment$enrichmentScore >  0, ]
      } else if(direction == "depleted") {
        df_enrichment = df_enrichment[df_enrichment$enrichmentScore <  0, ]
      }
      if(nrow(df_enrichment) == 0) next
      enrichmentResult@result = df_enrichment
      
      tiff_file = file.path(out_dir, paste0("Enrichment_", db, "_", contrast, ".tiff"))
      tiff(filename = tiff_file, width = widthValue, height = heightValue, res = 300)
      
      p = dotplot(enrichmentResult, showCategory = nrow(df_enrichment), x = "NES") +
        aes(fill = p.adjust) +
        ggtitle(plotTitle) +
        labs(title = "",
             x = "NES",
             fill = "Adj p") +
        geom_vline(xintercept = 0, color = "grey50", linewidth = 1.5, linetype = "dashed") +
        theme_void() +
        theme(plot.title = element_text(hjust = 0.5, margin = margin(b = 20)),
              axis.text.y = element_text(size = base_size, hjust = 1, margin = margin(r = 10)),
              axis.text.x = element_text(size = 26),
              axis.title.x = element_text(size = 26, face = "bold", margin = margin(t=5)),
              title = element_text(size = 34, face = "bold"),
              legend.text = element_text(size = 24),
              legend.title = element_text(size = 24, margin = margin(b = 20)),
              legend.position = "right",
              legend.box.just = "center",
              legend.title.align = 0.5,
              legend.spacing.y = unit(1.5, "cm"),
              legend.box.spacing = unit(1, "cm"),
              panel.border = element_rect(color = "black", fill = NA, linewidth = 2),
              axis.ticks = element_line(color = "black", linewidth = 1),
              plot.margin = unit(c(1, 1, 1, 1), "cm")) +
        guides(fill = guide_colorbar(
          title.hjust = 0.5,
          label.theme = element_text(hjust = 0.5),
          barheight = unit(6, "cm")
        ),
        size = guide_legend(
          title = "Counts",
          title.position = "top",
          title.hjust = 0.5,
          ncol = 1,
          label.theme = element_text(hjust = 0.5)
        )) +
        #annotate("rect", xmin = -Inf, xmax = 0, ymin = -Inf, ymax = Inf, fill = "dimgray", alpha = 0.1) +
        #annotate("rect", xmin = 0, xmax = Inf, ymin = -Inf, ymax = Inf, fill = "#870052", alpha = 0.1) +
        scale_x_continuous(expand = expansion(mult = c(0.15, 0.15)))+
        scale_size_continuous(range = c(10, 20),
                              breaks = function(x) unique(round(seq(min(x), max(x), length.out = 4)))) +
        scale_y_discrete(labels = function(x) str_wrap(x, width = 30)) +
        scale_fill_gradientn(labels = function(x) sprintf("%.3f", x),
                             colours =  colorRampPalette(c("#005C9A", "#A80022"))(100))
      print(p)
      dev.off()
      
      resultList[[contrast]][[db]] = list(
        enrichment = enrichmentResult,
        df = df_enrichment,
        plot = p,
        csv = csv_file,
        tiff = tiff_file
      )
    }
  } 
  return(resultList)
}


analyzeContrastFGSEA = function(ESScoreList, ESContrastSuffix, contrastList, gene_sets, organismName, padjval, gene_sets_Title, 
                                 out_dir, nPermSimple = 100000, eps = 0, direction = c("both", "enriched", "depleted"),
                                 widthValue = widthValue, heightValue = heightValue, base_size = base_size, rangeX = rangeX) {

  if(!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)
  
  direction = match.arg(direction)

  if(!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)
  
  contrast_results = list()
  
  for(z in contrastList) {
    
    contrast2analyzeString = paste0(z, ESContrastSuffix)
    contrast2analyze = ESScoreList[[contrast2analyzeString]]
    
    fgseaRes = fgseaMultilevel(pathways = gene_sets,
                                stats = contrast2analyze,
                                nPermSimple = nPermSimple,
                                eps = eps)
    
    res_fgseaRes = fgseaRes %>% unnest()
    
    fgseaRes_df = as.data.frame(fgseaRes)
    
    fgseaRes_save = fgseaRes_df
    
    fgseaRes_save$leadingEdge = sapply(fgseaRes_save$leadingEdge,
                                      function(x) paste(x, collapse = "/"))
    
    fgseaRes_sig = fgseaRes_df %>% filter(padj < padjval)
    
    csv_filename = paste0(out_dir, "gene_sets_Title", 
                           contrast2analyzeString, ".csv")
    write.csv(fgseaRes_save, csv_filename, row.names = FALSE)
    
    if(direction == "enriched"){
      fgseaRes_sig = fgseaRes_sig %>% filter(ES > 0)
    } else if(direction == "depleted"){
      fgseaRes_sig = fgseaRes_sig %>% filter(ES < 0)
    }
    
    gseaResult_obj = NA
    dotplot_obj = NA
    
    if(nrow(fgseaRes_sig) == 0) {
      message("No significant FGSEA results for ", contrast2analyzeString, ". Skipping plot.")
    } else {
      result_df = data.frame(
        ID = as.character(fgseaRes_sig$pathway),
        Description = as.character(fgseaRes_sig$pathway),
        setSize = fgseaRes_sig$size,
        enrichmentScore = fgseaRes_sig$ES,
        NES = fgseaRes_sig$NES,
        pvalue = fgseaRes_sig$pval,
        p.adjust = fgseaRes_sig$padj,
        qvalue = NA,
        rank = rank(fgseaRes_sig$pval, ties.method = "min"),
        leading_edge = sapply(fgseaRes_sig$leadingEdge, function(x) paste(x, collapse = "/")),
        core_enrichment = sapply(fgseaRes_sig$leadingEdge, function(x) paste(x, collapse = "/")),
        stringsAsFactors = FALSE
      )
      
      result_df = result_df[order(result_df$enrichmentScore), ]
      
      gseaResult_obj = new("gseaResult", result = result_df, geneSets = gene_sets,
                            organism = organismName, setType = "UserDefined", keytype = "ENTREZID")
      
      tiff_filename = paste0(out_dir, gene_sets_Title, "_",
                              contrast2analyzeString, ".tiff")
      tiff(filename = tiff_filename, width = widthValue, height = heightValue, res = 300)
      
      dotplot_obj = dotplot(gseaResult_obj,
                             showCategory = nrow(result_df),
                             x = "NES") +
        ggtitle("") +
        aes(fill = p.adjust) +
        labs(x = "NES",
             fill = "Adj p") +
        geom_vline(xintercept = 0, color = "grey50", linewidth = 1.5, linetype = "dashed") +
        theme_void() +
        theme(plot.title = element_text(hjust = 0.5, margin = margin(b = 20)),
              axis.text.y = element_text(size = base_size, hjust = 1, margin = margin(r = 15)),
              axis.text.x = element_text(size = 25),
              axis.title.x = element_text(size = 30, face = "bold", margin = margin(t=5)),
              title = element_text(size = 30, face = "bold"),
              legend.text = element_text(size = 18),
              legend.title = element_text(size = 20, margin = margin(b = 20), hjust = 0.5),
              legend.position = "right",
              legend.box.just = "center",
              legend.title.align = 0.5,
              legend.spacing.y = unit(1.5, "cm"),
              legend.box.spacing = unit(1, "cm"),
              panel.border = element_rect(color = "black", fill = NA, linewidth = 2),
              axis.ticks = element_line(color = "black", linewidth = 1),
              plot.margin = unit(c(1, 1, 1, 1), "cm")) +
        guides(color = guide_colorbar(
          title.hjust = 0.5,
          label.theme = element_text(hjust = 0.5)
        ),
        size = guide_legend(
          title = "Counts",
          title.position = "top",
          title.hjust = 0.5,
          ncol = 1,
          label.theme = element_text(hjust = 0.5)
        )) +
#        annotate("rect", xmin = -Inf, xmax = 0, ymin = -Inf, ymax = Inf,
#                 fill = "dimgray", alpha = 0.1) +
#        annotate("rect", xmin = 0, xmax = Inf, ymin = -Inf, ymax = Inf,
#                 fill = "#870052", alpha = 0.1) +
        scale_x_continuous(expand = expansion(mult = c(0.15, 0.15)),
                           breaks = scales::pretty_breaks(n = 4))+
        scale_size_continuous(range = c(10, 20),
                              breaks = function(x) unique(round(seq(min(x), max(x), length.out = 4)))) +
        scale_y_discrete(labels = function(x) str_wrap(x, width = 20)) +
        scale_fill_gradientn(labels = function(x) sprintf("%.3f", x),
          colours =  colorRampPalette(c("#005C9A", "#A80022"))(100))
      
      print(dotplot_obj)
      dev.off()
    } 
    
    contrast_results[[contrast2analyzeString]] = list(
      fgseaRes = fgseaRes,
      gseaResultObj = gseaResult_obj,
      dotplot = dotplot_obj
    )
  } 
  
  message("FGSEA analysis completed for all contrasts.")
  return(contrast_results)
}

analyzeCuratedAbundance = function(sigProt, gene_sets, translator, sigProtID, gene_setsID, gene_sets_Title,
                                   label_threshold, min_threshold = 0, term2Plot = NULL, out_dir, legendSize,
                                   is_log2 = TRUE, group_col = "Group") { 
  
  library(ggpubr)
  
  if(!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)
  
  expression_data = as.data.frame(assay(sigProt))
  
  if (is_log2) {
    expression_data = 2^expression_data
  }
  
  expression_data$MsSymbolName = rownames(expression_data)
  expression_data$MsENTREZID = translator[[gene_setsID]][match(tolower(expression_data$MsSymbolName), tolower(translator[[sigProtID]]))]
  
  annotation_df = tibble::enframe(gene_sets, name = "Category", value = "geneIDs") %>%
    dplyr::mutate(geneIDs = purrr::map(geneIDs, as.character)) %>%
    tidyr::unnest_longer(geneIDs) %>%
    dplyr::distinct() %>%
    dplyr::rename(MsENTREZID = geneIDs)
  
  combined_data = merge(expression_data, annotation_df, by = "MsENTREZID")
  
  long_data = combined_data %>%
    tidyr::pivot_longer(cols = -c(MsENTREZID, MsSymbolName, Category), 
                        names_to = "Sample",
                        values_to = "Expression")
  
  sample_info = as.data.frame(colData(sigProt))
  sample_info$Sample = rownames(sample_info)
  
  if (group_col %in% colnames(sample_info)) {
    sample_info$condition = as.character(sample_info[[group_col]])
  } else if (!"condition" %in% colnames(sample_info)) {
    stop(paste("Missing column", group_col))
  }
  
  plot_data = merge(long_data, sample_info, by = "Sample")
  
  plot_data = plot_data %>%
    dplyr::mutate(Expression = as.numeric(Expression))
  
  expression_sums = plot_data %>%
    dplyr::group_by(condition, Category) %>%
    dplyr::summarise(total_expression = sum(Expression, na.rm = TRUE), .groups = "drop")
  
  condition_totals = plot_data %>%
    dplyr::group_by(condition) %>%
    dplyr::summarise(condition_total = sum(Expression, na.rm = TRUE), .groups = "drop")
  
  abundance_data = expression_sums %>%
    dplyr::left_join(condition_totals, by = "condition") %>%
    dplyr::mutate(relative_abundance = total_expression / condition_total * 100)
  
  abundance_data_grouped = abundance_data %>%
    dplyr::group_by(condition) %>%
    dplyr::mutate(Category = ifelse(relative_abundance < min_threshold, "Others", as.character(Category))) %>%
    dplyr::group_by(condition, Category) %>%
    dplyr::summarise(total_expression = sum(total_expression, na.rm = TRUE),
                     condition_total = max(condition_total),
                     relative_abundance = sum(relative_abundance), .groups = "drop")
  
  abundance_data_grouped = abundance_data_grouped %>%
    dplyr::group_by(condition) %>%
    dplyr::arrange(dplyr::desc(relative_abundance)) %>%
    dplyr::mutate(Category = factor(Category, levels = unique(Category))) %>%
    dplyr::ungroup()
  
  replicate_abundance = plot_data %>%
    dplyr::group_by(condition, ID, Category) %>%
    dplyr::summarise(total_expression = sum(Expression, na.rm = TRUE), .groups = "drop")
  
  replicate_totals = plot_data %>%
    dplyr::group_by(condition, ID) %>%
    dplyr::summarise(rep_total = sum(Expression, na.rm = TRUE), .groups = "drop")
  
  replicate_abundance = replicate_abundance %>%
    dplyr::left_join(replicate_totals, by = c("condition", "ID")) %>%
    dplyr::mutate(relative_abundance = total_expression / rep_total * 100)
  
  localization_avg = replicate_abundance %>%
    dplyr::group_by(Category, condition) %>%
    dplyr::summarise(avg_abundance = mean(relative_abundance, na.rm = TRUE), .groups = "drop")
  
  categories_to_keep = localization_avg %>%
    dplyr::group_by(Category) %>%
    dplyr::summarise(max_avg = max(avg_abundance), .groups = "drop") %>%
    dplyr::filter(max_avg > min_threshold) %>%
    dplyr::pull(Category)
  
  replicate_abundance_filtered = replicate_abundance %>%
    dplyr::filter(Category %in% categories_to_keep)
  
  replicate_abundance_mean = replicate_abundance_filtered %>%
    dplyr::group_by(condition, Category) %>%
    dplyr::summarise(mean_total_expression = mean(relative_abundance, na.rm = TRUE), .groups = "drop")
  
  csv_filename = file.path(out_dir, "curated_annotation_abundance.csv")
  write.csv(abundance_data, csv_filename, row.names = FALSE)
  
  plots = list()
  
  unique_categories = unique(abundance_data_grouped$Category)
  num_categories = length(unique_categories)
  colorPalettePlot = colorRampPalette(brewer.pal(12, "Set3"))(num_categories)
  names(colorPalettePlot) = unique_categories
  
  unique_conditions = unique(abundance_data_grouped$condition)
  plot_list = list()
  
  for (cond in unique_conditions) {
    dat = abundance_data_grouped %>%
      dplyr::filter(condition == cond) %>%
      dplyr::arrange(dplyr::desc(relative_abundance))
    
    total = sum(dat$relative_abundance)
    
    dat = dat %>%
      dplyr::mutate(
        Category = factor(Category, levels = rev(unique(Category))),
        ymax = cumsum(relative_abundance),
        ymin = dplyr::lag(ymax, default = 0),
        labelPosition = (ymax + ymin) / 2,
        label = ifelse(relative_abundance < label_threshold, "", paste0(round(relative_abundance, 1), "%")),
        labelAngle = 90 - 360 * labelPosition / total,
        labelAngle = ifelse(labelAngle < -90, labelAngle + 180, labelAngle)
      )
    
    p = ggplot(dat, aes(x = 2, y = relative_abundance, fill = Category)) +
      geom_bar(stat = "identity", width = 1, color = "white") +
      coord_polar(theta = "y", start = 0, clip = "off") +
      xlim(c(0.7, 2.5)) +
      geom_text(aes(x = 2.2, y = labelPosition, label = label, angle = labelAngle),
                size = 12, fontface = "bold", hjust = 0.5)  +
      #annotate("text", x = 0.7, y = 0, label = cond,
      #         size = 18, fontface = "bold", angle = 0, hjust = 0.5) +
      scale_fill_manual(values = colorPalettePlot,
                        labels = function(x) stringr::str_wrap(x, width = 15),
                        guide = guide_legend(reverse = TRUE, ncol = 1,
                                             keyheight = unit(1.2, "cm"))) +
      theme_void() + 
      labs(fill = NULL) +
      theme(legend.position = "right",
            legend.text = element_text(size = legendSize, 
                                       margin = margin(t = 6, b = 6)), 
            legend.spacing.y = unit(0.6, "cm"),
            plot.title = element_text(size = 22),
            legend.box.margin = margin(-60, 0, 0, 0))
    
    plot_list[[cond]] = p
    
    individual_file = file.path(out_dir, paste0(cond, "_curated_", gene_sets_Title, "_annotation_abundance.tiff"))
    
    tiff(filename = individual_file, width = 4500, height = 4250, res = 300)
    print(p)
    dev.off()
  }
  
  plots$piechart = plot_list
  
  violin_file = file.path(out_dir, paste0("curated_", gene_sets_Title, "_annotation_abundance_violin.tiff"))
  tiff(filename = violin_file,  width = 4500, height = 4250, res = 300)
  
  violin_fig = ggplot(replicate_abundance_filtered, aes(x = Category, y = relative_abundance, fill = condition)) +
    geom_violin(trim = FALSE, position = position_dodge(width = 0.9)) +
    stat_compare_means(aes(group = condition), method = "anova", p.adjust.method = "BH", label = "p.signif") +
    labs(title = paste("Relative abundance (>", min_threshold, "%)", sep=""),
         x = NULL, y = "Relative abundance (%)") +
    theme(axis.text.x = element_text(angle = 90, hjust = 1, size = 20),
          axis.text.y = element_text(size = 20),
          axis.title.y = element_text(size = 15),
          plot.title = element_text(size = 20),
          legend.text = element_text(size = 20),
          legend.title = element_text(size = 20))
  
  print(violin_fig)
  dev.off()
  plots$violin = violin_fig
  
  if(term2Plot %in% unique(replicate_abundance_filtered$Category)) {
    mito_file = file.path(out_dir, paste0("curated_", gene_sets_Title, "_annotation_abundance_violin_", term2Plot, ".tiff"))
    tiff(filename = mito_file, width = 4500, height = 4250, res = 300)
    
    my_comparisons <- combn(unique_conditions, 2, simplify = FALSE)
    
    mito_fig = ggplot(dplyr::filter(replicate_abundance_filtered, Category == term2Plot),
                      aes(x = condition, y = relative_abundance, fill = condition)) + 
      geom_violin(trim = FALSE, alpha = 0.8) +
      geom_boxplot(width = 0.1, fill = "white", alpha = 0.5, outlier.shape = NA) + 
      stat_compare_means(comparisons = my_comparisons, method = "t.test", p.adjust.method = "BH", label = "p.signif", hide.ns = TRUE) + 
      stat_compare_means(method = "anova", label.y.npc = "bottom") + 
      labs(title = paste0("Relative abundance (", term2Plot, ")"),
           x = "Groups", y = "Relative abundance (%)") +
      theme(axis.text.x = element_text(hjust = 1, size = 20, angle = 45),
            axis.text.y = element_text(size = 20),
            axis.title.y = element_text(size = 15),
            plot.title = element_text(size = 20),
            legend.text = element_text(size = 20),
            legend.title = element_text(size = 20),
            legend.position = "none") 
    
    print(mito_fig)
    dev.off()
    plots$violin_mito = mito_fig
  }
  
  replicate_abundance_scaled = replicate_abundance_filtered %>%
    dplyr::group_by(Category) %>%
    dplyr::mutate(
      sd_val = sd(relative_abundance, na.rm = TRUE),
      zscore = ifelse(is.na(sd_val) | sd_val == 0, 0, 
                      (relative_abundance - mean(relative_abundance, na.rm = TRUE)) / sd_val)
    ) %>%
    dplyr::ungroup()
  
  replicate_matrix = replicate_abundance_scaled %>%
    dplyr::select(ID, Category, zscore) %>%
    tidyr::pivot_wider(names_from = Category, values_from = zscore) %>%
    tibble::column_to_rownames("ID") %>%
    as.matrix()
  
  heatmap_tile = Heatmap(
    replicate_matrix, 
    name = "Z-score",
    cluster_rows = TRUE,
    cluster_columns = TRUE,
    col = colorRamp2(
      c(-2, 0, 2), 
      c("#2E6F40", "white", "#4F0433")
    ),
    row_title = "Replicates",
    column_title = "Annotation Category",
    show_row_dend = TRUE,
    show_column_dend = TRUE
  )
  plots$heatmap_tile = heatmap_tile
  
  replicate_abundance_mean = replicate_abundance_mean %>%
    dplyr::group_by(Category) %>%
    dplyr::mutate(
      sd_val = sd(mean_total_expression, na.rm = TRUE),
      zscore = ifelse(is.na(sd_val) | sd_val == 0, 0, 
                      (mean_total_expression - mean(mean_total_expression, na.rm = TRUE)) / sd_val)
    ) %>%
    dplyr::ungroup()
  
  condition_matrix = replicate_abundance_mean %>%
    dplyr::select(condition, Category, zscore) %>%
    tidyr::pivot_wider(names_from = condition, values_from = zscore) %>%
    tibble::column_to_rownames("Category") %>%
    as.matrix()
  
  heatmap_mean = Heatmap(
    condition_matrix,
    name = "Z-score",
    cluster_rows = TRUE,
    cluster_columns = TRUE,
    col = colorRamp2(
      c(-2, 0, 2),
      c("#2E6F40", "white", "#4F0433")
    ),
    row_title = "Annotation Category",
    column_title = "Condition",
    show_row_dend = TRUE,
    show_column_dend = TRUE
  )
  plots$heatmap_mean = heatmap_mean
  
  return(list(
    abundance_data = abundance_data,
    replicate_abundance_mean = replicate_abundance_mean,
    plots = plots
  ))
}

runEnrichmentAnalysisORA = function(geneListData, 
                                     databases = c("GO_BP", "GO_MF", "GO_CC", "KEGG", "Reactome", "WikiPathways"),
                                     out_dir,
                                     universeSample,
                                     OrgDb, 
                                     organism,  
                                     keyType, 
                                     minGSSize, 
                                     maxGSSize, 
                                     pAdjustMethod, 
                                     pvalueCutoff,
                                     verbose = FALSE, 
                                     adjpval = TRUE,
                                     widthValue = widthValue,
                                     heightValue = heightValue,
                                     base_size = base_size) {
  
  if(!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)
  
  resultList = list()
  
  for (contrast in names(geneListData)) {
    geneList = geneListData[[contrast]]
    resultList[[contrast]] = list() 
    
    for (db in databases) {
      enrichmentResult = NULL
      basePlotTitle = NULL
      
      if (db == "KEGG") {
        enrichmentResult = enrichKEGG(gene = geneList, 
                                       organism = organism,
                                       universe = universeSample,
                                       keyType = keyType,
                                       minGSSize = minGSSize, 
                                       maxGSSize = maxGSSize,
                                       pvalueCutoff = pvalueCutoff, 
                                       pAdjustMethod = pAdjustMethod)
        basePlotTitle = "KEGG Terms"
      } else if (db %in% c("GO_BP", "GO_MF", "GO_CC")) {
        ont_value = switch(db,
                            "GO_BP" = "BP",
                            "GO_MF" = "MF",
                            "GO_CC" = "CC")
        enrichmentResult = enrichGO(gene = geneList, 
                                     ont = ont_value,
                                     OrgDb = OrgDb, 
                                     universe = universeSample,
                                     keyType = keyType,
                                     pvalueCutoff = pvalueCutoff, 
                                     pAdjustMethod = pAdjustMethod)
        basePlotTitle = paste0("GO-", ont_value, " Terms")
      } else if (db == "Reactome") {
        reactomeOrganism = if (organism == "mmu") "mouse" else organism
        enrichmentResult = enrichPathway(gene = geneList, 
                                          organism = reactomeOrganism,
                                          minGSSize = minGSSize, 
                                          universe = universeSample,
                                          maxGSSize = maxGSSize,
                                          pvalueCutoff = pvalueCutoff, 
                                          pAdjustMethod = pAdjustMethod)
        basePlotTitle = "Reactome Terms"
      } else if (db == "WikiPathways") {
        enrichmentResult = enrichWP(gene = geneList, 
                                     organism = organism,
                                     minGSSize = minGSSize,
                                     universe = universeSample,
                                     maxGSSize = maxGSSize,
                                     pvalueCutoff = pvalueCutoff, 
                                     pAdjustMethod = pAdjustMethod)
        basePlotTitle = "WikiPathways Terms"
      } else {
        warning("Database '", db, "' is not supported. Skipping this analysis...")
        next
      }
      
      df_enrichment = as.data.frame(enrichmentResult@result)
      if(nrow(df_enrichment) == 0) next
      
      df_top = df_enrichment %>% arrange(p.adjust, desc(Count))
      
      csv_file = file.path(out_dir, paste0("Enrichment_", db, "_", contrast, ".csv"))
      write.csv(df_enrichment, csv_file, row.names = FALSE)
      
      numTerms = sum(enrichmentResult@result$p.adjust < pvalueCutoff)
      if (numTerms > 40) {
        numTerms = 40
      }
      
      if (numTerms < 20) {
        labelSize = 8
      } else {
        labelSize = 6
      }
        
      if(numTerms > 0) {
        p = tryCatch({
          enrichmentResultCNET = pairwise_termsim(enrichmentResult)
          emapplot(enrichmentResultCNET, showCategory = numTerms,
                   layout = "kk") +
            ggtitle(paste0(basePlotTitle, " - ", contrast))
        }, error = function(e) {
          warning("Error in pairwise_termsim for contrast ", contrast, " in ", db, 
                  ": ", e$message, "\nCreating a placeholder plot.")
          ggplot(data = data.frame(x = 1, y = 1), aes(x, y)) +
            geom_point(alpha = 0) +
            scale_x_continuous(limits = c(0, 1)) +
            scale_y_continuous(limits = c(0, 1)) +
            ggtitle(paste0(basePlotTitle, " - ", contrast, "\n(No network plot available)")) +
            theme_minimal() +
            theme(legend.text = element_text(size = 32),
                  legend.title = element_text(size = 38),
                  text = element_text(size = 30),
                  plot.margin = grid::unit(c(1,1,1,1), "cm"))
        })
        
        tiff_file = file.path(out_dir, paste0("Enrichment_", db, "_", contrast, ".tiff"))
        tiff(filename = tiff_file, width = 4250, height = 4250, res = 300)
        tryCatch({
          p$layers[[3]] = NULL
          p = p +
            geom_text_repel(aes(label = label), size = labelSize) +
            scale_size_continuous(range = c(4, 10)) + 
            theme(legend.text = element_text(size = 15),
                  legend.title = element_text(size = 20),
                  text = element_text(size = 25))
          print(p)
        }, error = function(e) {
          warning("Error during printing plot for contrast ", contrast, " in ", db, ": ", e$message)
        })
        dev.off()
        
        tiff_file = file.path(out_dir, paste0("Enrichment_DOT_", db, "_", contrast, ".tiff"))
        tiff(filename = tiff_file, width = widthValue, height = heightValue, res = 300)
        d = NULL
        tryCatch({
          d = dotplot(enrichmentResult, showCategory = nrow(df_enrichment)) +
            aes(fill = p.adjust) +
            ggtitle(basePlotTitle) +
            labs(title = contrast,
                 x = "Gene Ratio",
                 fill = "Adj p") +
            geom_vline(xintercept = 0, linetype = "dashed", color = "gray50", linewidth = 1) + 
            theme_void() +
            theme(plot.title = element_text(hjust = 0.5, margin = margin(b = 20)),
                  axis.text.y = element_text(size = base_size, hjust = 1, margin = margin(r = 15)),
                  axis.text.x = element_text(size = 25),
                  axis.title.x = element_text(size = 30, face = "bold", margin = margin(t=5)),
                  title = element_text(size = 30, face = "bold"),
                  legend.text = element_text(size = 20),
                  legend.title = element_text(size = 24, margin = margin(b = 20), hjust = 0.5),
                  legend.position = "right",
                  legend.box.just = "center",
                  legend.title.align = 0.5,
                  legend.spacing.y = unit(1.5, "cm"),
                  legend.box.spacing = unit(1, "cm"),
                  panel.border = element_rect(color = "black", fill = NA, linewidth = 2),
                  axis.ticks = element_line(color = "black", linewidth = 1),
                  plot.margin = unit(c(1, 1, 1, 1), "cm")) +
            guides(fill = guide_colorbar(
              title.hjust = 0.5,
              label.theme = element_text(hjust = 0.5),
              barheight = unit(6, "cm")
            ),
            size = guide_legend(
              title = "Counts",
              title.position = "top",
              title.hjust = 0.5,
              ncol = 1,
              label.theme = element_text(hjust = 0.5)
            )) +
            #        annotate("rect", xmin = -Inf, xmax = 0, ymin = -Inf, ymax = Inf,
            #                 fill = "dimgray", alpha = 0.1) +
            #        annotate("rect", xmin = 0, xmax = Inf, ymin = -Inf, ymax = Inf,
            #                 fill = "#870052", alpha = 0.1) +
            scale_x_continuous(expand = expansion(mult = c(0.1, 0.2)), breaks = scales::pretty_breaks(n = 4), ) +
            scale_size_continuous(range = c(12, 25), ,
                                  breaks = function(x) unique(round(seq(min(x), max(x), length.out = 4)))) +
            scale_y_discrete(labels = function(x) str_wrap(x, width = 30)) +
            scale_fill_gradientn(labels = function(x) sprintf("%.3f", x),
                                 colours =  colorRampPalette(c("#005C9A", "#A80022"))(100))
          print(d)
        }, error = function(e) {
          warning("Error during printing plot for contrast ", contrast, " in ", db, ": ", e$message)
        })
        dev.off()
        
      resultList[[contrast]][[db]] = list(
        enrichObject = enrichmentResult,
        enrichment = df_top,
        df = df_enrichment,
        plot = d,
        dotplot = d,
        csv = csv_file,
        tiff = tiff_file
      )
      } else {
        resultList[[contrast]][[db]] = list(
          enrichObject = enrichmentResult,
          enrichment = df_top,
          df = df_enrichment,
          csv = csv_file)
      }
    }
  }
  return(resultList)
}


analyzeEnrichmentTermsclusters = function(enrich_df,
                                           p_adjust_thresh = 0.125,
                                           count_thresh = 3,
                                           minClusterSize = 40,
                                           plot_sim_thresh = 0.2,
                                           similarity_method = c("GO", "overlap"),
                                           overlap_metric = c("jaccard", "overlap_coef"),
                                           ont = "BP",
                                           output_dir = ".",
                                           max_keywords = 10,
                                           top_n_clusters = NULL,
                                           use_wordcloud = TRUE,
                                           use_leading_term = TRUE,
                                           titledotplot = "") {

  similarity_method = match.arg(similarity_method)
  overlap_metric = match.arg(overlap_metric)
  
  if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)
  
  ggrepel_available = requireNamespace("ggrepel", quietly = TRUE)
  
  required_cols = c("ID", "p.adjust", "Count", "Description", "geneID")
  missing_cols = setdiff(required_cols, colnames(enrich_df))
  if (length(missing_cols) > 0) stop("enrich_df is missing required columns: ", paste(missing_cols, collapse = ", "))
  
  enrich_filtered = enrich_df %>%
    dplyr::filter(p.adjust < p_adjust_thresh, Count > count_thresh)
  
  if (nrow(enrich_filtered) == 0) stop("No terms remain after filtering by p.adjust and Count. Adjust thresholds.")
  
  enrich_filtered$GeneRatio_num = sapply(enrich_filtered$GeneRatio, function(x) eval(parse(text = as.character(x))))
  
  go_ids = enrich_filtered$ID
  
  if (similarity_method == "GO") {
    mat = tryCatch({
      GO_similarity(go_ids, ont = ont)
    }, error = function(e) {
      message("GO_similarity failed with ont argument; trying without ont.")
      GO_similarity(go_ids)
    })
    mat[is.na(mat)] = 0
    diag(mat) = 1
    mat = (mat + t(mat)) / 2
  } else {
    gene_lists_raw = enrich_filtered$geneID
    names(gene_lists_raw) = enrich_filtered$ID
    gene_lists = lapply(gene_lists_raw, function(x) {
      if (is.na(x) || x == "") return(character(0))
      unique(unlist(strsplit(as.character(x), "/")))
    })
    ids_gl = names(gene_lists)
    n = length(gene_lists)
    mat = matrix(0, nrow = n, ncol = n, dimnames = list(ids_gl, ids_gl))
    for (i in seq_len(n)) {
      gi = gene_lists[[i]]
      for (j in i:n) {
        gj = gene_lists[[j]]
        if (length(gi) == 0 || length(gj) == 0) {
          s = 0
        } else {
          inter = length(intersect(gi, gj))
          if (overlap_metric == "jaccard") {
            union_len = length(union(gi, gj))
            s = if (union_len == 0) 0 else inter / union_len
          } else {
            s = inter / min(length(gi), length(gj))
          }
        }
        mat[i, j] = s
        mat[j, i] = s
      }
    }
    mat[is.na(mat)] = 0
    diag(mat) = 1
    mat = (mat + t(mat)) / 2
  }
  
  meanSim = rowMeans(mat)
  keep = meanSim >= 0.05
  if (sum(keep) < 2) stop("Too few terms pass the mean similarity >= 0.05 filter. Adjust thresholds.")
  mat_filt = mat[keep, keep]
  enrich_filt = enrich_filtered %>% filter(ID %in% rownames(mat_filt))
  
  cl_auto = cluster_by_dynamicTreeCut(mat_filt, minClusterSize = minClusterSize)
  res_clusters = data.frame(ID = rownames(mat_filt), cluster = as.integer(cl_auto), stringsAsFactors = FALSE)
  
  go_with_cluster = enrich_filt %>% left_join(res_clusters, by = "ID")
  go_with_cluster = go_with_cluster %>% mutate(gene_list = strsplit(as.character(geneID), "/"))
  
  gwc = go_with_cluster %>% filter(!is.na(cluster))
  if (nrow(gwc) == 0) stop("No terms have cluster assignments after join. Check IDs and clustering.")
  
  tab_n = tapply(gwc$ID, gwc$cluster, length)
  mean_padj = tapply(gwc$p.adjust, gwc$cluster, function(x) mean(-log10(x), na.rm = TRUE))
  mean_gene_ratio = tapply(gwc$GeneRatio_num, gwc$cluster, function(x) mean(x, na.rm = TRUE))
  enrich_by_cluster = data.frame(
    cluster = as.integer(names(tab_n)),
    n_terms = as.integer(as.vector(tab_n)),
    mean_minusLog10_padj = as.numeric(mean_padj[as.character(names(tab_n))]),
    mean_gene_ratio = as.numeric(mean_gene_ratio[as.character(names(tab_n))]),
    stringsAsFactors = FALSE
  )
  enrich_by_cluster = enrich_by_cluster[order(-enrich_by_cluster$n_terms), ]
  
  ids_mat = rownames(mat_filt)
  cluster_vec = res_clusters$cluster[match(ids_mat, res_clusters$ID)]
  cluster_mean_sim = tapply(seq_along(ids_mat), cluster_vec, function(idx) {
    idx = idx[!is.na(idx)]
    if (length(idx) <= 1) return(NA_real_)
    subm = mat_filt[idx, idx, drop = FALSE]
    mean(subm[lower.tri(subm)], na.rm = TRUE)
  })
  cluster_mean_sim = data.frame(cluster = as.integer(names(cluster_mean_sim)), mean_sim = as.numeric(cluster_mean_sim), stringsAsFactors = FALSE)
  cluster_mean_sim = cluster_mean_sim[!is.na(cluster_mean_sim$cluster), , drop = FALSE]
  cluster_mean_sim$cluster = as.integer(cluster_mean_sim$cluster)
  
  align_to = split(seq_along(ids_mat), cluster_vec)
  go_id_list = split(ids_mat, cluster_vec)
  non_empty = lengths(go_id_list) > 0
  align_to = align_to[non_empty]; go_id_list = go_id_list[non_empty]
  
  if (use_leading_term) {
    keywords_df = gwc %>%
      dplyr::group_by(cluster) %>%
      dplyr::arrange(p.adjust, desc(abs(Count))) %>% 
      dplyr::slice(1) %>% 
      dplyr::ungroup() %>%
      dplyr::select(cluster, Description) %>%
      dplyr::rename(keywords = Description) %>%
      dplyr::mutate(keywords = stringr::str_to_sentence(as.character(keywords)))
  } else {
    wc_list = NULL
    if (similarity_method == "GO" && isTRUE(use_wordcloud)) {
      wc_list = tryCatch({
        anno_word_cloud_from_GO(align_to = align_to, go_id = go_id_list, return_gbl = TRUE,
                               stat = "pvalue", min_stat = 0.01, max_words = 50)
      }, error = function(e) {
        message("anno_word_cloud_from_GO failed or was skipped: ", conditionMessage(e))
        NULL
      })
    }
  
    if (!is.null(wc_list) && !is.null(wc_list$word_table)) {
      wt = wc_list$word_table
      if (!"cluster" %in% colnames(wt) && "group" %in% colnames(wt)) wt = wt %>% rename(cluster = group)
      if (!"word" %in% colnames(wt) && "term" %in% colnames(wt)) wt = wt %>% rename(word = term)
      if (!"freq" %in% colnames(wt) && "count" %in% colnames(wt)) wt = wt %>% rename(freq = count)
      wt = wt %>% filter(!is.na(cluster))
      wt$cluster = as.integer(wt$cluster); wt$word = as.character(wt$word); wt$freq = as.numeric(wt$freq)
      clusters_present = sort(unique(wt$cluster))
      keywords_list = lapply(clusters_present, function(cl) {
        sub = wt[wt$cluster == cl, , drop = FALSE]
        sub = sub[order(-sub$freq, na.last = TRUE), , drop = FALSE]
        topw = unique(sub$word)[1:min(max_keywords, nrow(sub))]
        data.frame(cluster = as.integer(cl), keywords = paste(topw, collapse = " "), stringsAsFactors = FALSE)
      })
      if (length(keywords_list) == 0) {
        keywords_df = data.frame(cluster = integer(0), keywords = character(0), stringsAsFactors = FALSE)
      } else {
        keywords_df = do.call(rbind, keywords_list)
        rownames(keywords_df) = NULL
      }
    } else {
      data("stop_words", package = "tidytext")
      toks = gwc %>%
        dplyr::select(cluster, Description) %>%
        unnest_tokens(word, Description) %>%
        filter(!word %in% stop_words$word) %>%
        filter(nchar(word) > 2)
      clusters_present = sort(unique(toks$cluster))
      keywords_list = lapply(clusters_present, function(cl) {
        sub = toks[toks$cluster == cl, , drop = FALSE]
        tab = sort(table(sub$word), decreasing = TRUE)
        topw = names(tab)[1:min(max_keywords, length(tab))]
        data.frame(cluster = as.integer(cl), keywords = paste(topw, collapse = " "), stringsAsFactors = FALSE)
      })
        if (length(keywords_list) == 0) {
        keywords_df = data.frame(cluster = integer(0), keywords = character(0), stringsAsFactors = FALSE)
      } else {
        keywords_df = do.call(rbind, keywords_list)
        rownames(keywords_df) = NULL
      }
    }
  }
  
  if (!"cluster" %in% colnames(keywords_df)) stop("keywords_df does not contain column 'cluster' after processing.")
  
  plot_df = enrich_by_cluster %>%
    left_join(cluster_mean_sim, by = "cluster") %>%
    left_join(keywords_df, by = "cluster") %>%
    mutate(
      keywords = ifelse(is.na(keywords), paste0("cluster_", cluster), keywords),
      mean_sim = ifelse(is.na(mean_sim), 0, mean_sim)
    )
  
  plot_df = plot_df %>% filter(mean_sim > plot_sim_thresh)
  
  if (!is.null(top_n_clusters) && is.numeric(top_n_clusters) && top_n_clusters > 0) {
    top_ids = plot_df %>% arrange(desc(n_terms)) %>% slice_head(n = top_n_clusters) %>% pull(cluster)
    plot_df = plot_df %>% filter(cluster %in% top_ids)
  }
  
  if (nrow(plot_df) == 0) {
    warning("No clusters pass the plot similarity threshold. Returning computed objects without plotting.")
    return(list(plot_df = plot_df,
                enrich_by_cluster = enrich_by_cluster,
                keywords_df = keywords_df,
                cluster_mean_sim = cluster_mean_sim,
                res_clusters = res_clusters,
                go_with_cluster = go_with_cluster))
  }
  
  split_to_three_lines = function(text) {
    words = unlist(strsplit(text, "\\s+"))
    n = length(words)
    if (n <= 3) return(paste(words, collapse = " "))
    part1 = ceiling(n / 3)
    part2 = ceiling((n - part1) / 2)
    i1 = 1:part1
    i2 = (part1 + 1):(part1 + part2)
    i3 = (part1 + part2 + 1):n
    line1 = paste(words[i1], collapse = " ")
    line2 = if (length(i2) > 0) paste(words[i2], collapse = " ") else ""
    line3 = if (length(i3) > 0) paste(words[i3], collapse = " ") else ""
    lines = c(line1, line2, line3)
    lines = lines[lines != ""]
    paste(lines, collapse = "\n")
  }
  
  split_to_two_lines = function(text) {
    text = trimws(gsub("\\s+", " ", as.character(text)))
    if (nchar(text) == 0) return("")
    words = unlist(strsplit(text, "\\s+"))
    n = length(words)
    if (n <= 2) return(paste(words, collapse = " "))
    cut_point = ceiling(n / 2)
    if (cut_point == 1 && n > 2) cut_point = 2
    line1 = paste(words[1:cut_point], collapse = " ")
    if (cut_point + 1 <= n) {
      line2 = paste(words[(cut_point + 1):n], collapse = " ")
      return(paste(line1, line2, sep = "\n"))
    } else {
      return(line1)
    }
  }
  
  plot_df$label_y_wrapped = vapply(plot_df$keywords, split_to_two_lines, FUN.VALUE = character(1))
  plot_df = plot_df %>% arrange(desc(mean_sim))
  plot_df$label_y_wrapped = factor(plot_df$label_y_wrapped, levels = unique(plot_df$label_y_wrapped))
  
  common_ids = intersect(rownames(mat_filt), res_clusters$ID)
  mat_plot = mat_filt[common_ids, common_ids]
  ann = data.frame(cluster = factor(res_clusters$cluster[match(common_ids, res_clusters$ID)]))
  rownames(ann) = common_ids
  heatmap_file = file.path(output_dir, "Clusters_Heatmap.tiff")
  heatmapPlot = pheatmap(mat_plot,
                         annotation_row = ann,
                         show_rownames = FALSE,
                         show_colnames = FALSE,
                         clustering_method = "average",
                         main = "Similarity heatmap")
  ragg::agg_tiff(filename = heatmap_file, width = 20, height = 20, units = "in", res = 300)
  print(heatmapPlot)
  dev.off()
  
  plot_df = plot_df %>% arrange(desc(mean_sim))
  plot_df$label_y_wrapped = factor(plot_df$label_y_wrapped, levels = unique(plot_df$label_y_wrapped))
  plot_df$label_y_wrapped = factor(plot_df$label_y_wrapped, levels = rev(levels(plot_df$label_y_wrapped)))
  
  dotplot_file = file.path(output_dir, "Clusters_DotPlot.tiff")
  p = ggplot(plot_df, aes(x = mean_gene_ratio, y = label_y_wrapped)) +
    geom_point(aes(size = n_terms, fill = mean_minusLog10_padj), shape = 21, color = "black", alpha = 0.95) +
    labs(x = "Mean Gene Ratio",
         y = NULL,
         title = titledotplot) +
    theme_void() +
    theme(
      plot.title = element_text(hjust = 0.5, margin = margin(b = 20), size = 50, face = "bold"),
      axis.text.y = element_text(size = 30, hjust = 1, margin = margin(r = 15)),
      axis.text.x = element_text(size = 30),
      axis.title.x = element_text(size = 35, face = "bold", margin = margin(t = 5)),
      legend.text = element_text(size = 35),
      legend.title = element_text(size = 35, margin = margin(b = 20)),
      legend.position = "right",
      legend.box.just = "center",
      legend.title.align = 0.5,
      legend.spacing.y = unit(1.5, "cm"),
      legend.box.spacing = unit(1, "cm"),
      legend.key.height = unit(1.5, "cm"),
      panel.border = element_rect(color = "black", fill = NA, linewidth = 2),
      axis.ticks = element_line(color = "black", linewidth = 1),
      plot.margin = unit(c(1, 1, 1, 1), "cm")
    ) +
    scale_alpha_continuous(range = c(0.5, 1), name = "Similarity") +
    guides(
      fill = guide_colorbar(
        title = "-log10(Adj p)",
        title.hjust = 0.5,
        title.theme = element_text(size = 35, hjust = 0.5, face = "bold"),
        barwidth = unit(0.6, "cm"),
        barheight = unit(6, "cm")
      ),
      size = guide_legend(
        title = "Counts",
        title.position = "top",
        title.theme = element_text(size = 35, hjust = 0.5, face = "bold"),
        ncol = 1
      ),
      alpha = guide_legend(
        title = "Similarity",
        title.position = "top",
        title.hjust = 0.5,
        title.theme = element_text(size = 35, face = "bold", hjust = 0.5), 
        override.aes = list(size = 18)
      )
    ) +
    #annotate("rect", xmin = -Inf, xmax = 0, ymin = -Inf, ymax = Inf, fill = "dimgray", alpha = 0.1) +
    #annotate("rect", xmin = 0, xmax = Inf, ymin = -Inf, ymax = Inf, fill = "#870052", alpha = 0.1) +
    scale_x_continuous(expand = expansion(add = c(0.05, 0.1))) +
    scale_size(range = c(12, 25)) +
    scale_fill_gradientn(colours = colorRampPalette(c("#005C9A", "#A80022"))(100),
                         name = "-log10(Adj p)")
  
  #  if (ggrepel_available) {
  #    p = p + ggrepel::geom_text_repel(data = plot_df, aes(x = mean_sim, y = label_y_wrapped, label = n_terms),
  #                                      color = "black", size = 8, nudge_x = 0.02, direction = "y",
  #                                      segment.size = 0.2, max.overlaps = 20)
  #  } else {
  #    p = p + geom_text(data = plot_df, aes(x = mean_sim, y = label_y_wrapped, label = n_terms),
  #                       color = "black", size = 8, hjust = -0.3)
  #  }
  
  ragg::agg_tiff(filename = dotplot_file, width = 20, height = 20, units = "in", res = 300)
  print(p)
  dev.off()
  
  write.csv(plot_df, file.path(output_dir, "plot_df_clusters_summary.csv"), row.names = FALSE)
  write.csv(enrich_by_cluster, file.path(output_dir, "enrich_by_cluster.csv"), row.names = FALSE)
  write.csv(keywords_df, file.path(output_dir, "keywords_by_cluster.csv"), row.names = FALSE)
  
  return(list(
    plot_df = plot_df,
    enrich_by_cluster = enrich_by_cluster,
    keywords_df = keywords_df,
    cluster_mean_sim = cluster_mean_sim,
    res_clusters = res_clusters,
    go_with_cluster = go_with_cluster
  ))
}


analyzeEnrichmentTermsclustersGSEA = function(enrich_df,
                                          p_adjust_thresh = 0.125,
                                          count_thresh = 3,
                                          minClusterSize = 40,
                                          plot_sim_thresh = 0.2,
                                          similarity_method = c("GO", "overlap"),
                                          overlap_metric = c("jaccard", "overlap_coef"),
                                          ont = "BP",
                                          output_dir = ".",
                                          max_keywords = 10,
                                          top_n_clusters = NULL,
                                          use_wordcloud = TRUE,
                                          use_leading_term = TRUE,
                                          titledotplot = "") {
  
  similarity_method = match.arg(similarity_method)
  overlap_metric = match.arg(overlap_metric)
  
  if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)
  
  ggrepel_available = requireNamespace("ggrepel", quietly = TRUE)
  
  required_cols = c("ID", "p.adjust", "Count", "Description", "geneID")
  missing_cols = setdiff(required_cols, colnames(enrich_df))
  if (length(missing_cols) > 0) stop("enrich_df is missing required columns: ", paste(missing_cols, collapse = ", "))
  
  enrich_filtered = enrich_df %>%
    dplyr::filter(p.adjust < p_adjust_thresh, Count > count_thresh)
  
  if (nrow(enrich_filtered) == 0) stop("No terms remain after filtering by p.adjust and Count. Adjust thresholds.")
  
  go_ids = enrich_filtered$ID
  
  if (similarity_method == "GO") {
    mat = tryCatch({
      GO_similarity(go_ids, ont = ont)
    }, error = function(e) {
      message("GO_similarity failed with ont argument; trying without ont.")
      GO_similarity(go_ids)
    })
    mat[is.na(mat)] = 0
    diag(mat) = 1
    mat = (mat + t(mat)) / 2
  } else {
    gene_lists_raw = enrich_filtered$geneID
    names(gene_lists_raw) = enrich_filtered$ID
    gene_lists = lapply(gene_lists_raw, function(x) {
      if (is.na(x) || x == "") return(character(0))
      unique(unlist(strsplit(as.character(x), "/")))
    })
    ids_gl = names(gene_lists)
    n = length(gene_lists)
    mat = matrix(0, nrow = n, ncol = n, dimnames = list(ids_gl, ids_gl))
    for (i in seq_len(n)) {
      gi = gene_lists[[i]]
      for (j in i:n) {
        gj = gene_lists[[j]]
        if (length(gi) == 0 || length(gj) == 0) {
          s = 0
        } else {
          inter = length(intersect(gi, gj))
          if (overlap_metric == "jaccard") {
            union_len = length(union(gi, gj))
            s = if (union_len == 0) 0 else inter / union_len
          } else {
            s = inter / min(length(gi), length(gj))
          }
        }
        mat[i, j] = s
        mat[j, i] = s
      }
    }
    mat[is.na(mat)] = 0
    diag(mat) = 1
    mat = (mat + t(mat)) / 2
  }
  
  if ("NES" %in% colnames(enrich_filtered)) {
    nes_vec = enrich_filtered$NES[match(rownames(mat), enrich_filtered$ID)]
    same_sign_mat = outer(sign(nes_vec), sign(nes_vec), "==")
    mat[!same_sign_mat] = 0
  }
  
  meanSim = rowMeans(mat)
  keep = meanSim >= 0.05
  if (sum(keep) < 2) stop("Too few terms pass the mean similarity >= 0.05 filter. Adjust thresholds.")
  mat_filt = mat[keep, keep]
  enrich_filt = enrich_filtered %>% filter(ID %in% rownames(mat_filt))
  
  cl_auto = cluster_by_dynamicTreeCut(mat_filt, minClusterSize = minClusterSize)
  res_clusters = data.frame(ID = rownames(mat_filt), cluster = as.integer(cl_auto), stringsAsFactors = FALSE)
  
  go_with_cluster = enrich_filt %>% left_join(res_clusters, by = "ID")
  go_with_cluster = go_with_cluster %>% mutate(gene_list = strsplit(as.character(geneID), "/"))
  
  gwc = go_with_cluster %>% filter(!is.na(cluster))
  if (nrow(gwc) == 0) stop("No terms have cluster assignments after join. Check IDs and clustering.")
  
  tab_n = tapply(gwc$ID, gwc$cluster, length)
  mean_padj = tapply(gwc$p.adjust, gwc$cluster, function(x) mean(-log10(x), na.rm = TRUE))
  mean_NES = tapply(gwc$NES, gwc$cluster, function(x) mean(x, na.rm = TRUE))
  enrich_by_cluster = data.frame(
    cluster = as.integer(names(tab_n)),
    n_terms = as.integer(as.vector(tab_n)),
    mean_minusLog10_padj = as.numeric(mean_padj[as.character(names(tab_n))]),
    mean_NES = as.numeric(mean_NES[as.character(names(tab_n))]), # NUEVO
    stringsAsFactors = FALSE
  )
  enrich_by_cluster = enrich_by_cluster[order(-enrich_by_cluster$n_terms), ]
  
  ids_mat = rownames(mat_filt)
  cluster_vec = res_clusters$cluster[match(ids_mat, res_clusters$ID)]
  cluster_mean_sim = tapply(seq_along(ids_mat), cluster_vec, function(idx) {
    idx = idx[!is.na(idx)]
    if (length(idx) <= 1) return(NA_real_)
    subm = mat_filt[idx, idx, drop = FALSE]
    mean(subm[lower.tri(subm)], na.rm = TRUE)
  })
  cluster_mean_sim = data.frame(cluster = as.integer(names(cluster_mean_sim)), mean_sim = as.numeric(cluster_mean_sim), stringsAsFactors = FALSE)
  cluster_mean_sim = cluster_mean_sim[!is.na(cluster_mean_sim$cluster), , drop = FALSE]
  cluster_mean_sim$cluster = as.integer(cluster_mean_sim$cluster)
  
  align_to = split(seq_along(ids_mat), cluster_vec)
  go_id_list = split(ids_mat, cluster_vec)
  non_empty = lengths(go_id_list) > 0
  align_to = align_to[non_empty]; go_id_list = go_id_list[non_empty]
  
  if (use_leading_term) {
    keywords_df = gwc %>%
      dplyr::group_by(cluster) %>%
      dplyr::arrange(p.adjust, desc(abs(NES))) %>% 
      dplyr::slice(1) %>% 
      dplyr::ungroup() %>%
      dplyr::select(cluster, Description) %>%
      dplyr::rename(keywords = Description) %>%
      dplyr::mutate(keywords = stringr::str_to_sentence(as.character(keywords)))
  } else {
    wc_list = NULL
    if (similarity_method == "GO" && isTRUE(use_wordcloud)) {
      wc_list = tryCatch({
        anno_word_cloud_from_GO(align_to = align_to, go_id = go_id_list, return_gbl = TRUE,
                               stat = "pvalue", min_stat = 0.01, max_words = 50)
      }, error = function(e) {
        message("anno_word_cloud_from_GO failed or was skipped: ", conditionMessage(e))
        NULL
      })
    }
  
    if (!is.null(wc_list) && !is.null(wc_list$word_table)) {
      wt = wc_list$word_table
      if (!"cluster" %in% colnames(wt) && "group" %in% colnames(wt)) wt = wt %>% rename(cluster = group)
      if (!"word" %in% colnames(wt) && "term" %in% colnames(wt)) wt = wt %>% rename(word = term)
      if (!"freq" %in% colnames(wt) && "count" %in% colnames(wt)) wt = wt %>% rename(freq = count)
      wt = wt %>% filter(!is.na(cluster))
      wt$cluster = as.integer(wt$cluster); wt$word = as.character(wt$word); wt$freq = as.numeric(wt$freq)
      clusters_present = sort(unique(wt$cluster))
      keywords_list = lapply(clusters_present, function(cl) {
        sub = wt[wt$cluster == cl, , drop = FALSE]
        sub = sub[order(-sub$freq, na.last = TRUE), , drop = FALSE]
        topw = unique(sub$word)[1:min(max_keywords, nrow(sub))]
        data.frame(cluster = as.integer(cl), keywords = paste(topw, collapse = " "), stringsAsFactors = FALSE)
      })
      if (length(keywords_list) == 0) {
        keywords_df = data.frame(cluster = integer(0), keywords = character(0), stringsAsFactors = FALSE)
      } else {
        keywords_df = do.call(rbind, keywords_list)
        rownames(keywords_df) = NULL
      }
    } else {
      data("stop_words", package = "tidytext")
      toks = gwc %>%
        dplyr::select(cluster, Description) %>%
        unnest_tokens(word, Description) %>%
        filter(!word %in% stop_words$word) %>%
        filter(nchar(word) > 2)
      clusters_present = sort(unique(toks$cluster))
      keywords_list = lapply(clusters_present, function(cl) {
        sub = toks[toks$cluster == cl, , drop = FALSE]
        tab = sort(table(sub$word), decreasing = TRUE)
        topw = names(tab)[1:min(max_keywords, length(tab))]
        data.frame(cluster = as.integer(cl), keywords = paste(topw, collapse = " "), stringsAsFactors = FALSE)
      })
      if (length(keywords_list) == 0) {
        keywords_df = data.frame(cluster = integer(0), keywords = character(0), stringsAsFactors = FALSE)
      } else {
        keywords_df = do.call(rbind, keywords_list)
        rownames(keywords_df) = NULL
      }
    }
  }
  
  if (!"cluster" %in% colnames(keywords_df)) stop("keywords_df does not contain column 'cluster' after processing.")
  
  plot_df = enrich_by_cluster %>%
    left_join(cluster_mean_sim, by = "cluster") %>%
    left_join(keywords_df, by = "cluster") %>%
    mutate(
      keywords = ifelse(is.na(keywords), paste0("cluster_", cluster), keywords),
      mean_sim = ifelse(is.na(mean_sim), 0, mean_sim)
    )
  
  plot_df = plot_df %>% filter(mean_sim > plot_sim_thresh)
  
  if (!is.null(top_n_clusters) && is.numeric(top_n_clusters) && top_n_clusters > 0) {
    top_ids = plot_df %>% arrange(desc(n_terms)) %>% slice_head(n = top_n_clusters) %>% pull(cluster)
    plot_df = plot_df %>% filter(cluster %in% top_ids)
  }
  
  if (nrow(plot_df) == 0) {
    warning("No clusters pass the plot similarity threshold. Returning computed objects without plotting.")
    return(list(plot_df = plot_df,
                enrich_by_cluster = enrich_by_cluster,
                keywords_df = keywords_df,
                cluster_mean_sim = cluster_mean_sim,
                res_clusters = res_clusters,
                go_with_cluster = go_with_cluster))
  }
  
  split_to_three_lines = function(text) {
    words = unlist(strsplit(text, "\\s+"))
    n = length(words)
    if (n <= 3) return(paste(words, collapse = " "))
    part1 = ceiling(n / 3)
    part2 = ceiling((n - part1) / 2)
    i1 = 1:part1
    i2 = (part1 + 1):(part1 + part2)
    i3 = (part1 + part2 + 1):n
    line1 = paste(words[i1], collapse = " ")
    line2 = if (length(i2) > 0) paste(words[i2], collapse = " ") else ""
    line3 = if (length(i3) > 0) paste(words[i3], collapse = " ") else ""
    lines = c(line1, line2, line3)
    lines = lines[lines != ""]
    paste(lines, collapse = "\n")
  }
  
  split_to_two_lines = function(text) {
    text = trimws(gsub("\\s+", " ", as.character(text)))
    if (nchar(text) == 0) return("")
    words = unlist(strsplit(text, "\\s+"))
    n = length(words)
    if (n <= 2) return(paste(words, collapse = " "))
    cut_point = ceiling(n / 2)
    if (cut_point == 1 && n > 2) cut_point = 2
    line1 = paste(words[1:cut_point], collapse = " ")
    if (cut_point + 1 <= n) {
      line2 = paste(words[(cut_point + 1):n], collapse = " ")
      return(paste(line1, line2, sep = "\n"))
    } else {
      return(line1)
    }
  }
  
  plot_df$label_y_wrapped = vapply(plot_df$keywords, split_to_two_lines, FUN.VALUE = character(1))
  plot_df = plot_df %>% arrange(desc(mean_sim))
  plot_df$label_y_wrapped = factor(plot_df$label_y_wrapped, levels = unique(plot_df$label_y_wrapped))
  
  common_ids = intersect(rownames(mat_filt), res_clusters$ID)
  mat_plot = mat_filt[common_ids, common_ids]
  ann = data.frame(cluster = factor(res_clusters$cluster[match(common_ids, res_clusters$ID)]))
  rownames(ann) = common_ids
  heatmap_file = file.path(output_dir, "Clusters_Heatmap.tiff")
  heatmapPlot = pheatmap(mat_plot,
                         annotation_row = ann,
                         show_rownames = FALSE,
                         show_colnames = FALSE,
                         clustering_method = "average",
                         main = "Similarity heatmap")
  ragg::agg_tiff(filename = heatmap_file, width = 20, height = 20, units = "in", res = 300)
  print(heatmapPlot)
  dev.off()
  
  plot_df = plot_df %>% arrange(desc(mean_sim))
  plot_df$label_y_wrapped = factor(plot_df$label_y_wrapped, levels = unique(plot_df$label_y_wrapped))
  plot_df$label_y_wrapped = factor(plot_df$label_y_wrapped, levels = rev(levels(plot_df$label_y_wrapped)))
  
  dotplot_file = file.path(output_dir, "Clusters_DotPlot.tiff")
  p = ggplot(plot_df, aes(x = mean_NES, y = label_y_wrapped)) +
    geom_point(aes(size = n_terms, fill = mean_minusLog10_padj, alpha = mean_sim), 
               shape = 21, color = "black") +
    geom_vline(xintercept = 0, linetype = "dashed", color = "gray50", linewidth = 1) + 
    labs(x = "Mean NES",
         y = NULL,
         title = titledotplot) +
    theme_void() +
    theme(
      plot.title = element_text(hjust = 0.5, margin = margin(b = 20), size = 50, face = "bold"),
      panel.border = element_rect(color = "black", fill = NA, linewidth = 2),
      axis.text.y = element_text(size = 30, hjust = 1, margin = margin(r = 15)),
      axis.text.x = element_text(size = 30),
      axis.title.x = element_text(size = 35, face = "bold", margin = margin(t = 5)),
      legend.text = element_text(size = 35),
      legend.title = element_text(size = 35, margin = margin(b = 20)),
      legend.position = "right",
      legend.box.just = "center",
      legend.title.align = 0.5,
      legend.spacing.y = unit(1.5, "cm"),
      legend.box.spacing = unit(1, "cm"),
      legend.key.height = unit(1.5, "cm"),
      axis.ticks = element_line(color = "black", linewidth = 1),
      plot.margin = unit(c(1, 1, 1, 1), "cm")
    ) +
    scale_alpha_continuous(range = c(0.5, 1), name = "Similarity") +
    guides(
      fill = guide_colorbar(
        title = "-log10(Adj p)",
        title.hjust = 0.5,
        title.theme = element_text(size = 35, hjust = 0.5, face = "bold"),
        barwidth = unit(0.6, "cm"),
        barheight = unit(6, "cm")
      ),
      size = guide_legend(
        title = "Counts",
        title.position = "top",
        title.theme = element_text(size = 35, hjust = 0.5, face = "bold"),
        ncol = 1
      ),
      alpha = guide_legend(
        title = "Similarity",
        title.position = "top",
        title.hjust = 0.5,
        title.theme = element_text(size = 35, face = "bold", hjust = 0.5), 
        override.aes = list(size = 18)
      )
    ) +
    #annotate("rect", xmin = -Inf, xmax = 0, ymin = -Inf, ymax = Inf, fill = "dimgray", alpha = 0.1) +
    #annotate("rect", xmin = 0, xmax = Inf, ymin = -Inf, ymax = Inf, fill = "#870052", alpha = 0.1) +
    scale_x_continuous(breaks = c(-2.5, -2, -1.5, -1, -0.5, 0, 0.5, 1, 1.5, 2, 2.5),
                       labels = c("-2.5", "-2", "-1.5", "-1", "-0.5", "0", "0.5", "1", "1.5", "2", "2.5")) +
    scale_size(range = c(12, 25)) +
    scale_fill_gradientn(colours = colorRampPalette(c("#005C9A", "#A80022"))(100),
                         name = "-log10(Adj p)")
  
  #  if (ggrepel_available) {
  #    p = p + ggrepel::geom_text_repel(data = plot_df, aes(x = mean_sim, y = label_y_wrapped, label = n_terms),
  #                                      color = "black", size = 8, nudge_x = 0.02, direction = "y",
  #                                      segment.size = 0.2, max.overlaps = 20)
  #  } else {
  #    p = p + geom_text(data = plot_df, aes(x = mean_sim, y = label_y_wrapped, label = n_terms),
  #                       color = "black", size = 8, hjust = -0.3)
  #  }
  
  ragg::agg_tiff(filename = dotplot_file, width = 20, height = 20, units = "in", res = 300)
  print(p)
  dev.off()
  
  write.csv(plot_df, file.path(output_dir, "plot_df_clusters_summary.csv"), row.names = FALSE)
  write.csv(enrich_by_cluster, file.path(output_dir, "enrich_by_cluster.csv"), row.names = FALSE)
  write.csv(keywords_df, file.path(output_dir, "keywords_by_cluster.csv"), row.names = FALSE)
  
  return(list(
    plot_df = plot_df,
    enrich_by_cluster = enrich_by_cluster,
    keywords_df = keywords_df,
    cluster_mean_sim = cluster_mean_sim,
    res_clusters = res_clusters,
    go_with_cluster = go_with_cluster
  ))
}

# Generate Enriched/Depleted Heatmap for Enrichment Results
#
# generatePresenceHeatmap() builds and saves a heatmap in which each cell
# represents the enrichment direction of a term in a given contrast:
#   • +1 if the term is significantly enriched (NES > 0, p.adjust < cutoff)
#   • -1 if the term is significantly depleted (NES < 0, p.adjust < cutoff)
#   •  0 if non-significant or absent
# The function accepts a named list of enrichment results (each element can
# be a list with $df, an enrichResult/gseaResult, or a plain data.frame).
# Terms are extracted from “Description” (or “pathway”), filtered by
# adjusted p-value and NES sign, assembled into a matrix, and plotted via
# ComplexHeatmap. The output is written as a TIFF.
#
# Parameters:
# -resultList: Named list of contrasts. Each item is itself a list of
#  enrichment objects or data.frames for that contrast.
# -out_dir:     Character. Directory to save the output file
#                “divergent_heatmap.tiff” (created if missing).
# -contrasts:   Optional character vector of contrast names (and column
#                order). Defaults to names(resultList).
# -pAdjustCutoff: Numeric. Threshold on adjusted p-value (default 0.05).
#
# Output:
# Invisibly returns a list with:
# -matrix: Binary/ternary matrix (rows = terms, cols = contrasts) with
#           values -1, 0, or +1.
# -heatmap: The ComplexHeatmap Heatmap object (for further customization).
#
generatePresenceHeatmap <- function(resultList, out_dir  = ".", contrasts = NULL,
                                    pAdjustCutoff = 0.05,col_labels = NULL,
                                    bottom_group_label = NULL, 
                                    wrap_width = 30, widthVal = 3000,
                                    heightVal = 4250) {
  if (is.null(contrasts)) contrasts <- names(resultList)
  dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
  if (!requireNamespace("ComplexHeatmap", quietly = TRUE)) {
    stop("Please install the ComplexHeatmap package first.")
  }
  library(ComplexHeatmap)
  library(circlize)
  library(grid)
  
  # Extract enriched/depleted term names
  getTermStates <- function(obj) {
    if (is.list(obj) && !is.null(obj$fgseaRes) && is.data.frame(obj$fgseaRes)) {
      df <- as.data.frame(obj$fgseaRes)
    } else if (is.list(obj) && !is.null(obj$df) && is.data.frame(obj$df)) {
      df <- obj$df
    } else if (inherits(obj, c("enrichResult", "gseaResult"))) {
      df <- as.data.frame(obj)
    } else if (is.data.frame(obj)) {
      df <- obj
    } else {
      stop("Unsupported object type: ", class(obj)[1])
    }
    # unify padj → p.adjust
    if (!"p.adjust" %in% colnames(df) && "padj" %in% colnames(df)) {
      df$p.adjust <- df$padj
    }
    if (!"p.adjust" %in% colnames(df)) {
      stop("Data frame must contain 'p.adjust' or 'padj' column.")
    }
    df_sig <- df[df$p.adjust < pAdjustCutoff, , drop = FALSE]
    if (nrow(df_sig) == 0) return(list(enriched = character(0),
                                       depleted = character(0)))
    term_col <- if ("Description" %in% colnames(df_sig)) {
      "Description"
    } else if ("pathway" %in% colnames(df_sig)) {
      "pathway"
    } else {
      stop("Neither 'Description' nor 'pathway' present.")
    }
    enriched <- unique(df_sig[[term_col]][df_sig$NES > 0])
    depleted <- unique(df_sig[[term_col]][df_sig$NES < 0])
    list(enriched = enriched, depleted = depleted)
  }
  
  # Only keep objects getTermStates() can handle
  isValidEnrichObj <- function(obj) {
    (is.list(obj) && !is.null(obj$fgseaRes) && is.data.frame(obj$fgseaRes)) ||
      (is.list(obj) && !is.null(obj$df)      && is.data.frame(obj$df))      ||
      inherits(obj, c("enrichResult", "gseaResult")) ||
      is.data.frame(obj)
  }
  
  # Collect all unique terms from the specified contrasts
  states_list <- lapply(contrasts, function(ct) {
    lst       <- resultList[[ct]]
    valid_lst <- Filter(isValidEnrichObj, lst)
    if (length(valid_lst) == 0) return(character(0))
    st        <- lapply(valid_lst, getTermStates)
    c(unlist(lapply(st, "[[", "enriched"), use.names = FALSE),
      unlist(lapply(st, "[[", "depleted"),  use.names = FALSE))
  })
  all_terms <- sort(unique(unlist(states_list, use.names = FALSE)))
  
  if (!is.null(col_labels) && length(col_labels) == length(contrasts)) {
    c_names <- col_labels
  } else {
    c_names <- sub("_vs_.*$", "", contrasts)
  }
  
  # Build the presence/depletion matrix
  mat <- matrix(
    0L,
    nrow     = length(all_terms),
    ncol     = length(contrasts),
    dimnames = list(all_terms,c_names)
  )
  for (j in seq_along(contrasts)) {
    lst       <- resultList[[contrasts[j]]]
    valid_lst <- Filter(isValidEnrichObj, lst)
    if (length(valid_lst) == 0) next
    st        <- lapply(valid_lst, getTermStates)
    enriched  <- unique(unlist(lapply(st, "[[", "enriched"), use.names = FALSE))
    depleted  <- unique(unlist(lapply(st, "[[", "depleted"),  use.names = FALSE))
    mat[all_terms %in% enriched, j] <-  1L
    mat[all_terms %in% depleted,  j] <- -1L
  }
  mat <- mat[rowSums(abs(mat)) > 0L, , drop = FALSE]
  
  # Divergent color mapping & dynamic font sizing
  col_fun <- colorRamp2(
    breaks = c(-1, 0, +1),
    colors = c("#005C9A", "white", "#A80022")
  )
  
  if (nrow(mat)<20){
    fontRowSize = 20
  } else if (nrow(mat)<50){
    fontRowSize = 15
  } else if (nrow(mat)<100){
    fontRowSize = 10
  } else {
    fontRowSize = 5
  }
  
  # Draw heatmap
  if (!is.null(bottom_group_label)) {
    bot_anno <- HeatmapAnnotation(
      col_text = anno_text(c_names, rot = 90, just = "right", 
                           location = unit(1, "npc") - unit(2, "mm"), 
                           gp = gpar(fontsize = 20, fontface = "bold")),
      ref_label = anno_empty(border = FALSE, height = unit(1.5, "cm")),
      show_annotation_name = FALSE
    )
    show_c_names <- FALSE 
  } else {
    bot_anno <- NULL
    show_c_names <- TRUE
  }
  
  if (!is.null(wrap_width)) {
    wrapped_row_names <- stringr::str_wrap(rownames(mat), width = wrap_width)
  } else {
    wrapped_row_names <- rownames(mat)
  }
  
  ht <- Heatmap(
    mat,
    name = " ",
    row_labels         = wrapped_row_names,
    col                = col_fun,
    cluster_rows       = TRUE,
    cluster_columns    = FALSE,
    show_row_names     = TRUE,
    show_column_names  = show_c_names,
    row_names_side     = "left",
    row_dend_side      = "right",
    rect_gp            = gpar(col = "black", lwd = 1),
    border             = TRUE,
    row_names_gp       = gpar(fontsize =  fontRowSize, fontface = "bold"),
    column_names_gp    = gpar(fontsize = 2, fontface = "bold"),
    column_title = deparse(substitute(resultList)),
    show_heatmap_legend = FALSE,
    row_names_max_width = unit(50, "cm"),
    bottom_annotation  = bot_anno
  )
  
  # Save to TIFF
  tiff_file <- file.path(out_dir, paste0("EnrichmentHeatmap_", deparse(substitute(resultList)), ".tiff"))
  tiff(filename = tiff_file, width = widthVal, height = heightVal, res = 300)
  draw(ht)
  
  if (!is.null(bottom_group_label)) {
    decorate_annotation("ref_label", {
      grid.lines(x = c(0, 1), y = c(0.8, 0.8), gp = gpar(lwd = 2))
      grid.text(bottom_group_label, x = 0.5, y = 0.4, gp = gpar(fontsize = 18, fontface = "bold"))
    })
  }
  
  dev.off()
  
  invisible(list(matrix = mat, heatmap = ht))
}