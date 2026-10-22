createProteomicsLocalizationDb = function(file_subcell, file_proteomo, file_uniprot, translator,
                                           translatorHumanUniprot, translatorMouseEntrez, translatorHumanSymbol,
                                           translatorMouseUniprot) {
  
  `%||%` = function(a, b) if (is.null(a)) b else a
  
  subcell_df  = read.delim(file_subcell, stringsAsFactors = FALSE)
  proteomo_df = read_excel(file_proteomo, sheet = "subcell_annotation_summary")
  uniprot_df  = read.delim(file_uniprot, stringsAsFactors = FALSE)
  trim_df = function(df) df %>% mutate(across(where(is.character), ~ ifelse(is.na(.x), NA_character_, str_squish(.x))))
  subcell_df  = trim_df(subcell_df); proteomo_df = trim_df(proteomo_df); uniprot_df  = trim_df(uniprot_df)
  
  if (all(c("Main.location", "Additional.location") %in% names(subcell_df))) {
    subcell_df = subcell_df %>%
      mutate(Additional.location = ifelse(is.na(Additional.location) | Additional.location == "", NA_character_, Additional.location),
             Main.location = ifelse(is.na(Additional.location),
                                    str_squish(Main.location),
                                    str_c(str_squish(Main.location), str_squish(Additional.location), sep = "; "))) %>%
      mutate(Main.location = na_if(str_squish(Main.location), ""))
  } else subcell_df$Main.location = subcell_df$Main.location %||% NA_character_
  
  if ("Majority protein IDs" %in% names(proteomo_df)) {
    proteomo_expanded = proteomo_df %>%
      mutate(`Majority protein IDs` = as.character(`Majority protein IDs`)) %>%
      separate_rows(`Majority protein IDs`, sep = ";\\s*") %>%
      mutate(`Majority protein IDs` = str_squish(`Majority protein IDs`))
  } else proteomo_expanded = proteomo_df
  
  if ("Subcellular.location..CC." %in% names(uniprot_df)) {
    uniprot_df = uniprot_df %>%
      mutate(raw_loc = ifelse(is.na(Subcellular.location..CC.), "", Subcellular.location..CC.)) %>%
      mutate(raw_loc = raw_loc %>%
               str_replace_all("\\[Isoform.*?\\]:", "") %>%
               str_replace_all("\\[.*?\\]", "") %>%
               str_replace_all("(?i)Note=.*", "") %>%
               str_replace_all("(?i)SUBCELLULAR LOCATION:.", "") %>%
               str_squish()) %>%
      mutate(SubcellularLocation = ifelse(raw_loc == "", NA_character_, raw_loc))
  } else uniprot_df$SubcellularLocation = NA_character_
  
  protein_ids = tibble(
    mouseUniprotName  = as.character(translator[[translatorMouseUniprot]]),
    humanUniprotName  = as.character(translator[[translatorHumanUniprot]]),
    humanSymbolName   = as.character(translator[[translatorHumanSymbol]])
  ) %>% mutate(across(everything(), ~ ifelse(is.na(.x) | .x == "", NA_character_, str_squish(.x))))
  
  protein_ids = protein_ids %>%
    mutate(humanSymbol_trim = str_squish(humanSymbolName),
           humanUniprot_trim = str_squish(humanUniprotName)) %>%
    mutate(
      SubLocat_HPA = if ("Gene.name" %in% names(subcell_df)) subcell_df$Main.location[match(humanSymbol_trim, subcell_df$Gene.name)] else NA_character_,
      SubLocat_Proteomo = if ("Majority protein IDs" %in% names(proteomo_expanded) && "graph_localization_annotation" %in% names(proteomo_expanded))
        proteomo_expanded$graph_localization_annotation[match(humanUniprot_trim, proteomo_expanded$`Majority protein IDs`)]
      else NA_character_,
      SubLocat_UniProt = if ("Entry" %in% names(uniprot_df)) uniprot_df$SubcellularLocation[match(humanUniprot_trim, uniprot_df$Entry)] else NA_character_
    ) %>%
    mutate(
      SubLocat = coalesce(SubLocat_HPA, SubLocat_Proteomo, SubLocat_UniProt),
      SubLocat_source = case_when(
        !is.na(SubLocat_HPA) ~ "HPA",
        is.na(SubLocat_HPA) & !is.na(SubLocat_Proteomo) ~ "ProteomoDB",
        is.na(SubLocat_HPA) & is.na(SubLocat_Proteomo) & !is.na(SubLocat_UniProt) ~ "UniProt",
        TRUE ~ NA_character_
      )
    ) %>%
    mutate(SubLocat = ifelse(is.na(SubLocat) | SubLocat == "", NA_character_, as.character(SubLocat))) %>%
    select(-SubLocat_HPA, -SubLocat_Proteomo, -SubLocat_UniProt, -humanSymbol_trim, -humanUniprot_trim)
  
  protein_ids_split = protein_ids %>%
    mutate(SubLocat = ifelse(is.na(SubLocat), NA_character_, SubLocat)) %>%
    separate(SubLocat, into = paste0("Loc", 1:24), sep = ";\\s*|,\\s*|/\\s*|\\|\\s*", fill = "right", extra = "drop") %>%
    mutate(across(starts_with("Loc"), ~ na_if(str_trim(.x), "")))
  
  mouse_uniprot_vec = translator[[translatorMouseUniprot]]
  mouse_entrez_vec  = translator[[translatorMouseEntrez]]
  if (!is.null(names(mouse_uniprot_vec)) && length(names(mouse_uniprot_vec)) > 0) {
    trans_map = tibble(mouseUniprotName = str_squish(names(mouse_uniprot_vec)), mEntrez = as.character(mouse_entrez_vec))
  } else {
    trans_map = tibble(mouseUniprotName = str_squish(as.character(mouse_uniprot_vec)), mEntrez = as.character(mouse_entrez_vec))
  }
  trans_map = trans_map %>% distinct(mouseUniprotName, .keep_all = TRUE)
  protein_ids_split = protein_ids_split %>% mutate(mouseUniprotName = str_squish(mouseUniprotName)) %>% left_join(trans_map, by = "mouseUniprotName")
  
  clean_term_vec = function(x) {
    x = as.character(x)
    x[is.na(x)] = ""
    t = x
    t = str_replace_all(t, "[\\{\\}]", "")         
    t = str_replace_all(t, "\\[[^]]*\\]", "")     
    t = str_replace_all(t, "(?i)pubmed[:\\s0-9,|]*", "")
    t = str_replace_all(t, "(?i)uniprotkb[:\\s0-9,|]*", "")
    t = str_replace_all(t, "(?i)\\bref(erence)?\\b\\.?","") 
    t = str_replace_all(t, "^:\\s*", "")
    t = str_replace_all(t, "\\(.*?\\)", "")
    t = str_replace_all(t, "[\\.;]+$", "")
    t = str_replace_all(t, "[\\.;]", ",")
    t = str_replace_all(t, "_", " ")
    t = str_replace_all(t, "-isoform.*", "")
    t = str_replace_all(t, "-\\d+$", "")
    t = str_replace_all(t, "\\b(ec[0-9]+|eco:[0-9]+)\\b", "")
    t = str_replace_all(t, "\\s+,\\s*", ", ")
    t = str_replace_all(t, "\\s+", " ")
    t = str_trim(t)
    t[t == ""] = NA_character_
    return(t)
  }
  
  loc_long = protein_ids_split %>%
    pivot_longer(cols = starts_with("Loc"), names_to = "Loc_num", values_to = "raw_term") %>%
    filter(!is.na(raw_term)) %>%
    mutate(raw_term = str_trim(raw_term)) %>%
    mutate(term_clean = clean_term_vec(raw_term)) %>%
    filter(!is.na(term_clean)) %>%
    mutate(term_split = str_split(term_clean, ",\\s*")) %>%
    unnest_longer(term_split) %>%
    mutate(term_split = str_trim(term_split)) %>%
    filter(term_split != "") %>%
    select(mouseUniprotName, mEntrez, SubLocat_source, term_split)
  
  gene_df = loc_long %>%
    rename(Localization = term_split) %>%
    filter(!is.na(mEntrez)) %>%
    distinct()

  gene_sets_prelim = gene_df %>%
    group_by(Localization) %>%
    summarise(genes = list(unique(na.omit(mEntrez))), .groups = "drop") %>%
    deframe()
  
  raw_terms_df = gene_df %>% count(Localization, name = "n") %>% arrange(desc(n))
  
  is_bad_token_vec = function(x) {
    x = as.character(x)
    x[is.na(x)] = ""
    t = tolower(gsub("[\\{\\}\\[\\]]", "", x))
    pat_rule = "(hamap-rule|prosite-prorule|prosite|pdb:|pirnr:|pirnr|uniprotkb:|arba:|hamap|pfam|ipr)"
    is_rule = grepl(pat_rule, t, ignore.case = TRUE)
    is_uniprot_acc = grepl("^[OPQ][0-9][A-Z0-9]{3}[0-9]$", t, perl = TRUE)  
    is_a0a_acc = grepl("^a0a[0-9a-z]{3,}$", t)                            
    is_pdb = grepl("^pdb[:\\s]*[0-9a-z]+", t)
    is_generic_acc = grepl("^[a-z0-9]{5,}$", t) & grepl("[0-9]", t)         
    has_colon = grepl(":", t)
    only_digits = grepl("^\\d+$", t)
    letters_count = nchar(gsub("[^A-Za-z]", "", t))
    too_short = letters_count < 2
    res = is_rule | is_uniprot_acc | is_a0a_acc | is_pdb | is_generic_acc | has_colon | only_digits | too_short
    return(res)
  }
  
  tmp_df = tibble(original = names(gene_sets_prelim), genes = gene_sets_prelim) %>%
    tidyr::unnest_longer(genes) %>%
    mutate(cleaned = clean_term_vec(original)) %>%
    filter(!is.na(cleaned)) %>%
    mutate(term = strsplit(cleaned, ",\\s*")) %>%
    tidyr::unnest_longer(term) %>%
    mutate(term = trimws(term)) %>%
    filter(term != "") %>%
    filter(!is_bad_token_vec(term)) %>%
    filter(nchar(gsub("[^A-Za-z]", "", term)) >= 2)
  
  tmp_df = tmp_df %>%
    mutate(term_lower = tolower(term)) %>%
    filter(!term_lower %in% c("er", "er.", "ref")) %>%
    select(-term_lower)
  
  gene_sets_cleaned = tmp_df %>%
    group_by(term) %>%
    summarise(genes = list(unique(na.omit(genes))), .groups = "drop") %>%
    deframe()
  
  raw_terms_df = tmp_df %>% count(term, name = "n") %>% arrange(desc(n))
  
  canonical_map = c(
    "Nucleoplasm" = "Nucleus",
    "Cytosol" = "Cytoplasm",
    "Vesicles" = "Vesicles",
    "Plasma membrane" = "Plasma membrane",
    "Golgi apparatus" = "Golgi apparatus",
    "Mitochondria" = "Mitochondrion",
    "Nucleoli" = "Nucleus",
    "Nuclear bodies" = "Nucleus",
    "Endoplasmic reticulum" = "Endoplasmic reticulum",
    "Cell membrane" = "Plasma membrane",
    "Nuclear speckles" = "Nucleus",
    "Cytoplasm" = "Cytoplasm",
    "Multi-pass membrane protein" = "Others",
    "Centrosome" = "Microtubules",
    "Basal body" = "Primary cilium",
    "Primary cilium" = "Primary cilium",
    "Principal piece" = "Flagellum",
    "Secreted" = "Extracellular",
    "Cell Junctions" = "Cell Junctions",
    "Mid piece" = "Flagellum",
    "Microtubules" = "Microtubules",
    "Nucleoli fibrillar center" = "Nucleus",
    "Nuclear membrane" = "Nucleus",
    "Actin filaments" = "Actin cytoskeleton",
    "Cell projection" = "Others",
    "Centriolar satellite" = "Microtubules",
    "Cytokinetic bridge" = "Mitotic apparatus",
    "cytoskeleton" = "Others",
    "End piece" = "Flagellum",
    "Nucleoli rim" = "Nucleus",
    "Cytoplasmic vesicle" = "Vesicles",
    "Focal adhesion sites" = "Actin cytoskeleton",
    "Intermediate filaments" = "Intermediate filaments",
    "Mitotic spindle" = "Mitotic apparatus",
    "extracellular space" = "Extracellular",
    "Flagellar centriole" = "Flagellum",
    "Membrane" = "Plasma membrane",
    "Acrosome" = "Others",
    "Single-pass type I membrane protein" = "Others",
    "Nucleus" = "Nucleus",
    "Connecting piece" = "Flagellum",
    "Primary cilium transition zone" = "Primary cilium",
    "Cell junction" = "Cell Junctions",
    "microtubule organizing center" = "Microtubules",
    "Mitotic chromosome" = "Chromosome",
    "secretory vesicle" = "Vesicles",
    "Lipid-anchor" = "Extracellular",
    "Cytoplasmic bodies" = "Cytoplasmatic granules",
    "Equatorial segment" = "Mitotic apparatus",
    "Calyx" = "Others",
    "Perinuclear theca" = "Nucleus",
    "Primary cilium tip" = "Primary cilium",
    "Midbody" = "Mitotic apparatus",
    "Annulus" = "Flagellum",
    "cilium" = "Primary cilium",
    "mitochondrion" = "Mitochondrion",
    "cytosol" = "Cytoplasm",
    "extracellular matrix" = "Extracellular",
    "Lipid droplets" = "Lipid droplets",
    "plasma membrane" = "Plasma membrane",
    "Chromosome" = "Chromosome",
    "Peripheral membrane protein" = "Plasma membrane",
    "centrosome" = "Microtubules",
    "Peroxisomes" = "Peroxisomes",
    "Midbody ring" = "Mitotic apparatus",
    "nucleus" = "Nucleus",
    "Aggresome" = "Cytoplasmatic granules",
    "Endoplasmic reticulum membrane" = "Endoplasmic reticulum",
    "Single-pass type II membrane protein" = "Others",
    "unclassified" = "Others",
    "Lysosomes" = "Lysosomes",
    "Rods & Rings" = "Cytoplasmatic granules",
    "Postsynaptic cell membrane" = "Plasma membrane",
    "centromere" = "Chromosome",
    "Endosomes" = "Endosomes",
    "myofibril" = "Actin cytoskeleton",
    "GPI-anchor" = "Extracellular",
    "Synapse" = "Others",
    "lysosome" = "Lysosomes",
    "sarcomere" = "Actin cytoskeleton",
    "early endosome" = "Endosomes",
    "Recycling endosome membrane" = "Endosomes",
    "Cytoplasmic side" = "Cytoplasm",
    "Lysosome" = "Lysosomes",
    "Endosome" = "Endosomes",
    "Kinetochore" = "Chromosome",
    "Microtubule ends" = "Microtubules",
    "Single-pass membrane protein" = "Others",
    "actin cytoskeleton" = "Actin cytoskeleton",
    "Golgi apparatus membrane" = "Golgi apparatus",
    "Mitochondrion matrix" = "Mitochondrion",
    "kinetochore" = "Chromosome",
    "stress granule" = "Cytoplasmatic granules",
    "Early endosome" = "Endosomes",
    "Golgi" = "Golgi apparatus",
    "Mitochondrion inner membrane" = "Mitochondrion",
    "gap junction" = "Cell Junctions",
    "nucleolus" = "Nucleus",
    "14-3-3 scaffold" = "Others",
    "Cell surface" = "Plasma membrane",
    "Cytoplasmic granule" = "Cytoplasmatic granules",
    "Endosome membrane" = "Endosomes",
    "Golgi stack membrane" = "Golgi apparatus",
    "Microsome membrane" = "Endoplasmic reticulum",
    "Mitochondrion membrane" = "Mitochondrion",
    "sarcolemma" = "Plasma membrane",
    "spindle" = "Mitotic apparatus",
    "tight junction" = "Cell Junctions",
    "trans-Golgi" = "Golgi apparatus",
    "Cleavage furrow" = "Mitotic apparatus",
    "Lysosome membrane" = "Lysosomes",
    "Mitochondrion" = "Mitochondrion",
    "Mitochondrion outer membrane" = "Mitochondrion",
    "basement membrane" = "Extracellular",
    "recycling endosome" = "Endosomes",
    "surface film" = "Extracellular",
    "trans-Golgi network membrane" = "Golgi apparatus",
    "Apical cell membrane" = "Plasma membrane",
    "Cytoplasmic vesicle membrane" = "Vesicles",
    "Endomembrane system" = "Endosomes",
    "Extracellular side" = "Extracellular",
    "Intermembrane side" = "Others",
    "Matrix side" = "Extracellular",
    "Nucleus envelope" = "Nucleus",
    "Nucleus membrane" = "Nucleus",
    "P-body" = "Cytoplasmatic granules",
    "Recycling endosome" = "Endosomes",
    "Recycling endosome membrane" = "Endosomes",
    "Sarcoplasmic reticulum" = "Endoplasmic reticulum",
    "Single-pass type III membrane protein" = "Others",
    "perinuclear region" = "Nucleus",
    "spindle pole" = "Mitotic apparatus",
    "Basolateral cell membrane" = "Plasma membrane",
    "COPI-coated vesicle" = "Vesicles",
    "Cajal body" = "Nucleus",
    "Early endosome membrane" = "Endosomes",
    "Endoplasmic reticulum lumen" = "Endoplasmic reticulum",
    "Endoplasmic reticulum-Golgi intermediate compartment membrane" = "Endoplasmic reticulum",
    "Golgi stack" = "Golgi apparatus",
    "Late endosome" = "Endosomes",
    "Late endosome membrane" = "Endosomes",
    "M line" = "Actin cytoskeleton",
    "Melanosome" = "Others",
    "Nucleus speckle" = "Nucleus",
    "PML body" = "Nucleus",
    "Perikaryon" = "Others",
    "Peroxisome membrane" = "Peroxisomes",
    "Photoreceptor inner segment" = "Others",
    "Presynaptic cell membrane" = "Plasma membrane",
    "Sarcoplasmic reticulum membrane" = "Endoplasmic reticulum",
    "Single-pass type IV membrane protein" = "Others",
    "Z line" = "Actin cytoskeleton",
    "adherens junction" = "Cell Junctions",
    "autophagosome" = "Others",
    "caveola" = "Plasma membrane",
    "cell cortex" = "Actin cytoskeleton",
    "centriole" = "Microtubules",
    "cilium axoneme" = "Primary cilium",
    "cilium basal body" = "Primary cilium",
    "cilium membrane" = "Primary cilium",
    "cis-Golgi network membrane" = "Golgi apparatus",
    "clathrin-coated pit" = "Vesicles",
    "clathrin-coated vesicle" = "Vesicles",
    "focal adhesion" = "Actin cytoskeleton",
    "growth cone" = "Others",
    "interphotoreceptor matrix" = "Extracellular",
    "lamellipodium" = "Actin cytoskeleton",
    "nucleoplasm" = "Nucleus",
    "peroxisome" = "Peroxisomes",
    "phagocytic cup" = "Others",
    "photoreceptor outer segment" = "Others",
    "pseudopodium" = "Actin cytoskeleton",
    "ruffle membrane" = "Plasma membrane",
    "synaptic vesicle membrane" = "Vesicles",
    "telomere" = "Chromosome",
    "translation" = "Others"
  )
  
  originals = names(gene_sets_cleaned)
  
  map_keys = names(canonical_map)
  missing_keys = setdiff(map_keys, originals)
  if (length(missing_keys) > 0) {
    message("Term unification error",
            paste(head(missing_keys, 50), collapse = ", "),
            if (length(missing_keys) > 50) " ...")
  }
  
  merged = list()
  
  for (orig in originals) {
    key_match = which(trimws(map_keys) == trimws(orig))
    genes_orig = unique(na.omit(gene_sets_cleaned[[orig]]))
    if (length(key_match) == 1) {
      target = canonical_map[[ map_keys[key_match] ]]
      merged[[target]] = unique(c(merged[[target]] %||% character(0), genes_orig))
    } else {
      merged[[orig]] = unique(c(merged[[orig]] %||% character(0), genes_orig))
    }
  }
  
  idx_map = match(originals, map_keys)
  mapping_rows = tibble(
    original = originals,
    n_before = sapply(originals, function(x) length(gene_sets_cleaned[[x]])),
    canonical = vapply(seq_along(originals), function(i) {
      if (!is.na(idx_map[i])) return(canonical_map[[ map_keys[idx_map[i]] ]])
      return(originals[i])
    }, FUN.VALUE = character(1))
  ) %>%
    mutate(
      n_after = sapply(canonical, function(x) if (!is.null(merged[[x]])) length(merged[[x]]) else 0),
      example_genes = sapply(original, function(o) paste(head(gene_sets_cleaned[[o]] %||% character(0), 5), collapse = "; "))
    )
  
  gene_sets_final = merged
  
  return(list(
    gene_sets = gene_sets_final,
    raw_terms_df = raw_terms_df,
    protein_ids = protein_ids
  ))
}



translatorDb = function(geneList2translate){
  contrastMsUni_MsEntrez = translate_ID(geneList2translate, "UNIPROT", "ENTREZID", species_db = org.Mm.eg.db,
                                        manual_map = translatorMsUni2MsEn)
  
  NA_IDs = !is.na(contrastMsUni_MsEntrez$`10090_ENTREZID`)
  prev_original = contrastMsUni_MsEntrez$original[NA_IDs]
  prev_10090_ENTREZID = contrastMsUni_MsEntrez$`10090_ENTREZID`[NA_IDs]
  contrastMsUni_MsEntrez_noNA = data.frame(original = prev_original)
  contrastMsUni_MsEntrez_noNA$`10090_ENTREZID` = prev_10090_ENTREZID
  
  NULL_IDs = contrastMsUni_MsEntrez_noNA$`10090_ENTREZID` != ""
  prev_original2 = contrastMsUni_MsEntrez_noNA$original[NULL_IDs]
  prev_10090_ENTREZID2 = contrastMsUni_MsEntrez_noNA$`10090_ENTREZID`[NULL_IDs]
  contrastMsUni_MsEntrez_noNULL = data.frame(original = prev_original2)
  contrastMsUni_MsEntrez_noNULL$`10090_ENTREZID` = prev_10090_ENTREZID2
  
  contrastMsEntrez_HuEntrez = translate_ID(contrastMsUni_MsEntrez_noNULL$`10090_ENTREZID`, "ENTREZID", "ENTREZID",
                                           source_taxon = 10090, target_taxon = 9606, manual_map_inter = translatorMs2Hu)
  
  contrastHuEntrez_HuSymbol = translate_ID(contrastMsEntrez_HuEntrez$`9606_ENTREZID`, "ENTREZID", "SYMBOL", species_db = org.Hs.eg.db,
                                           manual_map = translatorMsUni2MsEn)
  
  contrastHuEntrez_HuUniprot = translate_ID(contrastMsEntrez_HuEntrez$`9606_ENTREZID`, "ENTREZID", "UNIPROT", species_db = org.Hs.eg.db)
  
  contrastMsEntrez_MsSymbol = translate_ID(contrastMsUni_MsEntrez_noNULL$`10090_ENTREZID`, "ENTREZID", "SYMBOL", species_db = org.Mm.eg.db,
                                           manual_map = translatorMsUni2MsEn)
  
  contrastsTranslator = data.frame(`10090_UNIPROT` = contrastMsUni_MsEntrez_noNULL$original, `10090_ENTREZID` =contrastMsUni_MsEntrez_noNULL$`10090_ENTREZID`,
                                   `10090_SYMBOL` = contrastMsEntrez_MsSymbol$`10090_SYMBOL`, `9606_ENTREZID` = contrastMsEntrez_HuEntrez$`9606_ENTREZID`,
                                   `9606_SYMBOL` = contrastHuEntrez_HuSymbol$`9606_SYMBOL`, `9606_UNIPROT` = contrastHuEntrez_HuUniprot$`9606_UNIPROT`, 
                                   check.names = FALSE)
  return(contrastsTranslator)
}


merge_terms = function(gene_sets, terms_to_merge, new_name) {
  merged_genes = unique(unlist(gene_sets[terms_to_merge]))
  gene_sets[terms_to_merge] = NULL
  gene_sets[[new_name]] = merged_genes
  return(gene_sets)
}


translatorHumanUni = function(geneList2translate){
  contrastMsUni_MsEntrez = translate_ID(geneList2translate, "UNIPROT", "ENTREZID", species_db = org.Hs.eg.db,
                                        manual_map = translatorMsUni2MsEn)
  
  NA_IDs = !is.na(contrastMsUni_MsEntrez$`9606_ENTREZID`)
  prev_original = contrastMsUni_MsEntrez$original[NA_IDs]
  prev_9606_ENTREZID = contrastMsUni_MsEntrez$`9606_ENTREZID`[NA_IDs]
  contrastMsUni_MsEntrez_noNA = data.frame(original = prev_original)
  contrastMsUni_MsEntrez_noNA$`9606_ENTREZID` = prev_9606_ENTREZID
  
  NULL_IDs = contrastMsUni_MsEntrez_noNA$`9606_ENTREZID` != ""
  prev_original2 = contrastMsUni_MsEntrez_noNA$original[NULL_IDs]
  prev_9606_ENTREZID2 = contrastMsUni_MsEntrez_noNA$`9606_ENTREZID`[NULL_IDs]
  contrastMsUni_MsEntrez_noNULL = data.frame(original = prev_original2)
  contrastMsUni_MsEntrez_noNULL$`9606_ENTREZID` = prev_9606_ENTREZID2
  
  contrastHuEntrez_MsEntrez = translate_ID(contrastMsUni_MsEntrez_noNULL$`9606_ENTREZID`, "ENTREZID", "ENTREZID",
                                           source_taxon = 9606, target_taxon = 10090)
  
  contrastMsEntrez_MsSymbol = translate_ID(contrastHuEntrez_MsEntrez$`10090_ENTREZID`, "ENTREZID", "SYMBOL", species_db = org.Mm.eg.db)
  
  contrastMsEntrez_MsuUniprot = translate_ID(contrastHuEntrez_MsEntrez$`10090_ENTREZID`, "ENTREZID", "UNIPROT", species_db = org.Mm.eg.db)
  
  contrastHuEntrez_HuSymbol = translate_ID(contrastMsUni_MsEntrez_noNULL$`9606_ENTREZID`, "ENTREZID", "SYMBOL", species_db = org.Hs.eg.db)
  
  contrastsTranslator = data.frame(`10090_UNIPROT` = contrastMsEntrez_MsuUniprot$`10090_UNIPROT`, `10090_ENTREZID` = contrastHuEntrez_MsEntrez$`10090_ENTREZID`,
                                   `10090_SYMBOL` = contrastMsEntrez_MsSymbol$`10090_SYMBOL`, `9606_ENTREZID` = contrastMsUni_MsEntrez_noNULL$`9606_ENTREZID`,
                                   `9606_SYMBOL` = contrastHuEntrez_HuSymbol$`9606_SYMBOL`, `9606_UNIPROT` = contrastMsUni_MsEntrez_noNULL$original, 
                                   check.names = FALSE)
  
  return(contrastsTranslator)
}


translatorHumanSym = function(geneList2translate){
  contrastMsUni_MsEntrez = translate_ID(geneList2translate, "SYMBOL", "ENTREZID", species_db = org.Hs.eg.db,
                                        manual_map = translatorMsUni2MsEn)
  
  NA_IDs = !is.na(contrastMsUni_MsEntrez$`9606_ENTREZID`)
  prev_original = contrastMsUni_MsEntrez$original[NA_IDs]
  prev_9606_ENTREZID = contrastMsUni_MsEntrez$`9606_ENTREZID`[NA_IDs]
  contrastMsUni_MsEntrez_noNA = data.frame(original = prev_original)
  contrastMsUni_MsEntrez_noNA$`9606_ENTREZID` = prev_9606_ENTREZID
  
  NULL_IDs = contrastMsUni_MsEntrez_noNA$`9606_ENTREZID` != ""
  prev_original2 = contrastMsUni_MsEntrez_noNA$original[NULL_IDs]
  prev_9606_ENTREZID2 = contrastMsUni_MsEntrez_noNA$`9606_ENTREZID`[NULL_IDs]
  contrastMsUni_MsEntrez_noNULL = data.frame(original = prev_original2)
  contrastMsUni_MsEntrez_noNULL$`9606_ENTREZID` = prev_9606_ENTREZID2
  
  contrastHuEntrez_MsEntrez = translate_ID(contrastMsUni_MsEntrez_noNULL$`9606_ENTREZID`, "ENTREZID", "ENTREZID",
                                           source_taxon = 9606, target_taxon = 10090)
  
  contrastMsEntrez_MsSymbol = translate_ID(contrastHuEntrez_MsEntrez$`10090_ENTREZID`, "ENTREZID", "SYMBOL", species_db = org.Mm.eg.db)
  
  contrastMsEntrez_MsuUniprot = translate_ID(contrastHuEntrez_MsEntrez$`10090_ENTREZID`, "ENTREZID", "UNIPROT", species_db = org.Mm.eg.db)
  
  contrastHuEntrez_HuSymbol = translate_ID(contrastMsUni_MsEntrez_noNULL$`9606_ENTREZID`, "ENTREZID", "UNIPROT", species_db = org.Hs.eg.db)
  
  contrastsTranslator = data.frame(`10090_UNIPROT` = contrastMsEntrez_MsuUniprot$`10090_UNIPROT`, `10090_ENTREZID` = contrastHuEntrez_MsEntrez$`10090_ENTREZID`,
                                   `10090_SYMBOL` = contrastMsEntrez_MsSymbol$`10090_SYMBOL`, `9606_ENTREZID` = contrastMsUni_MsEntrez_noNULL$`9606_ENTREZID`,
                                   `9606_SYMBOL` = contrastMsUni_MsEntrez_noNULL$original, `9606_UNIPROT` = contrastHuEntrez_HuSymbol$`9606_UNIPROT`, 
                                   check.names = FALSE)
  
  return(contrastsTranslator)
}


translatorMitoCarta = function(geneDf2translate){
  contrastHuEntrez_HuSym = translate_ID(geneDf2translate$`9606_ENTREZID`, "ENTREZID", "SYMBOL", species_db = org.Hs.eg.db,
                                        manual_map = translatorMs2Hu)
  
  contrastHuEntrez_HuUni = translate_ID(geneDf2translate$`9606_ENTREZID`, "ENTREZID", "UNIPROT", species_db = org.Hs.eg.db)
  
  contrastMsEntrez_MsUni = translate_ID(geneDf2translate$`10090_ENTREZID`, "ENTREZID", "UNIPROT", species_db = org.Mm.eg.db)
  
  contrastsTranslator = data.frame(`10090_UNIPROT` = contrastMsEntrez_MsUni$`10090_UNIPROT`, `10090_ENTREZID` = geneDf2translate$`10090_ENTREZID`,
                                   `10090_SYMBOL` = geneDf2translate$`10090_SYMBOL`, `9606_ENTREZID` = geneDf2translate$`9606_ENTREZID`,
                                   `9606_SYMBOL` = contrastHuEntrez_HuSym$`9606_SYMBOL`, `9606_UNIPROT` = contrastHuEntrez_HuUni$`9606_UNIPROT`, 
                                   check.names = FALSE)
  
  return(contrastsTranslator)
}

createGeneralCellOriginDb <- function(file_cell, 
                                      proteinList = NULL,
                                      contrastsTranslator, 
                                      nTPM_threshold = 2, 
                                      increase = 2.0,
                                      translatorHumanUniprot = "9606_UNIPROT",
                                      translatorMouseEntrez = "10090_ENTREZID",
                                      translatorMouseUniprot = "10090_UNIPROT",
                                      translatorHumanSymbol = "9606_SYMBOL") {
  
  df <- read.delim(file_cell, stringsAsFactors = FALSE)
  
  df <- df %>% dplyr::distinct(Gene.name, Cluster.type, .keep_all = TRUE)
  
  df_processed <- df %>%
    dplyr::group_by(Gene.name) %>%
    dplyr::mutate(mean_nTPM = mean(nTPM, na.rm = TRUE)) %>%
    dplyr::ungroup() %>%
    dplyr::filter(nTPM >= nTPM_threshold)
  
  if (!is.null(increase) && increase > 0) {
    df_processed <- df_processed %>%
      dplyr::filter(nTPM >= (mean_nTPM * increase))
  }
  
  map_df <- contrastsTranslator %>%
    dplyr::select(dplyr::all_of(c(translatorHumanSymbol, translatorMouseEntrez))) %>%
    dplyr::filter(!is.na(.data[[translatorHumanSymbol]]), 
                  !is.na(.data[[translatorMouseEntrez]]),
                  .data[[translatorMouseEntrez]] != "") %>%
    dplyr::distinct(.data[[translatorHumanSymbol]], .keep_all = TRUE) 
  
  colnames(map_df) <- c("Gene.name", "mEntrez")
  
  df_mapped <- df_processed %>%
    dplyr::inner_join(map_df, by = "Gene.name")
  
  if (!is.null(proteinList)) {
    valid_entrez <- contrastsTranslator %>%
      dplyr::filter(.data[[translatorMouseUniprot]] %in% proteinList) %>%
      dplyr::pull(dplyr::all_of(translatorMouseEntrez)) %>%
      na.omit() %>% unique()
    
    df_mapped <- df_mapped %>%
      dplyr::filter(mEntrez %in% valid_entrez)
  }
  
  gene_sets <- df_mapped %>%
    dplyr::group_by(Cluster.type) %>%
    dplyr::summarise(genes = list(unique(as.character(na.omit(mEntrez)))), .groups = "drop") %>%
    tibble::deframe()
  
  return(gene_sets)
}