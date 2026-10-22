# Modification based on the DEP repository  
# Original function from DEP by Arne Smits: https://github.com/arnesmits/DEP  
# Changes made in this version:  
# - Improved validation method for missing columns (`all(!indicate %in% columns)` →
#   `!any(indicate %in% columns)`) to make the condition more intuitive  
# - Optimized column selection using `dplyr::select(all_of(indicate))` to avoid warnings in `tidyselect`  
# - Enhanced annotation color assignment logic for better visualization  

get_annotation <- function(dep, indicate, custom_anno_cols = NULL, custom_anno_labels = NULL) {
  # Validate inputs
  assertthat::assert_that(
    inherits(dep, "SummarizedExperiment"),
    is.character(indicate)
  )
  
  # Extract column data as a data frame
  col_data <- colData(dep) %>%
    as.data.frame()
  columns <- colnames(col_data)
  
  # Change: Improved column existence validation using `!any()` for better readability and efficiency
  if (!any(indicate %in% columns)) {  
    stop("'",
         paste0(indicate, collapse = "' and/or '"),
         "' column(s) is/are not present in ",
         deparse(substitute(dep)),
         ".\nValid columns are: '",
         paste(columns, collapse = "', '"),
         "'.",
         call. = FALSE)
  }
  
  # Change: Warning refinement when some columns exist and others do not  
  if (any(!indicate %in% columns)) {  
    indicate <- indicate[indicate %in% columns]
    warning("Only used the following indicate column(s): '",
            paste0(indicate, collapse = "', '"),
            "'")
  }
  
  # Change: Optimized column selection using `dplyr::select(all_of(indicate))` to ensure compatibility with latest tidyselect standards
  anno <- dplyr::select(col_data, all_of(indicate))  
  
  # Generate annotation colors  
  names <- colnames(anno)
  anno_col <- vector(mode = "list", length = length(names))
  names(anno_col) <- names
  
  for (i in names) {
    if (!is.null(custom_anno_cols) && i %in% names(custom_anno_cols)) {
      anno_col[[i]] <- custom_anno_cols[[i]]
    } else {
      var <- anno[[i]] %>% unique() %>% sort()
      if (length(var) == 1) {
        cols <- c("black")
      } else if (length(var) == 2) {
        cols <- c("#2E6F40", "#4F0433")
      } else if (length(var) < 7) {
        cols <- RColorBrewer::brewer.pal(length(var), "Pastel1")
      } else {
        cols <- colorRampPalette(RColorBrewer::brewer.pal(12, "Set3"))(length(var))
      }
      names(cols) <- var
      anno_col[[i]] <- cols
    }
  }
  
  legend_params <- list()
  for (i in names) {
    l_param <- list(
      title_gp = grid::gpar(fontsize = 0, fontface = "bold"), 
      labels_gp = grid::gpar(fontsize = 18),                   
      grid_height = grid::unit(8, "mm"),                       
      grid_width  = grid::unit(8, "mm")                     
    )
    if (!is.null(custom_anno_labels) && i %in% names(custom_anno_labels)) {
      l_param$at <- names(custom_anno_labels[[i]])
      l_param$labels <- unname(custom_anno_labels[[i]])
    }
    legend_params[[i]] <- l_param
  }
  
  # Change: Return updated annotation for heatmap  
  HeatmapAnnotation(
    df = anno,
    col = anno_col,
    show_annotation_name = TRUE,
    annotation_legend_param = legend_params
  )
}

# Modification based on the DEP repository  
# Original function from DEP by Arne Smits: https://github.com/arnesmits/DEP  
# Changes made in this version:  
# - Explicit conversion of integer values to numeric (using as.numeric()) to prevent type-related issues  
# - Improved handling of missing values (obs_NA logic) to use 'gower' distance if any NAs are present  
# - Optimized k-means clustering by using the row_km parameter when kmeans = TRUE  
# - Enhanced integration of annotations using get_annotation() for smoother visualization  
# - Added option to limit the number of columns (col_limit) displayed in the heatmap  
# - Adjusted clustering logic by passing clustering_distance to both rows and columns  
# - Moved row and column font size settings into the Heatmap() call via row_names_gp and column_names_gp  
# - The 'plot' argument now controls whether the heatmap is drawn immediately (using draw()) or returned 
# - Add a set.sedd function
plot_heatmap <- function(dep, type = c("contrast", "centered"),
                         kmeans = FALSE, k = 6,
                         col_limit = 6, indicate = NULL,
                         clustering_distance = c("euclidean", "maximum", "manhattan", "canberra",
                                                 "binary", "minkowski", "pearson", "spearman", "kendall", "gower"),
                         row_font_size = 3, col_font_size = 5, plot = TRUE, seed,
                         custom_anno_cols = NULL, custom_anno_labels = NULL, ...) {
  
  # Explicitly convert integer values to numeric to avoid potential type-related issues  
  if (is.integer(k)) k <- as.numeric(k)
  if (is.integer(col_limit)) col_limit <- as.numeric(col_limit)
  if (is.integer(row_font_size)) row_font_size <- as.numeric(row_font_size)
  if (is.integer(col_font_size)) col_font_size <- as.numeric(col_font_size)
  
  # Change: Ensure valid input values
  assertthat::assert_that(inherits(dep, "SummarizedExperiment"),
                          is.character(type),
                          is.logical(kmeans),
                          is.numeric(k), length(k) == 1,
                          is.numeric(col_limit), length(col_limit) == 1,
                          is.numeric(row_font_size), length(row_font_size) == 1,
                          is.numeric(col_font_size), length(col_font_size) == 1,
                          is.logical(plot), length(plot) == 1)
  
  set.seed(seed)
  
  type <- match.arg(type)
  clustering_distance <- match.arg(clustering_distance)
  
  row_data <- rowData(dep, use.names = FALSE)
  col_data <- as.data.frame(colData(dep))
  
  # Change: Improved annotation integration using get_annotation() to ensure smooth visualization compatibility  
  if (!is.null(indicate) & type == "contrast") {
    warning("Heatmap annotation only applicable for type = 'centered'", call. = FALSE)
  }
  if (!is.null(indicate) & type == "centered") {
    ha1 <- get_annotation(dep, indicate, custom_anno_cols, custom_anno_labels)
  } else {
    ha1 <- NULL
  }
  
  # Change: Improved handling of missing values (obs_NA logic) to prevent clustering errors  
  filtered <- dep[row_data$significant, ]
  if (any(is.na(assay(filtered)))) {
    warning("Missing values in '", deparse(substitute(dep)), "'. Using clustering_distance = 'gower'", 
            call. = FALSE)
    clustering_distance <- "gower"
    obs_NA <- TRUE  
  } else {
    obs_NA <- FALSE
  }
  
  # Data transformation for heatmap representation
  if (type == "centered") {
    rowData(filtered)$mean <- rowMeans(assay(filtered), na.rm = TRUE)
    df <- assay(filtered) - rowData(filtered, use.names = FALSE)$mean
  }
  if (type == "contrast") {
    df <- rowData(filtered, use.names = FALSE) %>%
      data.frame() %>%
      column_to_rownames(var = "name") %>%
      dplyr::select(ends_with("_diff"))
    colnames(df) <- colnames(df) %>% 
      gsub("_diff", "", .) %>% 
      gsub("_vs_", " vs ", .)
    df <- as.matrix(df)
  }
  
  # Limit number of columns if the matrix is wider than col_limit  
  if (ncol(df) > col_limit) {
    df <- df[, seq_len(col_limit)]
  }
  
  # Format the column label names
  raw_cols <- colnames(df)
  math_cols <- paste0("'", gsub("_+", "'~'", raw_cols), "'")
  math_cols <- gsub("'APPKI'", "'APP'^'KI'", math_cols)
  parsed_cols <- parse(text = math_cols)
  # ---------------------------------------------------------
  
  # Build the heatmap using ComplexHeatmap. 
  # Pass the clustering distance for both rows and columns and use row_km if kmeans is enabled.
  # Set row and column font sizes via row_names_gp and column_names_gp.
  ht <- Heatmap(df,
                top_annotation = ha1,
                column_labels = parsed_cols,
                clustering_distance_rows = clustering_distance,
                clustering_distance_columns = clustering_distance,
                row_km = if (kmeans) k else NULL,
                row_names_gp = grid::gpar(fontsize = row_font_size),
                column_names_gp = grid::gpar(fontsize = col_font_size),
                heatmap_legend_param = list(title = ""),
                ...)
  
  # If plot is TRUE, draw the heatmap immediately; otherwise return the object
  if (plot) {
    draw(ht)
  }
  return(ht)
}


# Modification based on the DEP repository  
# Original function from DEP by Arne Smits: https://github.com/arnesmits/DEP  
# Changes made in this version:  
# - Improved validation method for required columns (`"name"` and `"ID"`)  
# - Optimized numeric column validation for better error handling  
# - Ensured compatibility with tibble objects by converting `expdesign` to a data frame  
# - Enhanced experimental condition creation using `paste()`  
# - Improved column name matching for assay data  

make_se <- function(proteins_unique, columns, expdesign, varNames) {
  # Change: Improved validation for required columns (`"name"` and `"ID"`) 
  assertthat::assert_that(is.data.frame(proteins_unique), is.integer(columns), is.data.frame(expdesign))
  
  if(any(!c("name", "ID") %in% colnames(proteins_unique))) {
    stop("'name' and/or 'ID' are not present in 'proteins_unique'. Run `make_unique()` to obtain the required columns", call. = FALSE)
  }
  
  # Change: Optimized numeric column validation for better error handling
  if(any(!apply(proteins_unique[, columns], 2, is.numeric))) {
    stop("The specified 'columns' must be numeric. Run `make_se_parse()` with the appropriate columns as an argument", call. = FALSE)
  }
  
  # Change: Ensured compatibility with tibble objects by converting `expdesign` to a data frame 
  if(tibble::is_tibble(expdesign)) expdesign <- as.data.frame(expdesign)
  
  # Change: Enhanced experimental condition creation using `paste()`
  varNames_condition <- setdiff(varNames, "replicate")
  if(length(varNames_condition) >= 1){
    expdesign$condition <- do.call(paste, c(expdesign[varNames_condition], sep = "_"))
  } else {
    message("Only 'replicate' provided in varNames; 'condition' remains unchanged.")
  }
  
  expdesign$ID <- expdesign$label
  rownames(expdesign) <- expdesign$ID
  
  rownames(proteins_unique) <- proteins_unique$name
  raw <- proteins_unique[, columns]
  
  # Change: Improved handling of zero values by converting them to NA before log2 transformation
  raw[raw == 0] <- NA
  raw <- log2(raw)
  
  # Change: Improved column name matching for assay data 
  matched <- match(make.names(expdesign$label), make.names(colnames(raw)))
  colnames(raw)[matched] <- expdesign$label
  raw <- raw[, !is.na(colnames(raw))]
  raw <- raw[, rownames(expdesign)]
  
  row_data <- proteins_unique[, -columns]
  rownames(row_data) <- row_data$name
  
  # Create SummarizedExperiment object  
  se <- SummarizedExperiment(assays = as.matrix(raw), colData = expdesign, rowData = row_data)
  return(se)
}

add_rejections <- function (diff, alpha = 0.05, lfc = 1) {
  if (is.integer(alpha)) 
    alpha <- as.numeric(alpha)
  if (is.integer(lfc)) 
    lfc <- as.numeric(lfc)
  assertthat::assert_that(inherits(diff, "SummarizedExperiment"), 
                          is.numeric(alpha), length(alpha) == 1, is.numeric(lfc), 
                          length(lfc) == 1)
  row_data <- rowData(diff, use.names = FALSE) %>% as.data.frame()
  if (any(!c("name", "ID") %in% colnames(row_data))) {
    stop("'name' and/or 'ID' columns are not present in '", 
         deparse(substitute(diff)), "'\nRun make_unique() and make_se() to obtain the required columns", 
         call. = FALSE)
  }
  if (length(grep("_p.adj|_diff", colnames(row_data))) < 1) {
    stop("'[contrast]_diff' and/or '[contrast]_p.adj' columns are not present in '", 
         deparse(substitute(diff)), "'\nRun test_diff() to obtain the required columns", 
         call. = FALSE)
  }
  cols_p <- grep("_p.adj", colnames(row_data))
  cols_diff <- grep("_diff", colnames(row_data))
  if (length(cols_p) == 1) {
    rowData(diff)$significant <- row_data[, cols_p] <= alpha & 
      abs(row_data[, cols_diff]) >= lfc
    rowData(diff)$contrast_significant <- rowData(diff, use.names = FALSE)$significant
    colnames(rowData(diff))[ncol(rowData(diff, use.names = FALSE))] <- gsub("p.adj", 
                                                                            "significant", colnames(row_data)[cols_p])
  }
  if (length(cols_p) > 1) {
    p_reject <- row_data[, cols_p, drop = FALSE] <= alpha
    diff_reject <- abs(row_data[, cols_diff, drop = FALSE]) >= lfc
    
    p_reject[is.na(p_reject)] <- FALSE
    diff_reject[is.na(diff_reject)] <- FALSE
    
    sign_df <- as.data.frame(p_reject & diff_reject)
    colnames(sign_df) <- gsub("_p.adj", "_significant", colnames(sign_df))
    
    sign_df$significant <- apply(sign_df, 1, any)
    
    for (col_name in colnames(sign_df)) {
      rowData(diff)[[col_name]] <- sign_df[[col_name]]
    }
  }
  return(diff)
}

# Modification based on the DEP repository  
# Original function from DEP by Arne Smits: https://github.com/arnesmits/DEP  
# Changes made in this version:
# - Added a new branch for type = "multiple" that incorporates a custom design formula, 
#  contrast specification and external differential analysis.
# - Constructs the design matrix using the externally defined formula and cleans its column names with make.names().
# - Preserves the original DEP behavior for type values "control", "all", and "manual".

test_diff <- function(se,
                      type = c("control", "all", "manual", "multiple"),
                      control = NULL, 
                      test = NULL,
                      design_formula = formula(~ 0 + condition),
                      contrastFormula = NULL,
                      contrastList = NULL) {
  
  type <- match.arg(type)
  
  # Validate that se is a SummarizedExperiment
  assertthat::assert_that(inherits(se, "SummarizedExperiment"))
  
  # Extract expression data, design information and variables
  col_data <- colData(se)
  raw <- assay(se)
  
  if(type == "multiple") {
    # Ensure that external design and contrast arguments are provided
    if(is.null(contrastFormula) || is.null(contrastList))
      stop("For type = 'multiple', 'contrastFormula', and 'contrastList' must be provided.")
    
    # Convert the externally defined 'design_formula' into a design formula and build the design matrix
    design_formula_ext <- as.formula(design_formula)
    design_matrix <- model.matrix(design_formula_ext, data = col_data)
    colnames(design_matrix) <- make.names(colnames(design_matrix))
    
    # Fit the linear model using limma
    fit <- limma::lmFit(raw, design_matrix)
    
    # Build the contrast matrix and fit contrasts
    contrast_matrix <- do.call(limma::makeContrasts, c(contrastFormula, list(levels = design_matrix)))
    fitCont <- limma::contrasts.fit(fit, contrast_matrix)
    fitStat <- limma::eBayes(fitCont)
    
    # Internal function to retrieve differential expression results per contrast
    retrieve_fun <- function(comp, fit = fitStat) {
      res <- limma::topTable(fit, sort.by = "t", coef = comp,
                             number = Inf, confint = TRUE)
      # Remove rows with missing t statistics
      res <- res[!is.na(res$t), , drop = FALSE]
      
      # Estimate FDR and local FDR using fdrtool
      fdr_res <- fdrtool::fdrtool(res$P.Value, statistic = "pvalue", plot = FALSE, verbose = FALSE)
      res$qval <- fdr_res$qval
      res$lfdr <- fdr_res$lfdr
      
      # Include the contrast name in the results
      res$comparison <- comp
      res <- tibble::rownames_to_column(as.data.frame(res), var = "rowname")
      return(res)
    }
    
    # Retrieve and combine results for each contrast in contrastList
    limma_res <- purrr::map_df(contrastList, retrieve_fun)
    
    # Reshape the results table
    table <- limma_res %>%
      dplyr::select(rowname, logFC, CI.L, CI.R, P.Value, qval, comparison) %>%
      tidyr::gather(variable, value, -c(rowname, comparison)) %>%
      dplyr::mutate(variable = dplyr::recode(variable,
                                             logFC = "diff",
                                             P.Value = "p.val",
                                             qval = "p.adj")) %>%
      tidyr::unite(temp, comparison, variable) %>%
      tidyr::spread(temp, value)
    
    # Merge the results table with rowData(se) using the "name" column as key
    rowData(se) <- merge(rowData(se, use.names = FALSE),
                         table,
                         by.x = "name", by.y = "rowname",
                         all.x = TRUE, sort = FALSE)
    
    return(se)
    
  } else {
    # Original DEP branch: using default contrast generation
    if(!is.null(control)) {
      # Show error if control input is not valid
      assertthat::assert_that(is.character(control),
                              length(control) == 1)
      if(!control %in% unique(col_data$condition)) {
        stop("run test_diff() with a valid control.\nValid controls are: '",
             paste0(unique(col_data$condition), collapse = "', '"), "'",
             call. = FALSE)
      }
    }
    
    col_data <- colData(se)
    raw <- assay(se)
    
    # variables in formula
    variables <- terms.formula(design_formula) %>%
      attr(., "variables") %>%
      as.character() %>%
      .[-1]
    
    # Throw error if variables are not col_data columns
    if(any(!variables %in% colnames(col_data))) {
      stop("run make_diff() with an appropriate 'design_formula'")
    }
    if(variables[1] != "condition") {
      stop("first factor of 'design_formula' should be 'condition'")
    }
    
    # Obtain variable factors
    for(var in variables) {
      temp <- factor(col_data[[var]])
      assign(var, temp)
    }
    
    # Make an appropriate design matrix
    design <- model.matrix(design_formula, data = environment())
    colnames(design) <- gsub("condition", "", colnames(design))
    
    # Generate contrasts to be tested
    # Either make all possible combinations ("all"),
    # only the contrasts versus the control sample ("control") or
    # use manual contrasts
    conditions <- as.character(unique(condition))
    if(type == "all") {
      # All possible combinations
      cntrst <- apply(utils::combn(conditions, 2), 2, paste, collapse = " - ")
      
      if(!is.null(control)) {
        # Make sure that contrast containing
        # the control sample have the control as denominator
        flip <- grep(paste("^", control, sep = ""), cntrst)
        if(length(flip) >= 1) {
          cntrst[flip] <- cntrst[flip] %>%
            gsub(paste(control, "- ", sep = " "), "", .) %>%
            paste(" - ", control, sep = "")
        }
      }
      
    }
    if(type == "control") {
      # Throw error if no control argument is present
      if(is.null(control))
        stop("run test_diff(type = 'control') with a 'control' argument")
      
      # Make contrasts
      cntrst <- paste(conditions[!conditions %in% control],
                      control,
                      sep = " - ")
    }
    if(type == "manual") {
      # Throw error if no test argument is present
      if(is.null(test)) {
        stop("run test_diff(type = 'manual') with a 'test' argument")
      }
      assertthat::assert_that(is.character(test))
      
      if(any(!unlist(strsplit(test, "_vs_")) %in% conditions)) {
        stop("run test_diff() with valid contrasts in 'test'",
             ".\nValid contrasts should contain combinations of: '",
             paste0(conditions, collapse = "', '"),
             "', for example '", paste0(conditions[1], "_vs_", conditions[2]),
             "'.", call. = FALSE)
      }
      
      cntrst <- gsub("_vs_", " - ", test)
      
    }
    # Print tested contrasts
    message("Tested contrasts: ",
            paste(gsub(" - ", "_vs_", cntrst), collapse = ", "))
    
    # Test for differential expression by empirical Bayes moderation
    # of a linear model on the predefined contrasts
    fit <- lmFit(raw, design = design)
    made_contrasts <- makeContrasts(contrasts = cntrst, levels = design)
    contrast_fit <- contrasts.fit(fit, made_contrasts)
    
    if(any(is.na(raw))) {
      for(i in cntrst) {
        covariates <- strsplit(i, " - ") %>% unlist
        single_contrast <- makeContrasts(contrasts = i, levels = design[, covariates])
        single_contrast_fit <- contrasts.fit(fit[, covariates], single_contrast)
        contrast_fit$coefficients[, i] <- single_contrast_fit$coefficients[, 1]
        contrast_fit$stdev.unscaled[, i] <- single_contrast_fit$stdev.unscaled[, 1]
      }
    }
    
    eB_fit <- eBayes(contrast_fit)
    
    # Function to retrieve the results of
    # the differential expression test using 'fdrtool'
    retrieve_fun <- function(comp, fit = eB_fit){
      res <- topTable(fit, sort.by = "t", coef = comp,
                      number = Inf, confint = TRUE)
      res <- res[!is.na(res$t),]
      fdr_res <- fdrtool::fdrtool(res$t, plot = FALSE, verbose = FALSE)
      res$qval <- fdr_res$qval
      res$lfdr <- fdr_res$lfdr
      res$comparison <- rep(comp, dim(res)[1])
      res <- rownames_to_column(res)
      return(res)
    }
    
    # Retrieve the differential expression test results
    limma_res <- map_df(cntrst, retrieve_fun)
    
    # Select the logFC, CI and qval variables
    table <- limma_res %>%
      dplyr::select(rowname, logFC, CI.L, CI.R, P.Value, qval, comparison) %>%
      mutate(comparison = gsub(" - ", "_vs_", comparison)) %>%
      gather(variable, value, -c(rowname,comparison)) %>%
      mutate(variable = recode(variable, logFC = "diff", P.Value = "p.val", qval = "p.adj")) %>%
      unite(temp, comparison, variable) %>%
      spread(temp, value)
    rowData(se) <- merge(rowData(se, use.names = FALSE), table,
                         by.x = "name", by.y = "rowname", all.x = TRUE, sort=FALSE)
    return(se)
  }
}

# Modification based on the DEP repository  
# Original function from DEP by Arne Smits: https://github.com/arnesmits/DEP  
# Changes made in this version:
# - Explicit conversion of label_size to numeric if provided as integer.
# - Validation for the presence of required columns ("name", "ID", "[contrast]_diff", 
#   "[contrast]_p.adj"/"[contrast]_p.val", and "[contrast]_significant") in rowData.
# - Enhanced error messages including a suggestion of valid contrasts when needed.
# - Construction of a tidy data frame for the volcano plot and conditional addition of protein names.
# - Dynamic y-axis labeling depending on whether p-values are adjusted.
plot_volcano <- function(dep, contrast, label_size = 3, add_names = TRUE, adjusted = FALSE,
                         plot = TRUE, proteins_to_label = NULL) {
  # Change: Convert label_size to numeric if provided as integer**  
  if (is.integer(label_size)) label_size <- as.numeric(label_size)
  
  # Validate input types and lengths
  assertthat::assert_that(
    inherits(dep, "SummarizedExperiment"),
    is.character(contrast),
    length(contrast) == 1,
    is.numeric(label_size), length(label_size) == 1,
    is.logical(add_names), length(add_names) == 1,
    is.logical(adjusted), length(adjusted) == 1,
    is.logical(plot), length(plot) == 1
  )
  
  # Extract rowData from the SummarizedExperiment
  row_data <- rowData(dep, use.names = FALSE)
  
  # Change: Validate that required columns ("name" and "ID") are present in rowData**
  if(any(!c("name", "ID") %in% colnames(row_data))) {
    stop(paste0("'name' and/or 'ID' columns are not present in '",
                deparse(substitute(dep)),
                "'.\nRun make_unique() to obtain required columns."),
         call. = FALSE)
  }
  
  # Change: Validate that differential analysis columns ([contrast]_diff and [contrast]_p.adj) exist**
  if (length(grep("_p.adj|_diff", colnames(row_data))) < 1) {
    stop(paste0("'[contrast]_diff' and '[contrast]_p.adj' columns are not present in '",
                deparse(substitute(dep)),
                "'.\nRun test_diff() to obtain the required columns."),
         call. = FALSE)
  }
  
  # Change: Validate that the [contrast]_significant column exists**
  if (length(grep("_significant", colnames(row_data))) < 1) {
    stop(paste0("'[contrast]_significant' columns are not present in '",
                deparse(substitute(dep)),
                "'.\nRun add_rejections() to obtain the required columns."),
         call. = FALSE)
  }
  
  # Change: Check that the provided contrast is valid (i.e. that the [contrast]_diff column exists)**
  if (length(grep(paste(contrast, "_diff", sep = ""), colnames(row_data))) == 0) {
    valid_cntrsts <- row_data %>%
      data.frame() %>%
      dplyr::select(ends_with("_diff")) %>%
      colnames() %>%
      gsub("_diff", "", .)
    valid_cntrsts_msg <- paste0("Valid contrasts are: '",
                                paste0(valid_cntrsts, collapse = "', '"),
                                "'")
    stop("Not a valid contrast, please run plot_volcano() with a valid contrast as argument\n",
         valid_cntrsts_msg,
         call. = FALSE)
  }
  
  # Identify columns corresponding to the contrast difference, p-value and significance
  diff <- grep(paste0("^", contrast, "_diff$"), colnames(row_data))
  if (adjusted) {
    p_values <- grep(paste0("^", contrast, "_p.adj$"), colnames(row_data))
  } else {
    p_values <- grep(paste0("^", contrast, "_p.val$"), colnames(row_data))
  }
  
  # Change: identify significance
  signif <- grep(paste0("^", contrast, "_significant$"), colnames(row_data))
  sig_type_col <- grep(paste0("^", contrast, "_sig_type$"), colnames(row_data))
  
  if (length(sig_type_col) == 1) {
    status_vals <- as.character(row_data[, sig_type_col])
  } else {
    status_vals <- ifelse(as.logical(row_data[, signif]), "DEP", "Not Significant")
  }
  
  # Change: construct the data frame for plotting the volcano plot:
  df <- data.frame(
    x = row_data[, diff],
    y = -log10(row_data[, p_values]),
    significant_base = status_vals,
    name = row_data$name,
    stringsAsFactors = FALSE
  ) %>%
    dplyr::mutate(significant = dplyr::case_when(
      significant_base == "Not Significant" ~ "Not Significant",
      significant_base == "DEP" & x > 0 ~ "DEP Up",
      significant_base == "DEP" & x < 0 ~ "DEP Down",
      significant_base == "PA" & x > 0 ~ "PA Up",
      significant_base == "PA" & x < 0 ~ "PA Down",
      TRUE ~ "Not Significant"
    )) %>%
    dplyr::filter(!is.na(significant)) %>%
    dplyr::mutate(significant = factor(significant,
                                       levels = c("Not Significant",
                                                  "DEP Down", "DEP Up", 
                                                  "PA Down", "PA Up"))) %>%
    dplyr::arrange(significant)

  
  # Change: Build the volcano plot using ggplot2
  p <- ggplot(df, aes(x = x, y = y)) +
    theme_void() +
    scale_x_continuous(breaks = scales::pretty_breaks(n = 10)) +
    scale_y_continuous(breaks = scales::pretty_breaks(n = 10)) +
    labs(title = "Volcano Plot") + 
    geom_vline(xintercept = 0) +
    geom_point(aes(col = significant), size=4) +
    geom_vline(xintercept = c(-log2(1.5), log2(1.5)), linetype = "dashed", color = "black", alpha = 0.8) + 
    geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "black", alpha = 0.8) +
    theme(plot.title   = element_text(size = 55, face = "bold", hjust = 0.5,
                                      margin = margin(b=10)),
          axis.title.x = element_text(size = 40, face = "bold", margin = margin(t = 10)),
          axis.title.y = element_text(size = 40,  face = "bold", angle = 90,
                                      margin = margin(r = 10)),
          legend.position = "none",
          axis.text.x = element_text(size = 30),
          axis.text.y = element_text(size = 30)) +
    scale_color_manual(values = c("Not Significant" = "darkgray", 
                                  "DEP Down" = "darkcyan", "DEP Up" = "orange", 
                                  "PA Down" = "purple", "PA Up" = "maroon"))
  
  # Change: Optionally add protein names via ggrepel if add_names is TRUE**
  if (add_names) {
    if (!is.null(proteins_to_label)) {
      df_label <- dplyr::filter(df, name %in% proteins_to_label)
    } else {
      df_label <- dplyr::filter(df, significant %in% c("DEP Down", "DEP Up", "PA Down", "PA Up"))
    }
    p = p + ggrepel::geom_text_repel(
      data = df_label,
      aes(label = name),
      size = label_size,
      box.padding = unit(0.4, "lines"),  
      point.padding = unit(0.1, "lines"),
      force = 5,
      max.overlaps = Inf,
      segment.color = "black",
      segment.size = 0.1,
      min.segment.length = 0
    )
  }
  
  # Change: Adjust y-axis labeling based on whether p-values are adjusted**
  if (adjusted) {
    p <- p + labs(y = expression(-log[10]~"Adjusted p-value"))
  } else {
    p <- p + labs(y = expression(-log[10]~"P-value"))
  }
  p <- p + labs(x = expression(log[2]~"Fold change"))
  
  # Return either the plot or the underlying data frame
  if (plot) {
    return(p)
  } else {
    df <- df %>%
      dplyr::select(name, x, y, significant) %>%
      dplyr::arrange(desc(x))
    colnames(df)[c(1, 2, 3)] <- c("protein", "log2_fold_change", "p_value_-log10")
    if (adjusted) {
      colnames(df)[3] <- "adjusted_p_value_-log10"
    }
    return(df)
  }
}