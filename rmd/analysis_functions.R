# ------------------------------------- Library used ---------------------------------------------#
library(splitstackshape)  # To split printed summary.  

# ------------------------------------- Plotting -------------------------------------------------#
# Basics. 
# Function for stacked bar chart.
stacked_bar <- function(df, 
                        group_var, 
                        groupvar_stack = NA, 
                        order_by_count = TRUE, 
                        string_wrap = NA){
  
  # To split the label text if it is too long. 
  if (!is.na(string_wrap)) {
    df[, group_var] <- str_wrap(pull(df[, group_var]), width = string_wrap)
  }
  
  temp_df <- df %>% group_by(!!sym(group_var), !!sym(groupvar_stack)) %>% summarise(Count = n()) 
  
  if (order_by_count) {
    p <- ggplot(data = temp_df, aes(x = reorder(!!sym(group_var), -Count, sum), y = Count, fill = !!sym(groupvar_stack))) 
  
  } else {
    p <- ggplot(data = temp_df, aes(x = !!sym(group_var), y = Count, fill = !!sym(groupvar_stack))) 
  
  }
  
  p + geom_bar(stat="identity", position = "stack") +
    scale_fill_jama() +
    labs(x = "", y = "Number of articles") + 
    geom_text(aes(label = Count), size = 3, color = "white", position = position_stack(vjust = 0.5)) + 
    theme_classic() +
    theme(legend.position = "top", axis.text.x = element_text(size = 9)) +
    guides(fill = guide_legend("")) 
}


# Function for horizontal bar chart.
horizontal_bar <- function(df, 
                           group_var, 
                           groupvar_stack = NA, 
                           restrict_count = FALSE, 
                           top_count = 999, 
                           ylim = NA, 
                           label_text = TRUE, 
                           title = NULL,
                           string_wrap = NA,
                           large_fig = FALSE,
                           legend_cap = "Unsupervised phenotype category"){
  
  # To split the label text if it is too long. 
  if (!is.na(string_wrap)) {
    df[, group_var] <- str_wrap(pull(df[, group_var]), width = string_wrap)
  }
  
  # Group by stacked variable if any. 
  if (is.na(groupvar_stack)) {
    temp_df <- df %>%
      group_by(!!sym(group_var)) %>%
      summarise(Count = n()) %>%
      mutate(n = Count)
  } else {
    temp_df <- df %>%
      group_by(!!sym(group_var)) %>%
      mutate(n = n()) %>%
      group_by(!!sym(group_var), !!sym(groupvar_stack)) %>%
      summarise(n = mean(n), Count = n()) 
  }
  
  # If restric the count of the major category of `group_var`.
  if (restrict_count) {
    temp_df <- temp_df %>% filter(n > 1) 
  }
  
  # Keep the top n bars. 
  temp_df <- temp_df %>% arrange(-n) %>% top_n(top_count, n)
  
  # Plotting. 
  if (is.na(groupvar_stack)) {
    p <- ggplot(data=temp_df, aes(x=reorder(!!sym(group_var), -Count, sum), y=Count)) + geom_bar(stat="identity") 
  } else {
    p <- ggplot(data=temp_df, aes(x=reorder(!!sym(group_var), -Count, sum), y=Count, fill = !!sym(groupvar_stack))) + geom_bar(stat="identity") 
  }
  
  # Label the number on top of the bars. 
  if (label_text) {
    p <- p + geom_text(aes(label = Count), position = position_dodge(width = 0.9), hjust = -0.5, size = 3) 
  }
  
  # Restric the y limit. 
  if (!is.na(ylim))  p <- p + ylim(0, ylim)
  
  p <- p + 
    scale_fill_jama() + 
    theme_classic() + 
    coord_flip() + 
    labs(title = title, y = "", x = "") +
    guides(fill = guide_legend(title = legend_cap)) +
    theme(legend.position = "right")
  
  
  if (large_fig) {
    p <- p + 
      theme(plot.title = element_text(face = "bold", size = 20),
            axis.text = element_text(size = 15),
            legend.title = element_text(size = 20), legend.text = element_text(size = 15)) 
  } else {
    p <- p + 
      theme(plot.title = element_text(face = "bold", size = 8),
            axis.text = element_text(size = 7),
            legend.title = element_text(size = 8), legend.text = element_text(size = 7)) 
  }
  p
}

# Figure to plot top phenotype. 
plot_top_pheno <- function(df, 
                           groupvar_stack = NA,
                           restrict_count = TRUE,
                           top_count = 5, 
                           title, 
                           ylim = 5, 
                           label_text = TRUE,
                           string_wrap = 25,
                           large_fig = FALSE) {
  
  phenotype_unnested <- unnest_string_var(df %>% filter(Competition_data_name == ''), "Phenotype")
  
  horizontal_bar(phenotype_unnested, 
                 "Phenotype_unnested", 
                 groupvar_stack = groupvar_stack,
                 restrict_count = restrict_count,
                 top_count = top_count, 
                 ylim = ylim,
                 label_text = label_text,
                 title = title, 
                 string_wrap = string_wrap,
                 large_fig = large_fig) 
}

# Figure to plot data types. 
plot_data_source <- function(df, 
                             group_var, 
                             order_by_count = FALSE, 
                             title = NULL,
                             facet = FALSE, 
                             facet_var = NULL,
                             large_fig = FALSE){
  
  temp_df <- df %>%
    group_by(!!sym(group_var), Category) %>%
    summarise(Count = n()) 
  
  p <- ggplot(data=temp_df, aes(x = !!sym(group_var), y = Count, fill = !!sym(group_var))) +
       geom_bar(stat="identity", position = position_dodge(width = 0.9)) + 
       scale_fill_jama() + 
       facet_wrap(~Category, nrow = 1) +
       labs(x = "", y = "", title = title) + 
       theme_bw() +
       theme(axis.title.x=element_blank(), axis.text.x=element_blank(), axis.ticks.x = element_blank()) +
       theme(legend.position = "bottom", legend.title=element_blank()) 
  
  if (large_fig) {
    p + geom_text(aes(label = Count), position = position_dodge(width = 0.9), vjust = -0.5, size = 5) +
        theme(axis.text.y = element_text(size = 15), legend.text = element_text(size = 15)) + 
        theme(plot.title = element_text(face = "bold", size = 15)) + 
        theme(strip.text.x = element_text(size = 15)) 
  } else {
    p + geom_text(aes(label = Count), position = position_dodge(width = 0.9), vjust = -0.5, size = 2) +
        theme(axis.text.y = element_text(size = 8), legend.text = element_text(size = 8)) + 
        theme(plot.title = element_text(face = "bold", size = 8)) + 
        theme(strip.text.x = element_text(size = 8)) 
  }
}

# Figure to plot validation metrics, two options: ML vs Rule and ML vs DL.
plot_validate_metrics <- function(df, comparator = "rule", large_fig = FALSE) {
  
  if (comparator == "rule") {
    # Replace the long phenotype description. 
    df$Phenotype[df$PMID == 31622801] <- "Obesity and multiple comorbidities"
    
    plot_rule_compare(df, large_fig = large_fig) 

  } else {
    # Replace the long phenotype description. 
    df$Phenotype[df$PMID == 34514351] <- "Other social determinants of health"
    df$Phenotype[df$PMID == 35007754] <- "Other social determinants of health"
    df$Phenotype[df$PMID == 34791302] <- "Aspects of frailty"
    
    df <- df %>% filter(Best_performing_model != "")
    df <- unnest_validate_string(df) 
    
    plot_deep_compare(df, large_fig = large_fig)
  }
}

# Figure to plot ML vs rule across all metrics. 
plot_rule_compare <- function(df, large_fig = FALSE) {
  p1 <- ml_rule_metrics(df, "Sensitivity", large_fig) 
  
  p2 <- ml_rule_metrics(df, "Specificity", large_fig)  +
    theme(axis.title.y=element_blank(), axis.text.y=element_blank(), axis.ticks.y=element_blank())
  
  p3 <- ml_rule_metrics(df, "PPV", large_fig)  +
    theme(axis.title.y=element_blank(), axis.text.y=element_blank(), axis.ticks.y=element_blank()) + 
    theme(legend.position = "bottom")
  
  p4 <- ml_rule_metrics(df, "NPV", large_fig)  +
    theme(axis.title.y=element_blank(), axis.text.y=element_blank(), axis.ticks.y=element_blank())
  
  p5 <- ml_rule_metrics(df, "AUROC", large_fig)  +
    theme(axis.title.y=element_blank(), axis.text.y=element_blank(), axis.ticks.y=element_blank())
  
  egg::ggarrange(p1,p2,p3,p4,p5, nrow = 1)
}

# Figure to plot ML vs DL across all metrics.  
plot_deep_compare <- function(df, large_fig = FALSE) {
  p1 <- ml_deep_metrics(df, "Sensitivity", large_fig)
  
  p2 <- ml_deep_metrics(df, "Specificity", large_fig) + 
    theme(axis.title.y=element_blank(), axis.text.y=element_blank(), axis.ticks.y=element_blank()) 
  
  p3 <- ml_deep_metrics(df, "Fscore", large_fig) + 
    theme(axis.title.y=element_blank(), axis.text.y=element_blank(), axis.ticks.y=element_blank()) + 
    theme(legend.position = "bottom")
  
  p4 <- ml_deep_metrics(df, "PPV", large_fig) + 
    theme(axis.title.y=element_blank(), axis.text.y=element_blank(), axis.ticks.y=element_blank()) 
  
  p5 <- ml_deep_metrics(df, "AUROC", large_fig) + 
    theme(axis.title.y=element_blank(), axis.text.y=element_blank(), axis.ticks.y=element_blank())
  
  egg::ggarrange(p1,p2,p3,p4,p5, nrow = 1)
}

# Utility function to plot ML vs rule for one metric. 
ml_rule_metrics <- function(df, metric = "Sensitivity", large_fig = FALSE) {
  
  metric_best <- paste0("Best_performing_", metric)
  metric_ml <- paste0("Best_comparator_traditional_", metric)
  metric_rule <- paste0("Best_comparator_Rule_", metric)
  
  df$Phenotype <- str_wrap(df$Phenotype, width = 25)
  
  p <- df %>%
    filter(Best_performing_model != "") %>%
    mutate(ML = case_when(str_detect(Best_performing_model, "Traditional ML") ~ as.numeric(!!sym(metric_best)),
                          TRUE ~ as.numeric(!!sym(metric_ml)))) %>%
    mutate(Rule = case_when(str_detect(Best_performing_model, "Rule-based") ~ as.numeric(!!sym(metric_best)),
                            TRUE ~ as.numeric(!!sym(metric_rule)))) %>%
    pivot_longer(cols = c(ML, Rule),
                 names_to = "Method", values_to = metric) %>%
    dplyr::select(Phenotype, !!sym(metric), Method) %>%
    mutate_if(~ all(. %in% c(0, NA)), ~ replace(., is.na(.), 0)) %>%
    ggplot(aes(x = Phenotype, y = as.numeric(!!sym(metric)), color = Method)) +
    scale_color_jama() + 
    coord_flip() + ylim(0, 1) +
    labs(x = "", y = "", title = metric) +
    theme_bw() + 
    theme(legend.position = "none", legend.title=element_blank()) 
  
  if (large_fig) {
    p + geom_point(size = 3) + 
      theme(plot.title = element_text(size=20, face="bold")) +
      theme(axis.title.x = element_blank(), axis.text.x = element_text(size = 10), 
            axis.text.y = element_text(size = 15), legend.text = element_text(size = 15))
  } else {
    p + geom_point(size = 2) +
      theme(plot.title = element_text(size=8, face="bold")) +
      theme(axis.title.x = element_blank(), axis.text.x = element_text(size = 5), 
            axis.text.y = element_text(size = 7), legend.text = element_text(size = 7))
  }
}

# Utility function to plot ML vs DL for one metric.
ml_deep_metrics <- function(df, metric = "Sensitivity", large_fig = FALSE) {
  
  metric_best <- paste0("Best_performing_", metric)
  metric_ml <- paste0("Best_comparator_traditional_", metric)
  metric_deep <- paste0("Best_comparator_DL_", metric)
  
  df$Phenotype <- str_wrap(df$Phenotype, width = 25)
  
  p <- df %>%
    filter(Best_performing_model != "") %>%
    mutate(DL = case_when(str_detect(Best_performing_model, "DL") ~ as.numeric(!!sym(metric_best)),
                          TRUE ~ as.numeric(!!sym(metric_deep)))) %>%
    mutate(ML = case_when(str_detect(Best_performing_model, "Traditional ML") ~ as.numeric(!!sym(metric_best)),
                          TRUE ~ as.numeric(!!sym(metric_ml)))) %>%
    pivot_longer(cols = c(ML, DL),names_to = "Method", values_to = metric) %>%
    dplyr::select(Phenotype, !!sym(metric), Method) %>%
    mutate_if(~ all(. %in% c(0, NA)), ~ replace(., is.na(.), 0)) %>%
    ggplot(aes(x = Phenotype, y = as.numeric(!!sym(metric)), color = Method)) +
    scale_color_jama() + 
    coord_flip() + ylim(0, 1) +
    labs(x = "", y = "", title = metric) +
    theme_bw() + 
    theme(legend.position = "none", legend.title=element_blank()) 
  
  
  if (large_fig) {
    p + 
      geom_point(size = 3) + 
      theme(plot.title = element_text(size=20, face="bold")) +
      theme(axis.title.x = element_blank(), axis.text.x = element_text(size = 10), 
            axis.text.y = element_text(size = 15), legend.text = element_text(size = 15))
  } else {
    p +
      geom_point(size = 2) + 
      theme(plot.title = element_text(size=8, face="bold")) +
      theme(axis.title.x = element_blank(), axis.text.x = element_text(size = 5), 
            axis.text.y = element_text(size = 7), legend.text = element_text(size = 7))
  }
}

# ------------------------------------- Tables -------------------------------------------------#
# Print a summary tables for grouped information. 
print_tables <- function(df, 
                         group_var, 
                         groupvar_stack = NA, 
                         groupvar_stack2 = NA,
                         title = NA, 
                         top_count = 5, 
                         restric_count = TRUE,
                         col_width = 5) {
  
  if (is.na(groupvar_stack)) {
    df <- df %>%
      group_by(!!sym(group_var)) %>%
      summarise(Count = n()) %>%
      slice_max(order_by = Count, n = top_count)
    
    if (restric_count) df <- df %>% filter(Count > 1)
    
  } else if (is.na(groupvar_stack2)) {
    df <- df %>%
      group_by(!!sym(group_var), !!sym(groupvar_stack)) %>%
      summarise(Count = n()) %>%
      slice_max(order_by = Count, n = top_count)
    
    if (restric_count) df <- df %>% filter(Count > 1)
    
  } else {
    df1 <- df %>%
      group_by(!!sym(group_var)) %>%
      summarise(Count = n()) %>%
      slice_max(order_by = Count, n = top_count) %>% 
      filter(Count > 1)
    
    df <- df %>%
      group_by(!!sym(group_var), !!sym(groupvar_stack), !!sym(groupvar_stack2)) %>%
      summarise(Count = n()) %>%
      pivot_wider(names_from = c("ML_type", "Traditional"), values_from = "Count", names_sep = " ", values_fill = 0) %>%
      inner_join(df1)
  }
  
  new_name <- strsplit(colnames(df)[1], "_")
  new_name <- new_name[[1]] 

  colnames(df)[1] <- str_c(new_name[-length(new_name)], collapse = " ")
  
  if (!is.na(groupvar_stack) & is.na(groupvar_stack2)) {
    new_name <- strsplit(colnames(df)[2], "_")
    new_name <- new_name[[1]] 
    colnames(df)[2] <- str_c(new_name[-length(new_name)], collapse = " ")
  }

  # To split the label text if it is too long. 
  if (!is.na(groupvar_stack2)) {
    
    num_col <- ncol(df) 
    
    df %>% 
      kbl(booktabs = T, caption = title) %>%
      kable_paper("striped") %>%
      kable_styling(position = "center", latex_options = "HOLD_position") %>%
      column_spec(c(1:num_col), width = paste0(col_width, "em"))
  } else {
    df %>% 
      kbl(booktabs = T, caption = title) %>%
      kable_paper("striped") %>%
      kable_styling(position = "center", latex_options = "HOLD_position") 
  }
}


# ------------------------------------- Summary -------------------------------------------------#
# Function for key summary stats.
print_summary_stats <- function(df, method_type = "traditional supervised learning"){
  
  print(paste('There are', 
              nrow(df), 
              'papers using', method_type))
  
  print(paste('There are', 
              nrow(df %>% filter(Unstructured == TRUE)), 
              'papers using', method_type, 'with unstructured data'))
  
  print(paste('There are', 
              nrow(df %>% filter(Unstructured == TRUE & NLP_software != "")), 
              'papers using', method_type, 'with NLP software'))
  
  print(paste('There are', 
              nrow(df %>% filter(Competition_data_name != '')), 
              'papers using', method_type, 'with competition data'))
  
  print(paste('There are', 
              nrow(df %>% filter(Multi_sites_data == 1)), 
              'papers using', method_type, 'with data from multiple sites'))
  
  print(paste('There are', 
              nrow(df %>% filter(Openly_available_data == 1)), 
              'papers using', method_type, 'with openly available data'))
  
  print(paste('There are', 
              nrow(df %>% filter(Multi_sites_data == 0 & Openly_available_data == 0)), 
              'papers using', method_type, 'with data from private single site'))
  
  print('-----------------------------')
  
  print(paste('There are', 
              nrow(df %>% filter(Compare_with_rule_based != "")), 
              'papers', method_type, 'compared with rule-based algorithms'))
  
  print(paste('There are', 
              nrow(df %>% filter(Compare_with_traditional_ML != "")), 
              'papers', method_type, 'compared with traditional ML algorithms'))
  
  print('-----------------------------')
  
  print(paste('There are', 
              nrow(df %>% filter(Reported_demographics == 1)), 
              'papers reported', method_type, 'demographics'))
  
  print(paste('There are', 
              nrow(df %>% filter(Open_code == 1)), 
              'papers released', method_type, 'source code'))
  
}

# Function for number of papers using multiple items.
print_multiple_items <- function(df, item) {
  print(paste('There are', 
              nrow(df %>% group_by(PMID) %>% summarise(n = n()) %>% filter(n > 1)), 
              'papers using multiple', item))
}


# ------------------------------------- Unnested functions ---------------------------------------#
# Function to unnest a column with strings separated by a semi-colon.
unnest_string_var <- function(df, var) {
  
  df %>% mutate(var_unnested = (strsplit(as.character(!!sym(var)), ";"))) %>%
    unnest(var_unnested) %>%
    mutate_if(is.character, trimws) %>%
    rename_with(stringr::str_replace, 
                pattern = "var_unnested",
                replacement = paste0(var, "_unnested"), 
                matches("var_unnested"))
}

# Function to unnest two column with strings separated by a semi-colon.
# Specifically for validation columns, they are recorded in paired columns by semi-colons. 
# Input:   col 1: phenotypes hypertension;obesity         col 2: specificity 0.99;0.22  
# Output:  hypertension_specificity: 0.99                 obesity_specificity: 0.22
unnest_validate_string <- function(df) {
  
  res <- df %>% select(PMID)
  
  for (metric in c("Sensitivity", "Specificity", "Fscore", "PPV", "AUROC")) {

    metric_best <- paste0("Best_performing_", metric)
    metric_ml <- paste0("Best_comparator_traditional_", metric)
    metric_deep <- paste0("Best_comparator_DL_", metric)
    
    best <- unnest_two_string(df, vars = c("Phenotype", metric_best))
    ml <- unnest_two_string(df, vars = c("Phenotype", metric_ml))
    deep <- unnest_two_string(df, vars = c("Phenotype", metric_deep))
    
    mfile <- merge(best, ml)
    mfile <- merge(mfile, deep, all = TRUE)
    res <- merge(res, mfile, all = TRUE)
  }
  
  res <- res[!is.na(res$Best_performing_model),]
  return(res)
}

# Function to separate phenotype and validation metrics by pair. 
unnest_two_string <- function(df, vars = c("Phenotype", "Best_performing_Sensitivity"), utility = "validation") {
  
  # Check if there is a paper with record of multiple values of a single phenotype. 
  # e.g. PMID 29447188 has 10 phenotypes and 10 sensitivity values. 
  pheno_view <- unnest_string_var(df, "Phenotype") %>%
    group_by(PMID) %>%
    summarise(n_pheno = n())
  
  metric_view <- unnest_string_var(df, vars[2]) %>%
    group_by(PMID) %>%
    summarise(n_metric = n())
  
  overview <- merge(pheno_view, metric_view, all = TRUE)
  
  # Split multiple values and matched them with the corresponding phenotype. 
  pmid_split <- overview %>% filter(n_metric > 1) %>% select(PMID) %>% pull()
  df_split <- df %>% filter(PMID %in% pmid_split)
  df_remain <- subset(df, !(PMID %in% pmid_split))
  
  if (utility == "validation") {
    # If there is no record with multiple values for a single metric, return itself.  
    if (nrow(df_split) == 0) {
      res <- df_remain %>% select(PMID, Best_performing_model, Phenotype, !!sym(vars[2])) %>% na.omit()
      
      # Otherwise, split the validation metrics by phenotypes.   
    } else {
      for (var in vars) {
        # From library splitstackshape.
        df_split <- cSplit(df_split, c(var), sep = ";")
      }
      
      n_pheno <- dim(df_split %>% select(starts_with(paste0(vars[1], "_"))))[2]
      n_metric <- dim(df_split %>% select(starts_with(paste0(vars[2], "_"))))[2]
      
      end_pname <- paste0("_0", c(1:9))
      end_pname <- c(end_pname, paste0("_", c(10:100)))
      end_mname <- paste0("_", c(1:100))
      
      res <- c()
      for (i in c(1:n_pheno)) {
        tmp <- df_split %>% select(PMID, Best_performing_model, ends_with(end_pname[i]), ends_with(end_mname[i])) 
        
        if (dim(tmp)[2] != 4) {
          break
        }
        colnames(tmp) <- c("PMID", "Best_performing_model", vars)
        res <- rbind(res, tmp)
      }
      
      df_remain[, vars[2]] <- as.numeric(df_remain[, vars[2]])
      
      res <- df_remain %>% select(PMID, Best_performing_model, Phenotype, !!sym(vars[2])) %>%
        bind_rows(res) %>%
        na.omit()
    } 
  } else {
    # If there is no record with multiple values for a single, return itself.  
    if (nrow(df_split) == 0) {
      res <- df_remain %>% select(PMID, Phenotype, !!sym(vars[2])) %>% na.omit()
      
      # Otherwise, split the by phenotypes.   
    } else {
      
      # From library splitstackshape.
      for (var in vars) {
        # From library splitstackshape.
        df_split <- cSplit(df_split, c(var), sep = ";")
      }
      
      n_pheno <- dim(df_split %>% select(starts_with(paste0(vars[1], "_"))))[2]
      n_metric <- dim(df_split %>% select(starts_with(paste0(vars[2], "_"))))[2]
      
      end_pname <- paste0("_0", c(1:9))
      end_pname <- c(end_pname, paste0("_", c(10:25)))
      end_mname <- paste0("_", c(1:25))
      
      res <- c()
      for (i in c(1:n_pheno)) {
        tmp <- df_split %>% select(PMID, ends_with(end_pname[i]), ends_with(end_mname[i])) 
        
        colnames(tmp) <- c("PMID", vars)
        res <- rbind(res, tmp)
      }

      res2 <- unnest_string_var(df_remain, "Phenotype") %>% 
        select(PMID, Phenotype_unnested, !!sym(vars[2]))
     
      colnames(res2) <- c("PMID", vars) 
      
      res <- res %>% bind_rows(res2) %>% na.omit()
    }
  }
  
  return(res)
}
