library(splitstackshape)

# Figure to plot data types. 
plot_data_source <- function(df, group_var, order_by_count = FALSE, title = NULL,
                          facet = FALSE, facet_var = NULL){
  
  temp_df <- df %>%
    group_by(!!sym(group_var), Category) %>%
    summarise(Count = n()) 
  
  p <- ggplot(data=temp_df, aes(x = !!sym(group_var), y = Count, fill = !!sym(group_var))) +
      geom_bar(stat="identity", position = position_dodge(width = 0.9)) + 
      scale_fill_jama() + facet_wrap(~Category, nrow = 1) +
      geom_text(aes(label = Count), position = position_dodge(width = 0.9), vjust = -0.5, size = 3) 
  
  p + labs(x = "", y = "", title = title) + theme(plot.title = element_text(face = "bold", size = 8)) + theme_bw() +
    theme(axis.title.x=element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank()) +
    theme(legend.position = "bottom", legend.title=element_blank()) +
    theme(strip.text.x = element_text(size = 8)) 
}

# Function for stacked bar chart
stacked_bar <- function(df, group_var, groupvar_stack, 
                        order_by_count = TRUE, string_wrap = NA){
  
  if (!is.na(string_wrap)) {
    df[, group_var] <- str_wrap(pull(df[, group_var]), width = string_wrap)
  }
  
  temp_df <- df %>%
    group_by(!!sym(group_var), !!sym(groupvar_stack)) %>%
    summarise(Count = n()) 
  
  if(order_by_count){
    p <- ggplot(data=temp_df, aes(x=reorder(!!sym(group_var), -Count, sum), y=Count, fill = !!sym(groupvar_stack))) 
  }else{
    p <- ggplot(data=temp_df, aes(x= !!sym(group_var), y=Count, fill = !!sym(groupvar_stack))) 
  }
  
  p + geom_bar(stat="identity", position = "stack") +
    scale_fill_jama() +
    labs(x = "", y = "Number of articles") + 
    geom_text(aes(label = Count), size = 3, color = "white", position = position_stack(vjust = 0.5)) + 
    theme_classic() +
    theme(legend.position = "top", axis.text.x = element_text(size = 9)) +
    guides(fill = guide_legend("")) 
}

# Function to unnest a column with strings separated by a semi-colon
unnest_string_var <- function(df, var) {
  
  df %>% mutate(var_unnested = (strsplit(as.character(!!sym(var)), ";"))) %>%
    unnest(var_unnested) %>%
    mutate_if(is.character, trimws) %>%
    rename_with(stringr::str_replace, pattern = "var_unnested",
                replacement = paste0(var, "_unnested"), 
                matches("var_unnested"))
}

# Function to unnest multiple column with strings separated by a semi-colon
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
unnest_two_string <- function(df, vars = c("Phenotype", "Best_performing_Sensitivity")) {
  
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
    end_pname <- c(end_pname, paste0("_", c(10:25)))
    end_mname <- paste0("_", c(1:25))
    
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
  
  return(res)
}

# Function for horizontal bar chart
horizontal_bar <- function(df, group_var, groupvar_stack = NA, 
                           string_wrap = NA, restrict_count = FALSE, 
                           ylim = NA, top_count = 999, 
                           label_text = TRUE, title = NULL){
  
  if (!is.na(string_wrap)) {
    df[, group_var] <- str_wrap(pull(df[, group_var]), width = string_wrap)
  }
  
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
  
  if (restrict_count) {
    temp_df <- temp_df %>%
      filter(n > 1) 
  }
  
  temp_df <- temp_df %>%
    arrange(-n) %>%
    top_n(top_count, n)
  
  if (is.na(groupvar_stack)) {
    p <- ggplot(data=temp_df, aes(x=reorder(!!sym(group_var), -Count, sum), y=Count)) +
          geom_bar(stat="identity") 
  } else {
    p <- ggplot(data=temp_df, aes(x=reorder(!!sym(group_var), -Count, sum), y=Count, fill = !!sym(groupvar_stack))) +
          geom_bar(stat="identity") 
  }
  
  if (label_text) {
    p <- p + geom_text(aes(label = Count), position = position_dodge(width = 0.9), hjust = -0.5, size = 3) 
  }
  
  p <- p + scale_fill_jama() + theme_classic() + 
    theme(plot.title = element_text(face = "bold", size = 8)) +
    theme(axis.text.y = element_text(size = 7)) +
    coord_flip() + labs(title = title, y = "", x = "") +
    guides(fill = guide_legend(title = "Unsupervised phenotype category")) + 
    theme(legend.position = "right", legend.title = element_text(size = 8),
          legend.text = element_text(size= 7)) 
  
  if (!is.na(ylim)) {
    p <- p + ylim(0, ylim)
  }
  p
}

# Figure to plot top phenotype 
plot_top_pheno <- function(df, top_count = 5, title, ylim = 5, 
                           restrict_count = TRUE, label_text = TRUE,
                           string_wrap = 25, groupvar_stack = NA) {
  
  phenotype_unnested <- unnest_string_var(df %>% filter(Competition_data_name == ''), "Phenotype")
  
  horizontal_bar(phenotype_unnested, 
                 "Phenotype_unnested", 
                 groupvar_stack = groupvar_stack,
                 restrict_count = restrict_count,
                 label_text = label_text,
                 string_wrap = string_wrap, top_count = top_count, 
                 title = title, ylim = ylim) 
}

# Function for key summary stats
print_summary_stats <- function(df, method_type = "traditional supervised learning"){
  
  print(paste('There are', 
              nrow(df), 
              'papers using', method_type))
  
  print(paste('There are', 
              nrow(df %>% filter(Unstructured == TRUE)), 
              'papers using', method_type, 'with unstructured data'))
  
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
  
  print(paste('There are', 
              nrow(df %>% filter(Reported_demographics == 1)), 
              'papers reported', method_type, 'demographics'))
  
  print(paste('There are', 
              nrow(df %>% filter(Open_code == 1)), 
              'papers released', method_type, 'source code'))

}

print_multiple_items <- function(df, item) {
  print(paste('There are', 
              nrow(df %>% group_by(PMID) %>% summarise(n = n()) %>% filter(n > 1)), 
              'papers using multiple', item))
}

print_tables <- function(df, group_var, groupvar_stack = NA, title = NA, 
                         top_count = 5, restric_count = TRUE) {
  
  if (is.na(groupvar_stack)) {
    df <- df %>%
      group_by(!!sym(group_var)) %>%
      summarise(Count = n()) 
  } else {
    df <- df %>%
      group_by(!!sym(group_var), !!sym(groupvar_stack)) %>%
      summarise(Count = n()) 
  }
  
  if (restric_count) {
    df <- df %>% filter(Count > 1)
  }
  
  df %>% 
    slice_max(order_by = Count, n = top_count) %>%
    kbl(booktabs = T, caption = title) %>%
    kable_paper("striped") %>%
    kable_styling(position = "center", latex_options = "HOLD_position") 
}


# Figure to plot validation metrics.
validate_metrics <- function(df, comparator = "rule") {
  
  if (comparator == "rule") {

    df$Phenotype[df$PMID == 31622801] <- "Obesity and multiple comorbidities"
    df$Phenotype <- str_wrap(df$Phenotype, width = 20)
    
    plot_rule_compare(df)
    
  } else {
    
    
    df$Phenotype[df$PMID == 34514351] <- "Other social determinants of health"
    df$Phenotype[df$PMID == 35007754] <- "Other social determinants of health"
    df$Phenotype[df$PMID == 34791302] <- "Aspects of frailty"
    
    df$Phenotype <- str_wrap(df$Phenotype, width = 20)
    df <- df %>% filter(Best_performing_model != "")
    df <- unnest_validate_string(df) 
    
    plot_deep_compare(df)
  }
}

plot_rule_compare <- function(df) {
  p1 <- ml_rule_metrics(df, "Sensitivity") 

  p2 <- ml_rule_metrics(df, "Specificity") +
        theme(axis.title.y=element_blank(), axis.text.y=element_blank(), axis.ticks.y=element_blank())

  p3 <- ml_rule_metrics(df, "PPV") +
    theme(axis.title.y=element_blank(), axis.text.y=element_blank(), axis.ticks.y=element_blank()) + 
    theme(legend.position = "bottom")

  p4 <- ml_rule_metrics(df, "NPV") +
    theme(axis.title.y=element_blank(), axis.text.y=element_blank(), axis.ticks.y=element_blank())

  p5 <- ml_rule_metrics(df, "AUROC") +
    theme(axis.title.y=element_blank(), axis.text.y=element_blank(), axis.ticks.y=element_blank())

  egg::ggarrange(p1,p2,p3,p4,p5, nrow = 1)
}

plot_deep_compare <- function(df) {
  p1 <- ml_deep_metrics(df, "Sensitivity") 
  
  p2 <- ml_deep_metrics(df, "Specificity") + 
    theme(axis.title.y=element_blank(), axis.text.y=element_blank(), axis.ticks.y=element_blank()) 
  
  p3 <- ml_deep_metrics(df, "Fscore") + 
    theme(axis.title.y=element_blank(), axis.text.y=element_blank(), axis.ticks.y=element_blank()) + 
    theme(legend.position = "bottom")
  
  p4 <- ml_deep_metrics(df, "PPV") + 
    theme(axis.title.y=element_blank(), axis.text.y=element_blank(), axis.ticks.y=element_blank()) 
  
  p5 <- ml_deep_metrics(df, "AUROC") + 
    theme(axis.title.y=element_blank(), axis.text.y=element_blank(), axis.ticks.y=element_blank())

  egg::ggarrange(p1,p2,p3,p4,p5, nrow = 1)

}

ml_rule_metrics <- function(df, metric = "Sensitivity") {
  
  metric_best <- paste0("Best_performing_", metric)
  metric_ml <- paste0("Best_comparator_traditional_", metric)
  metric_rule <- paste0("Best_comparator_Rule_", metric)

  df %>%
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
    coord_flip() +
    geom_point(size = 2) + ylim(0, 1) +
    labs(x = "", y = "", title = metric) +
    theme_bw() + 
    theme(plot.title = element_text(size=10, face="bold")) +
    theme(legend.position = "none", legend.title=element_blank()) +
    theme(axis.title.x = element_blank(), axis.text.x = element_text(size = 6), axis.text.y = element_text(size = 7))
}

ml_deep_metrics <- function(df, metric = "Sensitivity") {
  
  metric_best <- paste0("Best_performing_", metric)
  metric_ml <- paste0("Best_comparator_traditional_", metric)
  metric_deep <- paste0("Best_comparator_DL_", metric)
  
  df %>%
    filter(Best_performing_model != "") %>%
    mutate(DL = case_when(str_detect(Best_performing_model, "DL") ~ as.numeric(!!sym(metric_best)),
                          TRUE ~ as.numeric(!!sym(metric_deep)))) %>%
    mutate(ML = case_when(str_detect(Best_performing_model, "Traditional ML") ~ as.numeric(!!sym(metric_best)),
                          TRUE ~ as.numeric(!!sym(metric_ml)))) %>%
    pivot_longer(cols = c(ML, DL),names_to = "Method", values_to = metric) %>%
   # unite(Phenotype, c(Phenotype, PMID)) %>%
    dplyr::select(Phenotype, !!sym(metric), Method) %>%
    mutate_if(~ all(. %in% c(0, NA)), ~ replace(., is.na(.), 0)) %>%
    ggplot(aes(x = Phenotype, y = as.numeric(!!sym(metric)), color = Method)) +
    scale_color_jama() + 
    coord_flip() +
    geom_point(size = 2) + ylim(0, 1) +
    labs(x = "", y = "", title = metric) +
    theme_bw() + 
    theme(plot.title = element_text(size=10, face="bold")) +
    theme(legend.position = "none", legend.title=element_blank()) +
    theme(axis.title.x = element_blank(), axis.text.x = element_text(size = 6), axis.text.y = element_text(size = 7))
}

