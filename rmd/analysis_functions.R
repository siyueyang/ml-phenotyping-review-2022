# ------------------------------------- Library used ---------------------------------------------#
library(splitstackshape)  # To split printed summary.  
library(janitor)

# ------------------------------------- Plotting -------------------------------------------------#
# Basics. 
# Function for stacked bar chart.
stacked_bar <- function(df, 
                        group_var, 
                        groupvar_stack = NA, 
                        order_by_count = TRUE, 
                        string_wrap = NA,
                        large_fig = FALSE){
  
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
  
  p <- p + geom_bar(stat="identity", position = "stack") +
    scale_fill_jama() +
    labs(x = "", y = "") +
    theme_classic() +
    theme(legend.position = "bottom") +
    guides(fill = guide_legend("")) 
  
  if (large_fig) {
    p <- p +  
      geom_text(aes(label = Count), size = 8, color = "white", position = position_stack(vjust = 0.5)) + 
      theme(plot.title = element_text(face = "bold", size = 20),
            axis.text = element_text(size = 20), axis.title = element_text(size = 15),
            legend.title = element_text(size = 20), legend.text = element_text(size = 15)) 
    
  } else {
    p <- p + 
      geom_text(aes(label = Count), size = 3, color = "white", position = position_stack(vjust = 0.5)) 
  }
  
  p
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
                           color_grid = NA,
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
    
    df[, legend_cap] <- df[, groupvar_stack]
    
    temp_df <- df %>%
      group_by(!!sym(group_var)) %>%
      mutate(n = n()) %>%
      group_by(!!sym(group_var), !!sym(legend_cap)) %>%
      summarise(n = mean(n), Count = n()) 
  }
  
  # If restrict the count of the major category of `group_var`.
  if (restrict_count) {
    temp_df <- temp_df %>% filter(n > 1) 
  }
  
  # Keep the top n bars. 
  temp_df <- temp_df %>% arrange(-n) %>% top_n(top_count, n)
  
  # Plotting. 
  if (is.na(groupvar_stack)) {
    p <- ggplot(data=temp_df, aes(x=reorder(!!sym(group_var), -Count, sum), y=Count))
    
    if (!is.na(color_grid)) {
      p <- p + geom_bar(aes(fill = !!sym(group_var)), stat = "identity") 
    } else {
      p <- p + geom_bar(stat = "identity") 
    }
     
  } else {
    p <- ggplot(data=temp_df, aes(x=reorder(!!sym(group_var), -Count, sum), y=Count, alpha = !!sym(legend_cap))) + geom_bar(stat="identity") 
  }
  
  # Label the number on top of the bars. 
  if (label_text) {
    p <- p + geom_text(aes(label = Count), position = position_dodge(width = 0.9), hjust = -0.5, size = 3) 
  }
  
  # Restrict the y limit. 
  if (!is.na(ylim))  p <- p + ylim(0, ylim)
  
  if (!is.na(groupvar_stack)) {
    p <- p + 
      geom_col_pattern(aes(pattern = !!sym(legend_cap), pattern_angle = !!sym(legend_cap), pattern_spacing = !!sym(legend_cap)), 
                        fill            = 'white',
                        colour          = 'black', 
                        pattern_density = 0.1, 
                        pattern_fill    = 'black',
                        pattern_colour  = 'black') +
      theme_classic() + 
      coord_flip() + 
      labs(title = title, y = "", x = "") +
      theme(legend.position = "right")
  } else {
    p <- p + scale_fill_manual(values = color_grid) +
      theme_classic() + 
      coord_flip() + 
      labs(title = title, y = "", x = "") +
      theme(legend.position = "none")
  }
  
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

# Figure to plot the whole top phenotype.
plot_all_top_pheno <- function(traditional_supervised,
                               deep_supervised,
                               semi_supervised,
                               weakly_supervised,
                               un_supervised) {
  
  df_all <- rbind(traditional_supervised, deep_supervised, semi_supervised, weakly_supervised, un_supervised) 
  phenotype_unnested <- unnest_string_var(df_all %>% filter(Competition_data_name == ''), "Phenotype") 
  phenotype_common <- phenotype_unnested %>% group_by(Phenotype_unnested, ML_type, Traditional) %>% summarise(n = n()) %>% 
    filter(n > 1) %>% arrange(-n) %>% group_by(ML_type, Traditional) %>% top_n(5)
  phenotype_unique <- phenotype_common %>% group_by(Phenotype_unnested) %>% summarise(n = n()) %>% filter(n > 1) %>% 
    select(Phenotype_unnested) %>% pull()
  phenotype_common <- phenotype_common %>% select(Phenotype_unnested) %>% pull()
  
  color_grid <- rep("#000003", length(phenotype_common))
  color_panel <- pal_d3()(9)
  
  for (i in c(1:length(phenotype_unique))) {
    color_grid[which(phenotype_common == phenotype_unique[i])] <- color_panel[i]
  }
  
  names(color_grid) <- phenotype_common
              
  p1 <- plot_top_pheno(traditional_supervised, label_text = FALSE,
                       top_count = 5, title = "Traditional supervised learning", color_grid = color_grid) 
  
  p2 <- plot_top_pheno(deep_supervised, label_text = FALSE,       
                       top_count = 5, title = "Deep supervised learning", color_grid = color_grid) 
  
  p3 <- plot_top_pheno(semi_supervised, label_text = FALSE,        
                       top_count = 5, title = "Semi-supervised learning", color_grid = color_grid) 
  
  p4 <- plot_top_pheno(weakly_supervised, label_text = FALSE,     
                       top_count = 5, title = "Weakly-supervised learning", color_grid = color_grid) 
  
  p5 <- plot_top_pheno(un_supervised, label_text = FALSE,         
                       top_count = 5, title = "Unsupervised learning", 
                       string_wrap = 25, ylim = 5, restrict_count = TRUE, 
                       groupvar_stack = "Pheno_cluster_category")
  
  p5_unlegend <- p5 + theme(legend.position = "none", legend.title = element_blank())
  legend <- get_legend(p5)          

  
  plot_grid(p1, p2, p3, p4, p5_unlegend, legend, ncol = 2, nrow = 3, align = "v", axis = "l",
            labels = c('(a)', '(b)', '(c)', '(d)', '(e)', ''))
}


# Figure to plot top phenotypes. 
plot_top_pheno <- function(df, 
                           groupvar_stack = NA,
                           restrict_count = TRUE,
                           top_count = 5, 
                           title, 
                           ylim = 5, 
                           label_text = TRUE,
                           string_wrap = 25,
                           color_grid = NA,
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
                 color_grid = color_grid,
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
  
  p <- ggplot(data=temp_df, aes(x = Category, y = Count, fill = !!sym(group_var))) +
       geom_bar(stat="identity", position = position_stack()) + 
       scale_fill_jama() + 
       labs(y = "", title = title) + 
       theme_bw() +
       theme(axis.title.x=element_blank(), axis.text.x=element_blank(), axis.ticks.x = element_blank()) +
       theme(legend.position = "bottom", legend.title=element_blank()) 
  
  if (large_fig) {
    p + geom_text(aes(label = Count), position = position_stack(), vjust = -0.5, size = 5) +
        theme(axis.text.y = element_text(size = 15), legend.text = element_text(size = 15)) + 
        theme(plot.title = element_text(face = "bold", size = 15)) + 
        theme(strip.text.x = element_text(size = 15)) 
  } else {
    p + geom_text(aes(label = Count), position = position_stack(), vjust = -0.5, size = 2) +
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
    
    df <- df %>% 
      unite("Study1", c(Phenotype, Data_source, PMID), sep = "\n", remove = FALSE) %>% 
      unite("Study2", c(Phenotype, Data_source), sep = "\n", remove = FALSE)
    
    # Arrange by sensitivity. 
    df_sens <- ml_rule_metrics(df, "Sensitivity") %>% 
      select(Phenotype, diff, ML_better) %>% unique() %>% arrange(ML_better, -diff)
    
    df$Phenotype <- factor(df$Phenotype, levels = df_sens$Phenotype)
    ML_better <- c("Clinical trial eligibility for n2c2 2018 challenge",
                   "COVID-19", 
                   "Stroke")
    rule_better <- c("Obesity and multiple comorbidities", 
                     "Marital status")
    comparable <- c("Hepatorenal syndrome", "Alcohol abuse", 
                    "Metastatic prostate cancer",
                    "Acute respiratory distress syndrome",
                    "Rhabdomyolysis")
    
    df_ML_better <- df %>% filter(Phenotype %in% ML_better)
    df_rule_better <- df %>% filter(Phenotype %in% rule_better)
    df_comparable <- df %>% filter(Phenotype %in% comparable)
    
    g1 <- plot_rule_compare(df_ML_better, large_fig = large_fig)
    g3 <- plot_rule_compare(df_comparable, large_fig = large_fig, legend = TRUE)
    g2 <- plot_rule_compare(df_rule_better, large_fig = large_fig) 
    
    plot_grid(g1, g2, g3, nrow = 3, align = "v", labels = c('(a)', '(b)', '(c)'),
              rel_heights = c(3, 2, 5))

  } else if (comparator == "deep") {
   
    df <- df %>% filter(Compare_with_traditional_ML != "")
    df <- df %>% filter("Best_performing_AUROC" != "" & 
                        "Best_performing_PPV" != "" &
                        "Best_performing_Sensitivity" != "" &
                        "Best_performing_Specificity" != "")
    df <- unnest_validate_string(df) 
    df <- df[-which(rowMeans(is.na(df)) >= 0.75), ]
    
    # Replace the long phenotype description. 
    df$Phenotype[df$PMID == 34514351] <- "Social determinants of health"
    df$Phenotype[df$PMID == 35007754] <- "Social determinants of health"
    df$Phenotype[df$PMID == 34791302] <- "Aspects of frailty"
    df$Phenotype[df$PMID == 32548622] <- "Acute care conditions"
    
    df$Data_source[df$PMID == 34791302] <- "Penn Medicine; MIMIC-III database"
    df$Data_source[df$PMID == 34514351] <- "University of North Carolina Health System"
    df$Data_source[df$PMID == 32970677] <- "Optum Analytic"
    
    df <- df %>% 
      unite("Study1", c(Phenotype, PMID), sep = "\n", remove = FALSE) %>% 
      unite("Study2", c(Phenotype, Data_source), sep = "\n", remove = FALSE)
    
    # Arrange by sensitivity. 
    df_sens <- ml_deep_metrics(df, "Sensitivity") %>% 
      select(Study1, diff, DL_better) %>% unique() %>% arrange(DL_better, -diff)
    
    df$Study1 <- factor(df$Study1, levels = df_sens$Study1)
    
    ML_better <- c("Aspects of frailty\n34791302")
    
    DL_comparable <- c("Bleeding\n33936461",
                       "Diabetic retinopathy\n34423259",
                       "Social determinants of health\n35007754",
                       "COVID-19\n32449766",
                       "Chest injury\n33709067",
                       "Alcohol abuse\n29447188",
                       "Adverse drug event\n31197355",
                       "Chronic pain\n29447188",
                       "Metastatic cancer\n29447188")
    
    df_comparable <- df %>% filter(Study1 %in% DL_comparable)
    df_ML_better <- df %>% filter(Study1 %in% ML_better)
    df_DL_better <- df %>% filter(!(Study1 %in% c(DL_comparable, ML_better)))
    
    g1 <- plot_deep_compare(df_DL_better, large_fig = large_fig) 
    g2 <- plot_deep_compare(df_ML_better, large_fig = large_fig) 
    g3 <- plot_deep_compare(df_comparable, large_fig = large_fig, legend = TRUE)
    
    plot_grid(g1, g2, g3, nrow = 3, align = "v", labels = c('(a)', '(b)', '(c)'),
              rel_heights = c(15, 4, 8))
  
  } else if (comparator == "weakly") {
    
    df <- df %>% filter(Compare_with_rule_based != "")
    df <- unnest_validate_string(df, comapartor = "weakly") 
    
    df$Phenotype[df$PMID == 31613361] <- "Chronic diseases (MAP)"
    df$Phenotype[df$PMID == 29126253] <- paste0(df$Phenotype[df$PMID == 29126253], "\n (PheNorm)")
    df$Phenotype[df$PMID == 29788308] <- paste0(df$Phenotype[df$PMID == 29788308], "\n (PheProb)")
    df$Phenotype[df$PMID == 32974638] <- paste0(df$Phenotype[df$PMID == 32974638], "\n (PheMAP)")
    df$Phenotype[df$PMID == 32374408] <- paste0(df$Phenotype[df$PMID == 32374408], "\n (APHRODITE)")
    df$Phenotype[df$PMID == 30639392] <- paste0(df$Phenotype[df$PMID == 30639392], "\n (NimbleMiner)")

    df <- df %>% 
      unite("Study1", c(Phenotype, Data_source, PMID), sep = "\n", remove = FALSE) %>% 
      unite("Study2", c(Phenotype, Data_source), sep = "\n", remove = FALSE)
    
    # Arrange by sensitivity. 
    df_sens <- ml_rule_metrics(df, "Sensitivity") %>% 
      dplyr::select(Phenotype, rule, diff, ML_better) %>% unique() %>% arrange(-abs(diff))

    df$Phenotype <- factor(df$Phenotype, levels = df_sens$Phenotype)
    
    comparable <- c("Metastatic breast cancer",
                    "Rheumatoid arthritis\n (PheProb)",
                    "Obesity\n (APHRODITE)",
                    "Fall\n (NimbleMiner)",
                    "Glaucoma\n (APHRODITE)",
                    "Epilepsy\n (APHRODITE)",
                    "Type 2 diabetes mellitus\n (APHRODITE)",
                    "Cataracts\n (APHRODITE)",
                    "Venous thromboembolism\n (APHRODITE)",
                    "Heart failure\n (APHRODITE)",
                    "Peripheral arterial disease\n (APHRODITE)")
    
    df_ML_better <- df %>% filter(!(Phenotype %in% comparable))
    df_comparable <- df %>% filter(Phenotype %in% comparable)
    
    g1 <- plot_rule_compare(df_ML_better, comparison = "weakly", large_fig = large_fig)
    g3 <- plot_rule_compare(df_comparable, comparison = "weakly", large_fig = large_fig, legend = TRUE)
    
    plot_grid(g1, g3, nrow = 2, align = "v", labels = c('(a)', '(b)'),
              rel_heights = c(5, 5))
  
  } else {
    #weakly v.s. traditional
    df <- df %>% filter(Best_comparator_traditional != "")
    df <- unnest_validate_string(df, comapartor = "weakly") 
    
    df$Phenotype[df$PMID == 29126253] <- paste0(df$Phenotype[df$PMID == 29126253], "\n (PheNorm)")
    df$Phenotype[df$PMID == 33746080] <- paste0(df$Phenotype[df$PMID == 33746080], "\n (PheVis)")
    df$Phenotype[df$PMID == 32548637] <- paste0(df$Phenotype[df$PMID == 32548637], "\n (sureLDA)")
   
    df <- df %>% 
      unite("Study1", c(Phenotype, Data_source, PMID), sep = "\n", remove = FALSE) %>% 
      unite("Study2", c(Phenotype, Data_source), sep = "\n", remove = FALSE) %>%
      filter(Best_performing_AUROC != "")
  
    # Arrange by sensitivity. 
    df_auc <- weakly_rule_metrics(df, "AUROC") %>% 
      dplyr::select(Phenotype, diff) %>% unique() %>% arrange(-(diff)) 
    
    df$Phenotype <- factor(df$Phenotype, levels = df_auc$Phenotype)
    
    g1 <- plot_weakly_supervised_compare(df, comparison = "weakly deep", study = "Phenotype", large_fig = large_fig) 
    plot_grid(g1, nrow = 1)
  }
}

# Figure to plot weakly-supervised ML vs rule across all metrics. 
plot_weakly_supervised_compare <- function(df, comparison = "rule", study = "Phenotype", large_fig = FALSE) {
  
  p1 <- plot_metrics(df, comparison, study, metric = "Sensitivity", large_fig) 
  
  p2 <- plot_metrics(df, comparison, study, metric = "Specificity", large_fig)  +
    theme(axis.title.y=element_blank(), axis.text.y=element_blank(), axis.ticks.y=element_blank()) + 
    theme(legend.position = "bottom")
  
  p5 <- plot_metrics(df, comparison, study, metric = "AUROC", large_fig)  +
    theme(axis.title.y=element_blank(), axis.text.y=element_blank(), axis.ticks.y=element_blank())
  
  egg::ggarrange(p1,p2,p5, nrow = 1, draw = FALSE)
}

# Figure to plot ML vs rule across all metrics. 
plot_rule_compare <- function(df, comparison = "rule", study = "Phenotype", large_fig = FALSE, legend = FALSE) {

  p1 <- plot_metrics(df, comparison, study, metric = "Sensitivity", large_fig) 
  
  p2 <- plot_metrics(df, comparison, study, metric = "Specificity", large_fig)  +
    theme(axis.title.y=element_blank(), axis.text.y=element_blank(), axis.ticks.y=element_blank())
  
  p3 <- plot_metrics(df, comparison, study, metric = "PPV", large_fig)  +
    theme(axis.title.y=element_blank(), axis.text.y=element_blank(), axis.ticks.y=element_blank())
  
  if (legend == TRUE) {
    p3 <- p3 + theme(legend.position = "bottom")
  } else {
    p3 <- p3 + theme(legend.position = "none")
  }
  
  p4 <- plot_metrics(df, comparison, study, metric = "NPV", large_fig)  +
    theme(axis.title.y=element_blank(), axis.text.y=element_blank(), axis.ticks.y=element_blank())
  
  p5 <- plot_metrics(df, comparison, study, metric = "AUROC", large_fig)  +
    theme(axis.title.y=element_blank(), axis.text.y=element_blank(), axis.ticks.y=element_blank())
  
  egg::ggarrange(p1,p2,p3,p4,p5, nrow = 1, draw = FALSE)
}

# Figure to plot ML vs DL across all metrics.  
plot_deep_compare <- function(df, comparison = "deep", study = "Study1", large_fig = FALSE, 
                              legend = FALSE) {
  p1 <- plot_metrics(df, comparison, study, "Sensitivity", large_fig = large_fig)
  
  p2 <- plot_metrics(df, comparison, study, "Specificity", large_fig = large_fig) + 
    theme(axis.title.y=element_blank(), axis.text.y=element_blank(), axis.ticks.y=element_blank()) 
  
  if (legend == TRUE) {
    p2 <- p2 + 
      theme(legend.position = "bottom")
  } else {
    p2 <- p2 + 
      theme(legend.position = "none")
  }
  
  p4 <- plot_metrics(df, comparison, study, "PPV", large_fig = large_fig) + 
    theme(axis.title.y=element_blank(), axis.text.y=element_blank(), axis.ticks.y=element_blank()) 
  
  p5 <- plot_metrics(df, comparison, study, metric = "AUROC", large_fig = large_fig) + 
    theme(axis.title.y=element_blank(), axis.text.y=element_blank(), axis.ticks.y=element_blank())
  
  egg::ggarrange(p1,p2, p4,p5, nrow = 1, draw = FALSE)
}

plot_metrics <- function(df, comparison = "rule", study = "Phenotype", metric = "Sensitivity", large_fig = FALSE) {
  
  # Remove all the metrics are NAs. 
  df1 <- apply(df[, c(which(colnames(df) == "Best_performing_Sensitivity")):ncol(df)], 2, as.numeric)
  
  if (comparison == "deep" & typeof(df1) == "double") {
    df1 <- t(data.frame(df1))
  }
  
  missing_index <- rowMeans(is.na(df1))
  
  if (sum(missing_index == 1) > 0) {
    df <- df[-which(missing_index == 1), ]
  }
  
  if (comparison == "rule") {
    df <- ml_rule_metrics(df, metric)
  } else if (comparison == "weakly") {
    df <- weakly_metrics(df, metric)
  } else if (comparison == "deep") {
    df <- ml_deep_metrics(df, metric) 
  } else {
    df <- weakly_rule_metrics(df, metric)
  }
  
  if (study == "Phenotype") {
    p <- df %>%
      ggplot(aes(x = !!sym(study), y = as.numeric(!!sym(metric)), color = Method)) +
      scale_x_discrete(labels = function(x) str_wrap(x, width = 25)) 
    
  } else {
    p <- df %>%
      ggplot(aes(x = !!sym(study), y = as.numeric(!!sym(metric)), color = Method)) +
      scale_x_discrete(labels = function(x) str_extract(x, ".+\\n"))
  }
  
  p <- p + scale_y_continuous(limits = c(0.01, 1), labels = function(y) label_parsed(paste0(y*100))) +
    scale_color_jama() + 
    coord_flip() + 
    labs(x = "", y = "", title = metric) +
    theme_bw() + 
    theme(legend.position = "none", legend.title=element_blank()) 
  
  if (large_fig) {
    p + geom_point(size = 3) + 
      theme(plot.title = element_text(size=20, hjust = 0.5, face="bold")) +
      theme(axis.title.x = element_blank(), axis.text.x = element_text(size = 10), 
            axis.text.y = element_text(size = 15), legend.text = element_text(size = 15))
  } else {
    p + geom_point(size = 2) +
      theme(plot.title = element_text(size=8,  hjust = 0.5, face="bold")) +
      theme(axis.title.x = element_blank(), axis.text.x = element_text(size = 5), 
            axis.text.y = element_text(size = 7), legend.text = element_text(size = 7))
  }
}

# Utility function to plot supervised ML vs rule for one metric. 
ml_rule_metrics <- function(df, metric = "Sensitivity") {
  
  metric_best <- paste0("Best_performing_", metric)
  metric_ml <- paste0("Best_comparator_traditional_", metric)
  metric_rule <- paste0("Best_comparator_Rule_", metric)
  
  df %>%
    filter(Best_performing_model != "") %>%
    mutate(`Traditional supervised learning` = case_when(str_detect(Best_performing_model, "Traditional ML") ~ as.numeric(!!sym(metric_best)),
                          TRUE ~ as.numeric(!!sym(metric_ml)))) %>%
    mutate(`Rule-based algorithms` = case_when(str_detect(Best_performing_model, "Rule-based") ~ as.numeric(!!sym(metric_best)),
                            TRUE ~ as.numeric(!!sym(metric_rule)))) %>%
    mutate(diff = `Traditional supervised learning` - `Rule-based algorithms`) %>%
    mutate(ML_better = diff > 0) %>%
    mutate(rule = `Rule-based algorithms`) %>%
    pivot_longer(cols = c(`Traditional supervised learning`, `Rule-based algorithms`),
                 names_to = "Method", values_to = metric) %>%
    dplyr::select(Study1, Study2, PMID, Phenotype, rule, diff, ML_better, !!sym(metric), Method) %>%
    mutate_if(~ all(. %in% c(0, NA)), ~ replace(., is.na(.), 0)) 
}

# Utility function to plot supervised ML vs rule for one metric. 
weakly_metrics <- function(df, metric = "Sensitivity") {
  
  metric_best <- paste0("Best_performing_", metric)
  metric_ml <- paste0("Best_comparator_traditional_", metric)
  metric_rule <- paste0("Best_comparator_Rule_", metric)
  
  df %>%
    filter(Best_performing_model != "") %>%
    mutate(`Weakly-supervised learning` = case_when(str_detect(Best_performing_model, "Traditional ML") ~ as.numeric(!!sym(metric_best)),
                                                         TRUE ~ as.numeric(!!sym(metric_ml)))) %>%
    mutate(`Rule-based algorithms` = case_when(str_detect(Best_performing_model, "Rule-based") ~ as.numeric(!!sym(metric_best)),
                                               TRUE ~ as.numeric(!!sym(metric_rule)))) %>%
    mutate(diff = `Weakly-supervised learning` - `Rule-based algorithms`) %>%
    mutate(ML_better = diff > 0) %>%
    mutate(rule = `Rule-based algorithms`) %>%
    pivot_longer(cols = c(`Weakly-supervised learning`, `Rule-based algorithms`),
                 names_to = "Method", values_to = metric) %>%
    dplyr::select(Study1, Study2, PMID, Phenotype, rule, diff, ML_better, !!sym(metric), Method) %>%
    mutate_if(~ all(. %in% c(0, NA)), ~ replace(., is.na(.), 0)) 
}

# Utility function to plot supervised ML vs rule for one metric. 
weakly_rule_metrics <- function(df, metric = "Sensitivity") {
  
  metric_best <- paste0("Best_performing_", metric)
  metric_ml <- paste0("Best_comparator_traditional_", metric)
  
  df %>%
    mutate(`Weakly-supervised learning` = as.numeric(!!sym(metric_best))) %>%
    mutate(`Traditional supervised learning` = as.numeric(!!sym(metric_ml))) %>%
    mutate(diff = `Weakly-supervised learning` - `Traditional supervised learning`) %>%
    mutate(ML_better = diff > 0) %>%
    pivot_longer(cols = c(`Weakly-supervised learning`, `Traditional supervised learning`),
                 names_to = "Method", values_to = metric) %>%
    dplyr::select(Phenotype, diff, ML_better, !!sym(metric), Method) %>%
    mutate_if(~ all(. %in% c(0, NA)), ~ replace(., is.na(.), 0)) 
}

# Utility function to plot ML vs DL for one metric.
ml_deep_metrics <- function(df, metric = "Sensitivity") {
  
  metric_best <- paste0("Best_performing_", metric)
  metric_ml <- paste0("Best_comparator_traditional_", metric)
  metric_deep <- paste0("Best_comparator_DL_", metric)
  
  df %>%
    filter(Best_performing_model != "") %>%
    mutate(`Deep supervised learning` = case_when(str_detect(Best_performing_model, "DL") ~ as.numeric(!!sym(metric_best)),
                          TRUE ~ as.numeric(!!sym(metric_deep)))) %>%
    mutate(`Traditional supervised learning` = case_when(str_detect(Best_performing_model, "Traditional ML") ~ as.numeric(!!sym(metric_best)),
                          TRUE ~ as.numeric(!!sym(metric_ml)))) %>%
    mutate(diff = `Deep supervised learning` - `Traditional supervised learning`) %>%
    mutate(DL_better = diff > 0) %>%
    pivot_longer(cols = c(`Traditional supervised learning`, `Deep supervised learning`), 
                 names_to = "Method", values_to = metric) %>%
    dplyr::select(Study1, Study2, Phenotype, diff, DL_better, !!sym(metric), Method) %>%
    mutate_if(~ all(. %in% c(0, NA)), ~ replace(., is.na(.), 0)) 
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
                         col_width = 5,
                         hold_pos = FALSE) {
  
  if (is.na(groupvar_stack)) {
    df <- df %>%
      group_by(!!sym(group_var)) %>%
      summarise(Count = n()) %>%
      slice_max(order_by = Count, n = top_count) %>%
      arrange(-Count)
    
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
      slice_max(order_by = Count, n = top_count) 
    
    if (restric_count) df1 <- df1 %>% filter(Count > 1)
    
    df <- df %>%
      group_by(!!sym(group_var), !!sym(groupvar_stack), !!sym(groupvar_stack2)) %>%
      summarise(Count = n()) %>%
      pivot_wider(names_from = c("ML_type", "Traditional"), values_from = "Count", names_sep = " ", values_fill = 0) %>%
      inner_join(df1) %>%
      arrange(-Count)
    
    # %>%
    #   janitor::adorn_totals("row")
    
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
    
    if (hold_pos) {
      df %>% 
        kbl(booktabs = T, caption = title) %>%
        kable_paper("striped") %>%
        kable_styling(position = "center", latex_options = "HOLD_position") %>%
        column_spec(c(1:num_col), width = paste0(col_width, "em"))
      
    } else {
      df %>% 
        kbl(booktabs = T, caption = title) %>%
        kable_paper("striped") %>%
        kable_styling(position = "center", latex_options = "scale_down") %>%
        column_spec(c(1:num_col), width = paste0(col_width, "em"))
    }
    
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

print_summary_table <- function(df) { 
  
  res <- c()
  for (dataframe in df) {
  
    total <- nrow(dataframe)
    freetext <- nrow(dataframe %>% filter(Unstructured == TRUE & Unstructured_data_language != "Non-language" ))
    nlp <- nrow(dataframe %>% filter(Unstructured == TRUE & NLP_software != ""))
    comp_data <- nrow(dataframe %>% filter(Competition_data_name != ''))
    private_multisite <- nrow(dataframe %>% filter(Multi_sites_data == 1 & Openly_available_data == 0))
    open_data <-  nrow(dataframe %>% filter(Openly_available_data == 1 & Competition_data_name == ''))
    private_single <- nrow(dataframe %>% filter(Multi_sites_data == 0 & Openly_available_data == 0))
    compare_rule <- nrow(dataframe %>% filter(Compare_with_rule_based != ""))
    compare_ml <- nrow(dataframe %>% filter(Compare_with_traditional_ML != ""))
    demographics <- nrow(dataframe %>% filter(Reported_demographics == 1))
    open_code <- nrow(dataframe %>% filter(Open_code == 1))
    
    tmp <- c(total, freetext, nlp, comp_data, private_multisite, 
             open_data, private_single, compare_rule, 
             compare_ml, demographics, open_code)
    res <- rbind(res, tmp)
  }
  
  colnames(res) <- c("Total number of papers", "Used free-text", "Used NLP software", "Used competition data", "Used multisite data",
           "Used open data", "Used private single-site data", "Compared to rule-based algorithms",
           "Comapred to traditional ML", "Reported patient demographic", "Released open code")
  row_names <-c("TSL", "DSL", "SSL",
                "WSL", "USL", "Total")
  
  as.data.frame(res, row.names = row_names) %>% 
    kbl(booktabs = T) %>%
    kable_paper("striped") %>%
    kable_styling(position = "center", latex_options = "scale_down") %>%
    column_spec(c(1:12), width = paste0(4, "em"))
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
unnest_validate_string <- function(df, comapartor = "deep") {
  
  res <- df %>% select(PMID, Data_source)
  
  if (comapartor == "deep") {
    for (metric in c("Sensitivity", "Specificity", "PPV", "AUROC")) {
      
      metric_best <- paste0("Best_performing_", metric)
      metric_ml <- paste0("Best_comparator_traditional_", metric)
      metric_deep <- paste0("Best_comparator_DL_", metric)
      
      best <- unnest_two_string(df, vars = c("Phenotype", metric_best))
      ml <- unnest_two_string(df, vars = c("Phenotype", metric_ml))
      deep <- unnest_two_string(df, vars = c("Phenotype", metric_deep))
      
      mfile <- merge(best, ml, all = TRUE)
      mfile <- merge(mfile, deep, all = TRUE)
      res <- merge(res, mfile, all = TRUE)
    }
  } else {
    for (metric in c("Sensitivity", "Specificity", "PPV", "NPV", "AUROC")) {
      
      metric_best <- paste0("Best_performing_", metric)
      metric_ml <- paste0("Best_comparator_traditional_", metric)
      metric_rule <- paste0("Best_comparator_Rule_", metric)
      
      best <- unnest_two_string(df, vars = c("Phenotype", metric_best))
      ml <- unnest_two_string(df, vars = c("Phenotype", metric_ml))
      rule <- unnest_two_string(df, vars = c("Phenotype", metric_rule))
      
      mfile <- merge(best, ml, all = TRUE)
      mfile <- merge(mfile, rule, all = TRUE)
      res <- merge(res, mfile, all = TRUE)
    }
  }
  
  res <- res[!is.na(res$Best_performing_model),]
  #res <- res[!is.na(res$Data_source), ]
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
