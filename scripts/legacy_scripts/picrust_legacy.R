library(stats)
# Install the devtools package if not already installed
# install.packages("devtools")
# Install ggpicrust2 from GitHub
#devtools::install_github("cafferychen777/ggpicrust2")
library(ggpicrust2)
library(ggbiplot)
library(readr)
library(MicrobiomeStat)
library(tibble)
library(magrittr)
#BiocManager::install("DESeq2")
#library(DESEq2)
library(ggplot2)
library(ggpubr)
library(rstatix)
library(dplyr)
library(tidyr)
library(GGally)
library(cli)
library(rlang)
library(ggprism)
#BiocManager::install("ALDEx2")
library(ALDEx2)
#BiocManager::install("Maaslin2")
library(Maaslin2)

load(file = system.file("extdata", "kegg_reference.RData", package = "ggpicrust2"))
load(file = system.file("extdata", "KO_mest_reference.RData", package = "ggpicrust2"))

rownames(KO_mest_to_kegg_reference) <- KO_mest_to_kegg_reference$X1
l<-KO_mest_to_kegg_reference["KO_mest00010",]
l<-unlist(l)
l <- na.omit(l)
l <- l[2:length(l)]
l
length(l)

KO_mest_to_kegg_reference["KO_mest00010",]

rowSums(PWY_male[2:17])[order(rowSums(PWY_male[2:17]))]
#################################################################################################
# ****************************** Define pathway_errorbar3 function **************************** #
# ********************************************************************************************* #
# ************************************** Multiple groups ************************************** #

pathway_errorbar3 <-
  function(abundance,
           daa_results_df,
           Group,
           ko_to_kegg = FALSE,
           p_values_threshold = 0.05,
           order = "group",
           select = NULL,
           p_value_bar = TRUE,
           colors = NULL,
           x_lab = NULL) {
    ######
    #abundance = kegg_abund_F1_1year
    #daa_results_df = daa_annotated_final
    #Group = Group
    #p_values_threshold = 2
    #order = "pathway_class"
    #select = daa_annotated_final %>%
    #  arrange(p_adjust) %>%
    #  slice(1:15) %>%
    #  dplyr::select("feature") %>% pull()
    #ko_to_kegg = TRUE
    #p_value_bar = FALSE
    #colors =  c("#98F5FF","#B3EE3A","#CD2626","#4D4D4D")
    #x_lab = "pathway_name"
    ######
    # Identify pathways with missing annotation
    missing_pathways <- daa_results_df[is.na(daa_results_df$pathway_name), "feature"]
    
    # Inform the user about the missing annotations
    if(length(missing_pathways) > 0) {
      message("The following pathways are missing annotations and have been excluded: ",
              paste(missing_pathways, collapse = ", "))
      message("You can use the 'pathway_annotation' function to add annotations for these pathways.")
    }
    
    # Get the names of all columns in the data frame
    column_names <- colnames(daa_results_df)
    
    # Check if there are more than two 'group' columns
    group_columns <- grepl("^group", column_names)
    if (sum(group_columns) > 2) {
      p_value_bar <- FALSE
      message(
        "There are more than two 'group' columns in the 'daa_results_df' data frame. As a result, it is not possible to compute the log2 fold values. The 'p_value_bar' has been automatically set to FALSE."
      )
    }
    
    # Exclude rows with missing pathway annotation
    daa_results_df <- daa_results_df[!is.na(daa_results_df[,x_lab]),]
    
    if (is.null(x_lab)){
      if (ko_to_kegg == TRUE){
        x_lab <- "pathway_name"
      }else{
        x_lab <- "description"
      }
      
      if (is.null(daa_results_df$pathway_name)&is.null(daa_results_df$description)){
        message(
          "Please utilize the 'pathway_annotation' function to annotate the 'daa_results_df' data frame."
        )
      }
    }
    
    if (!(x_lab %in% colnames(daa_results_df))){
      message(
        "The 'x_lab' you defined does not exist as a column in the 'daa_results_df' data frame."
      )
    }
    
    if (nlevels(factor(daa_results_df$method)) != 1) {
      message(
        "The 'method' column in the 'daa_results_df' data frame contains more than one method. Please filter it to contain only one method."
      )
    }
    
    if (nlevels(factor(daa_results_df$group1)) != 1 || nlevels(factor(daa_results_df$group2)) != 1) {
      message(
        "The 'group1' or 'group2' column in the 'daa_results_df' data frame contains more than one group. Please filter each to contain only one group."
      )
    }
    
    if (is.null(colors)) {
      colors <- c("#d93c3e", "#3685bc", "#6faa3e", "#e8a825", "#c973e6", "#ee6b3d", "#2db0a7", "#f25292")[1:nlevels(as.factor(Group))]
    }
    errorbar_abundance_mat <- as.matrix(abundance)
    daa_results_filtered_df <-
      daa_results_df[daa_results_df$p_adjust < p_values_threshold,]
    if (!is.null(select)) {
      daa_results_filtered_sub_df <-
        daa_results_filtered_df[daa_results_filtered_df$feature %in% select, ]
    } else {
      daa_results_filtered_sub_df <- daa_results_filtered_df
    }
    if (nrow(daa_results_filtered_sub_df) > 30) {
      message(
        paste0(
          "The number of features with statistical significance exceeds 30, leading to suboptimal visualization. ",
          "Please use 'select' to reduce the number of features.\n",
          "Currently, you have these features: ",
          paste(paste0('"', daa_results_filtered_sub_df$feature, '"'), collapse = ", "), ".\n",
          "You can find the statistically significant features with the following command:\n",
          "daa_results_df %>% filter(p_adjust < 0.05) %>% select(c(\"feature\",\"p_adjust\"))"
        )
      )
      stop()
    }
    if (nrow(daa_results_filtered_sub_df) == 0){
      stop(
        "Visualization with 'pathway_errorbar' cannot be performed because there are no features with statistical significance. ",
        "For possible solutions, please check the FAQ section of the tutorial."
      )
    }
    # Convert to relative abundance
    relative_abundance_mat <- apply(t(errorbar_abundance_mat), 1, function(x)
      x / sum(x))
    
    # Subset to only include the features present in daa_results_filtered_sub_df$feature
    sub_relative_abundance_mat <- relative_abundance_mat[rownames(relative_abundance_mat) %in% daa_results_filtered_sub_df$feature,]
    
    # Create a matrix for the error bars
    error_bar_matrix <- cbind(
      sample = colnames(sub_relative_abundance_mat),
      group = Group,
      t(sub_relative_abundance_mat)
    )
    
    error_bar_df <- as.data.frame(error_bar_matrix)
    error_bar_df$group <- factor(Group,levels = levels(as.factor(Group)))
    error_bar_pivot_longer_df <- tidyr::pivot_longer(error_bar_df,-c(sample, group))
    error_bar_pivot_longer_tibble <-
      mutate(error_bar_pivot_longer_df, group = as.factor(group))
    error_bar_pivot_longer_tibble$sample <-
      factor(error_bar_pivot_longer_tibble$sample)
    error_bar_pivot_longer_tibble$name <-
      factor(error_bar_pivot_longer_tibble$name)
    error_bar_pivot_longer_tibble$value <-
      as.numeric(error_bar_pivot_longer_tibble$value)
    error_bar_pivot_longer_tibble_summarised <-
      error_bar_pivot_longer_tibble %>% group_by(name, group) %>%
      summarise(mean = mean(value), sd = stats::sd(value))
    error_bar_pivot_longer_tibble_summarised <-
      error_bar_pivot_longer_tibble_summarised %>% mutate(group2 = "nonsense")
    switch(
      order,
      "p_values" = {
        order <- order(daa_results_filtered_sub_df$p_adjust)
      },
      "name" = {
        order <- order(daa_results_filtered_sub_df$feature)
      },
      "group" = {
        daa_results_filtered_sub_df$pro <- 1
        for (i in levels(error_bar_pivot_longer_tibble_summarised$name)) {
          error_bar_pivot_longer_tibble_summarised_sub <-
            error_bar_pivot_longer_tibble_summarised[error_bar_pivot_longer_tibble_summarised$name ==
                                                       i,]
          pro_group <-
            error_bar_pivot_longer_tibble_summarised_sub[error_bar_pivot_longer_tibble_summarised_sub$mean ==
                                                           max(error_bar_pivot_longer_tibble_summarised_sub$mean),]$group
          pro_group <- as.vector(pro_group)
          daa_results_filtered_sub_df[daa_results_filtered_sub_df$feature ==
                                        i,]$pro <- pro_group
        }
        order <-
          order(daa_results_filtered_sub_df$pro,
                daa_results_filtered_sub_df$p_adjust)
      },
      "pathway_class" = {
        if (!"pathway_class" %in% colnames(daa_results_filtered_sub_df)) {
          stop(
            "The 'pathway_class' column is missing in the 'daa_results_filtered_sub_df' data frame. ",
            "Please use the 'pathway_annotation' function to annotate the 'pathway_daa' results."
          )
        }
        order <- order(
          daa_results_filtered_sub_df$pathway_class,
          daa_results_filtered_sub_df$p_adjust
        )
      },
      {
        order <- order
      }
    )
    daa_results_filtered_sub_df <-
      daa_results_filtered_sub_df[order,]
    error_bar_pivot_longer_tibble_summarised_ordered <-
      data.frame(
        name = NULL,
        group = NULL,
        mean = NULL,
        sd = NULL
      )
    for (i in daa_results_filtered_sub_df$feature) {
      error_bar_pivot_longer_tibble_summarised_ordered <-
        rbind(
          error_bar_pivot_longer_tibble_summarised_ordered,
          error_bar_pivot_longer_tibble_summarised[error_bar_pivot_longer_tibble_summarised$name ==
                                                     i,]
        )
    }
    if (ko_to_kegg == FALSE){
      error_bar_pivot_longer_tibble_summarised_ordered[, x_lab] <-
        rep(daa_results_filtered_sub_df[, x_lab], each = length(levels(
          factor(error_bar_pivot_longer_tibble_summarised_ordered$group)
        )))
    }
    
    if (ko_to_kegg == TRUE) {
      error_bar_pivot_longer_tibble_summarised_ordered$pathway_class <-
        rep(daa_results_filtered_sub_df$pathway_class,
            each = length(levels(
              factor(error_bar_pivot_longer_tibble_summarised_ordered$group)
            )))
    }
    error_bar_pivot_longer_tibble_summarised_ordered$name <- factor(error_bar_pivot_longer_tibble_summarised_ordered$name, levels = rev(unique(daa_results_filtered_sub_df$feature)))
    
    bar_errorbar <-
      ggplot2::ggplot(error_bar_pivot_longer_tibble_summarised_ordered, # nolint: object_usage_linter.
                      ggplot2::aes(mean, name, fill = group)) + # nolint
      ggplot2::geom_errorbar(
        ggplot2::aes(xmax = mean + sd, xmin = 0),
        position = ggplot2::position_dodge(width = 0.8),
        width = 0.5,
        size = 0.5,
        color = "black"
      ) +
      ggplot2::geom_bar(stat = "identity",
                        position = ggplot2::position_dodge(width = 0.8),
                        width = 0.8) +
      GGally::geom_stripped_cols(width = 10) +
      ggplot2::scale_fill_manual(values = colors) +
      ggplot2::scale_color_manual(values = colors) +
      ggprism::theme_prism() +
      ggplot2::scale_x_continuous(expand = c(0, 0),
                                  guide = "prism_offset_minor",) +
      ggplot2::scale_y_discrete(labels = rev(daa_results_filtered_sub_df[, x_lab])) +
      ggplot2::labs(x = "Relative Abundance", y = NULL) +
      ggplot2::theme(
        axis.ticks.y = ggplot2::element_blank(),
        axis.line.y = ggplot2::element_blank(),
        axis.line.x = ggplot2::element_line(size = 0.5),
        axis.ticks.x = ggplot2::element_line(size = 0.5),
        panel.grid.major.y = ggplot2::element_blank(),
        panel.grid.major.x = ggplot2::element_blank(),
        axis.text = ggplot2::element_text(size = 10, color = "black"), # nolint
        axis.text.x = ggplot2::element_text(margin = ggplot2::margin(r = 0)), # nolint
        axis.text.y = ggplot2::element_text(
          size = 10,
          color = "black",
          margin = ggplot2::margin(b = 6)
        ),
        axis.title.x =  ggplot2::element_text(
          size = 10,
          color = "black",
          hjust = 0.5
        ),
        legend.position = "top",
        legend.key.size = ggplot2::unit(0.1, "cm"),
        legend.direction = "vertical",
        legend.justification = "left",
        legend.text = ggplot2::element_text(size = 8, face = "bold"),
        legend.box.just = "right",
        plot.margin = ggplot2::margin(0, 0.5, 0.5, 0, unit = "cm")
      ) + ggplot2::coord_cartesian(clip = "off")
    
    if (ko_to_kegg == TRUE) {
      pathway_class_group_mat <-
        daa_results_filtered_sub_df$pathway_class %>%
        table() %>%
        data.frame() %>% column_to_rownames(".")
      pathway_class_group <- data.frame(.= unique(daa_results_filtered_sub_df$pathway_class),Freq = pathway_class_group_mat[unique(daa_results_filtered_sub_df$pathway_class),])
      start <-
        c(1, rev(pathway_class_group$Freq)[1:(length(pathway_class_group$Freq) - 1)]) %>%
        cumsum()
      end <- cumsum(rev(pathway_class_group$Freq))
      ymin <- start - 1 / 2
      ymax <- end + 1 / 2
      nPoints <- length(start)
      pCol <- c("#D51F26","#272E6A","#208A42","#89288F","#F47D2B",
                "#FEE500","#8A9FD1","#C06CAB","#E6C2DC","#90D5E4",
                "#89C75F","#F37B7D","#9983BD","#D24B27","#3BBCA8",
                "#6E4B9E","#0C727C", "#7E1416","#D8A767","#3D3D3D")[1:nPoints]
      pFill <- pCol
      for (i in 1:nPoints)  {
        bar_errorbar <- bar_errorbar +
          ggplot2::annotation_custom(
            grob = grid::rectGrob(
              gp = grid::gpar(
                col = pCol[i],
                fill = pFill[i],
                lty = NULL,
                lwd = NULL,
                alpha = 0.2
              )
            ),
            #xmin = ggplot2::unit(-2, 'native'),
            #xmax = ggplot2::unit(0, 'native'),
            #ymin = ggplot2::unit(ymin[i], 'native'),
            #ymax = ggplot2::unit(ymax[i], 'native')
            xmin = -2,
            xmax = 0,
            ymin = ymin[i],
            ymax = ymax[i]
          )
      }
    }
    daa_results_filtered_sub_df <-
      cbind(
        daa_results_filtered_df,
        #negative_log10_p = -log10(daa_results_filtered_sub_df$p_adjust),
        negative_log10_p = daa_results_filtered_sub_df$p_adjust,
        group_nonsense = "nonsense",
        log_2_fold_change = NA
      )
    
    for (i in daa_results_filtered_sub_df$feature){
      mean <- error_bar_pivot_longer_tibble_summarised_ordered[error_bar_pivot_longer_tibble_summarised_ordered$name %in% i,]$mean
      daa_results_filtered_sub_df[daa_results_filtered_sub_df$feature==i,]$log_2_fold_change <- log2(mean[1]/mean[2])
    }
    daa_results_filtered_sub_df$feature <- factor(daa_results_filtered_sub_df$feature,
                                                  levels = rev(unique(daa_results_filtered_sub_df$feature)))
    
    if (ko_to_kegg == TRUE) {
      pathway_class_y <- (ymax + ymin) / 2 - 0.5
      pathway_class_plot_df <-
        as.data.frame(
          cbind(
            nonsense = "nonsense",
            pathway_class_y = pathway_class_y,
            pathway_class = rev(unique(
              #daa_results_filtered_sub_df$pathway_class
              pathway_class_group$.
            ))
          )
        )
      pathway_class_plot_df$pathway_class_y <-
        as.numeric(pathway_class_plot_df$pathway_class_y)
      pathway_class_annotation <-
        pathway_class_plot_df %>% ggplot2::ggplot(ggplot2::aes(nonsense, pathway_class_y)) + ggplot2::geom_text(
          ggplot2::aes(nonsense, pathway_class_y, label = pathway_class),
          size = 3.5,
          color = "black",
          fontface = "bold",
          family = "sans"
        ) +
        ggplot2::scale_y_discrete(position = "right") +
        ggprism::theme_prism() +
        ggplot2::theme(
          axis.ticks = ggplot2::element_blank(),
          axis.line = ggplot2::element_blank(),
          panel.grid.major.y = ggplot2::element_blank(),
          panel.grid.major.x = ggplot2::element_blank(),
          panel.background = ggplot2::element_blank(),
          axis.text = ggplot2::element_blank(),
          plot.margin = ggplot2::unit(c(0, 0.2, 0, 0), "cm"),
          axis.title.y =  ggplot2::element_blank(),
          axis.title.x = ggplot2::element_blank(),
          legend.position = "non"
        )
    }
    daa_results_filtered_sub_df$negative_log10_p <-
      as.character(daa_results_filtered_sub_df$negative_log10_p)
    daa_results_filtered_sub_df$unique <-
      nrow(daa_results_filtered_sub_df) - seq_len(nrow(daa_results_filtered_sub_df)) + 1
    daa_results_filtered_sub_df$negative_log10_p <-
      substr(daa_results_filtered_sub_df$negative_log10_p, 1, 6)
    p_annotation <- daa_results_filtered_sub_df %>%
      ggplot2::ggplot(ggplot2::aes(group_nonsense, negative_log10_p)) +
      ggplot2::geom_text(
        ggplot2::aes(group_nonsense, unique, label = negative_log10_p),
        size = 3.5,
        color = "black",
        fontface = "bold",
        family = "sans"
      ) +
      #ggplot2::labs(y = "-log10(adj. p-value)") +
      ggplot2::labs(y = "adj. p-value") +
      ggplot2::scale_y_discrete(position = "right") +
      ggprism::theme_prism() +
      ggplot2::theme(
        axis.ticks = ggplot2::element_blank(),
        axis.line = ggplot2::element_blank(),
        panel.grid.major.y = ggplot2::element_blank(),
        panel.grid.major.x = ggplot2::element_blank(),
        panel.background = ggplot2::element_blank(),
        axis.text = ggplot2::element_blank(),
        plot.margin = ggplot2::unit(c(0, 0.2, 0, 0), "cm"),
        axis.title.y =  ggplot2::element_text(
          size = 11,
          color = "black",
          vjust = 0
        ),
        axis.title.x = ggplot2::element_blank(),
        legend.position = "non"
      )
    if (p_value_bar == TRUE) {
      if (ko_to_kegg == TRUE) {
        combination_bar_plot <-
          pathway_class_annotation + bar_errorbar + p_values_bar + p_annotation + patchwork::plot_layout(ncol = 4, widths =
                                                                                                           c(1, 1.2, 0.5, 0.1))
      }
      else{
        combination_bar_plot <-
          bar_errorbar + p_values_bar + p_annotation + patchwork::plot_layout(ncol = 3, widths = c(2.3, 0.7, 0.3))
      }
    }else{
      if (ko_to_kegg == TRUE) {
        combination_bar_plot <-
          pathway_class_annotation + bar_errorbar + p_annotation + patchwork::plot_layout(ncol = 3, widths =
                                                                                            c(1, 1.2, 0.1))
      }
      else{
        combination_bar_plot <-
          bar_errorbar + p_annotation + patchwork::plot_layout(ncol = 2, widths = c(2.5,  0.5))
      }
    }
    return(combination_bar_plot)
  }
# After debug pathway_errorbar(), run:
rm(ko_to_kegg,p_value_threshold,order,colors,x_lab,missing_pathways,column_names,group_columns,daa_results_df,errorbar_abundance_mat,daa_results_filtered_df,daa_results_filtered_sub_df,relative_abundance_mat,sub_relative_abundance_mat,error_bar_matrix,error_bar_df,error_bar_pivot_longer_df,error_bar_pivot_longer_tibble,error_bar_pivot_longer_tibble_summarised,error_bar_pivot_longer_tibble_summarised_ordered,bar_errorbar,pathway_class_group_mat,pathway_class_group,start,end,ymin,ymax,nPoints,pCol,pFill,packageurl,mean,pathway_class_y,pathway_class_plot_df,pathway_class_annotation,p_annotation,combination_bar_plot,grob,bar_errorbar2,test)
#######
p.adjust_bonferroni <- function (p, method = p.adjust.methods, n) 
{
  method <- match.arg(method)
  if (method == "fdr") 
    method <- "BH"
  nm <- names(p)
  p <- as.numeric(p)
  p0 <- setNames(p, nm)
  if (all(nna <- !is.na(p))) 
    nna <- TRUE
  else p <- p[nna]
  lp <- length(p)
  #stopifnot(n >= lp)
  if (n <= 1) 
    return(p0)
  if (n == 2 && method == "hommel") 
    method <- "hochberg"
  p0[nna] <- switch(method,
                    bonferroni = pmin(1, n * p),
                    none = p)
  return(p0)
}
#######
annotate_local <- function(daa_result){
  #daa_result <- daa_sub
  daa_result<-merge(x=daa_result,y=kegg_pathways,by.x = "feature",by.y = "KEGG")
  for(i in 1:length(daa_result[[1]])){
    buf<- unlist(strsplit(daa_result$Level3[i],split = ' '))
    strbuf<-""
    
    if (length(buf)==2) {
      strbuf <-buf[2]
    }else{
      for (j in 2:(length(buf)-1)) {
        strbuf <- paste(strbuf, buf[j])
      }
      strbuf <- substr(strbuf,2,nchar(strbuf))
      }
    daa_result$pathway_name[i] <- strbuf
  }
  for(i in 1:length(daa_result[[1]])){
    buf<- unlist(strsplit(daa_result$Level1[i],split = ' '))
    strbuf<-""
    if (length(buf)==2) {
      strbuf <-buf[2]
    }else{
      for (j in 2:length(buf)) {
        strbuf <- paste(strbuf, buf[j])
      }
      strbuf <- substr(strbuf,2,nchar(strbuf))
    }
    buf<- unlist(strsplit(daa_result$Level2[i],split = ' '))
    strbuf2<-""
    if (length(buf)==2) {
      strbuf2 <-buf[2]
    }else{
      for (j in 2:(length(buf))) {
        strbuf2 <- paste(strbuf2, buf[j])
      }
      strbuf2 <- substr(strbuf2,2,nchar(strbuf2))
    }
    
    strbuf<-paste(strbuf,'; ',sep = "")
    strbuf<-paste(strbuf,strbuf2,sep = " ")
    
    daa_result$pathway_class[i] <- strbuf
  }
  daa_result <- daa_result %>% dplyr::select(-c("Level1","Level2","Level3"))
  return(daa_result)
}
#######
calculate_ZScore_KO2kegg <- function(kegg_pathKO,
                                     KO_data,
                                     metadata,
                                     group_level){
  z_list <- list()
  for (m in 1:length(kegg_pathKO[[1]])) {
    l<- unlist(as.vector(kegg_pathKO[m,-1]))
    l <- na.omit(l)
    l <- as.vector(l)
    l <- KO_data %>% filter(ID %in% l)
    #l <- l$ID
    ID <- l$ID
    l <- l[,-1]
    rownames(l) <- ID
    l <- as.data.frame(t(l))
    l$sampleID <- rownames(l)
    metadata <- metadata %>% dplyr::select(c("sampleID","Grupo"))
    
    l <- metadata %>% dplyr::left_join(l, by="sampleID")
    ##############################
    #Build calculus buffer ZScore
    buf_ZScore <- data.frame(ID)
    buf_ZScore[,group_level] <- NA
    
    for(i in 1:length(group_level)){
      buf <- dplyr::select(l, -sampleID) %>% filter(Grupo == group_level[i])
      buf <- colMeans(buf[1:length(ID)+1])
      buf_ZScore[group_level[i]]<-buf
    }
    
    
    buf_ZScore[,"Mean"] <- apply(buf_ZScore[2:5],1,mean)
    buf_ZScore[,"sd"] <- apply(buf_ZScore[2:5],1,sd)
    ##############################
    # Build Z-Score dataframe
    ZScore <- dplyr::select(buf_ZScore, -Mean,-sd)
    for(i in 1:length(buf_ZScore$ID)){
      for (j in 1:length(group_level)) {
        ZScore[i,group_level[j]] <- (buf_ZScore[i,group_level[j]] - buf_ZScore[i,]$Mean)/buf_ZScore[i,]$sd
      }
    }
    z_list[[kegg_pathKO[[m,1]]]] <- ZScore
  }
  return(z_list)
}
#######
write_KOs_ZScore <- function(z_list,
                               path = "",
                               batch_id=""){
  for (i in 1:length(z_list)){
    #buf <- paste("images/picrust/ZScore/",batch_id,names(z_list)[i],"_Heatmap.tiff",sep = '')
    path_file <- paste(path,names(z_list)[i],"_",batch_id,".csv",sep = "")
    write.csv(z_list[[i]],path_file)
  }
}

# ********************************************************************************************* #
# ********************************************************************************************* #
#################################################################################################

######################################################################################################
# ************************************* Load Picrust files ***************************************** #
######## KOs ##############
path <- paste0(file.path("results", "legacy", "functional_prediction"), .Platform$file.sep)
# Mestrado: Dams e F1-PND21
KO_mest <- read_delim(paste(path,"Results_mestrado/KO2_PND21_Dams.tsv",sep = ""), delim = "\t")
KO_mest_ID <- KO_mest$ID
KO_mest <- KO_mest %>% dplyr::select(-c("ID","10F","22F"))
row.names(KO_mest)<-KO_mest_ID
# Alberto: F1 e F2-1year
KO_1year <- read_delim(paste(path,"Results_1year/KO2_1year_F1_F2.tsv",sep = ""), delim = "\t")
KO_1year_ID <- KO_1year$ID
KO_1year <- KO_1year %>% dplyr::select(-c("ID","t3"))
row.names(KO_1year)<-KO_1year_ID
######## ECs ##############
EC<- read_delim("Data_Metadata/Picrust/EC_pred.tsv", delim = "\t")
EC_ID <- EC$ID
EC <- EC %>% dplyr::select(-c("ID","10F","22F"))
row.names(EC)<-EC_ID

#########################################################################################################
# **************************************** Load metadata file ***************************************** #
path <- paste0(file.path("data", "metadata", "legacy"), .Platform$file.sep)
# Mestrado: Dams e F1-PND21
metadata_mest <- read_delim(paste(path,"Metadata_mehg.csv",sep = ""), delim = ";")
rownames(metadata_mest) <- metadata_mest$sampleID
metadata_mest <- metadata_mest[!(row.names(metadata_mest) %in% c("10F","22F")),]
rownames(metadata_mest) <- metadata_mest$sampleID
# Alberto: F1 e F2-1year
metadata_1year <- read_delim(paste(path,"metadata.csv",sep = ""), delim = ",")
rownames(metadata_1year) <- metadata_1year$sampleID
colnames(metadata_1year) <- c("sampleID","Geracao","Grupo","Sexo","obs")

group_level <- c("CN","VitA","MeHg","VitA+MeHg")
metadata_mest <- metadata_mest %>%  mutate(Grupo = factor(Grupo, levels = group_level))
metadata_1year <- metadata_1year %>%  mutate(Grupo = factor(Grupo, levels = group_level))
#########################################################################################################
# Subset KO data
# Dams and F1 PND21 (mestrado)
KO_mest_Dams <- KO_mest %>% dplyr::select(ID,
  as.vector(
    unlist(metadata_mest %>% filter(Geracao=="F0") %>% dplyr::select(sampleID))))
row.names(KO_mest_Dams)<-row.names(KO_mest)
KO_F1_PND21 <- KO_mest %>% dplyr::select(
  as.vector(
    unlist(metadata_mest %>% filter(Geracao=="F1") %>% dplyr::select(sampleID))))
row.names(KO_F1_PND21)<-row.names(KO_mest)
# Fa and F2 - 1year (alberto)
KO_F1_1year <- KO_1year %>% dplyr::select(
  as.vector(
    unlist(metadata_1year %>% filter(Geracao=="F1") %>% dplyr::select(sampleID))))
row.names(KO_F1_1year)<-row.names(KO_1year)
KO_F2_1year <- KO_1year %>% dplyr::select(
  as.vector(
    unlist(metadata_1year %>% filter(Geracao=="F2") %>% dplyr::select(sampleID))))
row.names(KO_F2_1year)<-row.names(KO_1year)
#########################################################################################################
####Principal Component Analysis PCA#####
#o que vai entre parenteses ? s? os grupos e sem diferenciar por amostra ex: CAF1, CAF2,CAF3 (ficaria tudo CAF)
# Realização de PCA

#pcdat<-prcomp(t(EC))

#Dams
pcdat <- prcomp(t(KO_mest_Dams))
metadata_dams<- metadata_mest %>% filter(Geracao=="F0")
pca_plot <- ggbiplot( pcdat, group=as.factor(metadata_dams$Grupo),
                      #labels = row.names(t(KO_mest)),     # Add labels as rownames
                      labels.size = 0.0004,
                      obs.scale = 1, var.scale = 1,
                      ellipse = TRUE, circle = FALSE,
                      var.axes = FALSE) + 
  scale_colour_discrete() +
  #  theme(legend.direction = 'horizontal', legend.position = 'top') +
  geom_point(aes(color = metadata_dams$Grupo),size = 2.5) +
  scale_fill_discrete(guide="none")+
  labs(color = "Group") +
  theme_bw()
pca_plot
tiff('images/pca_Picrust_KO_Dams.tiff', width = 6,height = 3, units = "in", res=300)
pca_plot
dev.off()

# F1 PND21 (mestrado)
pcdat<-prcomp(t(KO_F1_PND21))
metadata_F1_PND21<- metadata_mest %>% filter(Geracao=="F1")
pca_plot <- ggbiplot( pcdat, group=as.factor(metadata_F1_PND21$Grupo),
                      #labels = row.names(t(KO_mest)),     # Add labels as rownames
                      labels.size = 0.0004,
                      obs.scale = 1, var.scale = 1,
                      ellipse = TRUE, circle = FALSE,
                      var.axes = FALSE) + 
  scale_colour_discrete() +
  #  theme(legend.direction = 'horizontal', legend.position = 'top') +
  geom_point(aes(color = metadata_F1_PND21$Grupo,shape = metadata_F1_PND21$Sexo),size = 2.5) +
  scale_fill_discrete(guide="none")+
  labs(color = "Group",shape = "Sex") +
  theme_bw()
pca_plot
tiff('images/pca_Picrust_KO_F1_PND21.tiff', width = 6,height = 3, units = "in", res=300)
pca_plot
dev.off()

# F1 1year (alberto)
pcdat<-prcomp(t(KO_F1_1year))
metadata_F1_1year<- metadata_1year %>% filter(Geracao=="F1")
pca_plot <- ggbiplot( pcdat, group=as.factor(metadata_F1_1year$Grupo),
                      #labels = row.names(t(KO_mest)),     # Add labels as rownames
                      labels.size = 0.0004,
                      obs.scale = 1, var.scale = 1,
                      ellipse = TRUE, circle = FALSE,
                      var.axes = FALSE) + 
  scale_colour_discrete() +
  #  theme(legend.direction = 'horizontal', legend.position = 'top') +
  geom_point(aes(color = metadata_F1_1year$Grupo,shape = metadata_F1_1year$Sexo),size = 2.5) +
  scale_fill_discrete(guide="none")+
  labs(color = "Group",shape = "Sex") +
  theme_bw()
pca_plot
tiff('images/pca_Picrust_KO_F1_1year.tiff', width = 6,height = 3, units = "in", res=300)
pca_plot
dev.off()

# F2 1year (alberto)
pcdat<-prcomp(t(KO_F2_1year))
metadata_F2_1year<- metadata_1year %>% filter(Geracao=="F2")
pca_plot <- ggbiplot( pcdat, group=as.factor(metadata_F2_1year$Grupo),
                      #labels = row.names(t(KO_mest)),     # Add labels as rownames
                      labels.size = 0.0004,
                      obs.scale = 1, var.scale = 1,
                      ellipse = TRUE, circle = FALSE,
                      var.axes = FALSE) + 
  scale_colour_discrete() +
  #  theme(legend.direction = 'horizontal', legend.position = 'top') +
  geom_point(aes(color = metadata_F2_1year$Grupo,shape = metadata_F2_1year$Sexo),size = 2.5) +
  scale_fill_discrete(guide="none")+
  labs(color = "Group",shape = "Sex") +
  theme_bw()
pca_plot
tiff('images/pca_Picrust_KO_F2_1year.tiff', width = 6,height = 3, units = "in", res=300)
pca_plot
dev.off()

rm(pcdat)

#########################################################################################################
# ******************************************** Picrust ************************************************ #

######## KOs ##############
path <- paste0(file.path("results", "legacy", "functional_prediction"), .Platform$file.sep)
# Mestrado: Dams e F1-PND21
KO_mest <- read_delim(paste(path,"Results_mestrado/KO2_PND21_Dams.tsv",sep = ""), delim = "\t")
KO_mest_ID <- KO_mest$ID
KO_mest <- KO_mest %>% dplyr::select(-c("10F","22F"))
row.names(KO_mest)<-KO_mest_ID
# Alberto: F1 e F2-1year
KO_1year <- read_delim(paste(path,"Results_1year/KO2_1year_F1_F2.tsv",sep = ""), delim = "\t")
KO_1year_ID <- KO_1year$ID
KO_1year <- KO_1year %>% dplyr::select(-c("t3"))
row.names(KO_1year)<-KO_1year_ID
#############################################################################################
# Subset KO data
# Dams and F1 PND21 (mestrado)
KO_mest_Dams <- KO_mest %>% dplyr::select(ID,metadata_dams$sampleID)
KO_F1_PND21 <- KO_mest %>% dplyr::select(ID,metadata_F1_PND21$sampleID)
# F1 and F2 - 1year (alberto)
KO_F1_1year <- KO_1year %>% dplyr::select(ID,metadata_F1_1year$sampleID)
KO_F2_1year <- KO_1year %>% dplyr::select(ID,metadata_F2_1year$sampleID)
#############################################################################################
#Open Kegg pathways matrix from Metacyc
#PWY<- read_delim("Data_Metadata/Picrust/Pathway_abun_unstrat.tsv", delim = "\t")
#PWY <- dplyr::select(PWY,-'10F',-'22F')

#############################################################################################
# Transform KO_mest data to Kegg pathways
#KO
# mestrado
kegg_abund_dams <- ko2kegg_abundance(data = KO_mest_Dams)
kegg_abund_F1_PND21 <- ko2kegg_abundance(data = KO_F1_PND21)
# alberto
kegg_abund_F1_1year <- ko2kegg_abundance(data = KO_F1_1year)
kegg_abund_F2_1year <- ko2kegg_abundance(data = KO_F2_1year)
##################################################################################################
# Calculate differential abundant (DA) pathways
# The resulting kegg_abundance data frame can now be used for further analysis and visualization. The abundance table is better to be data.frame rather than tibble
#you can use ko2kegg_abundance output
#KO

#####################################

# OLHAR COLUNA NOME KEGG_ABUNDANCE_**
#abundance = kegg_abundance_** %>% column_to_rownames("ID")

#####################################

group <- "Grupo"
# mestrado
daa_KO_mest_dams <- pathway_daa(
  abundance = kegg_abund_dams,metadata = metadata_dams,group = group,
  daa_method = "ALDEx2",select = NULL,reference = NULL)
daa_KO_mest_dams$p_adjust <- p.adjust_bonferroni(daa_KO_mest_dams$p_values,"bonferroni",n=6)
#daa_KO_mest_dams$p_adjust <- 0.04
daa_KO_mest_dams$adj_method <- "bonferroni"

daa_KO_F1_PND21 <- pathway_daa(
  abundance = kegg_abund_F1_PND21,metadata = metadata_F1_PND21,group = group,
  daa_method = "ALDEx2",select = NULL,reference = NULL)
daa_KO_F1_PND21$p_adjust <- p.adjust_bonferroni(daa_KO_F1_PND21$p_values,"bonferroni",n=6)
daa_KO_F1_PND21$adj_method <- "bonferroni"

# alberto
daa_KO_F1_1year <- pathway_daa(
  abundance = kegg_abund_F1_1year,metadata = metadata_F1_1year,group = group,
  daa_method = "ALDEx2",select = NULL,reference = NULL)
daa_KO_F1_1year$p_adjust <- p.adjust_bonferroni(daa_KO_F1_1year$p_values,"bonferroni",n=6)
daa_KO_F1_1year$adj_method <- "bonferroni"

daa_KO_F2_1year <- pathway_daa(
  abundance = kegg_abund_F2_1year,metadata = metadata_F2_1year,group = group,
  daa_method = "ALDEx2",select = NULL,reference = NULL)
daa_KO_F2_1year$p_adjust <- p.adjust_bonferroni(daa_KO_F2_1year$p_values,"bonferroni",n=6)
daa_KO_F2_1year$adj_method <- "bonferroni"

# ***************************************************************************************** #
#KO
daa_sub <- daa_KO_F1_PND21 %>% filter(p_adjust <= 0.05)
daa_sub <- daa_KO_F1_1year %>% filter(p_adjust <= 0.1)
daa_sub <- daa_KO_F2_1year %>% filter(p_adjust <= 0.105)

daa_sub
#############################################################################################
# Subset to ensure the presence of only one method
# Please check the unique(daa_results_df$method) and choose one
#Metacyc and KO
daa_sub <- daa_sub %>% filter(method == "ALDEx2_glm test")
###########
# Annotate
daa_annotated <- annotate_local(daa_sub)

#daa_annotated <- pathway_annotation(
#  pathway = "KO",daa_results_df = daa_sub,ko_to_kegg = TRUE)
#daa_annotated$p_adjust <- p.adjust_bonferroni(daa_annotated$p_values,"bonferroni",n=6)

#############################################################################################
#Subset only significant DAs
#Metacyc and KO

#daa_sub_method_results_male_sig <- subset(daa_sub_method_results_male,p_adjust < 0.05)
#######################################################

# Please change Group to metadata$your_group_column if you are not using example dataset
#Group <- metadata_male$Grupo
#Group <- metadata_female$Grupo
Group <- metadata_F1_PND21$Grupo
Group <- metadata_F1_1year$Grupo
Group <- metadata_F2_1year$Grupo

daa_annotated_final <- dplyr::select(daa_annotated,-group3,-group4)
###########################
#KO

daa_annotated_filt <- daa_annotated_final %>% filter(feature != c("ko00514","ko00944","ko04011","ko04962","ko05219"))
daa_annotated_filt <- daa_annotated_filt %>% filter(feature != c("ko00514","ko00944"))

p <- pathway_errorbar3(
  abundance = kegg_abund_F1_PND21,      #change data
  daa_results_df = daa_annotated_filt,
  Group = Group,
  p_values_threshold = 2,
  order = "pathway_class",
  #select = daa_annotated_final %>%
  #  arrange(p_adjust) %>%
  #  slice(1:15) %>%
  #  dplyr::select("feature") %>% pull(),
  ko_to_kegg = TRUE,
  p_value_bar = FALSE,
  colors =  c("#98F5FF","#B3EE3A","#CD2626","#4D4D4D"),
  x_lab = "pathway_name")
p

tiff('images/picrust/Kegg_pathways_F2_1year_FILT.tiff', width = 15,height = 6, units = "in", res=150)
p
dev.off()
# ***
# ***
# ***
# ***
# ***
# ***
# ***
# ***
# ***
# ***
# ***
# ***
# ***
# ***
################################################################################################
################################################################################################
# GSEA Analysis
buf <- KO_F1_PND21 %>% dplyr::select(-c("ID"))
row.names(buf)<-KO_F1_PND21$ID
rownames(metadata_F1_PND21)<-metadata_F1_PND21$sampleID

gsea_F1_PND21 <- pathway_gsea(
  abundance = buf,
  metadata = metadata_F1_PND21,
  group = "Sexo",
  pathway_type = "KEGG",
  method = "fgsea",
  rank_method = "signal2noise"
)
head(gsea_F1_PND21)

# Annotate GSEA results
annotated_results <- gsea_pathway_annotation(
  gsea_results = gsea_F1_PND21,
  pathway_type = "KEGG"
)

# View the annotated results
head(annotated_results)
# Option 2: Use annotated results (automatically shows pathway names)
plot_gsea <- visualize_gsea(
  gsea_results = annotated_results,
  plot_type = "barplot",
  n_pathways = 10
)
# Option 3: Explicitly specify which column to use for labels
plot_gsea <- visualize_gsea(
  gsea_results = annotated_results,
  plot_type = "barplot",
  pathway_label_column = "pathway_name.y",
  #n_pathways = 10,
  sort_by = "p.adjust",
)
plot_gsea
#####################################
# Annotate
annotated_results <- annotated_results %>% dplyr::select(c("pathway_id","size",
                                                           "ES","NES","pvalue","p.adjust","method"))
colnames(annotated_results) <- c("feature","size",
                                 "ES","NES","p_values","p_adjust","method")
annotated_results$p_adjust <- p.adjust_bonferroni(annotated_results$p_values,"bonferroni",n=6)
annotated_results <- annotated_results %>% filter(p_adjust <= 0.01)

annotated_results <- annotate_local(annotated_results)
annotated_results$group1 <- "CN"
annotated_results$group2 <- "VitA"

p <- pathway_errorbar3(
  abundance = kegg_abund_F1_PND21,      #change data
  daa_results_df = annotated_results,
  Group = metadata_F1_PND21$Grupo,
  p_values_threshold = 2,
  order = "pathway_class",
  #select = annotated_results %>%
  #  arrange(p_adjust) %>%
  #  slice(1:15) %>%
  #  dplyr::select("feature") %>% pull(),
  ko_to_kegg = TRUE,
  p_value_bar = FALSE,
  colors =  c("#98F5FF","#B3EE3A","#CD2626","#4D4D4D"),
  x_lab = "pathway_name")
p
################################################################################################
################################################################################################
kegg_path_F1_PND21 <- daa_KO_F1_PND21 %>% filter(p_adjust <= 0.05 & method == "ALDEx2_glm test") %>% dplyr::select(feature) %>% pull()
kegg_path_F1_1year <- daa_KO_F1_1year %>% filter(p_adjust <= 0.103 & method == "ALDEx2_glm test") %>% dplyr::select(feature) %>% pull()
kegg_path_F2_1year <- daa_KO_F2_1year %>% filter(p_adjust <= 0.103 & method == "ALDEx2_glm test") %>% dplyr::select(feature) %>% pull()

kegg_pathKO_F1_PND21 <- ko_to_kegg_reference %>% filter(X1 %in% kegg_path_F1_PND21)
kegg_pathKO_F1_1year <- ko_to_kegg_reference %>% filter(X1 %in% kegg_path_F1_1year)
kegg_pathKO_F2_1year <- ko_to_kegg_reference %>% filter(X1 %in% kegg_path_F2_1year)

kegg_pathKO<-kegg_pathKO_F1_1year
KO_data<-KO_F1_1year
metadata<-metadata_F1_1year


##########
# Calcular ZScore de todas vias relevantes
# Para isso, metadata e kegg_pathKO precisam ambos ter uma coluna nomeada "sampleID"
# com os IDs das amostras
# Além disso, a função também realiza uma exclusão de colunas nomeadas como : 
# "Sexo", "Geracao", "Geracao_Sexo" e "Grupo_Geracao_Sexo"

# A função devolve uma lista de cada via linkada a um dataframe do ZScore de todos
# KOs presentes nos dados analisados
z_list <- calculate_ZScore_KO2kegg(kegg_pathKO = kegg_pathKO_F1_PND21,
                         KO_data = KO_F1_PND21,
                         metadata = metadata_F1_PND21,
                         group_level = group_level)
z_list <- calculate_ZScore_KO2kegg(kegg_pathKO = kegg_pathKO_F1_1year,
                                   KO_data = KO_F1_1year,
                                   metadata = metadata_F1_1year,
                                   group_level = group_level)
z_list <- calculate_ZScore_KO2kegg(kegg_pathKO = kegg_pathKO_F2_1year,
                                   KO_data = KO_F2_1year,
                                   metadata = metadata_F2_1year,
                                   group_level = group_level)
##############################
# Plot Z-Score
ZScore <- z_list[[1]]
rownames(ZScore) <- ZScore$ID
ZScore <- dplyr::select(ZScore,-ID)

y <- data.matrix(ZScore) #transformar o objeto numa matriz#
y <- na.omit(y)

library(gplots)
yb <-colorRampPalette(c("blue", "white", "red"))
#aqui pode transformar diretamente a planilha numa matriz, sem transform?-la num objeto antes#
heatmap.2(y) #pedir heatmap da matriz# #erro de margens muito grandes ? comum#

tiff('images/picrust/F1_PND21_ko00740_Heatmap.tiff', width = 10,height = 25, units = "in", res=150)
heatmap.2(y,key= FALSE, col=yb, Rowv = TRUE, Colv = FALSE,dendrogram = NULL,denscol = NULL,
          trace = 'none',cexRow = 3, cexCol = 3,density.info = 'none',key.xlab = NA, keysize = 0.1,
          margins = c(20,10))
dev.off()

ZScore$CN
buf<- data.frame(KO = rownames(ZScore),
           CN = ZScore$CN)
buf <- ZScore
buf$KO <- rownames(ZScore)
write.csv(buf,"images/picrust/color_ko00190.csv")
############
# For Loop
library(gplots)
yb <-colorRampPalette(c("blue", "white", "red"))
batch_id = "F1_PND21_"
batch_id = "F1_1year_"
batch_id = "F2_1year_"
for (i in 1:length(z_list)){
  ZScore <- z_list[[i]]
  rownames(ZScore) <- ZScore$ID
  ZScore <- dplyr::select(ZScore,-ID)
  
  y <- data.matrix(ZScore) #transformar o objeto numa matriz#
  y <- na.omit(y)
  if(length(y[,1]) == 1){
    next
  }else if(length(y[,1]) > 1 & length(y[,1]) < 5){
    buf <- paste("images/picrust/ZScore/",batch_id,names(z_list)[i],"_Heatmap.tiff",sep = '')
    tiff(buf, width = 4,height = 2, units = "in", res=150)
    heatmap.2(y,key= FALSE, col=yb, margins = c(6,6), Rowv = TRUE, Colv = FALSE,
              dendrogram = NULL, trace = 'none', denscol = NULL, cexRow = 1, cexCol = 1,
              key.title = "Z-Score",density.info = 'none',key.xlab = NA, keysize = 0.1)
    dev.off()
  }else if(length(y[,1]) >= 5 & length(y[,1]) < 12){
    buf <- paste("images/picrust/ZScore/",batch_id,names(z_list)[i],"_Heatmap.tiff",sep = '')
    tiff(buf, width = 4,height = 4, units = "in", res=150)
    heatmap.2(y,key= FALSE, col=yb, margins = c(8,8), Rowv = TRUE, Colv = FALSE,
              dendrogram = NULL, trace = 'none', denscol = NULL, cexRow = 1, cexCol = 1,
              key.title = "Z-Score",density.info = 'none',key.xlab = NA, keysize = 0.1)
    dev.off()
  }else{
    buf <- paste("images/picrust/ZScore/",batch_id,names(z_list)[i],"_Heatmap.tiff",sep = '')
    tiff(buf, width = 4,height = 12, units = "in", res=150)
    heatmap.2(y,key= FALSE, col=yb, margins = c(8,8), Rowv = TRUE, Colv = FALSE,
              dendrogram = NULL, trace = 'none', denscol = NULL, cexRow = 1, cexCol = 1,
              key.title = "Z-Score",density.info = 'none',key.xlab = NA, keysize = 0.1)
    dev.off()
  }
}
#######
write_KOs_ZScore(z_list = z_list,
                   path = "images/picrust/ZScore/",
                   batch_id = "_F1_PND21")
write_KOs_ZScore(z_list = z_list,
                   path = "images/picrust/ZScore/",
                   batch_id = "_F1_1year")
write_KOs_ZScore(z_list = z_list,
                   path = "images/picrust/ZScore/",
                   batch_id = "_F2_1year")
################################################################################################
################################################################################################
# KO Mercuric

mercuric_KOs <- read_delim("Data_Metadata/KEGG/KEGG_rest.kegg.jp_find_KO_mercuric.txt", delim = "\t")
for(i in 1:length(mercuric_KOs$KO)){
  mercuric_KOs$KO[i] <- strsplit(mercuric_KOs$KO[i],split=":")[[1]][2]
}
KO_mercuric <- KO
rownames(KO_mercuric) <- KO$ID
KO_mercuric <- KO_mercuric[rownames(KO_mercuric) %in% mercuric_KOs$KO,]
#rownames(KO_mercuric) <- KO_mercuric$ID
KO_ID <- KO_mercuric$ID
KO_mercuric <- dplyr::select(KO_mercuric,-ID)
KO_mercuric <- t(KO_mercuric)
colnames(KO_mercuric) <- KO_ID
KO_mercuric <- as.data.frame(KO_mercuric)
KO_mercuric$SampleID <- rownames(KO_mercuric)

KO_mercuric <- metadata %>% left_join(KO_mercuric, by = "SampleID")
KO_mercuric$Sum <- rowSums(KO_mercuric[7:9])
View(KO_mercuric)
rownames(mercuric_KOs) <- mercuric_KOs$KO

stat.test <- KO_mercuric %>%
  group_by(Geracao_Sexo) %>%
  t_test(Sum ~ Grupo,p.adjust.method="fdr")
stat.test

colnames(KO_mercuric[7:9])
stat.test <- KO_mercuric %>%
  group_by(Geracao_Sexo) %>%
  t_test(K08365 ~ Grupo,p.adjust.method="fdr",comparisons = list(c("CN","VitA"),c("CN","MeHg"),c("CN","VitA+MeHg")))
stat.test
stat.test <- stat.test[!stat.test$p.adj.signif=="ns",]
stat.test

#Anova geracao_sexo
p <- compare_means(K19058 ~ Geracao_Sexo, data = KO_mercuric, p.ajust.method = "fdr",method = "anova")
p$p.adj

# Box plots with p-values
stat.test <- stat.test %>% add_y_position()
#stat.test$y.position<-0.05
f <- ggboxplot(KO_mercuric, x = "Grupo", y = "K19058",palette = c("#98F5FF","#B3EE3A","#CD2626","#4D4D4D"),
               #outlier.shape = NA,
               fill = "Grupo", facet.by = "Geracao_Sexo") +
  #stat_pvalue_manual(stat.test, label = "p.adj.signif", tip.length = 0.01) +
  #coord_cartesian(ylim = c(0, 0.05))+
  theme(axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        legend.title = element_blank(),
        legend.direction='vertical',
        legend.position = "right",
        plot.subtitle = element_text(hjust = 1),
        panel.background = element_rect(fill = "gray93",colour = "grey", linetype = "solid"),
        panel.grid.major = element_line(size = 0.5, linetype = 'solid',colour = "white"), )+
  ggtitle(paste(mercuric_KOs["K19058",2],mercuric_KOs["K19058",3]),
          subtitle = paste("Anova =",p$p.adj))+
  xlab("Groups") +
  ylab("Picrust Score")

f

tiff("Figures/Picrust/MercuricKOs_K19058.tiff",width = 8,height = 5,units = "in", res=300)
f
dev.off()
################################################################################################

kegg_pathways <-load(system.file("extdata", "KO_reference.RData", package = "ggpicrust2"))
