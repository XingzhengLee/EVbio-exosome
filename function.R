
# load libraries ----------------------------------------------------------
library(tidyverse)
library(RColorBrewer)
library(cowplot)
library(pheatmap)
library(viridis)
library(ComplexHeatmap)
library(factoextra)
library(corrplot)
library(ggrepel)
library(leaps)
library(plotROC)
library(pROC)
library(sva)
library(survival)
library(survminer)

# input -------------------------------------------------------------------
# defaule parameters: quantification_file = 'Quantification.txt', protein_file = 'Protein.list.txt', sample_file = 'Sample.list.txt', fluorescence_label = '635'
file.in <- function(quantification_pattern = 'NULL', protein_file = 'Protein.list.txt', sample_file = 'Sample.list.txt', fluorescence_label = '635'){
  if (quantification_pattern == 'NULL') {
    cat("Please specify parameter 'quantification_pattern' to add quantification files.")
  } else {
    # protein list annotation
    protein_list <<- read_tsv(protein_file, col_names = F)%>%
      # add rows and cols (coordinate)
      dplyr::mutate(Row = row_number())%>%
      gather(Column, Protein, -Row)%>%
      dplyr::mutate(
        Column = str_sub(Column, 2),
        Column = as.integer(Column))
    # sample list annotation
    sample_list <<- read_tsv(
      sample_file,
      # returns the currently used default encoding
      locale = locale(
        encoding = stringi::stri_enc_get()
      )
    )
    # quantification raw data
    Qlist <- list.files(pattern = quantification_pattern)
    quantification <- data.frame()
    for (i in 1:length(Qlist)) {
      Q_single <- read_tsv(Qlist[i])%>%
        dplyr::select(
          Block,
          Column,
          Row,
          fluorescence_F = str_c('F', fluorescence_label, ' Median'),
          fluorescence_B = str_c('B', fluorescence_label),
          fluorescence_FT = str_c('F', fluorescence_label, ' Total Intensity')
        )%>%
        mutate(Batch = i)
      quantification <- bind_rows(quantification, Q_single)%>%
        tbl_df()
    }
    # raw data annotation (conbine data)
    raw_data <<- quantification%>%
      inner_join(
        # protein annotation
        protein_list, by = c('Row', 'Column')
      )%>%
      # sample annotation
      inner_join(sample_list, by = c('Batch', 'Block'))%>%
      dplyr::select(Batch, Block, Column, Row, SampleID, SampleName, Group, Protein, fluorescence_F, fluorescence_B, fluorescence_FT, everything())%>%
      mutate(Block = as.factor(Block))
    return(raw_data)
  }
}

# CV ----------------------------------------------------------------------
# CV calculation
# default parameters: data_input = raw_data
CV.cal <- function(data_input = raw_data){
  # raw data assessment (CV)
  CV <<- data_input%>%
    filter(
      # remove '+' and 'PBS'
      !(Protein %in% c('+', 'PBS')),
      Group != 'PBS'
    )%>%
    # CV calculation
    dplyr::group_by(Batch, Block, Protein)%>%
    dplyr::summarise(
      CV_F = sd(fluorescence_F)/mean(fluorescence_F),
      CV_B = sd(fluorescence_B)/mean(fluorescence_B),
      CV_FT = sd(fluorescence_FT)/mean(fluorescence_FT)
    )%>%
    ungroup()
  return(CV)
}
# CV violin
# default parameters: data_input = raw_data, type = 'CV_F|CV_B|CV_FT'
CV.violin <- function(type = 'NULL'){
  if (type == 'CV_F' | type == 'CV_B' | type == 'CV_FT'){
    col_Set2 <- brewer.pal(8, 'Set2')
    sample_num <- nrow(distinct(CV, Block))
    col_violin <- rep(col_Set2[1:3], sample_num)
    plot_suffix <- list(
      geom_violin(aes(fill = Block), trim = F),
      geom_boxplot(width = 0.1),
      facet_wrap(~Batch, ncol = 1),
      theme_minimal(),
      theme(
        legend.position = 'none',
        panel.grid.major.x = element_blank()
      ),
      scale_fill_manual(values = col_violin)
    )
    CV%>%
      dplyr::select(Batch, Block, value = type)%>%
      ggplot(aes(Block, value, group = Block)) +
      labs(y = type) +
      plot_suffix
  } else {cat("The value for the parameter 'type' is incorrect!\nPlease choose one of the following: 'CV_F, CV_B, CV_FT'.")}
}
# CV top
# default parameter: data_input = CV, topN = 20
CV.top <- function(data_input = CV, topN = 20){
  # unite Sample (block) and Protein name
  CVadd <- data_input%>%
    unite(Batch_Block_Protein, c('Batch', 'Block', 'Protein'), sep = '_', remove = F)
  # bar plot parameters
  plot_suffix <- list(
    geom_bar(stat = 'identity', fill = 'grey', color = 'grey50'),
    theme_minimal(),
    theme(legend.position = 'none'),
    coord_flip(),
    labs(x = 'Batch_Block_Protein')
  )
  # bar function
  Top_bar <- function(label){
    TopCV <- CVadd%>%
      dplyr::rename(value = label)%>%
      arrange(desc(value))%>%
      dplyr::slice(1:topN)
    TopCV%>%
      ggplot(aes(reorder(Batch_Block_Protein, value), value)) +
      labs(y = label) +
      plot_suffix
  }
  # plot bar
  CV_bar_1 <- Top_bar('CV_F')
  CV_bar_2 <- Top_bar('CV_B')
  CV_bar_3 <- Top_bar('CV_FT')
  plot_grid(CV_bar_1, CV_bar_2, CV_bar_3, nrow = 1)
}

# quality control ---------------------------------------------------------
# data for analysis (several options)
# default parameter: data_input = raw_data, option = '1|2|3'
data4analysis <- function(data_input = raw_data, option = 'NULL'){
  # normalization & reshape
  value_avg <- data_input%>%
    dplyr::group_by(Batch, Block, SampleID, Protein)%>%
    # mean value in two/three replicates
    dplyr::summarise(
      mean_F = mean(fluorescence_F),
      mean_B = mean(fluorescence_B),
      mean_FT = mean(fluorescence_FT)
    )%>%
    ungroup()
  if (option == 1) {
    Avalue <- value_avg%>%
      dplyr::mutate(avg = mean_F - mean_B)
  } else if (option == 2) {
    Avalue <- value_avg%>%
      dplyr::mutate(avg = mean_F)
  } else if (option == 3) {
    Avalue <- value_avg%>%
      dplyr::mutate(avg = mean_FT)
  } else {cat("The value for the parameter 'type' is incorrect!\nPlease choose one of the following: 1, 2, 3.")}
  Dvalue <<- Avalue%>%
    dplyr::select(-starts_with('mean'))
  return(Dvalue)
}
# negative control
# default parameter: data_input = Dvalue
negative.control <- function(data_input = Dvalue){
  # positive block (filter out 'PBS')
  block_positive <- data_input%>%
    dplyr::filter(SampleID != 'PBS')
  # negative block (protein correction)
  correction4protein <- data_input%>%
    dplyr::filter(SampleID == 'PBS')%>%
    dplyr::select(Batch, Protein, avg_negative4protein = avg)
  # negative sample in each block (block correction)
  correction4block <- block_positive%>%
    dplyr::filter(Protein == 'PBS')%>%
    dplyr::select(Batch, Block, avg_negative4block = avg)
  # data correction
  normalizated_data <- inner_join(block_positive, correction4protein, by = c('Batch', 'Protein'))%>%
    inner_join(correction4block, by = c('Batch', 'Block'))%>%
    # normalize for protein & block
    dplyr::mutate(norvalue = avg - avg_negative4protein - avg_negative4block)%>%
    dplyr::filter(!(Protein %in% c('PBS', '+')))%>%
    dplyr::select(-starts_with('avg'))
  # MartixQ
  MatrixQ <<- normalizated_data%>%
    dplyr::select(SampleID, Protein, norvalue)%>%
    spread(SampleID, norvalue)%>%
    tibble::column_to_rownames(var = 'Protein')
  return(MatrixQ)
}

# batch effect adjustment --------------------------------------------------
# defaule parameters: proein_expression = MatrixQ, phenotype = sample_list, type = 'box|density|data'
batch_adjust <- function(proein_expression = MatrixQ, phenotype = sample_list, type = 'NULL'){
  # load data
  pheno <- phenotype%>%
    filter(SampleID != 'PBS')%>%
    arrange(SampleID)
  edata <- as.matrix(proein_expression)
  # adjusting for batch effects with Combat
  batch <- pheno$Batch
  modcombat <-  model.matrix(~1, data = pheno)
  combat_edata <-  ComBat(dat = edata, batch = batch, mod = modcombat, par.prior = TRUE, prior.plots = FALSE)
  MatrixQ_adjusted <<- as.data.frame(combat_edata)
  # plot preparation
  anno <- pheno%>%
    select(Batch, SampleID)
  pre4plot <- function(data_input, label){
    as.data.frame(t(data_input))%>%
      tibble::rownames_to_column(var = 'SampleID')%>%
      tbl_df()%>%
      gather(protein, value, -SampleID)%>%
      left_join(anno, by = 'SampleID')%>%
      mutate(
        Batch = as.factor(Batch),
        type = label
      )
  }
  plot_none <- pre4plot(edata, 'Before adjustment')
  plot_ComBat <- pre4plot(combat_edata, 'ComBat adjustment')
  plot_bind <- bind_rows(plot_none, plot_ComBat)
  if (type == 'box') {
    ggplot(plot_bind, aes(Batch, value)) +
      geom_violin(aes(fill = Batch), trim = F) +
      geom_boxplot(width = 0.1) +
      theme_minimal() +
      theme(legend.position = 'none') +
      facet_wrap(~type) +
      scale_fill_brewer(palette = 'Spectral') +
      labs(
        x = NULL,
        y = 'Protein quantification'
      )
  } else if (type == 'density') {
    ggplot(plot_bind, aes(value)) +
      geom_density(aes(color = Batch), size = 1) +
      theme_minimal() +
      theme(legend.position = 'none') +
      facet_wrap(~type) +
      scale_color_brewer(palette = 'Spectral') +
      labs(x = NULL)
  } else if (type == 'data') {
    return(MatrixQ_adjusted)
  }
}

# error control -----------------------------------------------------------
# the percentage of positive values
# default parameter: data_input =  MatrixQ_adjusted
positive.percentage <- function(data_input =  MatrixQ_adjusted){
  # number of values above zero (percentage)
  Matrix_TF <- (data_input > 0)
  num_0 <- table(Matrix_TF)[2]
  num_all <- length(Matrix_TF)
  num_percentage <- round(num_0/num_all*100, 2)
  # print percentage
  cat(str_c(num_percentage, '% of the original values is positive.'))
}
# matrix translation and log transformation
# default parameter: data_input = MatrixQ_adjusted, drop_percentage
translation2log <- function(data_input =  MatrixQ_adjusted, drop_percentage){
  # find the threshold based on drop_percentage
  Matrix_vector <- sort(as.vector(as.matrix(data_input)))
  num_all <- length(Matrix_vector)
  num_vector <- round(num_all*drop_percentage)
  Matrix_threshold <- Matrix_vector[num_vector]
  data_input[data_input < Matrix_threshold] <- Matrix_threshold
  # matrix translation
  MatrixQ_positive <- 1 + data_input - min(data_input)
  MatrixQ_log <<- log2(MatrixQ_positive)
  return(MatrixQ_log)
}

# heatmap -----------------------------------------------------------------
# default parameters:data_input = MatrixQ_log, type = 'global|sample', sample_cluster = T/F
Matrix2heat <- function(data_input = MatrixQ_log, type = 'NULL', sample_cluster = 'NULL', protein_cluster = 'NULL'){
  data4heat <- data_input
  if (type == 'global') {
    if (sample_cluster != 'NULL') {
      # heatmap_global
      data4heat%>%
        pheatmap(
          scale = 'none',
          border_color = 'white',
          fontsize = 9, fontsize_row = 10, fontsize_col = 10,
          cluster_rows = protein_cluster, cluster_cols = sample_cluster
        )
    } else {cat("Please add parameter 'sample_cluster' (T/F).")}
  } else if (type == 'sample') {
    # sample correlation heatmap
    cor_sample <- cor(data_input, use = 'pairwise.complete.obs')
    cor_sample%>%
      pheatmap(
        color = rev(inferno(100)), border_color = NA,
        cluster_rows = T, cluster_cols = T
      )
  } else {cat("Please add parameter 'type' to display the data.")}
}

# k-means -----------------------------------------------------------------
# protein scale
# default parameters: data_input = MatrixQ_log
scale4prot <- function(data_input = MatrixQ_log){
  data_scale <<- as.data.frame(t(apply(data_input, 1, scale)))
  colnames(data_scale) <<- colnames(data_input)
  return(data_scale)
}
# kmeans 1:10
# default parameters: data_input = data_scale
kmeanstop <- function(data_input = data_scale) {
  wss <- NA
  for (i in 1:10) {
    wss[i] <- sum(kmeans(na.omit(data_input), centers = i, nstart = 24)$withinss)
  }
  as.data.frame(wss)%>%
    dplyr::mutate(cluster = row_number())%>%
    ggplot(aes(cluster, wss)) +
    geom_line(linetype = 'dashed') +
    geom_point(fill = 'white', color = 'black') +
    scale_x_continuous(breaks = c(1:8)) +
    theme_bw() +
    theme(panel.grid.major.x = element_blank(),
          panel.grid.minor.x = element_blank()) +
    labs(x = 'Number of Clusters',
         y = 'Within groups sum of squares')
}
# kmeans heatmap
# default parameter: data_intput = data_scale, kmeans_value
kmeansHeat <- function(data_intput = data_scale, kmeans_value){
  Heatmap(na.omit(data_intput), name = "heatmap",
          km = kmeans_value,
          column_names_side = "bottom",
          col = rev(brewer.pal(11, "Spectral")),
          cluster_columns = T,
          row_dend_side = "left",
          show_row_names = T)
}
# kmeans cluster
# default parameter: data_input = data_scale, kmeans_value, type = 'cluster|dend'
kmeansCluster <- function(data_input = data_scale, kmeans_value, type = 'NULL'){
  K <- kmeans(na.omit(data_input), kmeans_value, nstart = 24)
  clusterCol <- brewer.pal(9, 'Set1')[1:max(K$cluster)]
  if (type == 'cluster') {
    fviz_cluster(K, data = na.omit(data_input),
                 palette = clusterCol,
                 ellipse.type = "euclid",
                 star.plot = TRUE, 
                 repel = TRUE,
                 ggtheme = theme_minimal()
    )
  } else if (type == 'dend') {
    distRes <- dist(na.omit(data_input), method = "euclidean")
    hcRes <- hclust(distRes, method = "ward.D2")
    fviz_dend(hcRes, k = kmeans_value,
              k_colors = clusterCol,
              color_labels_by_k = T,
              rect = T)
  } else {cat("The value for the parameter 'type' is incorrect!\nPlease choose one of the following: 'cluster, dend'.")}
}

# principal components analysis (PCA) -------------------------------------
# PCA_normal+circle
# defaule parameters: data_input = MatrixQ_log, type = 'Group|Subgroup1|Subgroup2'
PCAplot <- function(data_input = MatrixQ_log, type = 'NULL'){
  PCA <- prcomp(na.omit(t(data_input)))
  PC1_importance <- round(summary(PCA)$importance[2,1]*100, 2)
  PC2_importance <- round(summary(PCA)$importance[2,2]*100, 2)
  pcaData <- data.frame(PCA$x[,1:2], check.names = F)%>%
    tibble::rownames_to_column(var = 'SampleID')%>%
    inner_join(sample_list, by = 'SampleID')%>%
    rename(condition = type)
  ggplot(pcaData, aes(x = PC1, y = PC2, color = condition)) +
    geom_point(size = 5, alpha = 2/3) +
    stat_ellipse(linetype = 'dashed', size = 1) +
    theme_bw() +
    scale_color_brewer(palette = 'Set1') +
    guides(color = guide_legend(title = type)) +
    labs(x = str_c('PC1 (', PC1_importance, '%)'),
         y = str_c('PC1 (', PC2_importance, '%)'))
}
# PCA_deep
# defaule parameters: data_input = MatrixQ_log, type = 'data|eigenvalue|variables|individuals|correlation|both', group = 'Group|Subgroup1|Subgroup2'
PCAdeep <- function(data_input = MatrixQ_log, type = 'NULL', group = 'NULL'){
  PCA <- prcomp(na.omit(t(data_input)))
  if (type == 'data') {
    get_eigenvalue(PCA)
  } else if (type == 'eigenvalue') {
    fviz_eig(PCA,
             addlabels = TRUE,
             barfill = "white",
             barcolor = "black",
             linecolor = "black")
  } else if (type == 'variables') {
    fviz_pca_var(PCA,
                 col.var = "contrib",
                 gradient.cols = brewer.pal(11, "Spectral"))
  } else if (type == 'correlation') {
    var <- get_pca_var(PCA)
    corrplot(var$cos2, is.corr = FALSE,
             method = 'shade',
             col = brewer.pal(9, "OrRd"),
             tl.cex = 0.7,
             tl.col = 'grey20',
             cl.pos = 'b',
             cl.cex = 0.5)
  } else if (type == 'individuals') {
    dt <- na.omit(t(data_input))
    dt_anno <- as.data.frame(dt)%>%
      tibble::rownames_to_column(var = 'SampleID')%>%
      left_join(sample_list, by = 'SampleID')%>%
      rename(condition = group)
    fviz_pca_ind(prcomp(dt),
                 geom.ind = "point",
                 col.ind = dt_anno$condition,
                 palette = brewer.pal(9, "Set1")[1:length(unique(dt_anno$condition))],
                 # whether to use ggrepel to avoid overplotting text labels
                 repel = T,
                 addEllipses = T,
                 legend.title = group)
  } else {cat("The value for the parameter 'type' is incorrect!\nPlease choose one of the following: 'data, eigenvalue, variables, individuals, correlation'.")}
}

# sub group ---------------------------------------------------------------
# defaule parameters: data_input = MatrixQ_log, annotation = sample_list, category_main = 'Group', object_main = c('SD', 'BL')
Matrix.sub <- function(data_input = MatrixQ_log, annotation = sample_list, category_main, object_main){
  annotation_sub <- annotation%>%
    dplyr::select(SampleID, condition = category_main)%>%
    dplyr::filter(condition %in% object_main)%>%
    distinct(SampleID, .keep_all = T)%>%
    filter(SampleID != 'PBS')
  # subset
  ID <- colnames(data_input)
  ID_sub <- annotation_sub$SampleID
  ID4join <- intersect(ID, ID_sub)
  MatrixQ_sub <<- data_input%>%
    dplyr::select(ID4join)
  return(MatrixQ_sub)
}

# differences between groups (single categroy) ----------------------------
# difference test
# default parameter: data_input = MatrixQ_log, object = c('A_BL', 'C_PD'), p_threshold = 0.05, fc_threshold = 0, test = 'nonpar|Ttest', type = 'data|plot'
diff.test <- function(data_input = MatrixQ_log, anno_input = raw_data, category, object, p_threshold = 0.05, fc_threshold = 0, test = 'NULL', type = 'NULL'){
  # add condition annotation
  annotation <- anno_input%>%
    dplyr::select(SampleID, condition = category)%>%
    distinct(SampleID, .keep_all = T)%>%
    filter(SampleID != 'PBS')
  Matrix_anno <- as.data.frame(t(data_input))%>%
    tibble::rownames_to_column(var = 'SampleID')%>%
    left_join(annotation, by = 'SampleID')
  Matrix_ID <- Matrix_anno%>%
    dplyr::select(-SampleID)
  # specify conditions for analysis
  DE4test <- Matrix_ID%>%
    filter(condition %in% object)
  # DE test
  diff_test <- data.frame(Protein = NA, p.value = NA, log2FC = NA)
  for (i in 1:(ncol(DE4test)-1)) {
    print(paste('Processing proteins:', i))
    # log2mean
    data4mean <- DE4test%>%
      dplyr::select(i, ncol(DE4test))
    colnames(data4mean) <- c('value', 'condition')
    mean <- group_by(data4mean, condition)%>%
      dplyr::summarise(mean = mean(value, na.rm = T))
    logvalue <- log2(mean$mean[2]/mean$mean[1])
    # creat result dataframe
    diff_test[i,1] <- colnames(DE4test)[i]
    diff_test[i,2] <- test(DE4test[,i] ~ DE4test[,ncol(DE4test)], data = DE4test, na.action = na.omit)$p.value
    diff_test[i,3] <- logvalue
  }
  DEP <<- arrange(diff_test, p.value)
  if (type == 'data') {
    return(DEP)
  } else if (type == 'plot') {
    # volcano plot
    FCcol <- brewer.pal(9, 'Set1')
    diff_test%>%
      mutate(
        logP = -log10(p.value),
        label = case_when(p.value < p_threshold & (log2FC > fc_threshold | log2FC < -fc_threshold) ~ Protein),
        FoldChange = ifelse(p.value < p_threshold & log2FC > fc_threshold, 'Up',
                            ifelse(p.value < p_threshold & log2FC < -fc_threshold, 'Down', 'Not'))
      )%>%
      ggplot(aes(log2FC, logP)) +
      geom_point(aes(color = FoldChange)) +
      geom_vline(xintercept = 0, linetype = 'dashed', color = 'grey50') +
      theme_test() +
      scale_color_manual(
        values = c(FCcol[1], 'grey', FCcol[2]),
        breaks = c('Down', 'Not', 'Up')
      ) +
      ggrepel::geom_text_repel(aes(label = label), box.padding = 0.1, label.padding = 0.1, point.padding = 0.1)
  } else {cat("The value for the parameter 'type' is incorrect!\nPlease choose one of the following: 'data, plot'.")}
}
# differential protein plot
# defaule parameters: data_input = MatrixQ_log, anno_input = raw_data, category, candiProtein: 'CD276', object: c('A_BL', 'C_PD'), type = 'boxplot|waterfall' 
diff.plot <- function(data_input = MatrixQ_log, anno_input = raw_data, category, candiProtein = 'NULL', object = 'NULL', type = 'NULL') {
  # add condition annotation
  annotation <- anno_input%>%
    dplyr::select(SampleID, condition = category)%>%
    distinct(SampleID, .keep_all = T)%>%
    filter(SampleID != 'PBS')
  Matrix_anno <- as.data.frame(t(data_input))%>%
    tibble::rownames_to_column(var = 'SampleID')%>%
    left_join(annotation, by = 'SampleID')
  Matrix_ID <- Matrix_anno%>%
    dplyr::select(-SampleID)
  if (type == 'boxplot') {
    # specify groups for analysis
    DE4test <- Matrix_ID%>%
      filter(condition %in% object)
    data4box <- DE4test%>%
      dplyr::select(candiProtein, condition)%>%
      rename(value = candiProtein)
    ggplot(data4box, aes(condition, value)) +
      geom_boxplot(aes(color = condition), outlier.color = NA) +
      geom_jitter(aes(color = condition), width = 0.3, height = 0) +
      theme_test() +
      theme(legend.position = 'none') +
      scale_color_brewer(palette = 'Set2') +
      labs(x = NULL,
           y = candiProtein)
  } else if (type == 'waterfall') {
    data4water <- Matrix_anno%>%
      dplyr::select(SampleID, condition, candiProtein)%>%
      rename(value = candiProtein)
    ggplot(data4water, aes(reorder(SampleID, -value), value)) +
      geom_bar(stat = 'identity', aes(fill = condition)) +
      theme_classic() +
      theme(legend.title = element_blank(),
            axis.text.x = element_text(angle = 45, hjust = 1)) +
      scale_fill_brewer(palette = 'Set2') +
      labs(x = NULL,
           y = candiProtein)
  } else {cat("Please input correct paremeters!\nIf 'type' is 'boxplot', please specify 'candiProtein' and 'object'.\nIf 'type' is 'waterfall', please specify 'candiProtein'.")}
}

# model construction and evaluation ---------------------------------------
# all subsets regression model: multivariable
# default parameter: data_input = MatrixQ_log, protein_personalized = c('glypican-3', 'CEA', 'VEGFR2', 'CD133'), category = 'Group|Subgroup1|Subgroup2|NULL'
reg.Model <- function(data_input = MatrixQ_log, protein_personalized = 'NULL', category = 'NULL'){
  if (protein_personalized == 'NULL') {
    MM <- data_input
  }
  if (protein_personalized != 'NULL') {
    MM <- data_input%>%
      tibble::rownames_to_column(var = 'protein')%>%
      filter(protein %in% protein_personalized)%>%
      tibble::column_to_rownames(var = 'protein')
  }
  anno <- sample_list%>%
    dplyr::select(SampleID, condition = category)
  data4model <- data.frame(t(MM))%>%
    tibble::rownames_to_column(var = 'SampleID')%>%
    inner_join(anno, by = 'SampleID')%>%
    dplyr::select(
      condition, everything(),
      -SampleID
    )%>%
    mutate(condition = as.factor(as.character(condition)))
  # time consuming
  leapsReg <<- regsubsets(condition ~ ., data4model, nvmax = 10, really.big = T)
  summary(leapsReg)
}
# leaps plot
# model_input = leapsReg
Mod.Plot <- function(model_input = leapsReg){
  model_input%>%
    plot(scale = "adjr2")
}
# model evaluate (ROC): variables from 2 to 10
# default parameter: data_input = MatrixQ_log, annotation = sample_list, predictor = 10, object = c('A_BL', 'B_PR'), type = 'data|plot|ROCpar', category = 'Group|Subgroup1|Subgroup2|NULL'
predict4ROC <- function(data_input = MatrixQ_log, annotation = sample_list, predictor, category = 'NULL', object, type = 'NULL'){
  # regression parameter
  coefficient <- coef(leapsReg, predictor)
  coef_protein <- names(coefficient)[-1]
  coef_value <- coefficient[-1]
  coef_inter <- coefficient[1]
  # regression calculation
  MM <- data_input
  anno <- annotation%>%
    dplyr::select(SampleID, category)
  data4model <- data.frame(t(MM))%>%
    tibble::rownames_to_column(var = 'SampleID')%>%
    inner_join(anno, by = 'SampleID')%>%
    rename(condition = category)%>%
    arrange(condition)
  data4cal <- data4model%>%
    dplyr::select(coef_protein)%>%
    as.matrix()
  data4ROC <- data.frame(
    condition = data4model$condition,
    value = data4cal%*%coef_value,
    check.names = F)%>%
    dplyr::mutate(value = value + coef_inter)
  # ROC by condition
  ROCbycondition <- data4ROC%>%
    filter(condition %in% object)%>%
    mutate(condition = as.character(condition))
  # ROC model
  ROC_model <- roc(ROCbycondition$condition ~ ROCbycondition$value, percent = T)
  if (type == 'data') {
    return(data4ROC)
  } else if (type == 'plot') {
    # predict by condition (ROC)
    AUC_value <- auc(ROC_model)%>%
      as.numeric()%>%
      round(2)
    ci_value <- ci(ROC_model)%>%
      as.numeric()%>%
      round(2)
    basicplot <- ggplot(ROCbycondition, aes(d = condition, m = value)) +
      geom_roc() +
      theme_classic()
    ROCplot <- basicplot +
      geom_abline(intercept = 0, slope = 1, linetype = 'dashed', color = 'grey50') +
      annotate('text',
               x = 0.65, y = 0.45,
               label = paste0('AUC: ', AUC_value, '%', ' (', ci_value[1], '%-', ci_value[3], '%)'))
    return(ROCplot)
  } else if (type == 'ROCpar') {
    # ROC
    roc_parameter <- coords(ROC_model, "best", ret = "all", transpose = F)
    return(roc_parameter)
  } else if (type == 'box') {
    ggplot(data4ROC, aes(condition, value, color = condition)) +
      geom_boxplot() +
      geom_jitter(height = 0, width = 0.3) +
      theme_minimal() +
      theme(legend.position = 'none') +
      scale_color_brewer(palette = 'Set1') +
      labs(x = NULL,
           y = 'ROC signal')
  } else {cat("The value for the parameter 'type' is incorrect!\nPlease choose one of the following: 'data, plot'.")}
}

# model construction and evaluation (personalized) ------------------------
# model evaluate (ROC): personalized
# default parameter: data_input = MatrixQ_log, annotation = sample_list, protein_candidate = c('glypican-3', 'CEA', 'VEGFR2', 'CD133'), class = 'Group|Subgroup1|Subgroup2|NULL', type = 'plot|data'
AUC.personalized.single <- function(data_input = MatrixQ_log, annotation = sample_list, protein_candidate = 'NULL', category, object, type = 'NULL'){
  anno <- annotation%>%
    dplyr::select(SampleID, category)%>%
    rename(condition = category)
  MM <- data_input
  df4ROC_personalized <- data.frame(t(MM), check.names = F)%>%
    select(one_of(protein_candidate))%>%
    tibble::rownames_to_column(var = 'SampleID')%>%
    inner_join(anno, by = 'SampleID')%>%
    filter(condition %in% object)%>%
    select(-SampleID)
  if (type == 'plot') {
    df4ROC_personalized%>%
      gather(Protein, value, -condition)%>%
      ggplot(aes(d = condition, m = value, color = Protein)) +
      geom_roc(n.cuts = 0) +
      geom_abline(intercept = 0, slope = 1, linetype = 'dashed', color = 'grey50') +
      theme_classic() +
      scale_color_brewer(palette = 'Set1')
  } else if (type == 'data') {
    AUC_pers <- data.frame(Protein = NA, AUC = NA)
    for (i in 1:(ncol(df4ROC_personalized)-1)) {
      ROC <- roc(df4ROC_personalized[, ncol(df4ROC_personalized)] ~ df4ROC_personalized[,i], percent = TRUE)
      AUC_value <- auc(ROC)%>%
        as.numeric()%>%
        round(2)
      AUC_name <- colnames(df4ROC_personalized)[i]
      AUC_pers[i, 1] <- AUC_name
      AUC_pers[i, 2] <- paste0(AUC_value, '%')
      AUC_pers <<- AUC_pers
      # best ROC
      roc_parameter <- coords(ROC, "best", ret = "all", transpose = F)
      print(AUC_name)
      print(roc_parameter)
    }
    return(AUC_pers)
  }
}
# model evaluate (ROC): combind plot
# default parameter: data_input = MatrixQ_log, annotation = sample_list, protein_candidate = c('glypican-3', 'CEA', 'VEGFR2', 'CD133'), category = 'Group|Subgroup1|Subgroup2|NULL', object, type = 'data|plot'
AUC.combind.plot <- function(data_input = MatrixQ_log, annotation = sample_list, protein_candidate = 'NULL', category = 'NULL', object, type = 'NULL'){
  # single protein
  anno <- annotation%>%
    dplyr::select(SampleID, category)%>%
    rename(condition = category)
  MM <- data_input
  df4ROC_personalized <- data.frame(t(MM), check.names = F)%>%
    select(one_of(protein_candidate))%>%
    tibble::rownames_to_column(var = 'SampleID')%>%
    inner_join(anno, by = 'SampleID')
  # multiple protein
  data4model <- data.frame(t(MM))%>%
    tibble::rownames_to_column(var = 'SampleID')%>%
    inner_join(anno, by = 'SampleID')%>%
    arrange(condition)
  data4ROC <- matrix(nrow = nrow(df4ROC_personalized))%>%
    as.data.frame()
  pre_num <- length(protein_candidate)
  for (i in 2:pre_num) {
    coefficient <- coef(leapsReg, i)
    coef_protein <- names(coefficient)[-1]
    coef_value <- coefficient[-1]
    coef_inter <- coefficient[1]
    # regression calculation
    data4cal <- data4model%>%
      dplyr::select(coef_protein)%>%
      as.matrix()
    data4ROC_single <- data.frame(
      value = data4cal%*%coef_value,
      check.names = F)%>%
      dplyr::mutate(value = value + coef_inter)
    names(data4ROC_single) <- paste(coef_protein, collapse = '+')
    data4ROC <- bind_cols(data4ROC, data4ROC_single)
  }
  data4ROC <- select(data4ROC, -V1)
  data4ROC$SampleID <- data4model$SampleID
  # combined plot
  df4ROC_combined <- inner_join(data4ROC, df4ROC_personalized, by = 'SampleID')%>%
    filter(condition %in% object)%>%
    select(-SampleID)
  if (type == 'plot') {
    df4ROC_combined%>%
      gather(Protein, value, -condition)%>%
      ggplot(aes(d = condition, m = value, color = Protein)) +
      geom_roc(n.cuts = 0) +
      geom_abline(intercept = 0, slope = 1, linetype = 'dashed', color = 'grey50') +
      theme_classic() +
      scale_color_brewer(palette = 'Set1')
  } else if (type == 'data') {
    AUC_coms <- data.frame(Protein = NA, AUC = NA)
    for (i in 1:(ncol(df4ROC_combined)-1)) {
      ROC <- roc(df4ROC_combined[, ncol(df4ROC_combined)] ~ df4ROC_combined[,i], percent = TRUE)
      AUC_value <- auc(ROC)%>%
        as.numeric()%>%
        round(2)
      AUC_name <- colnames(df4ROC_combined)[i]
      AUC_coms[i, 1] <- AUC_name
      AUC_coms[i, 2] <- paste0(AUC_value, '%')
      AUC_coms <<- AUC_coms
      # best ROC
      roc_parameter <- coords(ROC, "best", ret = "all", transpose = F)
      print(AUC_name)
      print(roc_parameter)
    }
    return(AUC_coms)
  }
}

# Survival analysis for proteins ------------------------------------------
# defaule parameters: data_input, PFS_info = sample_list, protein_name, plot_type = 'full|mini'
PFS4protein <- function(data_input, PFS_info = sample_list, protein_name, plot_type){
  # subset protein
  protein4pfs <- as.data.frame(t(data_input))%>%
    select(protein_name)%>%
    tibble::rownames_to_column(var = 'SampleID')%>%
    rename(protein_value = protein_name)
  # protein median value
  med <- median(protein4pfs$protein_value)
  # condition based on median value
  protein4condition <- protein4pfs%>%
    mutate(
      protein_condition = case_when(
        protein_value <= med ~ 1,
        protein_value > med ~ 2
      )
    )
  # join pfs information with protein value
  list4pfs <- PFS_info%>%
    inner_join(protein4condition, by = 'SampleID')%>%
    mutate(
      PFS = as.numeric(PFS),
      status = 1
    )
  # fit pfs module
  fit <- survminer::surv_fit(Surv(PFS, status) ~ protein_condition, data = list4pfs)
  # plot pfs module
  if (plot_type == 'full') {
    gg4pfs <- ggsurvplot(fit,
                         pval = TRUE, conf.int = TRUE,
                         risk.table = TRUE, # Add risk table
                         risk.table.col = "strata", # Change risk table color by groups
                         linetype = "strata", # Change line type by groups
                         surv.median.line = "hv", # Specify median survival
                         ggtheme = theme_bw(), # Change ggplot2 theme
                         palette = c("#E7B800", "#2E9FDF")
    )
  } else if (plot_type == 'mini') {
    gg4pfs <- ggsurvplot(
      fit,
      pval = TRUE, conf.int = TRUE,
      risk.table.col = "strata", # Change risk table color by groups
      linetype = "strata", # Change line type by groups
      surv.median.line = "hv", # Specify median survival
      ggtheme = theme_bw(), # Change ggplot2 theme
      palette = c("#E7B800", "#2E9FDF")
    )
  }
  # pfs results
  surv_diff <- survdiff(Surv(PFS, status) ~ protein_condition, data = list4pfs)
  return(gg4pfs)
}
