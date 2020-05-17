
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
      dplyr::select(Batch, Block, Column, Row, SampleID, ID, Group, Protein, fluorescence_F, fluorescence_B, fluorescence_FT, everything())%>%
      mutate(Block = as.factor(Block))
    return(raw_data)
  }
}

# CV ----------------------------------------------------------------------
# CV calculation
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
    col_violin <- rep(col_violin[1:3], sample_num)
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
      select(Batch, Block, value = type)%>%
      ggplot(aes(Block, value, group = Block)) +
      labs(y = type) +
      plot_suffix
  } else {cat("The value for the parameter 'type' is incorrect!\nPlease choose one of the following: 'CV_F, CV_B, CV_FT'.")}
}
# CV top
# default parameter: data_input = CV, topN = 20
CV.top <- function(data_input = CV, topN = 20){
  CVadd <- unite(
    CV,
    # unite Sample (block) and Protein name
    Batch_Block_Protein, c('Batch', 'Block', 'Protein'), sep = '_', remove = F)
  # CVF532 top20
  TopCVF532 <- arrange(CVadd, -CVF532)%>%
    dplyr::slice(1:topN)
  TopCVF532 <- dplyr::mutate(TopCVF532, Batch_Block_Protein = factor(Batch_Block_Protein, levels = TopCVF532$Batch_Block_Protein))
  # CVB532 top20
  TopCVB532 <- arrange(CVadd, -CVB532)%>%
    dplyr::slice(1:topN)
  TopCVB532 <- dplyr::mutate(TopCVB532, Batch_Block_Protein = factor(Batch_Block_Protein, levels = TopCVB532$Batch_Block_Protein))
  # CVF532T top20
  TopCVF532T <- arrange(CVadd, -CVF532T)%>%
    dplyr::slice(1:topN)
  TopCVF532T <- dplyr::mutate(TopCVF532T, Batch_Block_Protein = factor(Batch_Block_Protein, levels = TopCVF532T$Batch_Block_Protein))
  # plot
  bar_col <- brewer.pal(8, 'Set2')
  CV_bar_1 <- ggplot(
    TopCVF532,
    aes(Batch_Block_Protein, CVF532, fill = Batch_Block_Protein)) +
    geom_bar(stat = 'identity', color = 'black') +
    theme_minimal() +
    theme(legend.position = 'none') +
    scale_fill_manual(values = rep(bar_col[1:3], nrow(TopCVF532))) +
    coord_flip()
  CV_bar_2 <- ggplot(
    TopCVB532,
    aes(Batch_Block_Protein, CVB532, fill = Batch_Block_Protein)) +
    geom_bar(stat = 'identity', color = 'black') +
    theme_minimal() +
    theme(legend.position = 'none') +
    scale_fill_manual(values = rep(bar_col[1:3], nrow(TopCVB532))) +
    coord_flip()
  CV_bar_3 <- ggplot(
    TopCVF532T,
    aes(Batch_Block_Protein, CVF532T, fill = Batch_Block_Protein)) +
    geom_bar(stat = 'identity', color = 'black') +
    theme_minimal() +
    theme(legend.position = 'none') +
    scale_fill_manual(values = rep(bar_col[1:3], nrow(TopCVF532T))) +
    coord_flip()
  plot_grid(CV_bar_1, CV_bar_2, CV_bar_3, nrow = 1)
}

# quality control ---------------------------------------------------------
# data for analysis (several options)
data4analysis <- function(data_input = raw_data, option = 'NULL'){
  # normalization & reshape
  value_repeat <- dplyr::group_by(
    data_input,
    Block, Protein, 
    Batch, Block, ID, Group, Protein)%>%
    # mean value in two replicates
    dplyr::summarise(
      MF532 = mean(F532),
      MB532 = mean(B532),
      MF532T = mean(F532T)
    )%>%
    ungroup()
  if (option == 1) {
    Dvalue <<- dplyr::mutate(value_repeat, Value = MF532 - MB532)
  } else if (option == 2) {
    Dvalue <<- dplyr::mutate(value_repeat, Value = MF532)
  } else if (option == 3) {
    Dvalue <<- dplyr::mutate(value_repeat, Value = MF532T)
  } else {cat("The value for the parameter 'type' is incorrect!\nPlease choose one of the following: 1, 2, 3.")}
  return(Dvalue)
}
# negative control
# data_input = Dvalue, type = 'all|matrix'
negcontrol <- function(data_input = Dvalue, type = 'NULL'){
  # positive block
  block_positive <- dplyr::filter(data_input, ID != 'PBS')
  # negative panel (protein correction)
  correction4protein <- dplyr::filter(data_input, ID == 'PBS')%>%
    dplyr::select(Batch, Protein, Value_negative4protein = Value)
  # negative sample in each panel (block correction)
  correction4block <- dplyr::filter(block_positive, Protein == 'PBS')%>%
    dplyr::select(Batch, Block, Value_negative4block = Value)
  # data correction
  normalizated_data <<- inner_join(block_positive, correction4protein, by = c('Batch', 'Protein'))%>%
    inner_join(correction4block, by = c('Batch', 'Block'))%>%
    dplyr::mutate(Norvalue = Value - Value_negative4protein)%>%
    dplyr::filter(!(Protein %in% c('PBS', '+')))
  # normalizated_data <<- inner_join(block_positive, correction4protein, by = c('Batch', 'Protein'))%>%
  #   inner_join(correction4block, by = c('Batch', 'Block'))%>%
  #   dplyr::mutate(Norvalue = Value - Value_negative4protein - Value_negative4block)%>%
  #   dplyr::filter(!(Protein %in% c('PBS', '+')))
  # MartixQ
  MatrixQ <<- dplyr::select(normalizated_data, ID, Protein, Norvalue)%>%
    spread(ID, Norvalue)%>%
    tibble::column_to_rownames(var = 'Protein')
  if (type == 'all') {
    return(normalizated_data)
  } else if (type == 'matrix') {
    return(head(MatrixQ))
  }
}

# batch effect adjustment --------------------------------------------------
batch_adjust <- function(proein_expression = MatrixQ, phenotype = sample_list, type = 'NULL'){
  # load data
  pheno <- phenotype%>%
    filter(ID != 'PBS')%>%
    arrange(ID)
  edata <- as.matrix(proein_expression)
  # adjusting for batch effects with Combat
  batch <- pheno$Batch
  modcombat = model.matrix(~1, data = pheno)
  combat_edata = ComBat(dat = edata, batch = batch, mod = modcombat, par.prior = TRUE, prior.plots = FALSE)
  MatrixQ_ajusted <<- as.data.frame(combat_edata)
  # plot preparation
  anno <- select(pheno, Batch, ID)
  plot_none <- t(edata)%>%
    as.data.frame()%>%
    tibble::rownames_to_column(var = 'ID')%>%
    tbl_df()%>%
    gather(protein, value, -ID)%>%
    left_join(anno, by = 'ID')%>%
    mutate(Batch = as.factor(Batch),
           type = 'Before ajustment')
  plot_ComBat <- t(combat_edata)%>%
    as.data.frame()%>%
    tibble::rownames_to_column(var = 'ID')%>%
    tbl_df()%>%
    gather(protein, value, -ID)%>%
    left_join(anno, by = 'ID')%>%
    mutate(Batch = as.factor(Batch),
           type = 'ComBat ajustment')
  plot_bind <- bind_rows(plot_none, plot_ComBat)
  if (type == 'box') {
    ggplot(plot_bind, aes(Batch, value)) +
    geom_violin(aes(fill = Batch), trim = F) +
    geom_boxplot(width = 0.1) +
  theme_test() +
    theme(legend.position = 'none') +
    facet_wrap(~type) +
    scale_fill_brewer(palette = 'Set2') +
    labs(x = NULL,
         y = 'Protein quantification')
  } else if (type == 'density') {
    ggplot(plot_bind, aes(value)) +
    geom_density(aes(color = Batch), size = 1) +
    theme_test() +
    theme(legend.position = 'none') +
    facet_wrap(~type) +
    scale_color_brewer(palette = 'Set1') +
    labs(x = NULL)
  } else if (type == 'data') {
    return(MatrixQ_ajusted)
  }
}

# error control -----------------------------------------------------------
# initial distribution
# data_input = normalizated_data, xinter = 0
initialdist <- function(data_input = MatrixQ_ajusted, xinter = 0, xinter_min, xinter_max){
  # specify parameters
  xLow <<- xinter_min
  xHigh <<- xinter_max
  MatrixQ_ajusted%>%
    tibble::rownames_to_column(var = 'protein')%>%
    gather(ID, adjvalue, -protein)%>%
    ggplot(aes(adjvalue)) +
    geom_histogram(bins = 1000) +
    geom_vline(xintercept = xinter, color = 'red', linetype = 'dashed', size = 1) +
    geom_vline(xintercept = xinter_min, color = 'blue', linetype = 'dashed', size = 1) +
    geom_vline(xintercept = xinter_max, color = 'blue', linetype = 'dashed', size = 1) +
    theme_classic() +
    labs(x = NULL)
}
# double distribution
# xrange_min = -1000, x_max = 1500, x_inter1 = 0, x_inter2 = 180
doubledist <- function(data_input = MatrixQ_ajusted, xrange_min, x_max, x_inter1, x_inter2){
  # distribution
  MatrixQ_ajusted%>%
    tibble::rownames_to_column(var = 'protein')%>%
    gather(ID, adjvalue, -protein)%>%
    dplyr::filter(adjvalue < x_max)%>%
    ggplot(aes(adjvalue)) +
    geom_histogram(fill = 'white', color = 'black') +
    geom_vline(xintercept = x_inter1, color = 'red', linetype = 'dashed', size = 1) +
    geom_vline(xintercept = x_inter2, color = 'blue', linetype = 'dashed', size = 1) +
    theme_classic() +
    labs(x = NULL) +
    xlim(xrange_min, x_max)
}
# parameter estimation
# data_input = MatrixQ_ajusted, x_min = -1000, x_max = 1500, miu_error = 0, sigma_error = 150, miu_true = 100, sigma_true = 100
parestimation <- function(data_input = MatrixQ_ajusted, x_min, x_max, miu_error, sigma_error, miu_true, sigma_true){
  estimation <<- round(data_input[data_input < x_max & data_input > x_min])
  estimation <<- na.omit(estimation)
  mnf <- function(para, data){
    x <- dnorm(data, para[2], para[3])
    y <- dnorm(data, para[4], para[5])
    res = para[1]*x + (1 - para[1])*y
    l = sum(log(res))
    return(-l)
  }
  ML <- nlminb(c(0.5, miu_error, sigma_error, miu_true, sigma_true),
               mnf,
               data = estimation,
               lower = c(0.01, 0, 0, 0, 0),
               upper = c(0.99, 0, Inf, Inf, Inf))
  p <<- ML$par[1]
  miu1 <<- ML$par[2]
  sigma1 <<- ML$par[3]
  miu2 <<- ML$par[4]
  sigma2 <<- ML$par[5]
  cat("p: ", p, '\nmiu_error: ', miu1, '\nsigma_error: ', sigma1, '\nmiu_true: ', miu2, '\nsigma_true: ', sigma2,
      sep = '')
}
# regression fit
# xrange_min = -1500, xrange_max = 1500
regfit <- function(xrange_min, xrange_max){
  x = seq(xrange_min, xrange_max, length.out = 1000)
  density_estimation <- data.frame(
    x,
    f = p*dnorm(x,miu1,sigma1) + (1-p)*dnorm(x,miu2,sigma2))
  # regression fit
  as.data.frame(estimation)%>%
    ggplot() +
    geom_histogram(aes(x = estimation, y = ..density..), fill = 'white', color = 'black') +
    geom_line(data = density_estimation, aes(x, f)) +
    theme_classic() +
    xlim(xrange_min, xrange_max)
}
# FDR & true value & log
# xrange_max = 1500, fdr_threshold = 0.5, type = 'plot|data|fdr'
fdr2log <- function(xrange_max, fdr_threshold = 0.5, type = 'NULL'){
  x = seq(0, xrange_max, length.out = 1000)
  fdr <- data.frame(
    EX = x,
    FDR = 1 - (1 - p)*dnorm(x, miu2, sigma2)/(p*dnorm(x, miu1, sigma1) + (1 - p)*dnorm(x, miu2, sigma2)),
    check.names = F
  )
  fdrplot <- ggplot(fdr) +
    geom_point(aes(EX, FDR), shape = 21) +
    geom_hline(yintercept = fdr_threshold, linetype = 'dashed', color = 'grey50') +
    theme_classic()
  fdr_order <- filter(fdr, FDR > fdr_threshold)%>%
    tbl_df()%>%
    arrange(-EX)
  threshold <- fdr_order$EX[1]
  Ana <<- MatrixQ_ajusted
  Ana[Ana < threshold | Ana <= 0] <<- 1
  # is.na(Ana)
  Analog2 <<- log2(Ana)
  if (type == 'plot') {
    return(fdrplot)
  } else if (type == 'data') {
    return(Analog2)
  } else if (type == 'fdr') {
    print(paste('FDR threshold:', threshold))
  } else {cat("Please add parameter 'type' to display the data.")}
}

# heatmap -----------------------------------------------------------------
# data_input = Analog2, type = 'global|sample', sample_cluster = T/F
Analog2heat <- function(data_input = Analog2, type = 'NULL', sample_cluster = 'NULL', protein_cluster = 'NULL'){
  data4heat <- data_input
  data4heat[data4heat == 0] <- NA
  if (type == 'global') {
    if (sample_cluster != 'NULL') {
      # heatmap_global
      pheatmap(data4heat, border_color = 'white',
               fontsize = 9, fontsize_row = 10, fontsize_col = 10,
               cluster_rows = protein_cluster, cluster_cols = sample_cluster)
    } else {cat("Please add parameter 'sample_cluster' (T/F).")}
  } else if (type == 'sample') {
    # sample correlation heatmap
    cor_sample <- cor(data_input, use = 'pairwise.complete.obs')
    pheatmap(cor_sample,
             color = rev(inferno(100)), border_color = NA,
             cluster_rows = T, cluster_cols = T)
  } else {cat("Please add parameter 'type' to display the data.")}
}

# k-means -----------------------------------------------------------------
# protein scale
# default parameters: data_input = Analog2
scale4prot <- function(data_input = Analog2){
  data_scale <<- as.data.frame(t(apply(data_input, 1, scale)))
  colnames(data_scale) <<- colnames(data_input)
  return(data_scale)
}
# kmeans 1:10
# default parameters: data_input = data_scale
kmeanstop <- function(data_input = data_scale) {
  wss <- NA
  for (i in 1:10) {
    wss[i] <- sum(kmeans(na.omit(data_input), centers = i)$withinss)
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
# data_intput = data_scale, kmeans_value
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
# data_input = data_scale, kmeans_value, type = 'cluster|dend'
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
# data_input = Analog2, type = 'both|group|IAlias'
PCAplot <- function(data_input = Analog2, type = 'NULL'){
  PCA <- prcomp(na.omit(t(data_input)))
  PC1_importance <- summary(PCA)$importance[2,1]*100
  PC2_importance <- summary(PCA)$importance[2,2]*100
  pcaData <- data.frame(PCA$x[,1:2], check.names = F)%>%
    tibble::rownames_to_column(var = 'ID')%>%
    inner_join(sample_list, by = 'ID')
  if (type == 'Group') {
    ggplot(pcaData, aes(x = PC1, y = PC2, color = Group)) +
      geom_point(size = 5, alpha = 2/3) +
      stat_ellipse(linetype = 'dashed', size = 1) +
      theme_bw() +
      scale_color_brewer(palette = 'Set1')
  } else if (type == 'Subgroup1') {
    ggplot(pcaData, aes(x = PC1, y = PC2, color = Subgroup1)) +
      geom_point(size = 5, alpha = 2/3) +
      stat_ellipse(linetype = 'dashed', size = 1) +
      theme_bw() +
      scale_color_brewer(palette = 'Set1')
  } else if (type == 'Subgroup2') {
    ggplot(pcaData, aes(x = PC1, y = PC2, color = Subgroup2)) +
      geom_point(size = 5, alpha = 2/3) +
      stat_ellipse(linetype = 'dashed', size = 1) +
      theme_bw() +
      scale_color_brewer(palette = 'Set1')
  } else {cat("The value for the parameter 'type' is incorrect!\nPlease choose one of the following: 'both, group, IAlias'.")}
}
# PCA_deep
# data_input = Analog2, type = 'data|eigenvalue|variables|individuals|correlation|both'
PCAdeep <- function(data_input = Analog2, type = 'NULL', group = 'NULL'){
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
  } else if (type == 'individuals') {
    dt <- na.omit(t(data_input))
    dt_anno <- as.data.frame(dt)%>%
      tibble::rownames_to_column(var = 'ID')%>%
      left_join(sample_list, by = 'ID')
    if (group == 'Group') {
      fviz_pca_ind(prcomp(dt),
                   geom.ind = "point",
                   col.ind = dt_anno$Group,
                   palette = brewer.pal(9, "Set1")[1:length(unique(dt_anno$Group))],
                   # whether to use ggrepel to avoid overplotting text labels
                   repel = T,
                   addEllipses = T,
                   legend.title = "Groups")
    } else if (group == 'Subgroup1') {
      fviz_pca_ind(prcomp(dt),
                   geom.ind = "point",
                   col.ind = dt_anno$Subgroup1,
                   palette = brewer.pal(9, "Set1")[1:length(unique(dt_anno$Subgroup1))],
                   # whether to use ggrepel to avoid overplotting text labels
                   repel = T,
                   addEllipses = T,
                   legend.title = "Groups")
    } else if (group == 'Subgroup2') {
      fviz_pca_ind(prcomp(dt),
                   geom.ind = "point",
                   col.ind = dt_anno$Subgroup2,
                   palette = brewer.pal(9, "Set1")[1:length(unique(dt_anno$Subgroup2))],
                   # whether to use ggrepel to avoid overplotting text labels
                   repel = T,
                   addEllipses = T,
                   legend.title = "Groups")
    }
  } else if (type == 'correlation') {
    var <- get_pca_var(PCA)
    corrplot(var$cos2, is.corr = FALSE,
             method = 'shade',
             col = brewer.pal(9, "OrRd"),
             tl.cex = 0.7,
             tl.col = 'grey20',
             cl.pos = 'b',
             cl.cex = 0.5)
  } else {cat("The value for the parameter 'type' is incorrect!\nPlease choose one of the following: 'data, eigenvalue, variables, individuals, correlation'.")}
}

# differences between groups ----------------------------------------------
# wilcoxon test
# data_input = raw_data, object = c('A_BL', 'C_PD'), p_threshold = 0.05, fc_threshold = 0, test = 'nonpar|Ttest', type = 'data|plot'
diffTest <- function(data_input = Analog2, anno_input = raw_data, object, p_threshold = 0.05, fc_threshold = 0, test = 'nonpar', type = 'NULL'){
  # add group annotation
  annotation <- dplyr::select(anno_input, ID, Group)%>%
    distinct(ID, .keep_all = T)%>%
    filter(ID != 'PBS')
  DEgroup <<- t(data_input)%>%
    as.data.frame()%>%
    tibble::rownames_to_column(var = 'ID')%>%
    left_join(annotation, by = 'ID')
  DE <<- dplyr::select(DEgroup, -ID)
  # specify groups for analysis
  DE4wil <- filter(DE, Group %in% object)
  # wilcoxon test
  diff_test <- data.frame(Protein = NA, p.value = NA, log2FC = NA)
  for (i in 1:(ncol(DE4wil)-1)) {
    print(paste('Processing proteins:', i))
    # log2mean
    data4mean <- dplyr::select(DE4wil, i, ncol(DE4wil))
    colnames(data4mean) <- c('Value', 'Group')
    mean <- group_by(data4mean, Group)%>%
      dplyr::summarise(mean = mean(Value, na.rm = T))
    logvalue <- log2(mean$mean[2]/mean$mean[1])
    # creat result dataframe
    if (test == 'nonpar') {
      diff_test[i,1] <- colnames(DE4wil)[i]
      diff_test[i,2] <- wilcox.test(DE4wil[,i] ~ DE4wil[,ncol(DE4wil)], data = DE4wil, na.action = na.omit)$p.value
      diff_test[i,3] <- logvalue
    } else if (test == 'Ttest') {
      diff_test[i,1] <- colnames(DE4wil)[i]
      diff_test[i,2] <- t.test(DE4wil[,i] ~ DE4wil[,ncol(DE4wil)], data = DE4wil, na.action = na.omit)$p.value
      diff_test[i,3] <- logvalue
    }
  }
  DEP <<- arrange(diff_test, p.value)
  if (type == 'data') {
    return(diff_test)
  } else if (type == 'plot') {
    # volcano plot
    FCcol <- brewer.pal(9, 'Set1')
    mutate(
      diff_test,
      logP = -log10(p.value),
      label = case_when(p.value < p_threshold & (log2FC > fc_threshold | log2FC < -fc_threshold) ~ Protein),
      FoldChange = ifelse(p.value < p_threshold & log2FC > fc_threshold, 'Up',
                          ifelse(p.value < p_threshold & log2FC < -fc_threshold, 'Down', 'Not'))
    )%>%
      ggplot(aes(log2FC, logP)) +
      geom_point(aes(color = FoldChange)) +
      geom_vline(xintercept = 0, linetype = 'dashed', color = 'grey50') +
      theme_test() +
      scale_color_manual(values = c(FCcol[1], 'grey', FCcol[2])) +
      ggrepel::geom_text_repel(aes(label = label), box.padding = 0.1, label.padding = 0.1, point.padding = 0.1)
  } else {cat("The value for the parameter 'type' is incorrect!\nPlease choose one of the following: 'data, plot'.")}
}
# differential protein plot
# candiProtein: 'CD276', object: c('A_BL', 'C_PD'), type = 'boxplot|waterfall' 
diffplot <- function (candiProtein = 'NULL', object = 'NULL', type = 'NULL') {
  if (type == 'boxplot') {
    # specify groups for analysis
    DE4wil <- filter(DE, Group %in% object)
    data4box <- dplyr::select(DE4wil, candiProtein, Group)
    colnames(data4box) <- c('Value', 'Group')
    ggplot(data4box, aes(Group, Value)) +
      geom_boxplot(aes(color = Group), outlier.color = NA) +
      geom_jitter(aes(color = Group), width = 0.3, height = 0) +
      theme_test() +
      theme(legend.position = 'none') +
      scale_color_brewer(palette = 'Set1') +
      labs(x = NULL,
           y = candiProtein)
  } else if (type == 'waterfall') {
    data4water <- dplyr::select(DEgroup, ID, Group, candiProtein)
    colnames(data4water) <- c('ID', 'Group', 'Value')
    ggplot(data4water, aes(reorder(ID, -Value), Value)) +
      geom_bar(stat = 'identity', aes(fill = Group)) +
      theme_classic() +
      theme(legend.position = c(0.9, 0.8),
            legend.title = element_blank(),
            axis.text.x = element_text(angle = 45, hjust = 1)) +
      scale_fill_brewer(palette = 'Set2') +
      labs(x = NULL,
           y = candiProtein)
  } else {cat("Please input correct paremeters!\nIf 'type' is 'boxplot', please specify 'candiProtein' and 'object'.\nIf 'type' is 'waterfall', please specify 'candiProtein'.")}
}

# model construction and evaluation ---------------------------------------
# all subsets regression model: multivariable
# data_input = Analog2, protein_personalized = c('glypican-3', 'CEA', 'VEGFR2', 'CD133'), class = 'Group|Subgroup1|Subgroup2|NULL'
regModel <- function(data_input = Analog2, protein_personalized = 'NULL', class = 'NULL'){
  if (protein_personalized == 'NULL') {
    MM <- data_input
  }
  if (protein_personalized != 'NULL') {
    MM <- tibble::rownames_to_column(Analog2, var = 'protein')%>%
      filter(protein %in% protein_personalized)%>%
      tibble::column_to_rownames(var = 'protein')
  }
  anno <- dplyr::select(sample_list, ID, class)
  data4model <- data.frame(t(MM))%>%
    tibble::rownames_to_column(var = 'ID')%>%
    inner_join(anno, by = 'ID')%>%
    dplyr::select(class = last_col(), everything(),
                  -ID)%>%
    mutate(class = as.factor(as.character(class)))
  # time consuming
  leapsReg <<- regsubsets(class ~ ., data4model, nvmax = 10)
  summary(leapsReg)
}
# leaps plot
# model_input = leapsReg
ModPlot <- function(model_input = leapsReg){
  plot(model_input, scale = "adjr2")
}
# model evaluate (ROC): variables from 2 to 10
# datainput = Analog2, annotation = sample_list, predictor = 10, object = c('A_BL', 'B_PR'), type = 'data|plot|ROCpar', class = 'Group|Subgroup1|Subgroup2|NULL'
pred4ROC <- function(datainput = Analog2, annotation = sample_list, predictor, object, type = 'NULL', class = 'NULL'){
  # regression parameter
  coefficient <- coef(leapsReg, predictor)
  coef_protein <- names(coefficient)[-1]
  coef_value <- coefficient[-1]
  coef_inter <- coefficient[1]
  # regression calculation
  MM <- datainput
  anno <- dplyr::select(annotation, ID, class)
  data4model <- data.frame(t(MM))%>%
    tibble::rownames_to_column(var = 'ID')%>%
    inner_join(anno, by = 'ID')%>%
    select(class = last_col(), everything())%>%
    arrange(class)
  data4cal <- dplyr::select(data4model, coef_protein)%>%
    as.matrix()
  data4ROC <- data.frame(
    class = data4model$class,
    Value = data4cal%*%coef_value,
    check.names = F)%>%
    dplyr::mutate(Value = Value + coef_inter)
  if (type == 'data') {
    return(data4ROC)
  } else if (type == 'plot') {
    ROCbyGroup <- filter(data4ROC, class %in% object)%>%
      mutate(class = as.character(class))
    Roc <- roc(ROCbyGroup$class ~ ROCbyGroup$Value, percent = T)
    # predict by group (ROC)
    aucValue <- auc(Roc)%>%
      as.numeric()%>%
      round(2)
    ciValue <- ci(Roc)%>%
      as.numeric()%>%
      round(2)
    basicplot <- ggplot(ROCbyGroup, aes(d = class, m = Value)) +
      geom_roc() +
      theme_classic()
    ROCplot <- basicplot +
      geom_abline(intercept = 0, slope = 1, linetype = 'dashed', color = 'grey50') +
      annotate('text',
               x = 0.65, y = 0.45,
               label = paste0('AUC: ', aucValue, '%', ' (', ciValue[1], '%-', ciValue[3], '%)'))
    return(ROCplot)
  } else if (type == 'ROCpar') {
    ROCbyGroup <- filter(data4ROC, class %in% object)%>%
      mutate(class = as.character(class))
    Roc <- roc(ROCbyGroup$class ~ ROCbyGroup$Value, percent=TRUE)
    # ROC
    roc_parameter <- coords(Roc, "best", ret = "all", transpose = F)
    return(roc_parameter)
  } else if (type == 'box') {
    ggplot(data4ROC, aes(class, Value)) +
      geom_boxplot() +
      geom_jitter(height = 0, width = 0.3) +
      theme_classic() +
      labs(x = NULL,
           y = 'ROC signal')
  } else {cat("The value for the parameter 'type' is incorrect!\nPlease choose one of the following: 'data, plot'.")}
}
# model evaluate (ROC): personalized
# data_input = Analog2, annotation = sample_list, protein_candidate = c('glypican-3', 'CEA', 'VEGFR2', 'CD133'), class = 'Group|Subgroup1|Subgroup2|NULL', type = 'plot|data'
AUC_personalized <- function(data_input = Analog2, annotation = sample_list, object, protein_candidate = 'NULL', class = 'NULL', type = 'NULL'){
  anno <- dplyr::select(annotation, ID, class)%>%
    select(class = last_col(), everything())
  MM <- data_input
  df4ROC_personalized <- data.frame(t(MM), check.names = F)%>%
    select(one_of(protein_candidate))%>%
    tibble::rownames_to_column(var = 'ID')%>%
    inner_join(anno, by = 'ID')%>%
    filter(class %in% object)%>%
    select(-ID)
  if (type == 'plot') {
    gather(df4ROC_personalized, Protein, value, -class)%>%
      ggplot(aes(d = class, m = value, color = Protein)) +
      geom_roc(n.cuts = 0) +
      geom_abline(intercept = 0, slope = 1, linetype = 'dashed', color = 'grey50') +
      theme_classic() +
      scale_color_brewer(palette = 'Set1')
  } else if (type == 'data') {
    AUC_pers <- data.frame(Protein = NA, AUC = NA)
    for (i in 1:(ncol(df4ROC_personalized)-1)) {
      Roc <- roc(df4ROC_personalized[, ncol(df4ROC_personalized)] ~ df4ROC_personalized[,i], percent = TRUE)
      aucValue <- auc(Roc)%>%
        as.numeric()%>%
        round(2)
      aucName <- colnames(df4ROC_personalized)[i]
      AUC_pers[i, 1] <- aucName
      AUC_pers[i, 2] <- paste0(aucValue, '%')
      AUC_pers <<- AUC_pers
      # best ROC
      roc_parameter <- coords(Roc, "best", ret = "all", transpose = F)
      print(aucName)
      print(roc_parameter)
    }
    return(AUC_pers)
  }
}
# model evaluate (ROC): combind plot
# data_input = Analog2, annotation = sample_list, protein_candidate = c('glypican-3', 'CEA', 'VEGFR2', 'CD133'), class = 'Group|Subgroup1|Subgroup2|NULL', type = 'data|plot'
AUC_combind_plot <- function(data_input = Analog2, annotation = sample_list, object, protein_candidate = 'NULL', class = 'NULL', type = 'NULL'){
  # single protein
  anno <- dplyr::select(annotation, ID, class)%>%
    select(class = last_col(), everything())
  MM <- data_input
  df4ROC_personalized <- data.frame(t(MM), check.names = F)%>%
    select(one_of(protein_candidate))%>%
    tibble::rownames_to_column(var = 'ID')%>%
    inner_join(anno, by = 'ID')
  # multiple protein
  data4model <- data.frame(t(MM))%>%
    tibble::rownames_to_column(var = 'ID')%>%
    inner_join(anno, by = 'ID')%>%
    select(class = last_col(), everything())%>%
    arrange(class)
  data4ROC <- matrix(nrow = nrow(df4ROC_personalized))%>%
    as.data.frame()
  pre_num <- length(protein_candidate)
  for (i in 2:pre_num) {
    coefficient <- coef(leapsReg, i)
    coef_protein <- names(coefficient)[-1]
    coef_value <- coefficient[-1]
    coef_inter <- coefficient[1]
    # regression calculation
    data4cal <- dplyr::select(data4model, coef_protein)%>%
      as.matrix()
    data4ROC_single <- data.frame(
      Value = data4cal%*%coef_value,
      check.names = F)%>%
      dplyr::mutate(Value = Value + coef_inter)
    names(data4ROC_single) <- paste(coef_protein, collapse = '+')
    data4ROC <- bind_cols(data4ROC, data4ROC_single)
  }
  data4ROC <- select(data4ROC, -V1)
  data4ROC$ID <- data4model$ID
  # combined plot
  df4ROC_combined <- inner_join(data4ROC, df4ROC_personalized, by = 'ID')%>%
    filter(class %in% object)%>%
    select(-ID)
  if (type == 'plot') {
    gather(df4ROC_combined, Protein, value, -class)%>%
      ggplot(aes(d = class, m = value, color = Protein)) +
      geom_roc(n.cuts = 0) +
      geom_abline(intercept = 0, slope = 1, linetype = 'dashed', color = 'grey50') +
      theme_classic() +
      scale_color_brewer(palette = 'Set1')
  } else if (type == 'data') {
    AUC_coms <- data.frame(Protein = NA, AUC = NA)
    for (i in 1:(ncol(df4ROC_combined)-1)) {
      Roc <- roc(df4ROC_combined[, ncol(df4ROC_combined)] ~ df4ROC_combined[,i], percent = TRUE)
      aucValue <- auc(Roc)%>%
        as.numeric()%>%
        round(2)
      aucName <- colnames(df4ROC_combined)[i]
      AUC_coms[i, 1] <- aucName
      AUC_coms[i, 2] <- paste0(aucValue, '%')
      AUC_coms <<- AUC_coms
      # best ROC
      roc_parameter <- coords(Roc, "best", ret = "all", transpose = F)
      print(aucName)
      print(roc_parameter)
    }
    return(AUC_coms)
  }
}

# multiple regression -----------------------------------------------------
# multiple regression analysis
# group = c('Group', 'Subgroup1', 'Subgroup2')
multireg_analysis <- function(data_input = Analog2, annotation = sample_list, group = 'NULL'){
  group_name <- group
  anno <- annotation%>%
    select(ID, group_name)
  df_anno <- data_input%>%
    t()%>%
    as.data.frame()%>%
    tibble::rownames_to_column(var = 'ID')%>%
    tbl_df()%>%
    left_join(anno, by = 'ID')%>%
    select(group_name, everything(),
           -ID)
  group_number <- length(group_name)
  mulregressoin <<- matrix(ncol = 1 + group_number, nrow = 0)%>%
    as.data.frame()
  colnames(mulregressoin) <<- c('protein', group_name)
  for (i in (group_number+1):length(df_anno)) {
    df4avo <- select(df_anno, 1:3, i)
    protein4test <- colnames(df_anno)[i]
    fit <- aov(as.formula(paste(protein4test, '.', sep = " ~ ")), data = df4avo)
    p.table <- summary(fit)[[1]]%>%
      as.data.frame()%>%
      tibble::rownames_to_column(var = 'group')%>%
      mutate(group = str_trim(group))
    rowID <- i-group_number
    mulregressoin[rowID,1] <<- protein4test
    mulregressoin[rowID,2] <<- p.table%>%
      filter(group == group_name[1])%>%
      .[6]
    mulregressoin[rowID,3] <<- p.table%>%
      filter(group == group_name[2])%>%
      .[6]
    mulregressoin[rowID,4] <<- p.table%>%
      filter(group == group_name[3])%>%
      .[6]
  }
  return(mulregressoin)
}
# multiple regression pvalue plot
mulregressoin_pvalue <- function(data_input = mulregressoin, pvalue = 0.05){
  data_input%>%
    gather(group, value, -protein)%>%
    mutate(ln = -log(value),
           label = case_when(value < pvalue ~ protein))%>%
    ggplot(aes(group, ln, color = group)) +
    geom_jitter(height = 0) +
    geom_hline(yintercept = -log(pvalue), color = 'grey50', linetype = 'dashed') +
    theme_test() +
    theme(legend.position = 'none') +
    scale_color_brewer(palette = 'Set1') +
    labs(x = NULL,
         y = '-ln(protein)')
}
