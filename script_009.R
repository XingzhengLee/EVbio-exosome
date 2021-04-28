
# function file directory
script_directory <- 'D:/Leroy/GitHub/2020/EVbio-exosome/'
# working directory
working_directory <- 'D:/Leroy/Pea/EVbio/data4analysis/009'
setwd(working_directory)

# source script
source(paste0(script_directory, 'function.R'))

# input -------------------------------------------------------------------
# defaule parameters: quantification_file = 'Quantification.txt', protein_file = 'Protein.list.txt', sample_file = 'Sample.list.txt', fluorescence_label = '635'
file.in(quantification_pattern = 'L', protein_file = 'Protein.list.txt', sample_file = 'Sample.list.txt', fluorescence_label = '635')

# CV ----------------------------------------------------------------------
# CV calculation
# default parameters: data_input = raw_data
CV.cal(data_input = raw_data)
write_csv(CV, '01.coefficient_of_variation.csv')

# CV violin
# default parameters: data_input = raw_data, type = 'CV_F|CV_B|CV_FT'
CV.violin(type = 'CV_F')
CV.violin(type = 'CV_B')
CV.violin(type = 'CV_FT')

# CV top
# default parameter: data_input = CV, topN = 20
CV.top(data_input = CV, topN = 20)

# quality control ---------------------------------------------------------
# data for analysis (several options)
# defaule parameters: data_input = raw_data, option = '1|2|3'
data4analysis(data_input = raw_data, option = 1)

# negative control
# default parameters: data_input = Dvalue
negative.control(data_input = Dvalue)
write.csv(MatrixQ, '02.MatrixQ_raw.csv')

# batch effect adjustment --------------------------------------------------
# defaule parameters: proein_expression = MatrixQ, phenotype = sample_list, type = 'box|density|data'
batch_adjust(proein_expression = MatrixQ, phenotype = sample_list, type = 'box')
batch_adjust(proein_expression = MatrixQ, phenotype = sample_list, type = 'density')
batch_adjust(proein_expression = MatrixQ, phenotype = sample_list, type = 'data')
write.csv(MatrixQ_adjusted, '03.MatrixQ_batch_effect_adjusted.csv')

# error control -----------------------------------------------------------
# the percentage of positive values
# default parameter: data_input = MatrixQ_adjusted
positive.percentage(data_input = MatrixQ_adjusted)

# matrix translation and log transformation
# default parameter: data_input = MatrixQ_adjusted, drop_percentage
translation2log(data_input = MatrixQ_adjusted, drop_percentage = 0.1)
write.csv(MatrixQ_log, '04.MatrixQ_processed.csv')

# heatmap -----------------------------------------------------------------
# heatmap_global
# default parameters:data_input = MatrixQ_log, type = 'global|sample', sample_cluster = T/F
Matrix2heat(data_input = MatrixQ_log, type = 'global', sample_cluster = T, protein_cluster = T)
Matrix2heat(data_input = MatrixQ_log, type = 'global', sample_cluster = F, protein_cluster = F)

# sample correlation heatmap
# default parameters: data_input = MatrixQ_log, type = 'global|sample'
Matrix2heat(data_input = MatrixQ_log, type = 'sample')

# k-means -----------------------------------------------------------------
# protein scale
# default parameters: data_input = MatrixQ_log
scale4prot(data_input = MatrixQ_log)

# kmeans 1:10
# default parameters: data_input = data_scale
kmeanstop(data_input = data_scale)

k_means_value <- 5
# kmeans heatmap
# defaule parameters: data_intput = data_scale, kmeans_value
kmeansHeat(data_intput = data_scale, kmeans_value = k_means_value)

# kmeans cluster
# data_input = data_scale, kmeans_value, type = 'cluster|dend'
kmeansCluster(data_input = data_scale, kmeans_value = k_means_value, type = 'cluster')
kmeansCluster(data_input = data_scale, kmeans_value = k_means_value, type = 'dend')

# principal components analysis (PCA) -------------------------------------
# PCA_normal+circle
# defaule parameters: data_input = MatrixQ_log, type = 'Group|Subgroup1|Subgroup2'
PCAplot(data_input = MatrixQ_log, type = 'Group')
PCAplot(data_input = MatrixQ_log, type = 'subgroup1')
PCAplot(data_input = MatrixQ_log, type = 'subgroup2')
PCAplot(data_input = MatrixQ_log, type = 'subgroup3')
PCAplot(data_input = MatrixQ_log, type = 'subgroup4')

# PCA_deep
# defaule parameters: data_input = MatrixQ_log, type = 'data|eigenvalue|variables|individuals|correlation|both', group = 'Group|Subgroup1|Subgroup2'
PCAdeep(data_input = MatrixQ_log, type = 'data')
PCAdeep(data_input = MatrixQ_log, type = 'eigenvalue')
PCAdeep(data_input = MatrixQ_log, type = 'variables')
PCAdeep(data_input = MatrixQ_log, type = 'correlation')
PCAdeep(data_input = MatrixQ_log, type = 'individuals', group = 'Group')
PCAdeep(data_input = MatrixQ_log, type = 'individuals', group = 'subgroup1')
PCAdeep(data_input = MatrixQ_log, type = 'individuals', group = 'subgroup2')
PCAdeep(data_input = MatrixQ_log, type = 'individuals', group = 'subgroup3')
PCAdeep(data_input = MatrixQ_log, type = 'individuals', group = 'subgroup4')

# sub group ---------------------------------------------------------------
# defaule parameters: data_input = MatrixQ_log, annotation = sample_list, category_main = 'Group', object_main = c('SD', 'BL')
# MatrixQ_sub
Matrix.sub(data_input = MatrixQ_log, annotation = sample_list, category_main = 'subgroup4', object_main = c('CB2NCB'))
Matrix.sub(data_input = MatrixQ_log, annotation = sample_list, category_main = 'Group', object_main = c('SD', 'BL'))

# differences between groups (single categroy) ----------------------------
# difference test
# defaule parameters: data_input = MatrixQ_log, anno_input = raw_data, category = 'Group', object = c('A_BL', 'C_PD'), p_threshold = 0.05, fc_threshold = 0, test = 'wilcox.test|t.test', type = 'data|plot'
diff.test(data_input = MatrixQ_sub, anno_input = raw_data, category = 'Group', object = c('SD', 'BL'), p_threshold = 0.05, fc_threshold = 0, test = t.test, type = 'data')
diff.test(data_input = MatrixQ_sub, anno_input = raw_data, category = 'Group', object = c('SD', 'BL'), p_threshold = 0.05, fc_threshold = 0, test = t.test, type = 'plot')
write_csv(DEP, '05.differential_expressed_proteins.csv')

# differential protein plot
# defaule parameters: candiProtein: 'CD276', object: c('A_BL', 'C_PD'), type = 'boxplot|waterfall'
diff.plot(data_input = MatrixQ_sub, anno_input = raw_data, category = 'Group', candiProtein = 'CD73', object = c('SD', 'BL'), type = 'boxplot')
diff.plot(data_input = MatrixQ_sub, anno_input = raw_data, category = 'Group', candiProtein = 'CD73', type = 'waterfall')

# model construction and evaluation ---------------------------------------
# all subsets regression model: multivariable
# defaule parameters: data_input = MatrixQ_log, protein_personalized = c('glypican-3', 'CEA', 'VEGFR2', 'CD133')/'NULL', category = 'NULL'
reg.Model(data_input = MatrixQ_sub, protein_personalized = 'NULL', category = 'Group')

# leaps plot
# defaule parameters: model_input = leapsReg
Mod.Plot(model_input = leapsReg)

# model evaluate (ROC): variables from 2 to 10
# defaule parameters: data_input = MatrixQ_log, annotation = sample_list, predictor = 10, object = c('A_BL', 'B_PR'), type = 'data|plot|ROCpar|box', category = 'Group|Subgroup1|Subgroup2|NULL'
pre_num <- 4
predict4ROC(data_input = MatrixQ_log, annotation = sample_list, predictor = pre_num, category = 'Group', object = c('SD', 'BL'), type = 'data')
predict4ROC(data_input = MatrixQ_log, annotation = sample_list, predictor = pre_num, category = 'Group', object = c('SD', 'BL'), type = 'plot')
predict4ROC(data_input = MatrixQ_log, annotation = sample_list, predictor = pre_num, category = 'Group', object = c('SD', 'BL'), type = 'ROCpar')
predict4ROC(data_input = MatrixQ_log, annotation = sample_list, predictor = pre_num, category = 'Group', object = c('SD', 'BL'), type = 'box')

# model construction and evaluation (personalized) ------------------------
# all subsets regression model: multivariable
# defaule parameters: data_input = MatrixQ_log, protein_personalized = c('glypican-3', 'CEA', 'VEGFR2', 'CD133')/'NULL'
reg.Model(data_input = MatrixQ_log, protein_personalized = c('CD73', 'Tspan8', 'TSG101'), category = 'Group')

# leaps plot
# defaule parameters: model_input = leapsReg
Mod.Plot(model_input = leapsReg)

# model evaluate (ROC): personalized (single)
# defaule parameters: data_input = MatrixQ_log, annotation = sample_list, protein_candidate = c('glypican-3', 'CEA', 'VEGFR2', 'CD133'), category = 'Group|Subgroup1|Subgroup2|NULL', object = c('SD', 'BL'), type = 'plot|data'
AUC.personalized.single(data_input = MatrixQ_log, annotation = sample_list, protein_candidate = c('CD73', 'Tspan8', 'TSG101'), category = 'Group', object = c('SD', 'BL'), type = 'data')
AUC.personalized.single(data_input = MatrixQ_log, annotation = sample_list, protein_candidate = c('CD73', 'Tspan8', 'TSG101'), category = 'Group', object = c('SD', 'BL'), type = 'plot')

# model evaluate (ROC): combind plot
# defaule parameters: data_input = MatrixQ_log, annotation = sample_list, protein_candidate = c('glypican-3', 'CEA', 'VEGFR2', 'CD133'), category = 'Group|Subgroup1|Subgroup2|NULL', object = c('SD', 'BL'), type = 'plot|data
AUC.combind.plot(data_input = MatrixQ_log, annotation = sample_list, protein_candidate = c('CD73', 'Tspan8', 'TSG101'), category = 'Group', object = c('SD', 'BL'), type = 'data')
AUC.combind.plot(data_input = MatrixQ_log, annotation = sample_list, protein_candidate = c('CD73', 'Tspan8', 'TSG101'), category = 'Group', object = c('SD', 'BL'), type = 'plot')

# Survival analysis for proteins ------------------------------------------
# defaule parameters: data_input = MatrixQ_log, PFS_info = sample_list, protein_name
# run ROC first
AUC.personalized.single(data_input = MatrixQ_sub, annotation = sample_list, protein_candidate = c('Tspan8'), category = 'Group', object = c('SD', 'BL'), type = 'data')
# pfs analysis
# defaule parameters: data_input, PFS_info = sample_list, protein_name, method = 'median|ROC', plot_type = 'full|mini'
PFS4protein(data_input = MatrixQ_sub, PFS_info = sample_list, protein_name = 'Tspan8', method = 'ROC', plot_type = 'mini')
PFS4protein(data_input = MatrixQ_sub, PFS_info = sample_list, protein_name = 'Tspan8', method = 'ROC', plot_type = 'full')
p_value
hazard_ratio

# Dynamic curve -----------------------------------------------------------
# data_input = MatrixQ_sub, protein, xlab, ylab, object
dynamic_curve(data_input = MatrixQ_log, protein = 'IFN_gamma', xlab = 'Group', ylab = 'subgroup1', object = c('BL', 'PD', 'PR', 'SD'))
