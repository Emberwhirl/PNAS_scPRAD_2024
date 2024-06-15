library(dplyr)
library(Seurat)
library(SCENIC)
library(AUCell)
library(circlize)
library(ComplexHeatmap)
library(stringr)
library(ggplot2)
library(msigdbr)
library(RColorBrewer)
library('stats')
library(tidyHeatmap)
library(colorspace)
library(randomcoloR)
library(colorspace)
library(scCustomize)
library(dendextend)
library(cowplot)
library(fgsea)
library(Rphenograph)
library(clusterProfiler)
library(dendextend)
set.seed(1234)

cluster_between_groups_ward = function(mat, factor) {
  
  if (!is.factor(factor)) {
    factor = factor(factor, levels = unique(factor))
  }
  
  dend_list = list()
  order_list = list()
  for(le in unique(levels(factor))) {
    m = mat[, factor == le, drop = FALSE]
    if (ncol(m) == 1) {
      order_list[[le]] = which(factor == le)
      dend_list[[le]] = structure(which(factor == le), class = "dendrogram", leaf = TRUE,
                                  height = 0, label = 1, members = 1)
    } else if(ncol(m) > 1) {
      hc1 = hclust(dist(1:ncol(m)),method = "ward.D2")
      dend_list[[le]] = reorder(as.dendrogram(hc1), wts = 1:ncol(m), agglo.FUN = mean)
      order_list[[le]] = which(factor == le)[order.dendrogram(dend_list[[le]])]
      order.dendrogram(dend_list[[le]]) = order_list[[le]]
    }
    attr(dend_list[[le]], ".class_label") = le
  }
  
  parent = as.dendrogram(hclust(dist(t(sapply(order_list, function(x) rowMeans(mat[, x, drop = FALSE]))))))
  dend_list = lapply(dend_list, function(dend) dendrapply(dend, function(node) {
    attr(node, "height") = 0
    node
  }))
  dend = merge_dendrogram(parent, dend_list)
  order.dendrogram(dend) = unlist(order_list[order.dendrogram(parent)])
  return(dend)
}

# SCENIC run ####
print('Make input')
seura = readRDS('msk.integrated.remove.cellcycle.tumor.cells') 
seura = subset(seura, subset = subtype == 'CSPC', invert = T )
exprMat <- as.matrix(seura@assays$RNA@counts)
saveRDS(exprMat, "integrated_exprMat.rds")
seura$subtype = factor(seura$subtype, levels = c('CRPC','NEPC'))
Idents(seura) = seura$subtype
cellInfo <- data.frame(seuratCluster=Idents(seura))
rm(seura);gc()
colnames(cellInfo) <- "CellType"
saveRDS(cellInfo, "integrated_cellInfo.rds")
 
print('Initialize SCENIC')
exprMat = readRDS('integrated_exprMat.rds')
cellInfo = readRDS('integrated_cellInfo.rds')

org <- "hgnc"
dbDir = "cisTarget_databases"
myDatasetTitle <- "samir"
dbs <- c("hg38__refseq-r80__10kb_up_and_down_tss.mc9nr.feather", 'hg38__refseq-r80__500bp_up_and_100bp_down_tss.mc9nr.feather')
scenicOptions <- initializeScenic(org=org, dbDir=dbDir, dbs=dbs, datasetTitle=myDatasetTitle, nCores=10)
scenicOptions@inputDatasetInfo$cellInfo = cellInfo
saveRDS(scenicOptions, file="int/scenicOptions.Rds")
gc()
rm(dbDir)
rm(dbs)
rm(org)
 
print('Gene filtering')
genesKept <- geneFiltering(exprMat, scenicOptions)
exprMat_filtered <- exprMat[genesKept, ]
saveRDS(exprMat_filtered,'integrated_exprMat_filtered.rds')
gc()
rm(genesKept)
 
print('Run correlation')
exprMat_filtered = readRDS('integrated_exprMat_filtered_v2.rds')
runCorrelation(exprMat_filtered, scenicOptions)
gc()


print('Start Genie3')
exprMat_filtered_log <- log2(exprMat_filtered+1); rm(exprMat_filtered)
saveRDS(exprMat_filtered_log, 'integrated_exprMat_filtered_log.rds')
gc()

exprMat_filtered_log = readRDS('integrated_exprMat_filtered_log.rds')
scenicOptions = readRDS('int/scenicOptions.Rds')
gc()
runGenie3(exprMat_filtered_log, scenicOptions, resumePreviousRun = T)
rm(exprMat_filtered_log)
gc()

print('Almost done')
exprMat = readRDS('integrated_exprMat.rds')
exprMat_log <- log2(exprMat+1)
print('runSCENIC_1')
runSCENIC_1_coexNetwork2modules(scenicOptions)
saveRDS(scenicOptions,'int/scenicOptions.Rds')

print('runSCENIC_2')
scenicOptions = readRDS('int/scenicOptions.Rds')
runSCENIC_2_createRegulons(scenicOptions)
saveRDS(scenicOptions,'int/scenicOptions.Rds')

print('runSCENIC_3')
scenicOptions = readRDS('int/scenicOptions.Rds')
runSCENIC_3_scoreCells(scenicOptions, exprMat_log)
saveRDS(scenicOptions,'int/scenicOptions.Rds')


# Regulon activity matrix ####
cancer = readRDS('/Users/jooyoung/OneDrive - 고려대학교/samir/rds_jy/msk.integrated.remove.cellcycle.tumor.fastmnn.subtype.nohepato.rds')
cancer = subset(cancer, subtype == 'CSPC',invert =T)
scenicOptions = readRDS('int/scenicOptions.Rds')
regulonAUC <- loadInt(scenicOptions, "aucell_regulonAUC")
regulonAUC <- regulonAUC[onlyNonDuplicatedExtended(rownames(regulonAUC)),]
regulonAUC = regulonAUC %>% getAUC 
regulonAUC = regulonAUC[, colnames(cancer)]
regulonActivity_byCell_scaled <- scale(t(regulonAUC), center = T, scale=T) %>% t()
saveRDS(regulonActivity_byCell_scaled,'regulonActivity_byCell_scaled_raw.rds') # 300

ann = cancer@meta.data
ann$Cluster = cancer$pheno_cluster %>% as.numeric()
ann = arrange(ann ,subtype)
colnames(ann)
ann = ann[,c('subtype','Cluster','patient')]
saveRDS(ann,'ann_regulon.rds')


# Create dendrogram ####
regulonActivity_byCell_scaled = readRDS('regulonActivity_byCell_scaled_raw.rds')
regulonActivity_byCell_scaled = regulonActivity_byCell_scaled[,rownames(ann)]

## Cluster between phenocluster within each subtype ####
### CRPC ####
mat_crpc =regulonActivity_byCell_scaled[,rownames(ann)[ann$subtype == 'CRPC']]
ann_crpc = ann[ann$subtype == 'CRPC',]
mat_crpc = mat_crpc[,rownames(ann_crpc)]
# dend3 = cluster_between_groups_ward(mat_crpc,ann_crpc$Cluster)
dend3 = cluster_between_groups_ward(mat_crpc,ann_crpc$Cluster)
order_crpc = mat_crpc[,order.dendrogram(dend3)] %>% colnames()
saveRDS(order_crpc, 'pheno_cluster_CRPC_order_raw.rds')


### NEPC ####
mat_NEPC =regulonActivity_byCell_scaled[,rownames(ann)[ann$subtype == 'NEPC']]
ann_NEPC = ann[ann$subtype == 'NEPC',]
mat_NEPC = mat_NEPC[,rownames(ann_NEPC)]
dend4 = cluster_between_groups_ward(mat_NEPC,ann_NEPC$Cluster)
order_NEPC = mat_NEPC[,order.dendrogram(dend4)] %>% colnames()
saveRDS(order_NEPC, 'pheno_cluster_NEPC_order_raw.rds')

### Merge dendrogram ####
dend_p = as.dendrogram(hclust(dist(rbind(colMeans(t(mat_crpc)), colMeans(t(mat_NEPC))))))
dend_m = merge_dendrogram(dend_p, list(dend3, dend4))
order_col = c(colnames(mat_crpc),colnames(mat_NEPC))
saveRDS(dend_m,'dend_m_pheno_raw.rds')
saveRDS(order_col,'dend_m_pheno_ordercol_raw.rds')


## Make GRN with ARI Hclust ####
ann = readRDS('ann_regulon.rds')
regulonActivity_byCell_scaled = readRDS('regulonActivity_byCell_scaled_raw.rds')
regulonActivity_byCell_scaled = regulonActivity_byCell_scaled[,rownames(ann)]

#### CRPC ####
##### Hclust ####
mat_crpc =regulonActivity_byCell_scaled[,rownames(ann)[ann$subtype == 'CRPC']]
ann_crpc = ann[ann$subtype == 'CRPC',]
mat_crpc = mat_crpc[,rownames(ann_crpc)]
dend3 = cluster_between_groups_ward(mat_crpc,ann_crpc$Cluster)
dend3 = cluster_between_groups(mat_crpc,ann_crpc$Cluster)
hc1 = as.hclust(dend3)
saveRDS(dend3, '20231113_CRPC_dend.rds')

##### ARI ####
dend3 = readRDS('20231113_CRPC_dend.rds')
hc1 = as.hclust(dend3)
clusterings = list()
n = 1
for (k in seq(0,30,1)){
  print(k)
  cut_avg <- cutree(hc1, h = k) 
  clusterings[[n]] = cut_avg
  n = n+1
}
names(clusterings) = paste0('h_',seq(0,30,1))
combination = combn(names(clusterings) ,2) 

rand_score = list()
for (n in 1: dim(combination)[2]){
  print(n)
  comp1 = combination[1,n]
  comp2 = combination[2,n]
  score = genieclust::adjusted_rand_score(clusterings[[comp1]], clusterings[[comp2]])
  rand_score[[paste0(comp1,'-',comp2)]] = score
}
rand_score

rand_matrix = matrix(nrow = length(names(clusterings)) ,ncol = length(names(clusterings)) ) %>% as.data.frame()
rownames(rand_matrix) = colnames(rand_matrix) = paste0('h_',seq(0,30,1))
rand_matrix

for (i in 1:dim(combination)[2]){
  print(i)
  comp1 = combination[1,i]
  comp2 = combination[2,i]
  score = rand_score[[paste0(comp1,'-',comp2)]]
  rand_matrix[comp1,comp2] = rand_matrix[comp2,comp1] = score
}

for (i in 1:length(colnames(rand_matrix))){
  print(i)
  rand_matrix[i,i] = 1
}

col_fun = colorRamp2(c(0,0.65, 1), c( "blue", "white","red3"))

rand_matrix = read.csv('20231224_ARI_CRPC_step1_score.csv')
rownames(rand_matrix) = rand_matrix$X; rand_matrix$X = NULL
ComplexHeatmap::Heatmap(rand_matrix, cluster_rows = F, cluster_columns = F,name="Adjusted Rand Index", na_col = lighten('lightgray',0.3), col = col_fun)

write.csv(rand_matrix,'20231224_ARI_CRPC_step1_score.csv')

##### CRPC GRN ####
cut_avg <- cutree(hc1, h = 15) # h = 15
cut_avg %>% table
rect.hclust(hc1,  h = 15) # h = 15
crpc_reg = paste0('CRPC_',cut_avg)
crpc_reg = crpc_reg %>% as.data.frame()
crpc_reg$barcode = colnames(mat_crpc)
colnames(crpc_reg)[1] = 'group'
saveRDS(crpc_reg, '20231113_CRPC_regulon_h15.rds')

#### NEPC ####
##### Hclust ####
mat_crpc =regulonActivity_byCell_scaled[,rownames(ann)[ann$subtype == 'NEPC']]
ann_crpc = ann[ann$subtype == 'NEPC',]
mat_crpc = mat_crpc[,rownames(ann_crpc)]
dend3 = cluster_between_groups_ward(mat_crpc,ann_crpc$Cluster)
hc1 = as.hclust(dend4)
saveRDS(dend3, '20231113_NEPC_dend.rds')

##### NEPC GRN ####
cut_avg <- cutree(hc1,h = 15 )
rect.hclust(hc1, h = 15)
crpc_reg = paste0('NEPC_',cut_avg)
crpc_reg = crpc_reg %>% as.data.frame()
crpc_reg$barcode = colnames(mat_crpc)
colnames(crpc_reg)[1] = 'group'
saveRDS(crpc_reg, '20231113_NEPC_regulon_h15.rds')


# Add GRN information ####
crpc_reg = readRDS('20231113_CRPC_regulon_h15.rds')
nepc_reg = readRDS('20231113_NEPC_regulon_h15.rds')
colnames(crpc_reg)[1] = 'group'; colnames(nepc_reg)[1] = 'group'

reg = rbind(crpc_reg, nepc_reg)
rownames(reg) = reg$barcode
reg = reg[rownames(ann),]
which(rownames(reg) != rownames(ann))
ann$group = reg$group
saveRDS(ann,'ann_regulon.rds')

# RSS score for each group ####
scenicOptions = readRDS('int/scenicOptions.Rds')
regulonAUC <- loadInt(scenicOptions, "aucell_regulonAUC")
regulonAUC = regulonAUC[rownames(regulonActivity_byCell_scaled),rownames(meta)]
rss <- calcRSS(AUC=getAUC(regulonAUC), cellAnnotation=meta[colnames(regulonAUC), "Regulon_h15_modi_ano"])
rss
rssNorm <- scale(rss) 
write.csv(rssNorm, '20231129_Human_SCENIC_RSS.csv')

tf_bystsage = as.data.frame(matrix(ncol = 4)); colnames(tf_bystsage) = c('regulon','rank','rss','group')
for (i in meta$Regulon_h15_modi_ano_factor %>% unique){
  print(i)
  tfs = c()
  for (n in 1:dim(rssNorm)[1]){
    if(max(rssNorm[n,]) == rssNorm[n,i]){
      tfs = c(tfs, rownames(rssNorm)[n])
    }
  }
  if (length(tfs) > 0){
    tmp = rssNorm[tfs,]
    rssThisType <- sort(tmp [,i], decreasing=TRUE)
    thisRss <- data.frame(regulon=names(rssThisType), rank=seq_along(rssThisType), rss=rssThisType)
    thisRss$group = i
    tf_bystsage = rbind(tf_bystsage, thisRss)
  }
}

tf_bystsage = tf_bystsage[-1,]
saveRDS(tf_bystsage, '20231225_tf_bystsage_h15.rds')


# T-test ####
cancer = readRDS('../msk.integrated.remove.cellcycle.tumor.fastmnn.subtype.nohepato.rds')
tf_bystsage = readRDS('20231225_tf_bystsage_h15.rds')

regulonActivity_byCell_scaled = readRDS('regulonActivity_byCell_scaled_raw.rds')

ann = readRDS('ann_regulon.rds')
meta = cancer@meta.data
meta = meta[rownames(ann),]
ann$GRN = meta$Regulon_h15_modi_ano_factor

df = tf_bystsage %>% dplyr::select('regulon','group')
df = distinct(df)

t_test = as.data.frame(matrix(ncol = 5)); colnames(t_test) = c('GRN','Regulon','t_test_p','Mean','Mean_diff')

for (i in df$group %>% unique){
  print(i)
  
  cells = rownames(cancer@meta.data)[cancer$Regulon_h15_modi_ano_factor == i]
  cells = cells[!cells %>% is.na()]
  not = colnames(regulonActivity_byCell_scaled)[!colnames(regulonActivity_byCell_scaled) %in% cells]
  
  tmp = df %>% filter(group == i)
  reg_tmp = regulonActivity_byCell_scaled[tmp$regulon,] %>% t() %>% as.data.frame()
  colnames(reg_tmp) = colnames(reg_tmp) %>% str_remove(' ') %>% str_remove('[)]') %>% str_remove('[(]') %>% str_replace_all('-','_')
  tfs = colnames(reg_tmp) 
  
  reg_tmp$is_condi = ifelse(rownames(reg_tmp) %in% cells, 'Condi','Others')
  
  for (j in tfs){
    formula = as.formula( paste(j, 'is_condi', sep="~") )
    
    test = t.test(formula, data = reg_tmp)
    # if (test$p.value < 0.01 & test$estimate[1] > 1  & test$estimate[1] > test$estimate[2]){
    if (test$p.value < 0.01){
      t_test_tmp = as.data.frame(matrix(ncol = 5)); colnames(t_test_tmp) = c('GRN','Regulon','t_test_p','Mean','Mean_diff')
      t_test_tmp$GRN = i
      t_test_tmp$Regulon = j
      t_test_tmp$t_test_p = test$p.value
      t_test_tmp$Mean = test$estimate[1]
      t_test_tmp$Mean_diff = round(as.numeric(test$estimate[1] - test$estimate[2]),digits = 2)
      
      t_test = rbind(t_test, t_test_tmp)
    }
  }
}

t_test = t_test[-1,]
write.csv(t_test,'20240202_regulon_t_test_all.csv')


## Supp Figure 3A ####
reg = sub('[_][^_]+$', '', t_test$Regulon)
reg[reg %>% str_detect('NKX2_')] = 'NKX2-1_extended'

reg_fin = rownames(regulonActivity_byCell_scaled)[rownames(regulonActivity_byCell_scaled) %>% str_detect(paste(reg,collapse = '|'))]
reg_fin = c(reg_fin, "NKX2-1_extended (28g)")

idx = match( reg, (reg_fin %>% word(1, sep = ' ')) )

reg_fin = reg_fin[idx]
regulonActivity_byCell_scaled_t = regulonActivity_byCell_scaled[reg_fin,]
col_fun = colorRamp2(c(-2,0, 2), c( "blue", "white","red3"))

dend_m = readRDS('dend_m_pheno_raw.rds')

colours <- list( "Annotation" = setNames( c('skyblue','blue4','red') ,c('AR+','AR-','NEPC')),  'Patient' = setNames( patient_col ,cancer$patient %>% levels()), na.value =  setNames(patient_col ,cancer$patient %>% levels())["NA"],  "Site"  = setNames(c('#C7EA46','#738678','#4AFF00'), ann$Site %>% unique()) ,  "GRN" = setNames(col,lev) )

colnames(ann)[6] = 'Annotation'
ann$subtype = NULL
ann$patient = NULL
ann$Cluster = NULL
ann$group = NULL
ann$group_modi = NULL

ann$GRN = factor(ann$GRN , levels = lev)

ann = arrange(ann, 'GRN')

colAnn <- HeatmapAnnotation(df = ann,
                            which = 'col',
                            col = colours,
                            annotation_width = unit(c(1, 4), 'cm'),
                            gap = unit(1, 'mm'),annotation_legend_param = list(
                              Annotation = list(direction = "horizontal", show_annotation_name = F), GRN = list(direction = "horizontal", show_annotation_name = F)))


tiff("20240311_scenic_rss_h15_filt_ar-_ttest.tiff",units="cm", width=20, height=20, res=300)
ComplexHeatmap::Heatmap(regulonActivity_byCell_scaled_t[1: 102, rownames(ann)], name="Regulon activity", cluster_rows = F,cluster_columns = dend_m ,top_annotation = colAnn, show_column_names = F, row_names_gp = grid::gpar(fontsize = 5,fontface = 'bold'), column_split = 2,clustering_method_rows = 'ward.D2', column_gap = unit(c(3), "mm"),  use_raster = F, col = col_fun, column_title=NULL,row_title=NULL)
dev.off()

tiff("20240311_scenic_rss_h15_filt_ar+_ttest.tiff",units="cm", width=20, height=20, res=300)
ComplexHeatmap::Heatmap(regulonActivity_byCell_scaled_t[103:205, rownames(ann)], name="Regulon activity", cluster_rows = F,cluster_columns = dend_m ,top_annotation = colAnn, show_column_names = F, row_names_gp = grid::gpar(fontsize = 5,fontface = 'bold'), column_split = 2,clustering_method_rows = 'ward.D2', column_gap = unit(c(3), "mm"),  use_raster = F, col = col_fun, column_title=NULL,row_title=NULL)
dev.off()

tiff("20240311_scenic_rss_h15_filt_nepc_ttest.tiff",units="cm", width=20, height=20*(92/102), res=300)
ComplexHeatmap::Heatmap(regulonActivity_byCell_scaled_t[206:297, rownames(ann)], name="Regulon activity", cluster_rows = F,cluster_columns = dend_m ,top_annotation = colAnn, show_column_names = F, row_names_gp = grid::gpar(fontsize = 5,fontface = 'bold'), column_split = 2,clustering_method_rows = 'ward.D2', column_gap = unit(c(3), "mm"),  use_raster = F, col = col_fun, column_title=NULL,row_title=NULL)
dev.off()


# Figure 2D ####
regulonActivity_byCell_scaled = readRDS('regulonActivity_byCell_scaled_raw.rds')
ann = readRDS('ann_regulon.rds')

main = c('AR ','TCF7L2','FOXK1','BACH1','FOXP2_extended','IRF2','EGR2','NFATC1','NFATC2','SOX2','SOX4','FOXA2','TCF7L1','CTCF_extended','MAFG','NR1H4','FOSL1','FOSL2','BATF','GATA2','HOXA13','REST_extended','IRX4','IRF7','IRF9','STAT1','STAT2','HNF4G','RELA','NR5A2_extended','CREB3','ZNF189','FLI1','HOXB13','FOXA1','NR3C1','NEUROD1','NEUROD2','SOX11_extended','ONECUT1','ASCL1','NKX2-1_extended','HOXA7_extended','ONECUT2','SOX6','HOXD10_extended','HOXD11_extended')
idx = which(rownames(regulonActivity_byCell_scaled) %>% str_detect(paste(main, collapse = '|')))

reg_fin = sort(rownames(regulonActivity_byCell_scaled)[idx])
reg_fin = reg_fin[-c(7,8,5)]
main[1] = 'AR'

idx = match( main, (reg_fin %>% word(1, sep = ' ')) )

reg_fin = reg_fin[idx]
regulonActivity_byCell_scaled_t = regulonActivity_byCell_scaled[reg_fin,]
saveRDS(regulonActivity_byCell_scaled_t,'regulonActivity_byCell_scaled_h15_rss_filt_ttest.rds')

meta = cancer@meta.data
meta = meta[rownames(ann),]
ann$GRN = meta$Regulon_h15_modi_ano_factor

ann_row = rownames(regulonActivity_byCell_scaled_t) %>% as.data.frame()
colnames(ann_row) = 'TF'

ann_row$GRN = ann_row$TF
ann_row$GRN[2:5] = 'TCF7L2+ WNT'
ann_row$GRN[6:9] = 'IRF2+ Inflammatory' # 
ann_row$GRN[10:13] = 'SOX2/4+ Embryonic EMT'
ann_row$GRN[14:16] ="MAFG+"     
ann_row$GRN[17:19] = 'FOSL1+ AP-1'
ann_row$GRN[20:23] = 'AR+ HOXB13-'
ann_row$GRN[24:27] =  'AR+ IRF+ Inflammatory'
ann_row$GRN[28:30] =  'AR+ GI'
ann_row$GRN[31:33] =  'AR+ HOXB13+'
ann_row$GRN[34:36] =  'AR+ HOXB13+ FOXA1+'
ann_row$GRN[37:40] =  'NEPC-N'
ann_row$GRN[41:44] =  'NEPC-A'
ann_row$GRN[45:47] =  'NEPC-A/SOX6'

ann_row$GRN = factor(ann_row$GRN, levels = c("AR (12g)" ,lev))
saveRDS(ann_row,'ann_row_ttest.rds')

colours_row <- list( 'GRN' = setNames(c("white",col), ann_row$GRN %>% levels()))

ann_row$TF = NULL
ann_row = as.data.frame(ann_row)
colnames(ann_row) = 'GRN'
rowAnn <- HeatmapAnnotation(df = ann_row,
                            which = 'row',
                            col = colours_row,
                            annotation_width = unit(c(1, 4), 'cm'),
                            gap = unit(1, 'mm'),show_annotation_name = FALSE, show_legend = F)

col_fun = colorRamp2(c(-2,0, 2), c( "blue4", "white","red3"))

samir = readxl::read_xlsx('~/Downloads/20240207_Human_SCENIC_metadata.xlsx') %>% as.data.frame()
rownames(samir) = samir$...1
samir = samir[rownames(ann),]

ann$Site = samir$...9
ann$Patient = cancer@meta.data[rownames(ann), 'patient']

patient_col = c(c(c(darken('#EEDC82'),lighten('#F2A900',0.2),'#EEED09','#6F4E37','#CCA01D',darken('#CFFF04'),'gold','#CC7722')),c("#CDDEEB", "#89CFF0", "#0000FF", "#7FFFD4", "#367588", "#1ca9c9", "#0C2340","#A0B6AC" ,"#6F8078",'#CCCCFF' ), c( '#FF6EC7',darken('#FF2400', 0.1),darken('#FF2400', 0.4)))

scale_color_patient <- function(...){
  ggplot2:::manual_scale(
    'color', 
    values = setNames( patient_col ,cancer$patient %>% levels()), na.value =  setNames(patient_col ,cancer$patient %>% levels())["NA"],
    ...
  )
}

colours <- list( "Annotation" = setNames( c('skyblue','blue4','red') ,c('AR+','AR-','NEPC')),  'Patient' = setNames( patient_col ,cancer$patient %>% levels()), na.value =  setNames(patient_col ,cancer$patient %>% levels())["NA"],  "Site"  = setNames(c('#C7EA46','#738678','#4AFF00'), ann$Site %>% unique()) ,  "GRN" = setNames(col,lev) )

colnames(ann)[6] = 'Annotation'
ann$subtype = NULL
ann$Cluster = NULL
ann$group = NULL
ann$group_modi = NULL
ann$patient= NULL

ann = ann[,c(1,3,4,2)]

colAnn <- HeatmapAnnotation(df = ann,
                            which = 'col',
                            col = colours,
                            annotation_width = unit(c(1, 4), 'cm'),
                            gap = unit(1, 'mm'),annotation_legend_param = list(
                              Annotation = list(direction = "horizontal", show_annotation_name = F), GRN = list(direction = "horizontal", show_annotation_name = F)))


pdf('20240310_scenic_h15_filt_ttest_patient_site.pdf',width = 10, height = 10)
ComplexHeatmap::Heatmap(regulonActivity_byCell_scaled_t, name="Regulon activity", cluster_rows = F,cluster_columns = dend_m ,top_annotation = colAnn, show_column_names = F, row_names_gp = grid::gpar(fontsize = 10,fontface = 'bold'), right_annotation =rowAnn,column_split = 2,clustering_method_rows = 'ward.D2', column_gap = unit(c(3), "mm"),  use_raster = F, col = col_fun, column_title = '')
dev.off()   


# Make geneset for Human SCENIC ####
humtomusGeneList <- function(x){
  
  require("biomaRt")
  human = useMart(host="https://dec2021.archive.ensembl.org", "ENSEMBL_MART_ENSEMBL", dataset = "hsapiens_gene_ensembl")
  mouse = useMart(host = "https://dec2021.archive.ensembl.org", "ENSEMBL_MART_ENSEMBL", dataset = "mmusculus_gene_ensembl")
  
  genesV2 = getLDS(attributes = c("hgnc_symbol"), filters ="hgnc_symbol", values = musGenes , mart =human , attributesL = c( "mgi_symbol"), martL = mouse, uniqueRows=T)
  return(genesV2)
}

### Feed back : top 5 ####
#### Mouse gene ###############
regulon = readRDS('./int/2.6_regulons_asGeneSet.Rds')
tf_bystsage = readRDS('20231225_tf_bystsage_h15.rds')

df = as.data.frame(matrix(ncol = 2))
colnames(df) = c("gene_symbol", "gs_cat"     )

for (i in tf_bystsage$group  %>% unique){
  print(i)
  print('!')
  tmp = tf_bystsage[tf_bystsage$group == i,]
  
  if (dim(tmp)[1] >= 5){
    tmp_reg = tmp$regulon[1:5]
  } else {
    tmp_reg = tmp$regulon
  }
  tmp_reg = tmp_reg %>% unique()
  tmp_reg = tmp_reg %>% word(1,sep = ' ')
  genes = c()
  for (ii in tmp_reg){
    print(ii)
    musGenes = c(ii %>% str_remove('_extended'), regulon[[ii %>% word(start = 1, sep = ' ')]])
    covert = humtomusGeneList(musGenes)
    covert = covert$MGI.symbol %>% unique
    
    genes = c(genes, covert)
  }
  genes = genes %>% unique
  tmp = genes %>% as.data.frame()
  
  tmp$gs_cat = i
  colnames(tmp) = colnames(df)
  df = rbind(df,tmp)
}

df
df = df[-1,]
df$gene_symbol = df$gene_symbol %>% word(1, sep = ' ')
df = distinct(df)
saveRDS(df,'Human_h15_human_rss_top5.rds')

df = readRDS('./SCENIC/Human_h15_human_rss_top5.rds')
write.csv(df,'./SCENIC/Human_h15_human_rss_top5.csv')

####### Human gene ###############
regulon = readRDS('./int/2.6_regulons_asGeneSet.Rds')
tf_bystsage = readRDS('20231225_tf_bystsage_h15.rds')

df = as.data.frame(matrix(ncol = 2))
colnames(df) = c("gene_symbol", "gs_cat"     )

for (i in tf_bystsage$group  %>% unique){
  print(i)
  print('!')
  tmp = tf_bystsage[tf_bystsage$group == i,]
  
  if (dim(tmp)[1] >= 5){
    tmp_reg = tmp$regulon[1:5]
  } else {
    tmp_reg = tmp$regulon
  }
  
  tmp_reg = tmp_reg %>% unique()
  genes = c()
  for (ii in tmp_reg){
    print(ii)
    covert = c(ii %>% str_remove('_extended'), regulon[[ii %>% word(start = 1, sep = ' ')]])
    genes = c(genes, covert)
  }
  genes = genes %>% unique
  tmp = genes %>% as.data.frame()
  
  tmp$gs_cat = i
  colnames(tmp) = colnames(df)
  df = rbind(df,tmp)
}

df
df = df[-1,]
df$gs_cat %>% table
df$gene_symbol = df$gene_symbol %>% word(1,sep = ' ')
df = distinct(df)

saveRDS(df,'Human_h15_human_rss_top5_hum.rds')
