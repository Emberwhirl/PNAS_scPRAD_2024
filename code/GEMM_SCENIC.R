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
library(fgsea)
library(randomcoloR)
library(colorspace)
library(dendextend)
library(cowplot)
library(scCustomize)
library(export)
library(clusterProfiler)
set.seed(1234)

cluster_between_groups_mus = function(mat, factor) {
  
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
exprMat = readRDS('integrated_exprMat.rds')
cellInfo = readRDS('integrated_cellInfo.rds')

print('Initialize SCENIC')
org <- "mgi"
dbDir = "/home/jc2545/palmer_scratch/ss/myeloid/"
myDatasetTitle <- "samir"
dbs <- c("mm10__refseq-r80__10kb_up_and_down_tss.mc9nr.feather", 'mm10__refseq-r80__500bp_up_and_100bp_down_tss.mc9nr.feather')
scenicOptions <- initializeScenic(org=org, dbDir=dbDir, dbs=dbs, datasetTitle=myDatasetTitle, nCores=8) 
scenicOptions@inputDatasetInfo$cellInfo = cellInfo
saveRDS(scenicOptions, file="int/scenicOptions.Rds")

print('Gene filtering')
genesKept <- geneFiltering(exprMat, scenicOptions)
exprMat_filtered <- exprMat[genesKept, ]
saveRDS(exprMat_filtered,'integrated_exprMat_filtered.rds')

print('Run correlation')
exprMat_filtered = readRDS('integrated_exprMat_filtered.rds')
runCorrelation(exprMat_filtered, scenicOptions)

print('Start Genie3')
exprMat_filtered_log <- log2(exprMat_filtered+1)
saveRDS(exprMat_filtered_log, 'integrated_exprMat_filtered_log.rds')
scenicOptions = readRDS('int/scenicOptions.Rds')
exprMat_filtered_log = readRDS('integrated_exprMat_filtered_log.rds')
runGenie3(exprMat_filtered_log, scenicOptions, resumePreviousRun = T)

print('Almost done')
exprMat_log <- log2(exprMat+1)
runSCENIC_1_coexNetwork2modules(scenicOptions)
scenicOptions = readRDS('int/scenicOptions.Rds')
runSCENIC_2_createRegulons(scenicOptions)
scenicOptions = readRDS('int/scenicOptions.Rds')
runSCENIC_3_scoreCells(scenicOptions, exprMat_log)

saveRDS(scenicOptions, file="int/scenicOptions.Rds")


# Regulon activity matrix ####
scenicOptions = readRDS('int/scenicOptions.Rds')
regulonAUC <- loadInt(scenicOptions, "aucell_regulonAUC")
regulonAUC <- regulonAUC[onlyNonDuplicatedExtended(rownames(regulonAUC)),]
regulonActivity_byCell = regulonAUC %>% getAUC 
regulonActivity_byCell_Scaled <- scale(t(regulonActivity_byCell), center = T, scale=T) %>% t()
saveRDS(regulonActivity_byCell_Scaled,'regulonActivity_byCell_Scaled.rds')

# Create dendrogram ####
mus = readRDS('/Users/jooyoung/Dropbox/samir/rds/mus_cancer.rds')
regulonActivity_byCell_scaled = readRDS('~/Documents/samir_macbook/SCENIC_mus/regulonActivity_byCell_Scaled.rds')
DimPlot(mus, group.by = 'clusterings_phenoLV_tuned')

ann= mus@meta.data
regulonActivity_byCell_scaled = regulonActivity_byCell_scaled[,rownames(ann)]

dend3 = cluster_between_groups_mus(regulonActivity_byCell_scaled,ann$clusterings_phenoLV_tuned)
saveRDS(dend3,'dend_mus_pheno_cluster_tune.rds')

# Make GRN with ARI Hclust ####
dend3 = readRDS('dend_mus_pheno_cluster_tune.rds')
hc1 = as.hclust(dend3)
plot(hc1)

clusterings = list()
n = 1
for (k in seq(0,24,1)){
  print(k)
  cut_avg <- cutree(hc1, h = k) 
  clusterings[[n]] = cut_avg
  n = n+1
}
names(clusterings) = paste0('h_',seq(0,24,1))
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
rownames(rand_matrix) = colnames(rand_matrix) = paste0('h_',seq(0,24,1))
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

rand_matrix = read.csv('20231224_ARI_mus_pheno_tuned_step1.csv')
rownames(rand_matrix) = rand_matrix$X; rand_matrix$X = NULL
ComplexHeatmap::Heatmap(rand_matrix, cluster_rows = F, cluster_columns = F,name="Adjusted Rand Index", na_col = lighten('lightgray',0.3), col = col_fun)
write.csv(rand_matrix,'20231224_ARI_mus_pheno_tuned_step1.csv')

# Add GRN information ####
cut_avg <- cutree(hc1, h = 12) 
mus_reg = paste0('Mus_',cut_avg)
mus_reg = mus_reg %>% as.data.frame()
mus_reg$barcode = colnames(regulonActivity_byCell_scaled)
rownames(mus_reg) = mus_reg$barcode
mus_reg = mus_reg[colnames(mus),]

mus$test = mus_reg$.
dend3 = cluster_between_groups_mus(regulonActivity_byCell_scaled,mus$test)
hc1 = as.hclust(dend3)
plot(hc1)
saveRDS(dend3,'dend_mus_pheno_tune.rds')

# Figure 3A ####
## GRN annotation map ####
mus$GRN_pheno_tune = ifelse(mus$test %in% c('Mus_6','Mus_7'),'NEPC-A',ifelse(mus$test == 'Mus_9','Ar/Ascl1+', ifelse(mus$test %in% c('Mus_12','Mus_3'), 'Trp63+ Sox4/6',ifelse(mus$test %in% c('Mus_10','Mus_13'),'Twist1/2 EMT',ifelse(mus$test == 'Mus_5','Pou2f3+',ifelse(mus$test == 'Mus_1','Tff3+',ifelse(mus$test == 'Mus_11','Stat1/2 Inflam',ifelse(mus$test == 'Mus_2','Trp63 Basal','CK8 Luminal'))))))))
dend3 = cluster_between_groups_mus(regulonActivity_byCell_scaled,mus$GRN_pheno_tune)
hc1 = as.hclust(dend3)
plot(hc1)
saveRDS(dend3,'dend_mus_pheno_tune_ano.rds')

mus = readRDS('/Users/jooyoung/Dropbox/samir/rds/mus_cancer.rds')
regulonActivity_byCell_scaled = readRDS('regulonActivity_byCell_Scaled.rds')
colours <- list( 'GRN_pheno_tune' = setNames(c("#FF0000","orange","#8000FF" ,"#FF00BF","#80461B","#0040FF" , "#00FFFF","gold" ,"#00FF40"   ), unique(mus$GRN_pheno_tune)))


ann_col = mus@meta.data %>% dplyr::select('GRN_pheno_tune')
colAnn <- HeatmapAnnotation(df = ann_col,
                            which = 'col',
                            col = colours,
                            annotation_width = unit(c(1, 4), 'cm'),
                            gap = unit(1, 'mm'),show_annotation_name =T,show_legend = T)

dend3 = readRDS('dend_mus_pheno_tune_ano.rds')

ann_row = readRDS('ann_row_ano.rds') # row annotation: data frame with TF annotated as specific GRN group
col_fun = colorRamp2(c(-2,0, 2), c( "blue4", "white","red3"))
dim(ann_row)
ann_row$reg = rownames(ann_row)
ann_row$GRN[rownames(ann_row) %in% c('Atoh1_extended (51g)','Tff3 (62g)',"Spdef (34g)"  )] = 'GEMM GRN 5_1'
ann_row = ann_row[! rownames(ann_row) %>% str_detect('Gata2|Onecut2|Hnf4a|Cdx|Irf1'),]
ann_row[c('Ar_extended (14g)','Ascl1 (150g)'),'GRN'] = 'NAN'
ann_row$GRN = factor(ann_row$GRN, levels = c('NAN', "GEMM GRN 1" ,  "GEMM GRN 2" ,  "GEMM GRN 3" ,"GEMM GRN 8","GEMM GRN 5_1", "GEMM GRN 4" ,"GEMM GRN 5" ,"GEMM GRN 6", "GEMM GRN 7" ))
ann_row = arrange(ann_row,GRN)
ann_row$reg = NULL


colours_row <- list( 'GRN' = setNames(c("white","#00FFFF","#0040FF", "gold",'#80461B',"#FF0000", "orange", "#FF00BF", "#00FF40" , "#8000FF" ), unique(ann_row$GRN)))

rowAnn <- HeatmapAnnotation(df = ann_row,
                            which = 'row',
                            col = colours_row,
                            annotation_width = unit(c(1, 4), 'cm'),
                            gap = unit(1, 'mm'),show_annotation_name = FALSE)

ann_col = mus@meta.data
m_sel  = regulonActivity_byCell_scaled[rownames(ann_row),rownames(ann_col)]

hmap = ComplexHeatmap::Heatmap(m_sel,cluster_columns = dend3, name="Regulon activity", cluster_rows = F,col = col_fun,top_annotation = colAnn, show_column_names = F, row_names_gp = grid::gpar(fontsize = 20,fontface = 'bold'), show_column_dend = T, column_title = " ",heatmap_legend_param = list(label_gp = gpar(fontsize = 10)), right_annotation =rowAnn, column_split = 9,  column_gap = unit(c(5), "mm")) # 

# pdf('20231224_GEMM_pheno_tuned_ano_main.pdf', width = 17, height = 20)
pdf('20240313_GEMM_pheno_tuned_ano_main_rev_v1.pdf', width = 17, height = 20)
draw(hmap, heatmap_legend_side="right", annotation_legend_side="right",background = "transparent",legend_grouping = "original")
dev.off()


### RSS ####
scenicOptions = readRDS('int/scenicOptions.Rds')
regulonAUC <- loadInt(scenicOptions, "aucell_regulonAUC")
regulonAUC = regulonAUC[rownames(regulonActivity_byCell_scaled),rownames(ann_col)]
rss <- calcRSS(AUC=getAUC(regulonAUC), cellAnnotation=ann_col[colnames(regulonAUC), "GRN_pheno_tune"])
rssNorm <- scale(rss) 
write.csv(rssNorm,'20231224_tf_rss_pheno_tuned_ano.csv')


tf_bystsage = as.data.frame(matrix(ncol = 4)); colnames(tf_bystsage) = c('regulon','rank','rss','group')
# for (i in ann$group %>% unique){
for (i in ann_col$GRN_pheno_tune %>% unique){
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
tf_bystsage$group = factor(tf_bystsage$group, levels = c('NEPC-A','Ar/Ascl1+','Trp63+ Sox4/6','Twist1/2 EMT','Pou2f3+','Tff3+','Stat1/2 Inflam','Trp63 Basal','CK8 Luminal'))
tf_bystsage = arrange(tf_bystsage,group, rank)


mat = regulonActivity_byCell_scaled[tf_bystsage$regulon,]
dend3 = readRDS('dend_mus_pheno_tune_ano.rds')
col_fun = colorRamp2(c(-2,0, 2), c( "blue", "white","red3"))
tf_bystsage = readRDS('20231224_tf_rss_pheno_tuned_ano.rds')
colAnn = readRDS('colAnn_pheno_tune_ano.rds')

tf_bystsage$group  = tf_bystsage$group %>% as.character()
tf_bystsage$group[tf_bystsage$group %>% is.na()] = 'NEPC-A'
tf_bystsage$group = factor(tf_bystsage$group, levels = c('NEPC-A','Ar/Ascl1+','Trp63+ Sox4/6','Twist1/2 EMT','Pou2f3+','Tff3+','Stat1/2 Inflam','Trp63 Basal','CK8 Luminal'))
tf_bystsage = arrange(tf_bystsage,group, rank)
table(tf_bystsage$group)

pdf('20231224_GEMM_pheno_tuned_ano_all_rss_1.pdf', width = 10, height =16*75/117)
ComplexHeatmap::Heatmap(mat[1:75,],cluster_columns = dend3, name="Regulon activity", cluster_rows = F,top_annotation = colAnn, show_column_names = F,  show_column_dend = T, column_title = " ",heatmap_legend_param = list(label_gp = gpar(fontsize = 15,fontface = 'bold')),  column_split = 13,  column_gap = unit(c(2), "mm"), col = col_fun, row_names_gp = grid::gpar(fontsize = 10,fontface = 'bold'))
dev.off()

pdf('20231224_GEMM_pheno_tuned_ano_all_rss_2.pdf', width = 10, height =16*106/117)
ComplexHeatmap::Heatmap(mat[76:181,],cluster_columns = dend3, name="Regulon activity", cluster_rows = F,top_annotation = colAnn, show_column_names = F,  show_column_dend = T, column_title = " ",heatmap_legend_param = list(label_gp = gpar(fontsize = 15,fontface = 'bold')),  column_split = 13,  column_gap = unit(c(2), "mm"), col = col_fun, row_names_gp = grid::gpar(fontsize = 10,fontface = 'bold'))
dev.off()


pdf('20231224_GEMM_pheno_tuned_ano_all_rss_3.pdf', width = 10, height =16*137/117)
ComplexHeatmap::Heatmap(mat[182:319,],cluster_columns = dend3, name="Regulon activity", cluster_rows = F,top_annotation = colAnn, show_column_names = F,  show_column_dend = T, column_title = " ",heatmap_legend_param = list(label_gp = gpar(fontsize = 15,fontface = 'bold')),  column_split = 13,  column_gap = unit(c(2), "mm"), col = col_fun, row_names_gp = grid::gpar(fontsize = 10,fontface = 'bold'))
dev.off()


# Make geneset for GEMM SCENIC ####
mustohumGeneList <- function(x){
  
  require("biomaRt")
  human = useMart(host="https://dec2021.archive.ensembl.org", "ENSEMBL_MART_ENSEMBL", dataset = "hsapiens_gene_ensembl")
  mouse = useMart(host = "https://dec2021.archive.ensembl.org", "ENSEMBL_MART_ENSEMBL", dataset = "mmusculus_gene_ensembl")
  
  genesV2 = getLDS(attributes = c("mgi_symbol"), filters ="mgi_symbol", values = musGenes , mart =mouse , attributesL = c( "hgnc_symbol"), martL = human, uniqueRows=T)
  return(genesV2)
}

### Feed back : top 5 ####
#### Human gene ###############
tf_bystsage = readRDS('20231224_tf_rss_pheno_tuned_ano.rds')
tf_bystsage = tf_bystsage[tf_bystsage$rss > 0,]
tf_bystsage$group = tf_bystsage$group  %>% as.character()
tf_bystsage$group[tf_bystsage$group %>% is.na()] = 'NEPC-A'

mus_reg = readRDS('int/2.6_regulons_asGeneSet.Rds')
top5 = list()
for ( i in unique(tf_bystsage$group)){
  dff = tf_bystsage[tf_bystsage$group == i,]
  if (max(dff$rank) >= 5){
    top5[[i]] =  dff$regulon[1:5]
  } 
}

df = matrix(ncol = 3) %>% as.data.frame(); colnames(df) = c('GRN','reg','gene')
for (i in names(top5)){
  print(i)
  
  dff = matrix(ncol = 3) %>% as.data.frame(); colnames(dff) = c('GRN','reg','gene')
  
  for (j in top5[[i]]){
    g = c()
    g = c(g,mus_reg[[j %>% word(1,sep = ' ')]])
    musGenes = c(g,(j %>% word(1,sep = ' ') %>% str_remove('_extended')))
    reg_tohum = mustohumGeneList(musGenes) 
    
    tmp = matrix(ncol = 3, nrow = length(reg_tohum)) %>% as.data.frame(); colnames(tmp) = c('GRN','reg','gene')
    tmp$gene = reg_tohum
    tmp$reg = j
    dff = rbind(dff, tmp)
  }
  
  dff = dff[-1,]; dff$GRN = i
  df = rbind(df,dff)
}

df = df[-1,]
df$reg = NULL
df = distinct(df)
saveRDS(df,'GEMM_GRN_geneset_rss_top5.rds')
