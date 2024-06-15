setwd('/Users/jooyoung/Dropbox/samir/rds_jy/SCENIC/robust/')

library(dplyr)
library(stringr)
set.seed(1234)

# Input prepare - By cluster ####
#### 10 cells ####
cancer = readRDS('/Users/jooyoung/Dropbox/samir/rds_jy/msk.integrated.remove.cellcycle.tumor.fastmnn.subtype.nohepato.figure4.rds')
meta = cancer@meta.data

cells = c()
for (i in cancer$pheno_cluster %>% unique){
  print(i)
  df = meta[meta$pheno_cluster == i,]
  rand = sample(rownames(df), 10)
  cells = c(cells, rand)
}

cellinfo = readRDS('/Users/jooyoung/Dropbox/samir/rds_jy/SCENIC/integrated_cellInfo.rds')
cellinfo = cellinfo[cells,]
dir.create('/Users/jooyoung/Dropbox/samir/rds_jy/SCENIC/robust')
dir.create('/Users/jooyoung/Dropbox/samir/rds_jy/SCENIC/robust/sample_10')
saveRDS(cellinfo , './integrated_cellInfo.rds')

exp = readRDS('/Users/jooyoung/Dropbox/samir/rds_jy/SCENIC/integrated_exprMat.rds')
idx = colnames(exp) %in% cells
exp = exp[,idx]
saveRDS(exp , './integrated_exprMat.rds') # 310


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

# Analysis start - Get data ####
msk = readRDS('../../../msk.integrated.remove.cellcycle.tumor.fastmnn.subtype.nohepato.figure4.rds')
scenicOptions = readRDS('int/scenicOptions.Rds')
regulonAUC <- loadInt(scenicOptions, "aucell_regulonAUC")
regulonAUC <- regulonAUC[onlyNonDuplicatedExtended(rownames(regulonAUC)),]
regulonAUC = regulonAUC %>% getAUC 
dim(regulonAUC)
regulonActivity_byCell_scaled <- scale(t(regulonAUC), center = T, scale=T) %>% t()
saveRDS(regulonActivity_byCell_scaled,'regulonActivity_byCell_scaled_raw.rds') # 416 310

## Compare regulon ####
prev = readRDS('../../regulonActivity_byCell_scaled_raw.rds')
prev = rownames(prev)  %>% word(1,sep = ' ')
now = rownames(regulonActivity_byCell_scaled) %>% word(1, sep = ' ')

idx = intersect(prev, now)

length(prev)
length(now)

# Hclust ####
ann = subset(msk, cells = colnames(regulonActivity_byCell_scaled))
ann = ann@meta.data

regulonActivity_byCell_scaled = readRDS('regulonActivity_byCell_scaled_raw.rds')
regulonActivity_byCell_scaled = regulonActivity_byCell_scaled[,rownames(ann)]

## CRPC ####
mat_crpc =regulonActivity_byCell_scaled[,rownames(ann)[ann$subtype == 'CRPC']]
ann_crpc = ann[ann$subtype == 'CRPC',]
ann_crpc$pheno_cluster = ann_crpc$pheno_cluster %>% as.character()
mat_crpc = mat_crpc[,rownames(ann_crpc)]
dend3 = cluster_between_groups_ward(mat_crpc,ann_crpc$pheno_cluster)
dend3
hc1 = as.hclust(dend3)
plot(hc1)
saveRDS(dend3, '20231228_CRPC_dend.rds')

##### Adjusted Rand Index ####
dend3 = readRDS('20231228_CRPC_dend.rds')
hc1 = as.hclust(dend3)

plot(hc1)
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
write.csv(rand_matrix,'20231228_ARI_CRPC_step1_score.csv')

col_fun = colorRamp2(c(0,0.65, 1), c( "blue", "white","red3"))

# rand_matrix = read.csv('20231227_ARI_CRPC_step1_score_wrong.csv')
# rownames(rand_matrix) = rand_matrix$X; rand_matrix$X = NULL
ComplexHeatmap::Heatmap(rand_matrix, cluster_rows = F, cluster_columns = F,name="Adjusted Rand Index", na_col = lighten('lightgray',0.3), col = col_fun) # 14

cut_avg <- cutree(hc1, h = 12) # 12
cut_avg %>% table

crpc_reg = paste0('CRPC_',cut_avg)
crpc_reg = crpc_reg %>% as.data.frame()
crpc_reg$barcode = colnames(mat_crpc)
colnames(crpc_reg)[1] = 'group'
saveRDS(crpc_reg, '20231228_CRPC_regulon_h12.rds')

## NEPC ####
mat_nepc =regulonActivity_byCell_scaled[,rownames(ann)[ann$subtype == 'NEPC']]
ann_nepc = ann[ann$subtype == 'NEPC',]
ann_nepc$pheno_cluster = ann_nepc$pheno_cluster %>% as.character()
mat_nepc = mat_nepc[,rownames(ann_nepc)]
dend3 = cluster_between_groups_ward(mat_nepc,ann_nepc$pheno_cluster)
dend3
hc1 = as.hclust(dend3)
plot(hc1)
saveRDS(dend3, '20231228_NEPC_dend.rds')

plot(hc1)

##### Adjusted Rand Index ####
dend3 = readRDS('20231228_NEPC_dend.rds')
hc1 = as.hclust(dend3)

plot(hc1)
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
  # score = adjustedRandIndex(clusterings[[comp1]], clusterings[[comp2]])
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
write.csv(rand_matrix,'20231228_ARI_NEPC_step1_score.csv')

col_fun = colorRamp2(c(0,0.65, 1), c( "blue", "white","red3"))
ComplexHeatmap::Heatmap(rand_matrix, cluster_rows = F, cluster_columns = F,name="Adjusted Rand Index", na_col = lighten('lightgray',0.3), col = col_fun) # 19

# rand_matrix = read.csv('20231227_ARI_NEPC_step1_score_wrong.csv')
# rownames(rand_matrix) = rand_matrix$X; rand_matrix$X = NULL

cut_avg <- cutree(hc1, h = 12)  
cut_avg %>% table

nepc_reg = paste0('NEPC_',cut_avg)
nepc_reg = nepc_reg %>% as.data.frame()
nepc_reg$barcode = colnames(mat_nepc)
colnames(nepc_reg)[1] = 'group'
saveRDS(nepc_reg, '20231228_NEPC_regulon_h12.rds')

## Add GRN information ####
crpc_reg = readRDS('20231228_CRPC_regulon_h12.rds')
nepc_reg = readRDS('20231228_NEPC_regulon_h12.rds')
colnames(crpc_reg)[1] = 'group'; colnames(nepc_reg)[1] = 'group'

reg = rbind(crpc_reg, nepc_reg)
rownames(reg) = reg$barcode
reg = reg[rownames(ann),]
which(rownames(reg) != rownames(ann))
ann$group = reg$group
saveRDS(ann,'20231228_ann_regulon.rds')

## Merge dendrogram ####
dend3 = readRDS('20231228_CRPC_dend.rds')
dend4 = readRDS('20231228_NEPC_dend.rds')
dend_p = as.dendrogram(hclust(dist(rbind(colMeans(t(mat_crpc)), colMeans(t(mat_nepc))))))
dend_m = merge_dendrogram(dend_p, list(dend3, dend4))
order_col = c(colnames(mat_crpc),colnames(mat_nepc))
ann[order_col[order.dendrogram(dend_m)],'group'] %>% unique
ann = ann[order_col,]
saveRDS(dend_m,'20231228_dend_m.rds')

## Calculate ARI ####
# score = genieclust::adjusted_rand_score(ann$group, ann$Regulon_h15_modi_ano) # needs to be numeric
library(mclust)
table(ann$group, ann$subtype)
adjustedRandIndex(ann$group, ann$Regulon_h15_modi_ano) # 0.69

# Heatmap ####
## Row clustering ####
all = readRDS('../../../msk.integrated.remove.cellcycle.tumor.fastmnn.subtype.nohepato.rds')

ann$Regulon_h15_modi_ano
ann_new = ann %>% dplyr::select('group', 'subtype','Regulon_h15_modi_ano' )
ann_new$Regulon_h15_modi_ano[ann_new$Regulon_h15_modi_ano  == 'TCFL2+ WNT'] = 'TCF7L2+ WNT'
col  = readRDS('../../../regulon_col.rds')

colours <- list(
  'subtype' = c( 'CRPC' = '#2171B5','NEPC' ='#CB181D'),
  "Regulon_h15_modi_ano" = setNames( col , lev),
  "group" = setNames(sample(rainbow(length(unique(ann$group)))), ann$group %>% unique())
)

colAnn <- HeatmapAnnotation(df = ann_new,
                            which = 'col',
                            col = colours,
                            annotation_width = unit(c(1, 4), 'cm'),
                            gap = unit(1, 'mm'),annotation_legend_param = list(
                              subtype = list(direction = "horizontal")))

col_fun = colorRamp2(c(-2,0, 2), c( "blue", "white","red3"))

regulonActivity_byCell_scaled = regulonActivity_byCell_scaled[,rownames(ann_new)]

png('test.png',width = 20*100, height = 80*100)
ComplexHeatmap::Heatmap(regulonActivity_byCell_scaled, name="Regulon activity", cluster_rows = T,cluster_columns = dend_m ,top_annotation = colAnn, show_column_names = F, row_names_gp = grid::gpar(fontsize = 25,fontface = 'bold'), column_split = 2,clustering_method_rows = 'ward.D2', column_gap = unit(c(10), "mm"),  use_raster = F, col = col_fun)
dev.off()

## RSS ####
scenicOptions = readRDS('int/scenicOptions.Rds')
regulonAUC <- loadInt(scenicOptions, "aucell_regulonAUC")
regulonAUC = regulonAUC[rownames(regulonActivity_byCell_scaled),]
rss <- calcRSS(AUC=getAUC(regulonAUC), cellAnnotation=ann[colnames(regulonAUC), "group"])
rssNorm <- scale(rss) 

tf_bystsage = as.data.frame(matrix(ncol = 4)); colnames(tf_bystsage) = c('regulon','rank','rss','group')
# for (i in ann$group %>% unique){
for (i in ann$group %>% unique){
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
tf_bystsage$group = factor(tf_bystsage$group, levels = ann_new$group[order.dendrogram(dend_m)] %>% unique)
tf_bystsage = arrange(tf_bystsage,group, rank)
tf_bystsage$group %>% table

rowsplit = c()
n = 1
for (i in tf_bystsage$group %>% levels()){
  n = n + 1
  print(i)
  tmp = tf_bystsage[tf_bystsage$group == i,]
  rowsplit = c(rowsplit, rep(paste0(LETTERS[n],'_',i),tmp$regulon %>% length() ) )
}

pdf('test_rss.pdf',width = 20, height = 50)
ComplexHeatmap::Heatmap(regulonActivity_byCell_scaled[tf_bystsage$regulon,], name="Regulon activity", cluster_rows = F,cluster_columns = dend_m ,top_annotation = colAnn, show_column_names = F, row_names_gp = grid::gpar(fontsize = 10,fontface = 'bold'), column_split = 2,clustering_method_rows = 'ward.D2', column_gap = unit(c(10), "mm"),  use_raster = F, col = col_fun,row_split = rowsplit)
dev.off()

pdf("20231228_rss_1.pdf",  width=15, height=10*135/96)
ComplexHeatmap::Heatmap(regulonActivity_byCell_scaled[tf_bystsage$regulon[1: 135], ], name="Regulon activity", cluster_rows = F,cluster_columns = dend_m ,top_annotation = colAnn, show_column_names = F, row_names_gp = grid::gpar(fontsize = 5,fontface = 'bold'), column_split = 2,clustering_method_rows = 'ward.D2', column_gap = unit(c(3), "mm"),  use_raster = F, col = col_fun, column_title=NULL,row_title=NULL)
dev.off()

pdf("20231228_rss_2.pdf",width=15, height= 10*151/96)
ComplexHeatmap::Heatmap(regulonActivity_byCell_scaled[tf_bystsage$regulon[136:285] , ], name="Regulon activity", cluster_rows = F,cluster_columns = dend_m ,top_annotation = colAnn, show_column_names = F, row_names_gp = grid::gpar(fontsize = 5,fontface = 'bold'), column_split = 2,clustering_method_rows = 'ward.D2', column_gap = unit(c(3), "mm"),  use_raster = F, col = col_fun, column_title=NULL,row_title=NULL)
dev.off()

pdf("20231228_rss_3.pdf",  width=15, height= 10*129/96)
ComplexHeatmap::Heatmap(regulonActivity_byCell_scaled[tf_bystsage$regulon[286:length(tf_bystsage$group)], ], name="Regulon activity", cluster_rows = F,cluster_columns = dend_m ,top_annotation = colAnn, show_column_names = F, row_names_gp = grid::gpar(fontsize = 5,fontface = 'bold'), column_split = 2,clustering_method_rows = 'ward.D2', column_gap = unit(c(3), "mm"),  use_raster = F, col = col_fun, column_title=NULL,row_title=NULL)
dev.off()





# Input prepare - By patient 100 random cells ####
## 6 patients ####
cancer = readRDS('/Users/jooyoung/Dropbox/samir/rds_jy/msk.integrated.remove.cellcycle.tumor.fastmnn.subtype.nohepato.figure4.rds')
cellinfo = readRDS('/Users/jooyoung/Dropbox/samir/rds_jy/SCENIC/integrated_cellInfo.rds')
exp = readRDS('/Users/jooyoung/Dropbox/samir/rds_jy/SCENIC/integrated_exprMat.rds')

meta = cancer@meta.data

meta$patient %>% table
patients_over100 = table(meta$patient)[table(meta$patient) > 100] # 11 patients

for (numb in 1:10){
  print(numb)
  dir.create('./patient_6/samir', showWarnings = F)
  patient = names(patients_over100) %>% sample(6) %>% as.character()
  
  cells = c()
  for (i in patient){
    print(i)
    df = meta[meta$patient == i,]
    rand = sample(rownames(df), 100)
    cells = c(cells, rand)
  }
  
  cellinfo$barcode = rownames(cellinfo)
  cellinfo_sub = cellinfo[cells,]
  cellinfo_sub$barcode = NULL
  
  dir.create(paste0('./patient_6/samir/',numb), showWarnings = F)
  saveRDS(cellinfo_sub , paste0('./patient_6/samir/',numb,'/integrated_cellInfo.rds') )
  
  idx = colnames(exp) %in% cells
  exp_sub = exp[,idx]
  
  dir.create(paste0('./patient_6/samir/',numb), showWarnings = F)
  saveRDS(exp_sub , paste0('./patient_6/samir/',numb,'/integrated_exprMat.rds') )
}

cellinfo = readRDS('/Users/jooyoung/Dropbox/samir/rds_jy/SCENIC/integrated_cellInfo.rds')
exp = readRDS('/Users/jooyoung/Dropbox/samir/rds_jy/SCENIC/integrated_exprMat.rds')
for (numb in c(1,6)){
  print(numb)
  dir.create('./patient_6/samir', showWarnings = F)
  patient = names(patients_over100) %>% sample(6) %>% as.character()
  
  cells = c()
  for (i in patient){
    print(i)
    df = meta[meta$patient == i,]
    rand = sample(rownames(df), 100)
    cells = c(cells, rand)
  }
  
  cellinfo$barcode = rownames(cellinfo)
  cellinfo_sub = cellinfo[cells,]
  cellinfo_sub$barcode = NULL
  
  dir.create(paste0('./patient_6/samir/',numb), showWarnings = F)
  saveRDS(cellinfo_sub , paste0('./patient_6/samir/',numb,'/integrated_cellInfo.rds') )
  
  idx = colnames(exp) %in% cells
  exp_sub = exp[,idx]
  
  dir.create(paste0('./patient_6/samir/',numb), showWarnings = F)
  saveRDS(exp_sub , paste0('./patient_6/samir/',numb,'/integrated_exprMat.rds') )
}


# Analysis start - Get data ####
msk = readRDS('../../../../msk.integrated.remove.cellcycle.tumor.fastmnn.subtype.nohepato.figure4.rds')

# Loop ####
dir = list.dirs( recursive = F) %>% str_remove_all('./')

dir = c('1','6')
for (dd in dir){
  print(dd)
  
  setwd(dd)
  tryCatch( {
    ## Regulon matrix ####
    print('Regulon matrix')
    instep = "Regulon matrix"
    
    scenicOptions = readRDS('int/scenicOptions.Rds')
    regulonAUC <- loadInt(scenicOptions, "aucell_regulonAUC")
    regulonAUC <- regulonAUC[onlyNonDuplicatedExtended(rownames(regulonAUC)),]
    regulonAUC = regulonAUC %>% getAUC 
    regulonActivity_byCell_scaled <- scale(t(regulonAUC), center = T, scale=T) %>% t()
    saveRDS(regulonActivity_byCell_scaled,'regulonActivity_byCell_scaled_raw.rds')
    
    ## Compare regulon ####
    print('Compare regulon')
    instep = "Compare regulon"
    
    prev = readRDS('../../../../regulonActivity_byCell_scaled_raw.rds')
    prev = rownames(prev)  %>% word(1,sep = ' ')
    now = rownames(regulonActivity_byCell_scaled) %>% word(1, sep = ' ')
    idx = intersect(prev, now)
    
    n = length(idx)
    
    fisher = data.frame(c(n,300-n),c( length(now) - n ,1841-(300 +length(now) - n)))
    test <- fisher.test(fisher) 
    write.csv(c(test$p.value, n %>% as.integer()),'fisher.csv')
    
    # Hclust ####
    print('Hclust _ CRPC')
    instep = "Hclust _ CRPC"
    
    ann = subset(msk, cells = colnames(regulonActivity_byCell_scaled))
    ann = ann@meta.data
    
    regulonActivity_byCell_scaled = readRDS('regulonActivity_byCell_scaled_raw.rds')
    regulonActivity_byCell_scaled = regulonActivity_byCell_scaled[,rownames(ann)]
    
    
    ## CRPC ####
    mat_crpc =regulonActivity_byCell_scaled[,rownames(ann)[ann$subtype == 'CRPC']]
    ann_crpc = ann[ann$subtype == 'CRPC',]
    ann_crpc$pheno_cluster = ann_crpc$pheno_cluster %>% as.character()
    mat_crpc = mat_crpc[,rownames(ann_crpc)]
    dend3 = cluster_between_groups_ward(mat_crpc,ann_crpc$pheno_cluster)
    dend3
    hc1 = as.hclust(dend3)
    plot(hc1)
    saveRDS(dend3, '20240226_CRPC_dend.rds')
    
    ##### Feedback - Adjusted Rand Index ####
    dend3 = readRDS('20240226_CRPC_dend.rds')
    hc1 = as.hclust(dend3)
    
    plot(hc1)
    clusterings = list()
    n = 1
    for (k in seq(0,40,1)){
      print(k)
      cut_avg <- cutree(hc1, h = k) 
      clusterings[[n]] = cut_avg
      n = n+1
    }
    names(clusterings) = paste0('h_',seq(0,40,1))
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
    rownames(rand_matrix) = colnames(rand_matrix) = paste0('h_',seq(0,40,1))
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
    write.csv(rand_matrix,'20240226_ARI_CRPC_step1_score.csv')
    
    hcutoff = as.numeric(colnames(rand_matrix)[which(rand_matrix[1,] < 0.65 ) [1]] %>% word(start = 2, sep = '_')) - 1
    
    cut_avg <- cutree(hc1, h = hcutoff)
    crpc_reg = paste0('CRPC_',cut_avg)
    crpc_reg = crpc_reg %>% as.data.frame()
    crpc_reg$barcode = colnames(mat_crpc)
    colnames(crpc_reg)[1] = 'group'
    saveRDS(crpc_reg, paste0('20240226_CRPC_regulon_',hcutoff,'.rds'))
    
    
    ## NEPC ####
    print('Hclust _ NEPC')
    instep = "Hclust _ NEPC"
    
    mat_nepc =regulonActivity_byCell_scaled[,rownames(ann)[ann$subtype == 'NEPC']]
    ann_nepc = ann[ann$subtype == 'NEPC',]
    ann_nepc$pheno_cluster = ann_nepc$pheno_cluster %>% as.character()
    mat_nepc = mat_nepc[,rownames(ann_nepc)]
    dend3 = cluster_between_groups_ward(mat_nepc,ann_nepc$pheno_cluster)
    dend3
    hc1 = as.hclust(dend3)
    plot(hc1)
    saveRDS(dend3, '20240226_NEPC_dend.rds')
    
    plot(hc1)
    
    ##### Feedback - Adjusted Rand Index ####
    dend3 = readRDS('20240226_NEPC_dend.rds')
    hc1 = as.hclust(dend3)
    
    cut_avg <- cutree(hc1, h = hcutoff)  
    cut_avg %>% table
    
    nepc_reg = paste0('NEPC_',cut_avg)
    nepc_reg = nepc_reg %>% as.data.frame()
    nepc_reg$barcode = colnames(mat_nepc)
    colnames(nepc_reg)[1] = 'group'
    saveRDS(nepc_reg, paste0('20240226_NEPC_regulon_',hcutoff,'.rds'))
    
    ## Add GRN information ####
    print('Add GRN information')
    instep = "Add GRN information"
    colnames(crpc_reg)[1] = 'group'; colnames(nepc_reg)[1] = 'group'
    
    reg = rbind(crpc_reg, nepc_reg)
    rownames(reg) = reg$barcode
    reg = reg[rownames(ann),]
    which(rownames(reg) != rownames(ann))
    ann$group = reg$group
    saveRDS(ann,'20240226_ann_regulon.rds')
    
    
    ## Merge dendrogram ####
    print('Merge dendrogram')
    instep = "Merge dendrogram"
    
    dend3 = readRDS('20240226_CRPC_dend.rds')
    dend4 = readRDS('20240226_NEPC_dend.rds')
    dend_p = as.dendrogram(hclust(dist(rbind(colMeans(t(mat_crpc)), colMeans(t(mat_nepc))))))
    dend_m = merge_dendrogram(dend_p, list(dend3, dend4))
    order_col = c(colnames(mat_crpc),colnames(mat_nepc))
    ann[order_col[order.dendrogram(dend_m)],'group'] %>% unique
    ann = ann[order_col,]
    saveRDS(dend_m,'20240226_dend_m.rds')
    
    
    ## colAnn ####
    print('colAnn')
    instep = "colAnn"
    
    ann_new = ann %>% dplyr::select('group', 'Regulon_h15_modi_ano' )
    ann_new$Regulon_h15_modi_ano[ann_new$Regulon_h15_modi_ano  == 'TCFL2+ WNT'] = 'TCF7L2+ WNT'
    col  = readRDS('../../../../../regulon_col.rds')
    lev = readRDS('../../../../../regulon_lev.rds')
    
    colours <- list(
      "Regulon_h15_modi_ano" = setNames( col , lev),
      "group" = setNames(sample(rainbow(length(unique(ann$group)))), ann$group %>% unique())
    )
    
    colAnn <- HeatmapAnnotation(df = ann_new,
                                which = 'col',
                                col = colours,
                                annotation_width = unit(c(1, 4), 'cm'),
                                gap = unit(1, 'mm'))
    
    ## RSS ####
    print('RSS')
    instep = "RSS"
    scenicOptions = readRDS('int/scenicOptions.Rds')
    regulonAUC <- loadInt(scenicOptions, "aucell_regulonAUC")
    regulonAUC = regulonAUC[rownames(regulonActivity_byCell_scaled),]
    rss <- calcRSS(AUC=getAUC(regulonAUC), cellAnnotation=ann[colnames(regulonAUC), "group"])
    rssNorm <- scale(rss) 
    
    tf_bystsage = as.data.frame(matrix(ncol = 4)); colnames(tf_bystsage) = c('regulon','rank','rss','group')
    # for (i in ann$group %>% unique){
    for (i in ann$group %>% unique){
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
    tf_bystsage$group = factor(tf_bystsage$group, levels = ann_new$group[order.dendrogram(dend_m)] %>% unique)
    tf_bystsage = arrange(tf_bystsage,group, rank)
    tf_bystsage$group %>% table
    
    rowsplit = c()
    n = 1
    for (i in tf_bystsage$group %>% levels()){
      n = n + 1
      print(i)
      tmp = tf_bystsage[tf_bystsage$group == i,]
      rowsplit = c(rowsplit, rep(paste0(LETTERS[n],'_',i),tmp$regulon %>% length() ) )
    }
    
    saveRDS(tf_bystsage,'20240226_tf_bystsage.rds')
    
    col_fun = colorRamp2(c(-2,0, 2), c( "blue", "white","red3"))
    pdf('test_rss.pdf',width = 20*100, height = 80*100)
    print(ComplexHeatmap::Heatmap(regulonActivity_byCell_scaled[tf_bystsage$regulon, rownames(ann_new)], name="Regulon activity", cluster_rows = F,cluster_columns = dend_m ,top_annotation = colAnn, show_column_names = F, row_names_gp = grid::gpar(fontsize = 10,fontface = 'bold'), column_split = 2,clustering_method_rows = 'ward.D2', column_gap = unit(c(10), "mm"),  use_raster = F, col = col_fun,row_split = rowsplit))
    dev.off()
    
    
    ## Calculate ARI - regulon based ####
    print('ARI regulon')
    instep = "ARI regulon"
    tf_bystsage_prev = readRDS('../../../../20231225_tf_bystsage_h15.rds')
    tf_bystsage_prev = tf_bystsage_prev[-c(138,213,268),] 
    rownames(tf_bystsage_prev) = tf_bystsage_prev$regulon %>% word(sep = ' ', start = 1)
    
    rownames(tf_bystsage) = tf_bystsage$regulon %>% word(sep = ' ', start = 1)
    
    idx = intersect(rownames(tf_bystsage_prev), rownames(tf_bystsage))
    ari = adjustedRandIndex(tf_bystsage[idx,'group'], tf_bystsage_prev[idx,'group'])
    write.csv(ari,'ari_regulon.csv')
    
    tf_bystsage_inter = tf_bystsage[idx,]
    tf_bystsage_prev_inter = tf_bystsage_prev[idx,]
    
    tf_bystsage_inter$prev = tf_bystsage_prev_inter$group
    
    df = table(tf_bystsage_inter$group, tf_bystsage_inter$prev) %>% as.data.frame.array()
    df_normalized <- apply(df, 2, prop.table)
    
    test_new = table(tf_bystsage_inter$group, tf_bystsage_inter$prev)  %>% as.data.frame()
    test = as.table(df_normalized) %>% as.data.frame()
    colnames(test) = c('New','Prev','Freq')
    test$Num = test_new$Freq
    test = arrange(test,New)
    test$New = factor(test$New, levels = rev(levels(test$New)))
    
    
    
    pdf('heatmap_overlap_pct.pdf', width = 10,height = 10)
    print(ComplexHeatmap::Heatmap(df_normalized ,cell_fun = function(j, i, x, y, width, height, fill) {
      grid.text(sprintf("%.2f", df_normalized[i, j]), x, y, gp = gpar(fontsize = 10))
    }, cluster_rows = F, cluster_columns = F))
    dev.off()
    
    pdf('heatmap_overlap_num.pdf', width = 10,height = 10)
    print(ComplexHeatmap::Heatmap(df ,cell_fun = function(j, i, x, y, width, height, fill) {
      grid.text(sprintf("%1.f", df[i, j]), x, y, gp = gpar(fontsize = 10))
    }, cluster_rows = F, cluster_columns = F))
    dev.off()
    
    p = ggplot(test, aes(x = Prev , y = New , color = Freq, size = Num)) + geom_point()  + scale_size_continuous(breaks = seq(min(test$Num), max(test$Num), by = 2)) + theme_classic() + scale_colour_gradient2(low = lighten("lightgray",0.4),mid =lighten("lightgray",0.8),high = 'red3',midpoint = 0.1, name = 'Overlap Freq',guide=guide_colorbar(reverse=F) ) + theme(axis.text.x = element_text(angle = 90, hjust = 1))
    
    pdf('heatmap_overlap_dotplot.pdf', width = 9/1.1,height = 8/1.1)
    print(p)
    dev.off()
    
    
    
    setwd('../') }, error = function(e){
      write.csv(c(dd,instep ),'error.csv')
      setwd('/Users/jooyoung/Dropbox/samir/rds_jy/SCENIC/robust/patient_6/samir/')}
  )}



# Loop- Fisher regulon ####
dir = list.dirs( recursive = F) %>% str_remove_all('./')

for (dd in dir){
  print(dd)
  
  setwd(dd)
  tryCatch( {
    ## Fisher's exact test ####
    print('Fishers regulon')
    instep = "Fishers regulon"
    tf_bystsage_prev = readRDS('../../../../20231225_tf_bystsage_h15.rds')
    tf_bystsage_prev = tf_bystsage_prev[-c(138,213,268),] 
    rownames(tf_bystsage_prev) = tf_bystsage_prev$regulon %>% word(sep = ' ', start = 1)
    
    tf_bystsage = readRDS('20240226_tf_bystsage.rds')
    rownames(tf_bystsage) = tf_bystsage$regulon %>% word(sep = ' ', start = 1)
    
    idx = intersect(rownames(tf_bystsage_prev), rownames(tf_bystsage))
    
    tf_bystsage_inter = tf_bystsage[idx,]
    tf_bystsage_prev_inter = tf_bystsage_prev[idx,]
    tf_bystsage_inter$prev = tf_bystsage_prev_inter$group
    
    df = table(tf_bystsage_inter$group, tf_bystsage_inter$prev) %>% as.data.frame.array()
    df_normalized <- apply(df, 2, prop.table)
    
    df_fish = as.data.frame(matrix(ncol = 7)); colnames(df_fish ) = c('Pair','a','b','c','d','pval','TFs')
    
    for ( old in colnames(df)){
      print(old)
      for (new in rownames(df)){
        print(new)
        
        if (df_normalized[new, old] >= 0.5){
          
          a = df[new, old]
          b = sum(df[,old]) - a
          c = sum(df[new, ]) - a
          d = (dim(tf_bystsage_inter)[1]) - (a+b+c)
          
          
          fisher = data.frame(c(a,b),c( c ,d))
          test <- fisher.test(fisher) 
          
          if (test$p.value < 0.05){
            tmp = as.data.frame(matrix(ncol = 7)); colnames(tmp ) = c('Pair','a','b','c','d','pval','TFs')
            tmp$Pair = paste0(old,'_',new)
            tmp$a = a
            tmp$b = b
            tmp$c = c
            tmp$d = d
            tmp$pval = test$p.value
            tmp$TFs = tf_bystsage_inter[(tf_bystsage_inter$prev == old) & (tf_bystsage_inter$group == new),'regulon'] %>% word(sep = ' ', start = 1) %>% paste(collapse = ' | ')
            
            
            df_fish = rbind(df_fish, tmp)
          }
        }
      }
    }
    df_fish = df_fish[-1,]
    # write.csv(df_fish,'fisher_regulon.csv') # no TF info
    write.csv(df_fish,'fisher_regulon_tfs.csv')
    
    setwd('../') }, error = function(e){
      write.csv(c(dd,instep ),'error_tfs.csv')
      setwd('/Users/jooyoung/Dropbox/samir/rds_jy/SCENIC/robust/patient_6/samir/')}
  )}


# Loop- Fisher regulon Heatmap ####
lev = readRDS('../../../../regulon_lev.rds')
heat_df = as.data.frame(matrix( ncol =length(lev), nrow = 10 )); colnames(heat_df) = lev; rownames(heat_df) = seq(1,10)

dir = list.dirs( recursive = F) %>% str_remove_all('./')

for (dd in dir){
  print(dd)
  
  setwd(dd)
  tryCatch( {
    ## Fisher's exact test ####
    print('Fishers Heat')
    instep = "Fishers Heat"
    
    df_fish = read.csv('fisher_regulon.csv', row.names = 1)
    
    old = df_fish$Pair %>% word(sep = '_', start = 1) %>% unique
    
    for ( j in old){
      heat_df[dd,j] = 'Yes'
    }
    
    setwd('../') }, error = function(e){
      write.csv(c(dd,instep ),'error.csv')
      setwd('/Users/jooyoung/Dropbox/samir/rds_jy/SCENIC/robust/patient_10/samir/')}
  )}

my_color = 'red3'
names(my_color) = 'Yes'

pdf("20240301_recurrent_grn_heat.pdf",  width=11, height=9)
ComplexHeatmap::Heatmap(t(heat_df),  cluster_rows = F,  use_raster = F, column_title=NULL,row_title=NULL, row_names_side = 'left',column_names_side = 'top', show_heatmap_legend = F, col = my_color , na_col =  '#F5F5F5',column_names_rot = 0 )
dev.off()

