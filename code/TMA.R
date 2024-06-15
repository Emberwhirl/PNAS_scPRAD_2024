set.seed(1234)
library(dplyr)
library(Seurat)
library(SCENIC)
library(AUCell)
library(circlize)
library(ComplexHeatmap)
library(stringr)
library(ggplot2)
library(dendextend)
library(msigdbr)
library(RColorBrewer)
library('stats')
library(tidyHeatmap)
library(colorspace)
library(randomcoloR)
library(colorspace)
library(dendextend)
library(cowplot)
library(scCustomize)
#library(export)
library(ggpubr)
library(readxl)


# Figure 1A ####
## Input prepare ####
excel = read_xlsx('Supplementary_table1.xlsx', sheet = 1) %>% as.data.frame()
df <- excel[,colSums(is.na(excel))<nrow(excel)]
df = distinct(df)
rownames(df) = paste0('random_',c(1:131))

meta = df %>% dplyr::select('Anatomic Site - large','Morphology','Identifier') # Morphology == Histology
counts = df[,c(4:21)]; rownames(counts) = rownames(df)
counts_nomki = counts[,-12]
counts_mki = counts[,12] %>% as.data.frame(); rownames(counts_mki) = rownames(df)
counts_nomki = t(counts_nomki)
counts = t(counts)

colnames(meta)[1] = 'Site'
meta$Site = factor(meta$Site, levels = c('Bone',"Liver/Lung" , "Prostate" ,   "Lymph Node",    "Other viscera"))
colnames(meta)[2] = 'Histology'
meta$Histology %>% unique %>% rev()
colnames(meta)[3] = 'Patient'

matrix = matrix(0, nrow = 2, ncol = dim(counts_nomki)[2]) %>% as.data.frame()

counts_nomki = as.data.frame(counts_nomki)
counts_nomki = counts_nomki[c( "AR"   ,   "NKX3.1"  ,"CK8"   ,  "TP63" ,"SYP","INSM1",'ASCL1','NEUROD1','FOXA2',"YAP1",'POU2F3','MYC','SOX2','EZH2','TFF3','TROP2','DLL3'),]


## Make dendrogram for each patient ####
#### p1 ####
m_p1 = counts[,  rownames(meta)[meta$Patient == "Patient 1"], drop = FALSE]
hc1 = hclust(dist(t(m_p1)), method = "ward.D2")
dend_p1 = as.dendrogram(hc1)
order_p1 = m_p1[,order.dendrogram(dend_p1)] %>% colnames()
saveRDS(dend_p1,'20231206_dend_p1_new.rds')

#### p2 ####
m_p2 = counts[,  rownames(meta)[meta$Patient == "Patient 2"], drop = FALSE]
hc2 = hclust(dist(t(m_p2)), method = "ward.D2")
dend_p2 = as.dendrogram(hc2)
order_p2 = m_p2[,order.dendrogram(dend_p2)] %>% colnames()
saveRDS(dend_p2,'20231206_dend_p2_new.rds')

#### p3 ####
m_p3 = counts[,  rownames(meta)[meta$Patient == "Patient 3"], drop = FALSE]
hc3 = hclust(dist(t(m_p3)), method = "ward.D2")
dend_p3 = as.dendrogram(hc3)
order_p3 = m_p3[,order.dendrogram(dend_p3)] %>% colnames()
saveRDS(dend_p3,'20231206_dend_p3_new.rds')

#### p4 ####
m_p4 = counts[,  rownames(meta)[meta$Patient == "Patient 4"], drop = FALSE]
hc4 = hclust(dist(t(m_p4)), method = "ward.D2")
dend_p4 = as.dendrogram(hc4)
order_p4 = m_p4[,order.dendrogram(dend_p4)] %>% colnames()
saveRDS(dend_p4,'20231206_dend_p4_new.rds')

#### p5 ####
m_p5 = counts[,  rownames(meta)[meta$Patient == "Patient 5"], drop = FALSE]
hc5 = hclust(dist(t(m_p5)), method = "ward.D2")
dend_p5 = as.dendrogram(hc5)
order_p5 = m_p5[,order.dendrogram(dend_p5)] %>% colnames()
saveRDS(dend_p5,'20231206_dend_p5_new.rds')

#### p6 ####
m_p6 = counts[,  rownames(meta)[meta$Patient == "Patient 6"], drop = FALSE]
hc6 = hclust(dist(t(m_p6)), method = "ward.D2")
dend_p6 = as.dendrogram(hc6)
order_p6 = m_p6[,order.dendrogram(dend_p6)] %>% colnames()
saveRDS(dend_p6,'20231206_dend_p6_new.rds')

#### p7 ####
m_p7 = counts[,  rownames(meta)[meta$Patient == "Patient 7"], drop = FALSE]
hc7 = hclust(dist(t(m_p7)), method = "ward.D2")
dend_p7 = as.dendrogram(hc7)
order_p7 = m_p7[,order.dendrogram(dend_p7)] %>% colnames()
saveRDS(dend_p6,'20231206_dend_p7_new.rds')

#### p8 ####
m_p8 = counts[,  rownames(meta)[meta$Patient == "Patient 8"], drop = FALSE]
hc8 = hclust(dist(t(m_p8)), method = "ward.D2")
dend_p8 = as.dendrogram(hc8)
order_p8 = m_p8[,order.dendrogram(dend_p8)] %>% colnames()
saveRDS(dend_p8,'20231206_dend_p8_new.rds')

#### p9 ####
m_p9 = counts[,  rownames(meta)[meta$Patient == "Patient 9"], drop = FALSE]
hc9 = hclust(dist(t(m_p9)), method = "ward.D2")
dend_p9 = as.dendrogram(hc9)
order_p9 = m_p9[,order.dendrogram(dend_p9)] %>% colnames()
saveRDS(dend_p9,'20231206_dend_p9_new.rds')

#### p10 ####
m_p10 = counts[,  rownames(meta)[meta$Patient == "Patient 10"], drop = FALSE]
hc10 = hclust(dist(t(m_p10)), method = "ward.D2")
dend_p10 = as.dendrogram(hc10)
order_p10 = m_p10[,order.dendrogram(dend_p10)] %>% colnames()
saveRDS(dend_p10,'20231206_dend_p10_new.rds')

#### p11 ####
m_p11 = counts[,  rownames(meta)[meta$Patient == "Patient 11"], drop = FALSE]
hc11 = hclust(dist(t(m_p11)), method = "ward.D2")
dend_p11 = as.dendrogram(hc11)
order_p11 = m_p11[,order.dendrogram(dend_p11)] %>% colnames()
saveRDS(dend_p11,'20231206_dend_p11_new.rds')

#### p12 ####
m_p12 = counts[,  rownames(meta)[meta$Patient == "Patient 12"], drop = FALSE]
hc12 = hclust(dist(t(m_p12)), method = "ward.D2")
dend_p12 = as.dendrogram(hc12)
order_p12 = m_p12[,order.dendrogram(dend_p12)] %>% colnames()
saveRDS(dend_p12,'20231206_dend_p12_new.rds')

#### p13 ####
m_p13 = counts[,  rownames(meta)[meta$Patient == "Patient 13"], drop = FALSE]
hc13 = hclust(dist(t(m_p13)), method = "ward.D2")
dend_p13 = as.dendrogram(hc13)
order_p13 = m_p13[,order.dendrogram(dend_p13)] %>% colnames()
saveRDS(dend_p13,'20231206_dend_p13_new.rds')

#### p14 ####
m_p14 = counts[,  rownames(meta)[meta$Patient == "Patient 14"], drop = FALSE]
hc14 = hclust(dist(t(m_p14)), method = "ward.D2")
dend_p14 = as.dendrogram(hc14)
order_p14 = m_p14[,order.dendrogram(dend_p14)] %>% colnames()
saveRDS(dend_p14,'20231206_dend_p14_new.rds')

#### p15 ####
m_p15 = counts[,  rownames(meta)[meta$Patient == "Patient 15"], drop = FALSE]
hc15 = hclust(dist(t(m_p15)), method = "ward.D2")
dend_p15 = as.dendrogram(hc15)
order_p15 = m_p15[,order.dendrogram(dend_p15)] %>% colnames()
saveRDS(dend_p15,'20231206_dend_p15_new.rds')

#### p16 ####
m_p16 = counts[,  rownames(meta)[meta$Patient == "Patient 16"], drop = FALSE]
hc16 = hclust(dist(t(m_p16)), method = "ward.D2")
dend_p16 = as.dendrogram(hc16)
order_p16 = m_p16[,order.dendrogram(dend_p16)] %>% colnames()
saveRDS(dend_p16,'20231206_dend_p16_new.rds')


orders = ls()[ls() %>% str_detect('order_p')]
order_list = c()
for (i in orders){
  order_list= c(order_list, get(i))
}
meta = meta[order_list,]
m = counts
# meta = arrange(meta,Patient)
m = m[,rownames(meta)]
hc_patients = cluster_between_groups_ward(m , meta$Patient)
dend_test = as.dendrogram(hc_patients )
plot(dend_test )
meta[order.dendrogram(hc_patients),]
order = meta[m  %>% .[,order.dendrogram(hc_patients)] %>% colnames(),'Patient'] %>% unique



## Make dendrogram across patient ####
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

orders = ls()[ls() %>% str_detect('order_p')]
order_list = c()
for (i in orders){
  order_list= c(order_list, get(i))
}
meta = meta[order_list,]
m = counts
m = m[,rownames(meta)]
hc_patients = cluster_between_groups_ward(m , meta$Patient)
dend_test = as.dendrogram(hc_patients )
plot(dend_test )
meta[order.dendrogram(hc_patients),]
order = meta[m  %>% .[,order.dendrogram(hc_patients)] %>% colnames(),'Patient'] %>% unique # 7 5 6 16 1 10 12 13 14 9 2 11 15 8 3
saveRDS(dend_test,'20231206_cluster_between_pateints_dend_new.rds')

order_col = c(colnames(m_p7),colnames(m_p5), colnames(m_p6), colnames(m_p16),colnames(m_p1),colnames(m_p10),colnames(m_p12),colnames(m_p13),colnames(m_p14),colnames(m_p9),colnames(m_p2),colnames(m_p11),colnames(m_p15),colnames(m_p8),colnames(m_p3),colnames(m_p4))
dend_p = as.dendrogram(hclust(dist(rbind(colMeans(t(m_p7)), rbind(colMeans(t(m_p5)),rbind(colMeans(t(m_p6)), rbind(colMeans(t(m_p16))),rbind(colMeans(t(m_p1)), rbind(colMeans(t(m_p10)), rbind(colMeans(t(m_p12)), rbind(colMeans(t(m_p13)), rbind(colMeans(t(m_p14)), rbind(colMeans(t(m_p9)), rbind(colMeans(t(m_p2)), rbind(colMeans(t(m_p11)), rbind(colMeans(t(m_p15)), rbind(colMeans(t(m_p8)),rbind(colMeans(t(m_p3)), colMeans(t(m_p4)))))))))))))))))))
dend_m = merge_dendrogram(dend_p, list(dend_p7, dend_p5, dend_p6, dend_p16 ,dend_p1,dend_p10,dend_p12,dend_p13,dend_p14,dend_p9,dend_p2,dend_p11,dend_p15,dend_p8,dend_p3,dend_p4))

saveRDS(dend_m,'20231206_dend_patients_new.rds')
saveRDS(order_col ,'20231206_order_col_new.rds')


## Figure 1A heatmap ####
order_col = readRDS('20231206_order_col_new.rds')
dend_m = readRDS('20231206_dend_patients_new.rds')

meta = meta[order_col,]
counts_nomki = counts_nomki[,rownames(meta)]

col_fun = colorRamp2(c(0, 200), c(  "#F2F3F4","red3"))
colours <- list( 
  'Site' = setNames(c('gold','#D6CADD','#4B0082','#DF00FF',"#2E8B57"), levels(meta$Site )),
  'Histology' =setNames(c('#73C2FB','orange','red3'), c('PRAD','HGC','NEPC'))) 
meta$Patient = NULL

colAnn <- HeatmapAnnotation(df = meta,
                            which = 'col',
                            col = colours,
                            annotation_width = unit(c(1, 4), 'cm'),
                            gap = unit(1, 'mm'),show_annotation_name = FALSE)


hmap = ComplexHeatmap::Heatmap(counts_nomki, name="H-score", cluster_rows = F,cluster_columns = dend_m, row_names_gp = grid::gpar(fontsize = 15,fontface = 'bold'), column_split = 16 , na_col = "#91A3B0", col = col_fun, row_split = c(rep('1',4), rep('2',1),rep('3',1),rep('5',9), rep('6',2)),top_annotation = colAnn,row_gap = unit(c(5), "mm"), show_column_names = F, column_title = " ", row_title = " ")
draw(hmap)

counts_mki = counts_mki[rownames(meta),]
counts_mki = counts_mki %>% t() %>% as.data.frame() 
colnames(counts_mki) = rownames(meta)
rownames(counts_mki) = 'KI67'

column_split = c(rep('1',length(which(metaa$Patient == 'Patient 7'))), rep('2',length(which(metaa$Patient == 'Patient 5'))), rep('3',length(which(metaa$Patient == 'Patient 6'))), rep('4',length(which(metaa$Patient == 'Patient 16'))), rep('5',length(which(metaa$Patient == 'Patient 1'))), rep('51',length(which(metaa$Patient == 'Patient 10'))), rep('6',length(which(metaa$Patient == 'Patient 12'))), rep('7',length(which(metaa$Patient == 'Patient 13'))), rep('8',length(which(metaa$Patient == 'Patient 14'))),rep('9',length(which(metaa$Patient == 'Patient 9'))), rep('10',length(which(metaa$Patient == 'Patient 2'))), rep('11',length(which(metaa$Patient == 'Patient 11'))), rep('12',length(which(metaa$Patient == 'Patient 15'))), rep('13',length(which(metaa$Patient == 'Patient 8'))),rep('14',length(which(metaa$Patient == 'Patient 3'))), rep('15',length(which(metaa$Patient == 'Patient 4')))  )

col_fun = colorRamp2(c(0, 100), c(  "#F2F3F4","gray4"))
hmapp = ComplexHeatmap::Heatmap(counts_mki, name="KI67", cluster_rows = F,cluster_columns = F, row_names_gp = grid::gpar(fontsize = 15,fontface = 'bold'), column_split = column_split  , na_col = "#91A3B0", col = col_fun,row_gap = unit(c(5), "mm"), show_column_names = F, column_title = " ", row_title = " ")

ht_list = hmap %v% hmapp

pdf('Figure 1A.pdf', width = 17, height =6)
draw(ht_list)
dev.off()



# Figure 1B-I and SuppFig1 ####
## Dotplot ####
excel = read_xlsx('Supplementary_table1.xlsx', sheet = 1) %>% as.data.frame()
df <- excel[,colSums(is.na(excel))<nrow(excel)]
df = distinct(df)
rownames(df) = paste0('random_',c(1:131))

meta = df %>% dplyr::select('Anatomic Site - large','Morphology','Identifier') # Morphology == Histology
counts = df[,c(4:21)]; rownames(counts) = rownames(df)
counts_nomki = counts[,-12]
counts_nomki = counts_nomki[,c('ASCL1','DLL3','YAP1','INSM1','TROP2','SYP','AR')]

colnames(meta)[1] = 'Site'
meta$Site = factor(meta$Site, levels = c('Bone',"Liver/Lung" , "Prostate" ,   "Lymph Node",    "Other viscera"))
colnames(meta)[2] = 'Histology'
meta$Histology %>% unique %>% rev()
colnames(meta)[3] = 'Patient'

dff = cbind(counts_nomki, meta)
dff$NEPC_score_mean = rowMeans(dff[c("ASCL1", "INSM1","FOXA2")], na.rm = TRUE)


pdf('dotplot.pdf')
ggplot(dff, aes(x = NEPC_score_mean,y=DLL3)) + geom_point(size= 3.5 , aes( col = Histology, fill =Histology))  + scale_color_manual(values = c('#73C2FB','orange','red3'))  + theme_classic()  + guides(col=guide_legend(title="Histology")) + xlab('ASCL1 H-score')+ ylab('DLL3 H-score')+  scale_fill_manual(values = c('#73C2FB','orange','red3')) + theme(text = element_text(family = 'Helvetica'))
dev.off()

### CCC ####
ccc = CCC(counts_nomki$NEPC_score_mean, counts_nomki$DLL3)


## Boxplot ####
for (y in colnames(df)[4:21] ){
  if ( !  y %in% c('TP63','NEUROD1','POU2F3')){
    formula = as.formula( paste(y, 'Morphology', sep="~") )
    
    ylab = y
    if (y == 'Ki67'){
      ylab = 'Ki-67'
    }
    
    p =ggplot(df, aes(x= Morphology, y = df[,y], col =Morphology )) + geom_boxplot() + geom_jitter()  + scale_color_manual(values = c('#73C2FB','orange','red3'))  + theme_classic()  + guides(col=guide_legend(title="Histology")) + theme(axis.title.x = element_blank()) + ylab(paste0(ylab, ' H-Score'))
    
    pdf(paste0('/Users/jooyoung/Google Drive/My Drive/PNAS_fin/20240131_rev_v1/SF1/SF1_jitter_',y,'.pdf'), width = 4, height = 4)
    print(p)
    dev.off()
  }
}

### zero-inflated modified wilcox test ####
library(rstatix)
library(ZIR)

for (y in colnames(df)[4:21] ){
  if ( !  y %in% c('TP63','NEUROD1','POU2F3')){
    formula = as.formula( paste(y, 'Morphology', sep="~") )
    
    tmp = df[,c(y,'Morphology')]
    
    PRAD = tmp[tmp$Morphology == "PRAD",y][!tmp[tmp$Morphology == "PRAD",y] %>% is.na()]
    HGC = tmp[tmp$Morphology == "HGC",y][!tmp[tmp$Morphology == "HGC",y] %>% is.na()]
    NEPC = tmp[tmp$Morphology == "NEPC",y][!tmp[tmp$Morphology == "NEPC",y] %>% is.na()]
    
    stat.test = data.frame(matrix(ncol = 13, nrow = 3)); colnames(stat.test) = c(".y." ,"group1","group2","n1","n2","statistic","p","p.adj","p.adj.signif","y.position","groups","xmin","xmax")
    
    stat.test[,1] = y
    stat.test$group1 = c('PRAD','PRAD','HGC')
    stat.test$group2 = c('HGC','NEPC','NEPC')
    
    stat.test$n1[1] = length(PRAD)
    stat.test$n1[2] = length(PRAD)
    stat.test$n1[3] = length(HGC)
    stat.test$n2[1] = length(HGC)
    stat.test$n2[2] = stat.test$n2[3] = length(NEPC)
    
    stat.test$statistic[1] = ziw(PRAD, HGC, perm = T)$statistics
    stat.test$statistic[2] = ziw(PRAD, NEPC, perm = T)$statistics
    stat.test$statistic[3] = ziw(HGC, NEPC, perm = T)$statistics
    
    stat.test$p[1] = ziw(PRAD, HGC, perm = T)$p.value
    stat.test$p[2] = ziw(PRAD, NEPC, perm = T)$p.value
    stat.test$p[3] = ziw(HGC, NEPC, perm = T)$p.value
    
    stat.test$p.adj[1] = p.adjust(ziw(PRAD, HGC, perm = T)$p.value, method = "bonferroni")
    stat.test$p.adj[2] = p.adjust(ziw(PRAD, NEPC, perm = T)$p.value, method = "bonferroni")
    stat.test$p.adj[3] = p.adjust(ziw(HGC, NEPC, perm = T)$p.value, method = "bonferroni")
    
    stat.test = stat.test %>% add_significance() 
    maxval = df[,y][! df[,y] %>% is.na()] %>% max()
    
    stat.test$manual_position <- maxval * c(1.05,1.2,1.11)
    stat.test$label <- stat.test$p.adj.signif
    
    ylab = y
    if (y == 'Ki67'){
      ylab = 'Ki-67'
    }
    
    p = ggplot(df, aes(x= Morphology, y = df[,y], col =Morphology )) + geom_boxplot()   + scale_color_manual(values = c('#73C2FB','orange','red3'))  + theme_classic()  + guides(col=guide_legend(title="Histology")) + theme(axis.title.x = element_blank()) +  ggsignif::geom_signif(data=stat.test,aes(xmin=group1,xmax=group2,annotations=label,y_position=manual_position),manual=TRUE, inherit.aes=FALSE) + ylab(paste0(y, ' H-Score'))
    
    pdf(paste0('/Users/jooyoung/Google Drive/My Drive/PNAS_fin/20240131_rev_v1/Supp_Figure1_ziw_',y,'.pdf'), width = 4, height = 4)
    print(p)
    dev.off()
  }
}


