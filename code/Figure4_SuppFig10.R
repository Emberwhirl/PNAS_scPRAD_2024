library(Seurat)
library(dplyr)
library(ggplot2)
library(ggpubr)
library(ggrepel)
library(cowplot)
library(patchwork)
library(stringr)
library(grid)
library(SummarizedExperiment)
library(SeuratWrappers)
library(ggplotify)
library(ggpubr)
library(ggalluvial)
library(randomcoloR)
library(Rphenograph)
library(rstatix)
library(colorspace)
library(ggh4x)
library(scCustomize)
library(colorspace)
library(purrr)
library(Matrix)
library(DoubletFinder)
library(ggrepel)
library(readxl)
library(ggpubr)
library(RColorBrewer)
library(ComplexHeatmap)
set.seed(1234)
library(DescTools)


# Exclude CSPC ####
msk = readRDS('msk.integrated.remove.cellcycle.tumor.cells.rds')
msk  = subset(msk, subset = subtype == 'CSPC', invert = T )
saveRDS(msk, 'msk.integrated.remove.cellcycle.tumor.fastmnn.subtype.nohepato.figure4.rds')

# Figure 4A ####
cancer = readRDS('/Users/jooyoung/Dropbox/samir/rds_jy/msk.integrated.remove.cellcycle.tumor.fastmnn.subtype.nohepato.figure4.rds')

geneset = readxl::read_xlsx('/Users/jooyoung/OneDrive - 고려대학교/samir/rds/Table S18.JAK_STAT_FGFR_Misc_Signatures.xlsx')
AR = geneset$AR [!geneset$AR %>% is.na()]
AR = AR[AR %in% rownames(cancer)]
AR = list(AR)
cancer <- AddModuleScore(
  object = cancer,
  features  = AR,
  name = "AR_modulescore", seed = 1234
)

cancer$STEAP1_exp = (cancer@assays$RNA@data['STEAP1',] - min(cancer@assays$RNA@data['STEAP1',])) / (max(cancer@assays$RNA@data['STEAP1',]) - min(cancer@assays$RNA@data['STEAP1',]))
cancer$STEAP2_exp = (cancer@assays$RNA@data['STEAP2',] - min(cancer@assays$RNA@data['STEAP2',])) / (max(cancer@assays$RNA@data['STEAP2',]) - min(cancer@assays$RNA@data['STEAP2',]))
cancer$TROP2_exp = (cancer@assays$RNA@data['TACSTD2',] - min(cancer@assays$RNA@data['TACSTD2',])) / (max(cancer@assays$RNA@data['TACSTD2',]) - min(cancer@assays$RNA@data['TACSTD2',]))
cancer$EZH2_exp = (cancer@assays$RNA@data['EZH2',] - min(cancer@assays$RNA@data['EZH2',])) / (max(cancer@assays$RNA@data['EZH2',]) - min(cancer@assays$RNA@data['EZH2',]))
cancer$EPCAM_exp = (cancer@assays$RNA@data['EPCAM',] - min(cancer@assays$RNA@data['EPCAM',])) / (max(cancer@assays$RNA@data['EPCAM',]) - min(cancer@assays$RNA@data['EPCAM',]))
cancer$YAP1_exp = (cancer@assays$RNA@data['YAP1',] - min(cancer@assays$RNA@data['YAP1',])) / (max(cancer@assays$RNA@data['YAP1',]) - min(cancer@assays$RNA@data['YAP1',]))
cancer$CEACAM5_exp = (cancer@assays$RNA@data['CEACAM5',] - min(cancer@assays$RNA@data['CEACAM5',])) / (max(cancer@assays$RNA@data['CEACAM5',]) - min(cancer@assays$RNA@data['CEACAM5',]))
cancer$PSMA_exp = (cancer@assays$RNA@data['FOLH1',] - min(cancer@assays$RNA@data['FOLH1',])) / (max(cancer@assays$RNA@data['FOLH1',]) - min(cancer@assays$RNA@data['FOLH1',]))
cancer$DLL3_exp = (cancer@assays$RNA@data['DLL3',] - min(cancer@assays$RNA@data['DLL3',])) / (max(cancer@assays$RNA@data['DLL3',]) - min(cancer@assays$RNA@data['DLL3',]))
cancer$ASCL1_exp = (cancer@assays$RNA@data['ASCL1',] - min(cancer@assays$RNA@data['ASCL1',])) / (max(cancer@assays$RNA@data['ASCL1',]) - min(cancer@assays$RNA@data['ASCL1',]))
saveRDS(cancer,'/Users/jooyoung/OneDrive - 고려대학교/samir/rds_jy/msk.integrated.remove.cellcycle.tumor.fastmnn.subtype.nohepato.figure4.rds')


obj = SplitObject(cancer, split.by = "Regulon_h15_modi_ano") 
df = matrix(1, ncol = 15) %>% as.data.frame()
colnames(df) =c('Subtype','AR_modulescore','Patient','Cluster','STEAP1_Zscore','STEAP2_Zscore','TROP2_Zscore','EZH2_Zscore','EPCAM_Zscore','YAP1_Zscore','CEACAM5_Zscore','PSMA_Zscore','DLL3_Zscore','ASCL1_Zscore','Label')

for (n in 1:length(names(obj))){
  sample = names(obj)[n]
  print(sample)
  
  ob = obj[[n]]
  steap1 = ob$STEAP1_exp
  steap1_avg = mean(steap1)
  
  steap2 = ob$STEAP2_exp
  steap2_avg = mean(steap2)
  
  trop2 = ob$TROP2_exp
  trop2_avg = mean(trop2)
  
  exh2 = ob$EZH2_exp
  exh2_avg = mean(exh2)
  
  ar = ob$AR_modulescore1
  ar_avg = mean(ar)
  
  epcam = ob$EPCAM_exp
  epcam_avg = mean(epcam)
  
  ceacam = ob$CEACAM5_exp
  ceacam_avg = mean(ceacam)
  
  yap1 = ob$YAP1_exp
  yap1_avg = mean(yap1)
  
  psma = ob$PSMA_exp
  psma_avg = mean(psma)
  
  dll3 = ob$DLL3_exp
  dll3_avg = mean(dll3)
  
  ascl1 = ob$ASCL1_exp
  ascl1_avg = mean(ascl1)
  
  cl = ob$Regulon_h15_modi_ano %>% unique() %>% as.character()
  
  type = names(table(ob$subtype) )[table(ob$subtype)  == (table(ob$subtype) %>% max() )]
  sample = ob$patient %>% unique() %>% paste(collapse = ',')
  
  tmp = matrix(c(as.character(type),as.numeric(ar_avg),sample,cl,as.numeric(steap1_avg),as.numeric(steap2_avg),as.numeric(trop2_avg),as.numeric(exh2_avg),as.numeric(epcam_avg),as.numeric(yap1_avg),as.numeric(ceacam_avg),as.numeric(psma_avg),as.numeric(dll3_avg),as.numeric(ascl1_avg),""),ncol=15) %>% as.data.frame()
  colnames(tmp) =c('Subtype','AR_modulescore','Patient','Cluster','STEAP1_Zscore','STEAP2_Zscore','TROP2_Zscore','EZH2_Zscore','EPCAM_Zscore','YAP1_Zscore','CEACAM5_Zscore','PSMA_Zscore','DLL3_Zscore','ASCL1_Zscore','Label')
  
  df = rbind(df,tmp)
}

df = df[-1,]

df$AR_modulescore = as.numeric(df$AR_modulescore)
df$STEAP1_Zscore = as.numeric(df$STEAP1_Zscore)
df$STEAP2_Zscore = as.numeric(df$STEAP2_Zscore)
df$TROP2_Zscore = as.numeric(df$TROP2_Zscore)
df$EZH2_Zscore = as.numeric(df$EZH2_Zscore)
df$EPCAM_Zscore = as.numeric(df$EPCAM_Zscore)
df$YAP1_Zscore = as.numeric(df$YAP1_Zscore)
df$CEACAM5_Zscore = as.numeric(df$CEACAM5_Zscore)
df$PSMA_Zscore = as.numeric(df$PSMA_Zscore)
df$Subtype = factor(df$Subtype, levels = c('CSPC','CRPC','NEPC'))
df$DLL3_Zscore = as.numeric(df$DLL3_Zscore)
df$ASCL1_Zscore = as.numeric(df$ASCL1_Zscore)
df$Cluster[df$Cluster ==  "TCFL2+ WNT"] = 'TCF7L2+ WNT'
df$Cluster = factor(df$Cluster, levels = c('TCF7L2+ WNT','MAFG+','IRF2+ Inflammatory','SOX2/4+ Embryonic EMT','FOSL1+ AP-1','AR+ HOXB13+','AR+ HOXB13-','AR+ IRF+ Inflammatory','AR+ GI', 'AR+ HOXB13+ FOXA1+', 'NEPC-N', 'NEPC-A', 'NEPC-A/SOX6'))

col = readRDS('../regulon_col.rds')
ggplot(df, aes(AR_modulescore,PSMA_Zscore)) + geom_point(size= 2.5 ,aes(col = Cluster)) + theme_classic()  + scale_color_manual(values = col)


# Figure 4D ####
for (x in c('AR_modulescore')){ # 
  for (y in c('CEACAM5_Zscore','STEAP1_Zscore','STEAP2_Zscore','TROP2_Zscore','EZH2_Zscore','PSMA_Zscore')){
    print(paste0(x,' and ',y))
    
    p =ggplot(df, aes(df[,x],df[,y])) + geom_point(size= 2.5, aes( col = Cluster) )  + scale_color_manual(values = col) + theme_classic()  + xlab(x) + ylab(paste0('Scaled ',y %>% word(start = 1, sep='_'),' expression') )  + geom_smooth(aes(x = filter(df, Subtype == "CRPC")[,x], y =filter(df, Subtype == "CRPC")[,y]),color = '#2171B5',method = "lm",data = filter(df, Subtype == "CRPC"),size = 0.5,alpha  = 0.1, fill = alpha(darken( '#2171B5',0.1)))  + geom_smooth(aes(x = filter(df, Subtype == "NEPC")[,x], y =filter(df, Subtype == "NEPC")[,y]),color = '#CB181D',method = "lm",data = filter(df, Subtype == "NEPC"),size = 0,alpha  = 0.1, fill = alpha(darken( '#CB181D',0.1))) + stat_cor(label.x.npc = 0.65, label.y = max(df[,y])+ 0.9, aes(x = filter(df, Subtype == "NEPC")[,x], y =filter(df, Subtype == "NEPC")[,y]), color = '#CB181D', data = filter(df, Subtype == "NEPC") ) + stat_cor(label.x.npc = 0.65, label.y = max(df[,y]) + 0.8, aes(x = filter(df, Subtype == "CRPC")[,x], y =filter(df, Subtype == "CRPC")[,y]), color = '#2171B5', data = filter(df, Subtype == "CRPC")  )+ guides(col=guide_legend(title="GRN Annotation")) 
    pdf(paste0('/Users/jooyoung/Google Drive/My Drive/Main_Figure_v3/20231115_updated/20231115_Figure4D_',x,'_',y,'_updated_cor.pdf'), width = 5.5*1.2, height = 3.9*1.2)
    print(p)
    dev.off()
    
    #### For CEACAM5 #####
    p = ggplot(df, aes(df[,x],df[,y])) + geom_point(size= 2.5, aes( col = Cluster) )  + scale_color_manual(values = col) + theme_classic()  + xlab(x) + ylab(paste0('Scaled ',y %>% word(start = 1, sep='_'),' expression') )  + geom_smooth(aes(x = filter(df, Subtype == "CRPC")[,x], y =filter(df, Subtype == "CRPC")[,y]),color = '#2171B5',method = "lm",data = filter(df, Subtype == "CRPC"),size = 0.5,alpha  = 0.1, fill = alpha(darken( '#2171B5',0.1)))  + geom_smooth(aes(x = filter(df, Subtype == "NEPC")[,x], y =filter(df, Subtype == "NEPC")[,y]),color = '#CB181D',method = "lm",data = filter(df, Subtype == "NEPC"),size = 0.5,alpha  = 0.1, fill = alpha(darken( '#CB181D',0.1)),level=0.70) + stat_cor(label.x.npc = 0.65, label.y = max(df[,y])+ 0.01, aes(x = filter(df, Subtype == "NEPC")[,x], y =filter(df, Subtype == "NEPC")[,y]), color = '#CB181D', data = filter(df, Subtype == "NEPC") ) + stat_cor(label.x.npc = 0.65, label.y = max(df[,y]) , aes(x = filter(df, Subtype == "CRPC")[,x], y =filter(df, Subtype == "CRPC")[,y]), color = '#2171B5', data = filter(df, Subtype == "CRPC")  )+ guides(col=guide_legend(title="GRN Annotation")) + ylim(c(-0.05,0.08)) # for CEACAM5
    pdf(paste0('/Users/jooyoung/Google Drive/My Drive/Main_Figure_v3/20231223_Figure4_AR_CEAcAM5_cor.pdf'), width = 5.5*1.2, height = 3.9*1.2)
    print(p)
    dev.off()
    
    #### For TROP2 #####
    p = ggplot(df, aes(df[,x],df[,y])) + geom_point(size= 2.5, aes( col = Cluster) )  + scale_color_manual(values = col) + theme_classic()  + xlab(x) + ylab(paste0('Scaled ',y %>% word(start = 1, sep='_'),' expression') )  + geom_smooth(aes(x = filter(df, Subtype == "CRPC")[,x], y =filter(df, Subtype == "CRPC")[,y]),color = '#2171B5',method = "lm",data = filter(df, Subtype == "CRPC"),size = 0.5,alpha  = 0.1, fill = alpha(darken( '#2171B5',0.1)))  + geom_smooth(aes(x = filter(df, Subtype == "NEPC")[,x], y =filter(df, Subtype == "NEPC")[,y]),color = '#CB181D',method = "lm",data = filter(df, Subtype == "NEPC"),size = 0.5,alpha  = 0.1, fill = alpha(darken( '#CB181D',0.1)),level=0.70) + stat_cor(label.x.npc = 0, label.y = max(df[,y])+ 0.02, aes(x = filter(df, Subtype == "NEPC")[,x], y =filter(df, Subtype == "NEPC")[,y]), color = '#CB181D', data = filter(df, Subtype == "NEPC") ) + stat_cor(label.x.npc = 0, label.y = max(df[,y]) , aes(x = filter(df, Subtype == "CRPC")[,x], y =filter(df, Subtype == "CRPC")[,y]), color = '#2171B5', data = filter(df, Subtype == "CRPC")  )+ guides(col=guide_legend(title="GRN Annotation"))  
    
    pdf(paste0('/Users/jooyoung/Google Drive/My Drive/PNAS_fin/Main/Main_updated/Figure4_AR_Trop2_cor.pdf'), width = 5.5*1.2, height = 3.9*1.2)
    print(p)
    dev.off()
    
  }
}

x = 'ASCL1_Zscore'
y = 'DLL3_Zscore'

pdf(paste0('/Users/jooyoung/Google Drive/My Drive/Main_Figure_v3/20231115_updated/20231115_Figure4D_ASCL1_DLL3_updated_cor.pdf'),width = 5.5*1.2, height = 3.9*1.2)
ggplot(df, aes(df[,x],df[,y])) + geom_point(size= 2.5, aes( col = Cluster) ) + scale_color_manual(values = col) + theme_classic()  + xlab(x) + ylab(paste0('Scaled ',y %>% word(start = 1, sep='_'),' expression') )  + geom_smooth(aes(x = filter(df, Subtype == "CRPC")[,x], y =filter(df, Subtype == "CRPC")[,y]),color = '#2171B5',method = "lm",data = filter(df, Subtype == "CRPC"),size = 0.5,alpha  = 0.1, fill = alpha(darken( '#2171B5',0.1)))  + geom_smooth(aes(x = filter(df, Subtype == "NEPC")[,x], y =filter(df, Subtype == "NEPC")[,y]),color = '#CB181D',method = "lm",data = filter(df, Subtype == "NEPC"),size = 0,alpha  = 0.1, fill = alpha(darken( '#CB181D',0.1))) + stat_cor(label.x.npc = 0, label.y = max(df[,y]) + 0.28, aes(x = filter(df, Subtype == "NEPC")[,x], y =filter(df, Subtype == "NEPC")[,y]), color = '#CB181D', data = filter(df, Subtype == "NEPC") ) + stat_cor(label.x.npc = 0, label.y = max(df[,y]) + 0.33, aes(x = filter(df, Subtype == "CRPC")[,x], y =filter(df, Subtype == "CRPC")[,y]), color = '#2171B5', data = filter(df, Subtype == "CRPC")  )+ guides(col=guide_legend(title="GRN Annotation")) 
dev.off()


# Figure 4B ####
### Selected top 10 ####
cancer = readRDS('msk.integrated.remove.cellcycle.tumor.fastmnn.subtype.nohepato.figure4.rds')
regulonActivity_byCell_scaled = readRDS('./SCENIC/regulonActivity_byCell_scaled_h15_rss_filt.rds')
regulonActivity_byCell_scaled %>% dim()
ann_row = readRDS('./SCENIC/ann_row_h15.rds')

cancer$PSMA_type = ifelse(cancer$Regulon_h15_modi_ano %in% c('AR+ HOXB13+ FOXA1+','AR+ IRF+ Inflammatory'), 'PSMA high', ifelse(cancer$Regulon_h15_modi_ano == 'AR+ HOXB13-','PSMA low', ifelse( cancer$Regulon_h15_modi_ano %in% c('NEPC-A/SOX6',                'NEPC-N' ,'NEPC-A' ), 'NEPC', ifelse(cancer$Regulon_h15_modi_ano %in% c('AR+ HOXB13+','AR+ GI'), 'Else','AR low'))))

ann_row$PSMA_type = ifelse(ann_row$GRN %in% c('AR+ HOXB13+ FOXA1+','AR+ IRF+ Inflammatory'), 'PSMA high', ifelse(ann_row$GRN == 'AR+ HOXB13-','PSMA low', ifelse( ann_row$GRN %in% c('NEPC-A/SOX6',                'NEPC-N' ,'NEPC-A' ), 'NEPC', ifelse(ann_row$GRN %in% c('AR+ HOXB13+','AR+ GI'), 'Else','AR low'))))

ann = cancer@meta.data %>% dplyr::select('barcode','PSMA_type')
ann = ann[ann$PSMA_type != 'Else',]

ann_row = ann_row[ann_row$GRN %in% unique(cancer$Regulon_h15_modi_ano ),]
ann_row$PSMA_type %>% unique
ann_row = ann_row[ann_row$PSMA_type != 'Else',]

regulonActivity_byCell_scaled = regulonActivity_byCell_scaled[,rownames(ann)]
tf = readRDS('./SCENIC/20231114_tf_rss_h15.rds')
tf$group = as.character(tf$group)
tf$group_ano = ifelse(tf$group == 'CRPC_10','TCFL2+ WNT', ifelse(tf$group == 'CRPC_1', 'IRF2+ Inflammatory', ifelse(tf$group == 'CRPC_2','SOX2/4+ Embryonic EMT', ifelse(tf$group == 'CRPC_4','MAFG+', ifelse(tf$group == 'CRPC_9','FOSL1+ AP-1', ifelse(tf$group == 'CRPC_3','AR+ HOXB13+', ifelse(tf$group == 'CRPC_7', 'AR+ HOXB13-',ifelse(tf$group == 'CRPC_8','AR+ IRF+ Inflammatory',ifelse(tf$group == 'CRPC_5','AR+ GI',ifelse(tf$group == 'CRPC_6','AR+ HOXB13+ FOXA1+', ifelse(tf$group == 'NEPC_6','NEPC-N',ifelse(tf$group %in% c('NEPC_1','NEPC_2'), 'NEPC-A/SOX6', ifelse(tf$group %in% c('NEPC_3','NEPC_4'), 'NEPC-A','CSPC')))))))))))))

tf$PSMA_type = ifelse(tf$group_ano %in% c('AR+ HOXB13+ FOXA1+','AR+ IRF+ Inflammatory'), 'PSMA high', ifelse(tf$group_ano == 'AR+ HOXB13-','PSMA low', ifelse( tf$group_ano %in% c('NEPC-A/SOX6',                'NEPC-N' ,'NEPC-A' ), 'NEPC', ifelse(tf$group_ano %in% c('AR+ HOXB13+','AR+ GI'), 'Else','AR low'))))
tf$PSMA_type %>% unique
tf= tf[tf$PSMA_type != 'Else',]

saveRDS(tf, '20231126_tf_rss_df_for_fig4c.rds')

tf = readRDS('20231126_tf_rss_df_for_fig4c.rds')
setwd('./SCENIC')
scenicOptions = readRDS('int/scenicOptions.Rds')
regulonAUC <- loadInt(scenicOptions, "aucell_regulonAUC")
regulonAUC = regulonAUC[tf$regulon, rownames(ann)]

rss <- calcRSS(AUC=getAUC(regulonAUC), cellAnnotation=ann[colnames(regulonAUC), "PSMA_type"])
rssNorm <- scale(rss) 
rssNorm[rssNorm < 0] 

tf_bystsage = as.data.frame(matrix(ncol = 4)); colnames(tf_bystsage) = c('regulon','rank','rss','group')
for (i in ann$PSMA_type %>% unique){
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
tf_bystsage = arrange(tf_bystsage,group, rank)
saveRDS(tf_bystsage,'../20231126_tf_rss_psma.rds')

tf = readRDS('./20231126_tf_rss_psma.rds')
tf_use = tf %>%  group_by(group) %>%  slice_head(n = 10) %>%ungroup()
tf_use = rbind(tf[c('SOX2 (384g)','SOX4 (379g)',"NEUROD1 (59g)","FOXA1 (77g)","IRF7 (50g)"  ,"IRF2 (94g)" ),], tf_use)
tf_use = tf_use[!tf_use$regulon %in% c("SREBF1_extended (78g)", "SOX13_extended (31g)",'SETBP1 (98g)'),]
tf_use$group = factor(tf_use$group, levels =c('PSMA high','PSMA low','AR low','NEPC'))
tf_use[tf_use$regulon == "IRF2 (94g)",'group'] = 'AR low'
tf_use = arrange(tf_use, group, desc(rss))


### Heatmap ####
regulonActivity_byCell_scaled = readRDS('./SCENIC/regulonActivity_byCell_scaled_h15_rss_filt.rds')
regulonActivity_byCell_scaled %>% dim()
cancer = readRDS('msk.integrated.remove.cellcycle.tumor.fastmnn.subtype.nohepato.rds')
cancer = subset(cancer, subtype == 'CSPC', invert= T)

cancer$PSMA_type = ifelse(cancer$Regulon_h15_modi_ano %in% c('AR+ HOXB13+ FOXA1+','AR+ IRF+ Inflammatory'), 'PSMA high', ifelse(cancer$Regulon_h15_modi_ano == 'AR+ HOXB13-','PSMA low', ifelse( cancer$Regulon_h15_modi_ano %in% c('NEPC-A/SOX6',                'NEPC-N' ,'NEPC-A' ), 'NEPC', ifelse(cancer$Regulon_h15_modi_ano %in% c('AR+ HOXB13+','AR+ GI'), 'Else','AR low'))))

ann = cancer@meta.data %>% dplyr::select('barcode','PSMA_type')
ann = ann[ann$PSMA_type != 'Else',]

count = regulonActivity_byCell_scaled[tf_use$regulon,]
count = as.data.frame(count)
which(!rownames(ann) %in% colnames(count))
count = count[,rownames(ann)]

meta = ann
meta$subtype = cancer@meta.data[meta$barcode, 'subtype']
meta$PSMA_type = factor(meta$PSMA_type, levels = c('PSMA high','PSMA low','AR low','NEPC'))
meta = arrange(meta, PSMA_type,subtype)
meta$subtype = NULL
meta$barcode = NULL
colnames(meta) = 'PSMA/AR Group'
count = count[,rownames(meta)]

colours <- list( 
  'PSMA/AR Group' = setNames(c('pink','green3','lightgray','orange1'), levels(meta$`PSMA/AR Group` ))) # 'Patient' = setNames(distinctColorPalette(meta$Patient %>% unique %>% length), sort(unique(meta$Patient ))), 
colAnn <- HeatmapAnnotation(df = meta,
                            which = 'col',
                            col = colours,
                            annotation_width = unit(c(1, 4), 'cm'),
                            gap = unit(1, 'mm'),show_annotation_name = FALSE,show_legend = F
) 


col_fun = colorRamp2(c(-2,0, 2), c(  "blue",'white',"red"))
meta$`PSMA/AR Group` %>% table
column_split = c(rep('1', 3992),rep('2', 5223      ),rep('3',8505     +11602  ))

tf_use$group %>% table
row_split = c(rep('A',12), rep('B',10), rep('C',11), rep('D',10))

rownames(count)[1:13]
count = count[c(1:12,c(29,30,21:23,28,24:27),31:40),]
col_fun = colorRamp2::colorRamp2(c(-2, 0, 2), c(  "blue",'white',"red3"))
hmap = ComplexHeatmap::Heatmap(count, cluster_rows = F,cluster_columns = F,  name="Gene Expression", top_annotation = colAnn, row_names_gp = grid::gpar(fontsize = 10) ,column_split = column_split,show_column_names = F, column_title = " ", row_title = " ", col = col_fun,  show_heatmap_legend = T, row_split = row_split)
draw(hmap, heatmap_legend_side="right", annotation_legend_side="right",background = "transparent")


pdf('20231208_Figure4C_regulon_heapmap_rss_updated.pdf', width = 5*1.7,  height =7*1.3)
draw(hmap, heatmap_legend_side="right", annotation_legend_side="right",background = "transparent")
dev.off()


# Supp Fig 10D ####
cancer =  readRDS('msk.integrated.remove.cellcycle.tumor.fastmnn.subtype.nohepato.figure4.rds')
count = cancer@assays$RNA@data %>% as.data.frame()

msk = readRDS('../msk.integrated.remove.cellcycle.tumor.cells.rds')

meta = cancer@meta.data
meta$Regulon_h15_modi_ano[meta$Regulon_h15_modi_ano == 'TCFL2+ WNT']  = 'TCF7L2+ WNT'
meta$Regulon_h15_modi_ano_factor = factor(meta$Regulon_h15_modi_ano, levels = levels(msk$Regulon_h15_modi_ano_factor))
meta = meta %>% dplyr::select('subtype','Regulon_h15_modi_ano_ar','barcode','Regulon_h15_modi_ano_factor')
meta$Regulon_h15_modi_ano_ar = factor(meta$Regulon_h15_modi_ano_ar , levels = c('AR-','AR+','NEPC'))
meta = arrange(meta, subtype,  Regulon_h15_modi_ano_factor,Regulon_h15_modi_ano_ar )
meta$barcode = NULL
meta$subtype = NULL

DEG = readRDS('../SCENIC/Human_SCENIC_DEG.rds') 
DEG = DEG[DEG$avg_log2FC > 0.4 & DEG$p_val_adj < 0.05,]
DEG$cluster %>% table

cellsurface = read_xlsx('../Cell_Surface_Markers.xlsx')
cellsurface = cellsurface$`Cell Surface Markers`
cellsurface = c(cellsurface, 'FOLH1','STEAP1','STEAP2','DLL3') %>% unique

idx = DEG$gene %in% cellsurface
DEG_cellsurf = DEG[idx,]

idx = DEG_cellsurf$gene[DEG_cellsurf$gene %>% duplicated()]
idxx = ! DEG_cellsurf$gene %in% idx
DEG_cellsurf = DEG_cellsurf[idxx,]

n = 1
column_split = c()
for (i in meta$Regulon_h15_modi_ano_factor %>% levels()){
  n = n + 1
  print(i)
  column_split = c(column_split, rep(paste0(LETTERS[n],'_',i),which(meta$Regulon_h15_modi_ano_factor == i) %>% length() ) )
}

col  = readRDS('../regulon_col.rds')
colours <- list( 'Regulon_h15_modi_ano_ar' = setNames(c('blue4','skyblue','red'), meta$Regulon_h15_modi_ano_ar %>% levels()), 'Regulon_h15_modi_ano_factor' = setNames(col, levels(meta$Regulon_h15_modi_ano_factor)))

colAnn <- HeatmapAnnotation(df = meta,
                            which = 'col',
                            col = colours,
                            annotation_width = unit(c(1, 4), 'cm'),
                            gap = unit(1, 'mm'),show_annotation_name = FALSE,show_legend = T
)

deg = DEG_cellsurf %>%  group_by(cluster) %>% slice_head(n = 10) %>% ungroup() %>% .[,'gene']
deg = deg$gene

count_top10 = count[deg,]
count_top10 = count_top10[,rownames(meta)]

row_split = rep(c('a','b','c'),10) %>% sort()
col_fun = colorRamp2::colorRamp2(c(0, 1), c( 'lightgray',"red3"))

png('./ARNEPC/20231227_top10_heatmap.png', width = 10*90,  height =13*90)
ComplexHeatmap::Heatmap(count_top10, cluster_rows = F,cluster_columns = F,  name="Gene Expression",top_annotation = colAnn,  row_names_gp = grid::gpar(fontsize = 10,fontface = 'bold') ,show_column_names = F, column_title = " ", row_title = " ", column_split = column_split, row_split = row_split,column_gap = unit(c(0.2), "mm"), col = col_fun)
dev.off()


# Figure 4C, Supp Fig 10B ####
patient = 'HMP14'

for (reg in patient){ 
  
  seura = subset(cancer, patient == reg)
  DefaultAssay(seura)
  seura = DietSeurat(seura)
  
  seura = NormalizeData(seura, assay = 'RNA')
  s.genes <- cc.genes$s.genes ; g2m.genes <- cc.genes$g2m.genes; seura <- CellCycleScoring(seura, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE); seura <- ScaleData(seura, features = rownames(seura),vars.to.regress = c("S.Score", "G2M.Score"),assay = 'RNA') 
  
  markers = c('FOLH1','DLL3','STEAP1','STEAP2','TACSTD2') # drug gene
  seura <- RunPCA(seura, features = markers,assay = 'RNA')  
  
  data <- as.matrix(seura@reductions$pca@cell.embeddings[,1:4]) 
  
  Rphenograph_out <- Rphenograph(data, k = 30) 
  pheno = membership(Rphenograph_out[[2]]) %>% as.factor()
  seura$pheno_cluster = pheno
  names(pheno) = rownames(seura@reductions$pca@cell.embeddings)
  seura$pheno_cluster = pheno
  seura$barcode <- colnames(seura)
  
  write.csv(seura@meta.data, file=paste0(reg,'/meta_drug.csv'), quote=F, row.names=F)
  write.csv(seura@reductions$pca@cell.embeddings, file=paste0(reg,'/pca_drug.csv'), quote=F, row.names=F) 
  saveRDS(seura,paste0(reg,'/seura_drug.rds'))
  
  umap = read.csv(paste0(h,'/umap_drug.csv'))
  rownames(umap) = umap$barcode; umap$barcode = NULL
  colnames(umap) = c('X1','X2')
  seura[['umap']] = CreateDimReducObject(embeddings = umap %>% as.matrix(), key = 'UMAP_', assay = 'RNA')
  saveRDS(seura,paste0(h,'/seura_drug_rpheno.rds'))
}



## Supp Fig 10B ####
hmp14 = readRDS('HMP14/seura_drug_rpheno.rds')
Idents(hmp14) = hmp14$pheno_cluster
regulon = readRDS('../SCENIC/regulonActivity_byCell_scaled_raw.rds')

regulon = regulon[,colnames(hmp14)]
rownames(regulon)[regulon %>% rownames() %>% str_detect('HOXB13')]

hmp14$HOXB13_regulon = regulon["HOXB13 (14g)",]
FeaturePlot(hmp14,'HOXB13_regulon', reduction = 'umap', label = T)

deg = c('AR', 'FOLH1',  'TACSTD2','DLL3') # ,'HOXB13','KLF15','KFIL3','IRF9','CEBPD','LEF1','SOX5','NFKB1','ZBTB16','NR3C1','HOXA6','IRX4','ZNF274','GATA2','HOXA13','AR','HES4','ELK4','RELB','ATF6'
count = hmp14@assays$RNA@data
deg = deg[deg %in% rownames(count)]
count = count[deg,]
count = as.data.frame(count)
count['HOXB13 (14g)',] = t(hmp14$HOXB13_regulon)

meta = hmp14@meta.data
meta = meta %>% dplyr::select('pheno_cluster')
meta = arrange(meta, pheno_cluster)
colnames(meta) = 'pheno_cluster'

count = count[,rownames(meta)]

meta$PSMA = ifelse(meta$pheno_cluster == 1,'PSMA high','PSMA low')


ddd = meta$pheno_cluster %>% table %>% as.data.frame.array()

column_split = c()
for (i in rownames(ddd)){
  if (as.numeric(i) < 10){
    column_split = c(column_split, rep(as.numeric(i),ddd[i,]))
  } else{
    n = as.numeric(i)+80
    column_split = c(column_split, rep(as.numeric(n),ddd[i,]))
  }
}
meta$pheno_cluster = NULL

col = DiscretePalette_scCustomize(num_colors = meta$pheno_cluster %>% unique %>% length(), palette = "polychrome")
col = c('#5A5156FF',colorspace::darken('lavender'))
# colours <- list( 'pheno_cluster' = setNames(col, meta$pheno_cluster %>% unique)) # 'Patient' = setNames(distinctColorPalette(meta$Patient %>% unique %>% length), sort(unique(meta$Patient ))), 
colours <- list( 'PSMA' = setNames(col, meta$PSMA %>% unique))

colAnn <- HeatmapAnnotation(df = meta,
                            which = 'col',
                            col = colours,
                            annotation_width = unit(c(1, 4), 'cm'),
                            gap = unit(1, 'mm'),show_annotation_name = FALSE,show_legend = T
)
col_fun = colorRamp2::colorRamp2(c(0, 1), c( 'lightgray',"red3"))

hmap = ComplexHeatmap::Heatmap(count[1:4,], cluster_rows = F,cluster_columns = F,  name="Gene Expression", col = col_fun,top_annotation = colAnn, row_names_gp = grid::gpar(fontsize = 10,fontface = 'bold') ,show_column_names = F, column_title = " ", row_title = " ", column_split = column_split,column_gap = unit(c(0.2), "mm")) 
hmapp = ComplexHeatmap::Heatmap(count[5,] %>% as.matrix(), cluster_rows = F,cluster_columns = F,  name="Regulon Activity",  row_names_gp = grid::gpar(fontsize = 10,fontface = 'bold') ,show_column_names = F, column_title = " ", row_title = " ", column_split = column_split,column_gap = unit(c(0.2), "mm")) 

ht_list = hmap %v% hmapp

pdf('/Users/jooyoung/Google Drive/My Drive/PNAS_fin/Main/Main_updated/Figure4_HMP14_heatmap.pdf', width = 7,  height =4.5)
# draw(hmap, heatmap_legend_side="right", annotation_legend_side="right",background ="transparent")
draw(ht_list)
dev.off()