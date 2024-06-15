library(Seurat)
library(dplyr)
library(stringr)
library(cowplot)
library(SeuratDisk)
library(reticulate)
library(ggplot2)
library(readxl)
library(patchwork)
library(dplyr)
library(stringr)
library(showtext)
library(biomaRt)
library(anndata)
library(Matrix)
library(SoupX)
library(MAST)
library(colorspace)
library(scCustomize)
library(fgsea)
library(clusterProfiler)
library(readxl)
library(rio)
set.seed(1234)

# GEMM GRN with Human SCENIC Geneset ####
mus = readRDS('/Users/jooyoung/Dropbox/samir/rds/mus_cancer.rds')
mus$GRN_pheno_tune [mus$GRN_pheno_tune == 'Ascl1+'] = 'NEPC-A'
Idents(mus) = mus$GRN_pheno_tune

DEG = FindAllMarkers(mus, only.pos = F, logfc.threshold = 0, min.pct = 0.3, test.use = 'MAST')
saveRDS(DEG,'~/Documents/samir_macbook/GEMM/samir_mouse/DEG_GEMM_pheno_tune.rds')

## ORA ####
DEG = readRDS('~/Documents/samir_macbook/GEMM/samir_mouse/DEG_GEMM_pheno_tune.rds')
df = readRDS('/Users/jooyoung/Dropbox/samir/rds_jy/SCENIC/Human_h15_human_rss_top5.rds')
df = df[,c(2,1)]
colnames(df) = c("gs_name", "human_gene_symbol")

for (ii in DEG$cluster %>% unique ){
  
  DEG_filt = DEG[DEG$cluster == ii,]
  DEG_filt = arrange(DEG_filt,desc(avg_log2FC))
  DEG_filt = DEG_filt[DEG_filt$p_val_adj < 0.05,]
  DEG_filt = DEG_filt[DEG_filt$avg_log2FC > 0.25,]
  genes = DEG_filt$gene
  
  i =  ii  %>% str_remove('[/]')
  print(i)
  
  dir.create('./ORA/H15_top5', showWarnings = F)
  go_gene_sets = df
  em <- enricher(genes, TERM2GENE=go_gene_sets, pvalueCutoff =1, qvalueCutoff = 1, minGSSize = 1, maxGSSize = 5000)
  em = em@result
  
  for (j in 1:length(em$Description)){
    em$`Obs-Exp`[j] = ((em$GeneRatio[j] %>% str_split('[/]') %>% unlist %>% as.numeric() %>% .[1]) /  (em$GeneRatio[j] %>% str_split('[/]') %>% unlist %>% as.numeric() %>% .[2])) /((em$BgRatio[j] %>% str_split('[/]') %>% unlist %>% as.numeric() %>% .[1]) /  (em$BgRatio[j] %>% str_split('[/]') %>% unlist %>% as.numeric() %>% .[2]))
  }
  em = arrange(em, desc(`Obs-Exp`))
  
  if ( length(em[,1]) != 0){
    write.csv(em, paste0('./ORA/H15_top5/20231225_', i,'.csv'))
  }
}

## GSEA ####
DEG = readRDS('DEG_GEMM_pheno_tune.rds')
df = readRDS('/Users/jooyoung/Dropbox/samir/rds_jy/SCENIC/Human_h15_human_rss_top5.rds')

for ( ii in DEG$cluster %>% unique ){
  dir.create('./GSEA/H15_top5/')
  
  deg = DEG[DEG$cluster == ii, ]
   
  i =  ii %>% str_remove('[/]')
  print(i)
  deg <- deg %>% mutate(statistics = qnorm(p_val/2, lower.tail=F) * sign(avg_log2FC),
                         max_statistics_second = max(statistics[statistics != Inf]),
                         min_statistics_second = min(statistics[statistics != -Inf]),
                         statistics_corrected = ifelse(statistics == Inf, max_statistics_second + 1,
                                                       ifelse(statistics == -Inf,  min_statistics_second - 1,
                                                              statistics)))
  statdf =  deg %>% dplyr::select('gene','statistics_corrected')
  colnames(statdf) = c('symbol','stat')
  ranks <- tibble::deframe(statdf)
  
  pathwaysH = split(x = df$gene_symbol, f = df$gs_cat)
  fgseaRes <- fgsea(pathways=pathwaysH, stats=ranks, nproc = 1,nPermSimple = 100000)
  fgseaResTidy_p <- fgseaRes %>%
    as_tibble() %>%
    arrange(desc(NES))
  
  fgseaResTidy_p$enriched = ifelse(fgseaResTidy_p$NES < 0 , 'Others',i)
  fgseaResTidy_p$enriched = factor(fgseaResTidy_p$enriched, levels = rev(c('Others',i )))
  
  for (j in 1:length(fgseaResTidy_p$pathway)){
    print(j)
    fgseaResTidy_p$leadingEdge_genes[j] = fgseaResTidy_p$leadingEdge[j] %>% unlist() %>% paste(collapse = ', ')
  }
  fgseaResTidy_p$leadingEdge = NULL
  write.csv(fgseaResTidy_p,paste0('./GSEA/H15_top5/20231225_',i,'_regulon.csv'))
}

## Supple figure ####
jy = import_list('/Users/jooyoung/Dropbox/samir/rds_jy/SCENIC/20231225_GEMM_Human_h15_rss_overlap_pheno_tuned_top5_v7.xlsx')

gsea = as.data.frame(matrix(ncol = 7)); colnames(gsea) = c("pathway"   ,        "pval"          ,    "padj"            ,  "NES"              , "enriched"       ,   "leadingEdge_genes", "Group")
ora =  as.data.frame(matrix(ncol = 10)); colnames(ora) = c("Description", "GeneRatio" ,  "BgRatio"   ,  "pvalue"   ,   "p.adjust" ,   "qvalue"    ,  "geneID"   , "Group" ,'GeneRatio_old' ,'Obs-Exp')

for (i in names(jy)){
  print(i)
  tmp = jy[[i]] %>% as.data.frame()
  colnames(tmp) = tmp[1,] ; tmp = tmp[-1,]
  
  tmp_gsea = tmp[,1:6]
  tmp_gsea$Group= i
  gsea = rbind(gsea, tmp_gsea)
  
  tmp_ora = tmp[,7:13]
  tmp_ora = tmp_ora[! tmp_ora$Description %>% is.na(),]
  tmp_ora$Group = i
  tmp_ora$GeneRatio_old = tmp_ora$GeneRatio
  tmp_ora$`Obs-Exp` = NA
  for (j in 1:length(tmp_ora$Description)){
    tmp_ora$GeneRatio[j] = (tmp_ora$GeneRatio_old[j] %>% str_split('[/]') %>% unlist %>% as.numeric() %>% .[1]) /  (tmp_ora$GeneRatio_old[j] %>% str_split('[/]') %>% unlist %>% as.numeric() %>% .[2])
    tmp_ora$`Obs-Exp`[j] = ((tmp_ora$GeneRatio_old[j] %>% str_split('[/]') %>% unlist %>% as.numeric() %>% .[1]) /  (tmp_ora$GeneRatio_old[j] %>% str_split('[/]') %>% unlist %>% as.numeric() %>% .[2])) /((tmp_ora$BgRatio[j] %>% str_split('[/]') %>% unlist %>% as.numeric() %>% .[1]) /  (tmp_ora$BgRatio[j] %>% str_split('[/]') %>% unlist %>% as.numeric() %>% .[2]))
  }
  
  ora = rbind(ora, tmp_ora)
}

##### GSEA #####
gsea = gsea[-1,]
gsea$Group %>% unique

gsea$Group[gsea$Group == "Ck8+" ] = "CK8 Luminal"; gsea$Group[gsea$Group == "Ar+Ascl1+"  ] = "Ar/Ascl1+";gsea$Group[gsea$Group =="STAT1|2" ] = "Stat1/2 Inflam";gsea$Group[gsea$Group == "Trp63+"  ] = "Trp63 Basal"; gsea$Group[gsea$Group == "Trp63+ Sox2|4+"   ] = "Trp63+ Sox4/6";gsea$Group[gsea$Group == "Pou2f3+" ] = "Pou2f3+" ;gsea$Group[gsea$Group == "Twist1|2 EMT"    ] = "Twist1/2 EMT"; gsea$Group[gsea$Group == "Ascl1+"  ] = "NEPC-A"

gsea$Group = factor(gsea$Group, levels = c("Ar/Ascl1+"    ,  "NEPC-A" ,"Twist1/2 EMT", "Pou2f3+" , 'Tff3+',  "Trp63 Basal", "CK8 Luminal" , "Stat1/2 Inflam", "Trp63+ Sox4/6" ))

gsea$pathway[gsea$pathway == 'TCFL2+ WNT'] = "TCF7L2+ WNT"
gsea$pathway = factor(gsea$pathway, levels = rev(c("TCF7L2+ WNT","MAFG+","IRF2+ Inflammatory","SOX2/4+ Embryonic EMT" ,"FOSL1+ AP-1", "AR+ HOXB13-","AR+ IRF+ Inflammatory" ,"AR+ GI","AR+ HOXB13+","AR+ HOXB13+ FOXA1+","NEPC-A","NEPC-A/SOX6" ,"NEPC-N")))

gsea$padj = as.numeric(gsea$padj)
gsea$NES = as.numeric(gsea$NES)

gsea$padj_size_log = -log(gsea$padj)
gsea$padj_size_log[gsea$padj >= 0.05] = 0.01

pdf('/Users/jooyoung/Google Drive/My Drive/PNAS_fin/20231225_SF5_GEMM_Human_ref_GSEA_dotplpt_top5.pdf', width = 6.8, height = 6.7)
ggplot(gsea, aes(x=Group, y=pathway,color=NES, size = padj_size_log)) +geom_point() + scale_colour_gradient2(low = lighten("lightgray",0.4),mid = lighten("lightgray",0.8),high = 'red3',midpoint = 0, name = 'NES',guide=guide_colorbar(reverse=F) ) +xlab(NULL) + ylab(NULL) + scale_size(range=c(1, 10)) + theme_classic() + theme(text = element_text(size = 15), axis.text.x = element_text(size = 10, angle = 20, hjust = 1) ) + guides(size=guide_legend(title="-log(p.adj)")) # , limits = c(-2.5, 2.5),oob = scales::squish
dev.off()

##### ORA #####
ora = ora[-1,]

ora$Group %>% unique

ora$Group[ora$Group == "Ck8+" ] = "CK8 Luminal"; ora$Group[ora$Group == "Ar+Ascl1+"  ] = "Ar/Ascl1+";ora$Group[ora$Group =="STAT1|2" ] = "Stat1/2 Inflam";ora$Group[ora$Group == "Trp63+"  ] = "Trp63 Basal"; ora$Group[ora$Group == "Trp63+ Sox2|4+"   ] = "Trp63+ Sox4/6";ora$Group[ora$Group == "Pou2f3+" ] = "Pou2f3+" ;ora$Group[ora$Group == "Twist1|2 EMT"    ] = "Twist1/2 EMT"; ora$Group[ora$Group == "Ascl1+"  ] = "NEPC-A"

ora$Group = factor(ora$Group, levels = c("Ar/Ascl1+"    ,  "NEPC-A" ,"Twist1/2 EMT", "Pou2f3+" , 'Tff3+',  "Trp63 Basal", "CK8 Luminal" , "Stat1/2 Inflam", "Trp63+ Sox4/6" ))

ora$Description[ora$Description == 'TCFL2+ WNT'] = "TCF7L2+ WNT"
ora$Description = factor(ora$Description, levels = rev(c("TCF7L2+ WNT","MAFG+","IRF2+ Inflammatory","SOX2/4+ Embryonic EMT" ,"FOSL1+ AP-1", "AR+ HOXB13-","AR+ IRF+ Inflammatory" ,"AR+ GI","AR+ HOXB13+","AR+ HOXB13+ FOXA1+","NEPC-A","NEPC-A/SOX6" ,"NEPC-N")))

ora$GeneRatio = as.numeric(ora$GeneRatio)
ora$qvalue = as.numeric(ora$qvalue) ; ora$qvalue_size_log = -log(ora$qvalue)
ora$qvalue_size_log [ora$qvalue >= 0.05 ] = 0.01

pdf('/Users/jooyoung/Google Drive/My Drive/PNAS_fin/Suppl/20231225_SF5_GEMM_Human_ref_ORA_dotplpt_top5.pdf', width = 6.7, height = 6.7)
ggplot(ora, aes(x=Group, y=Description,color=`Obs-Exp`)) +geom_point(aes(size = qvalue_size_log)) + scale_color_continuous(low=lighten("lightgray",0.8), high="red3", name = 'Obs/Exp',guide=guide_colorbar(reverse=F) ) +xlab(NULL) + ylab(NULL) + scale_size(range=c(1, 10)) + theme_classic() + theme(text = element_text(size = 15) , axis.text.x = element_text(size = 10, angle = 20, hjust = 1)  ) + guides(size=guide_legend(title="-log(q-value)"))
dev.off()



# Human GRN with GEMM SCENIC Geneset ####
cancer = readRDS('/Users/jooyoung/OneDrive - 고려대학교/samir/rds_jy/msk.integrated.remove.cellcycle.tumor.fastmnn.subtype.nohepato.figure4.rds') # object without CSPC cells

Idents(cancer) = cancer$Regulon_h15_modi_ano # GRN annotation
DEG = FindAllMarkers(cancer,  min.pct = 0.3, min.diff.pct = 0,logfc.threshold = 0, test.use = 'MAST')
saveRDS(DEG,'Human_SCENIC_DEG.rds')

## ORA ####
df = readRDS('~/Documents/samir_macbook/SCENIC_mus/GEMM_GRN_geneset_rss_top5.rds')
colnames(df) = c("gs_name", "human_gene_symbol")

DEG = readRDS('Human_SCENIC_DEG.rds')

for (ii in DEG$cluster %>% unique){
  print(ii)
  DEG_filt = DEG[DEG$cluster == ii,]
  DEG_filt = arrange(DEG_filt,desc(avg_log2FC))
  DEG_filt = DEG_filt[DEG_filt$p_val_adj < 0.05,]
  DEG_filt = DEG_filt[DEG_filt$avg_log2FC > 0.25,]
  genes = DEG_filt$gene
  
  i =  ii %>% str_remove('[/]')
  print(i)
  
  dir.create('./ORA/GEMM_top5', showWarnings = F)
  go_gene_sets = df
  em <- clusterProfiler::enricher(genes, TERM2GENE=go_gene_sets, pvalueCutoff =1, qvalueCutoff = 1, minGSSize = 1, maxGSSize = 5000)
  
  em = em@result
  
  for (j in 1:length(em$Description)){
    em$`Obs-Exp`[j] = ((em$GeneRatio[j] %>% str_split('[/]') %>% unlist %>% as.numeric() %>% .[1]) /  (em$GeneRatio[j] %>% str_split('[/]') %>% unlist %>% as.numeric() %>% .[2])) /((em$BgRatio[j] %>% str_split('[/]') %>% unlist %>% as.numeric() %>% .[1]) /  (em$BgRatio[j] %>% str_split('[/]') %>% unlist %>% as.numeric() %>% .[2]))
  }
  em = arrange(em, desc(`Obs-Exp`))
  
  if ( length(em[,1]) != 0){
    write.csv(em, paste0('./ORA/GEMM_top5/', i,'_regulon.csv'))
  }
}

## GSEA ####
df = readRDS('~/Documents/samir_macbook/SCENIC_mus/GEMM_GRN_geneset_rss_top5.rds')
df = df[,c(2,1)]
colnames(df) = c("gene_symbol", "gs_cat"    )

list = list.files(path = './GSEA/GEMM/', pattern = '^20231123_.+rank.rds')
list = list[-1]
list = list %>% str_remove('20231123_') %>% str_remove('_rank.rds')

for ( i in DEG$cluster %>% unique()){
  dir.create('./GSEA/GEMM_top5/',showWarnings = F)
  
  ii = i %>% str_remove('[+]') %>% str_remove(" ") %>% str_remove('[/]')
  
  deg <- deg %>% mutate(statistics = qnorm(p_val/2, lower.tail=F) * sign(avg_log2FC),
                        max_statistics_second = max(statistics[statistics != Inf]),
                        min_statistics_second = min(statistics[statistics != -Inf]),
                        statistics_corrected = ifelse(statistics == Inf, max_statistics_second + 1,
                                                      ifelse(statistics == -Inf,  min_statistics_second - 1,
                                                             statistics)))
  statdf =  deg %>% dplyr::select('gene','statistics_corrected')
  colnames(statdf) = c('symbol','stat') 
  ranks <- tibble::deframe(statdf)
  
  iii = i %>% str_remove('[/]')
  pathwaysH = split(x = df$gene_symbol, f = df$gs_cat)
  fgseaRes <- fgsea(pathways=pathwaysH, stats=ranks, nproc = 1,nPermSimple = 100000)
  fgseaResTidy_p <- fgseaRes %>%
    as_tibble() %>%
    arrange(desc(NES))
  
  fgseaResTidy_p$enriched = ifelse(fgseaResTidy_p$NES < 0 , 'Others',i)
  fgseaResTidy_p$enriched = factor(fgseaResTidy_p$enriched, levels = rev(c('Others',i )))
  
  for (j in 1:length(fgseaResTidy_p$pathway)){
    print(j)
    fgseaResTidy_p$leadingEdge_genes[j] = fgseaResTidy_p$leadingEdge[j] %>% unlist() %>% paste(collapse = ', ')
  }
  fgseaResTidy_p$leadingEdge = NULL
  write.csv(fgseaResTidy_p,paste0('./GSEA/GEMM_top5/',ii,'_regulons.csv'))
}


## Supple figure ####
jy = import_list('SCENIC/20231225_Human_GEMM_pheno_overlap_top5_v4.xlsx')

gsea = as.data.frame(matrix(ncol = 7)); colnames(gsea) = c("pathway"   ,        "pval"          ,    "padj"            ,  "NES"              , "enriched"       ,   "leadingEdge_genes", "Group")
ora =  as.data.frame(matrix(ncol = 10)); colnames(ora) = c("Description", "GeneRatio" ,  "BgRatio"   ,  "pvalue"   ,   "p.adjust" ,   "qvalue"    ,  "geneID"   , "Group" ,'GeneRatio_old' ,'Obs-Exp')

for (i in names(jy)){
  print(i)
  tmp = jy[[i]] %>% as.data.frame()
  colnames(tmp) = tmp[1,] ; tmp = tmp[-1,]
  
  tmp_gsea = tmp[,1:6]
  tmp_gsea$Group= i
  gsea = rbind(gsea, tmp_gsea)
  
  tmp_ora = tmp[,7:13]
  tmp_ora = tmp_ora[! tmp_ora$Description %>% is.na(),]
  tmp_ora$Group = i
  tmp_ora$GeneRatio_old = tmp_ora$GeneRatio
  tmp_ora$`Obs-Exp` = NA
  for (j in 1:length(tmp_ora$Description)){
    tmp_ora$GeneRatio[j] = (tmp_ora$GeneRatio_old[j] %>% str_split('[/]') %>% unlist %>% as.numeric() %>% .[1]) /  (tmp_ora$GeneRatio_old[j] %>% str_split('[/]') %>% unlist %>% as.numeric() %>% .[2])
    tmp_ora$`Obs-Exp`[j] = ((tmp_ora$GeneRatio_old[j] %>% str_split('[/]') %>% unlist %>% as.numeric() %>% .[1]) /  (tmp_ora$GeneRatio_old[j] %>% str_split('[/]') %>% unlist %>% as.numeric() %>% .[2])) /((tmp_ora$BgRatio[j] %>% str_split('[/]') %>% unlist %>% as.numeric() %>% .[1]) /  (tmp_ora$BgRatio[j] %>% str_split('[/]') %>% unlist %>% as.numeric() %>% .[2]))
  }
  
  ora = rbind(ora, tmp_ora)
}

##### GSEA #####
gsea = gsea[-1,]

c("Ar+ Ascl1+" , "Ascl1+" ,"Twist+ EMT", "Pou2f3+ Tuft cell" , 'Tff3+',  "Trp63+ basal like", "Ck8+ luminal like" , "Ck8+ Irf+ Hnf4a Inflammatory, GI", "Trp63+ Sox4/6+" )
gsea$pathway = factor(gsea$pathway, levels = rev(c("Ar/Ascl1+" , "NEPC-A" ,"Twist1/2 EMT",'Pou2f3+','Tff3+',  "Trp63 Basal", "CK8 Luminal" , "Stat1/2 Inflam", "Trp63+ Sox4/6"  )))

gsea$Group %>% unique
gsea$Group [gsea$Group  ==  "NEPC-A|SOX6" ] = "NEPC-A/SOX6"
gsea$Group [gsea$Group  ==  "SOX2|4+ Embyonic" ] ="SOX2/4+ Embryonic EMT"
gsea$Group [gsea$Group  ==  "AR+ IRF+ Inflam"  ] ="AR+ IRF+ Inflammatory" 

gsea$Group = factor(gsea$Group, levels = c("TCF7L2+ WNT","MAFG+","IRF2+ Inflammatory","SOX2/4+ Embryonic EMT" ,"FOSL1+ AP-1",'AR+ HOXB13+', "AR+ HOXB13-","AR+ IRF+ Inflammatory" ,"AR+ GI","AR+ HOXB13+ FOXA1+","NEPC-A","NEPC-A/SOX6" ,"NEPC-N"))

gsea$padj = as.numeric(gsea$padj)
gsea$NES = as.numeric(gsea$NES)

gsea$padj_size_log = -log(gsea$padj)
gsea$padj_size_log[gsea$padj >= 0.05] = 0.01

pdf('/Users/jooyoung/Google Drive/My Drive/PNAS_fin/20231225_SF5_Human_GEMM_ref_GSEA_dotplpt_top5.pdf', width = 8, height = 6.7)
ggplot(gsea, aes(x=Group, y=pathway,color=NES, size = padj_size_log)) +geom_point() + scale_colour_gradient2(low = lighten("lightgray",0.4),mid = lighten("lightgray",0.8),high = 'red3',midpoint = 0, name = 'NES',guide=guide_colorbar(reverse=F) ) +xlab(NULL) + ylab(NULL) + scale_size(range=c(1, 10)) + theme_classic() + theme(text = element_text(size = 15), axis.text.x = element_text(size = 10, angle = 30, hjust = 1) ) + guides(size=guide_legend(title="-log(p.adj)")) # , limits = c(-2.5, 2.5),oob = scales::squish
dev.off()

##### ORA #####
ora = ora[-1,]

ora$Group %>% unique
ora$Group [ora$Group  ==  "NEPC-A|SOX6" ] = "NEPC-A/SOX6"
ora$Group [ora$Group  ==  "SOX2|4+ Embyonic" ] ="SOX2/4+ Embryonic EMT"
ora$Group [ora$Group  ==  "AR+ IRF+ Inflam"  ] ="AR+ IRF+ Inflammatory" 

ora$Group = factor(ora$Group, levels = c("TCF7L2+ WNT","MAFG+","IRF2+ Inflammatory","SOX2/4+ Embryonic EMT" ,"FOSL1+ AP-1",'AR+ HOXB13+', "AR+ HOXB13-","AR+ IRF+ Inflammatory" ,"AR+ GI","AR+ HOXB13+ FOXA1+","NEPC-A","NEPC-A/SOX6" ,"NEPC-N"))

ora$Description = factor(ora$Description, levels =  rev(c("Ar/Ascl1+" , "NEPC-A" ,"Twist1/2 EMT",'Pou2f3+','Tff3+',  "Trp63 Basal", "CK8 Luminal" , "Stat1/2 Inflam", "Trp63+ Sox4/6"  )))

ora$GeneRatio = as.numeric(ora$GeneRatio)
ora$qvalue = as.numeric(ora$qvalue) ; ora$qvalue_size_log = -log(ora$qvalue)
ora$qvalue_size_log [ora$qvalue >= 0.05 ] = 0.01

pdf('/Users/jooyoung/Google Drive/My Drive/PNAS_fin/20231225_SF5_Human_GEMM_ref_ORA_dotplpt_top5.pdf', width = 8, height = 6.7)
ggplot(ora, aes(x=Group, y=Description,color=`Obs-Exp`)) +geom_point(aes(size = qvalue_size_log)) + scale_color_continuous(low=lighten("lightgray",0.8), high="red3", name = 'Obs/Exp',guide=guide_colorbar(reverse=F),limits = c(0, 3),oob = scales::squish ) +xlab(NULL) + ylab(NULL) + scale_size(range=c(1, 10)) + theme_classic() + theme(text = element_text(size = 15) , axis.text.x = element_text(size = 10, angle = 20, hjust = 1)  ) + guides(size=guide_legend(title="-log(q-value)"))
dev.off()

