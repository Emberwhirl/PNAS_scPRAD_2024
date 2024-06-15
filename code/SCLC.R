library(Seurat)
library(dplyr)
library(stringr)
library(SeuratDisk)
library(anndata)
library(ggplot2)
library(cowplot)
library(msigdbr)
library(fgsea)
library(ComplexHeatmap)
library(randomcoloR)
library(circlize)
library(SoupX)
library(readxl)
library(scCustomize)
library(clusterProfiler)
library(ggrepel)
library(ggnewscale)
set.seed(1234)


# Figure 3C-D ####
## SCLC-A / N ####
sclc = readRDS('sclc.rds')
DEG = FindAllMarkers(sclc, only.pos = F, logfc.threshold = 0, min.pct = 0.3, test.use = 'MAST')
write.csv(DEG,'SCLC_DEG.csv')

df = readRDS('/Users/jooyoung/Dropbox/samir/rds_jy/SCENIC/Human_h15_human_rss_top5_hum.rds')
df$gene_symbol = df$gene_symbol %>% word(1,sep = ' ')

for ( i in c('SCLC-A','SCLC-N')){
  
  ranks = readRDS(paste0(i,'_rank.rds'))
  
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
  write.csv(fgseaResTidy_p,paste0('20231225_h15_top5_',i,'.csv'))
  
  fgseaResTidy_p = fgseaResTidy_p %>% filter(padj<0.05) %>% filter(!is.na(log2err))
  cl <-  c('red3','blue3')
  names(cl) <- paste( levels(fgseaResTidy_p$enriched  ))
  
  pdf(paste0('20231225_h15_top5_',i,'.pdf'),width = 15, height = 15)
  print(ggplot(fgseaResTidy_p, aes(reorder(pathway, NES), NES)) +
          geom_col(aes(fill=enriched)) +
          coord_flip() +
          labs(x="Pathway", y="Normalized Enrichment Score",
               title=i) +  theme_minimal() + theme(legend.title = element_blank(), axis.text = element_text(size = 15, face = 'bold')) +  scale_fill_manual(values = cl, na.value = cl["NA"])) 
  dev.off()
  
}

### Dotplot ####
df1 = read.csv('20231225_h15_top5_SCLC-A.csv')
df1$Group = 'SCLC-A'
df2 = read.csv('20231225_h15_top5_SCLC-N.csv')
df2$Group = 'SCLC-N'
df = rbind(df1, df2)
df %>% colnames()

for (i in 1:length(df$leadingEdge_genes)){
  df$N_leadingEdge_genes[i] = df$leadingEdge_genes[i] %>% str_split(',') %>% unlist %>% length() %>% as.numeric()
}

path = df$pathway %>% unique
path[!path %in% c("NEPC-A" ,"NEPC-A/SOX6"  ,"NEPC-N" ,"TCF7L2+ WNT","MAFG+","IRF2+ Inflammatory","SOX2/4+ Embryonic EMT" ,"FOSL1+ AP-1", "AR+ HOXB13-","AR+ IRF+ Inflammatory" ,"AR+ GI","AR+ HOXB13+","AR+ HOXB13+ FOXA1+"  ) ]
c("NEPC-A" ,"NEPC-A/SOX6"  ,"NEPC-N" ,"TCF7L2+ WNT","MAFG+","IRF2+ Inflammatory","SOX2/4+ Embryonic EMT" ,"FOSL1+ AP-1", "AR+ HOXB13-","AR+ IRF+ Inflammatory" ,"AR+ GI","AR+ HOXB13+","AR+ HOXB13+ FOXA1+"  ) [!c("NEPC-A" ,"NEPC-A/SOX6"  ,"NEPC-N" ,"TCF7L2+ WNT","MAFG+","IRF2+ Inflammatory","SOX2/4+ Embryonic EMT" ,"FOSL1+ AP-1", "AR+ HOXB13-","AR+ IRF+ Inflammatory" ,"AR+ GI","AR+ HOXB13+","AR+ HOXB13+ FOXA1+"  )  %in% path]

df$pathway = factor(df$pathway,rev(c("NEPC-A" ,"NEPC-A/SOX6"  ,"NEPC-N" ,"TCF7L2+ WNT","MAFG+","IRF2+ Inflammatory","SOX2/4+ Embryonic EMT" ,"FOSL1+ AP-1", "AR+ HOXB13-","AR+ IRF+ Inflammatory" ,"AR+ GI","AR+ HOXB13+","AR+ HOXB13+ FOXA1+"  ) ))

df$padj = as.numeric(df$padj)
df$NES = as.numeric(df$NES)

df$padj_size_log = -log(df$padj)
df$padj_size_log[df$padj >= 0.05] = 0.01

pdf('/Users/jooyoung/Google Drive/My Drive/PNAS_fin/Suppl/Suppl_updated/SF6_SCLC_AN_GSEA_dot_top5.pdf', width = 6,height = 8/1.3)
ggplot(df, aes(x=Group, y=pathway,color=NES, size = padj_size_log)) +geom_point() +scale_colour_gradient2(low = lighten("lightgray",0.4),mid = lighten("lightgray",0.8),high = 'red3',midpoint = 0, name = 'NES',guide=guide_colorbar(reverse=F) ) +xlab(NULL) + ylab(NULL) + scale_size(range=c(1,10)) + theme_classic() + theme(text = element_text(size = 15))  + guides(size=guide_legend(title="-log(p.adj)"))
dev.off()

### DEG overlap ####
humantf = read.table('../DEG_gsea/Homo_sapiens_TF.gene_list.txt')
sclc_deg = read.csv('DEG_SCLC-AN.csv', row.names = 1) 
nepc_deg = readRDS('DEG_nepc_an.rds')

genes = intersect(rownames(sclc_deg), rownames(nepc_deg))

sclc_deg = sclc_deg[genes,]; colnames(sclc_deg) = paste0('SCLC_',colnames(sclc_deg))
nepc_deg = nepc_deg[genes,]; colnames(nepc_deg) = paste0('NEPC_',colnames(nepc_deg))

df = cbind(sclc_deg,nepc_deg)
df %>% colnames()


df$gene = rownames(df)
df$Type = NA
df$Type[df$SCLC_avg_log2FC > 0.4 & df$NEPC_avg_log2FC > 0.4] = 'ASCL1'
df$Type[df$SCLC_avg_log2FC <  -0.4 & df$NEPC_avg_log2FC <  -0.4] = 'NEUROD1'

df$tf = NA
df[rownames(df) %in% humantf$V1,'tf'] = rownames(df) [rownames(df) %in% humantf$V1]

df$Label = NA
df$Label[df$SCLC_avg_log2FC > 0.4 & df$NEPC_avg_log2FC > 0.4] = df$gene[df$SCLC_avg_log2FC > 0.4 & df$NEPC_avg_log2FC > 0.4]
df$Label[df$SCLC_avg_log2FC < -0.4 & df$NEPC_avg_log2FC < -0.4] = df$gene[df$SCLC_avg_log2FC < -0.4 & df$NEPC_avg_log2FC< -0.4]

df$Label_show = NA
df$Label_show[( (!df$Label %>% is.na()) & (!df$tf %>% is.na() ))] = df$Label[( (!df$Label %>% is.na()) & (!df$tf %>% is.na() ))]
df$Label_show[!df$Label_show %>% is.na()] = 'Transcription factor'
df$dotsize = ifelse(df$Type  %in% c('ASCL1','NEUROD1'), 4, 1)

for (i in df$gene){
  if ( ! (df[df$gene == i,'Label_show'] %>% is.na())){
    df[df$gene == i, 'Label'] = paste0("bold('",i,"')")
  }
}

pdf('/Users/jooyoung/Google Drive/My Drive/Main_Figure_v3/20231122_updated/Figure3D_SCLC_Human_A_N_dotplot_allgenes.pdf', width = 9/1.3,height = 8/1.3)
ggplot(df, aes(x = NEPC_avg_log2FC, y = SCLC_avg_log2FC, label = Label, col = Type)) + geom_point(aes(size=dotsize)) + scale_color_manual(values = c(col[12],col[11])) + new_scale_color()  + geom_text_repel(aes(col = Label_show), show.legend = F, parse = TRUE ,box.padding = 0.3, max.overlaps = Inf) + scale_color_manual(values = c('black'))+ theme_classic()  + ggtitle('A group vs N group',subtitle = 'ASCL1: log2FC > 0.4 for both \nNEUROD1: log2FC < -0.4 for both') + xlab('NEPC-N vs NEPC-A') + ylab('SCLC-N vs SCLC-A') + scale_size_continuous(range = c(0.5, 2.5)) +guides(size = "none")
dev.off()


## SCLC-P ####
df = readRDS('~/Documents/samir_macbook/SCENIC_mus/GEMM_GRN_geneset_rss_top5.rds')
colnames(df) = c("gs_cat", "gene_symbol")
DEG = read.csv('SCLC_DEG.csv')

for ( i in "SCLC-P" ){ # DEG$cluster %>% unique
  print(i)
  ranks = readRDS(paste0(i,'_rank.rds'))
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
  write.csv(fgseaResTidy_p,paste0('20231225_',i,'_GEMMGRN_top5.csv'))
  
  fgseaResTidy_p = fgseaResTidy_p %>% filter(padj<0.05) %>% filter(!is.na(log2err))
  cl <-  c('red3','blue3')
  names(cl) <- paste( levels(fgseaResTidy_p$enriched  ))
  pdf('20231225_GEMM_top5_SCLC-P.pdf',width = 15, height = 15)
  print(ggplot(fgseaResTidy_p, aes(reorder(pathway, NES), NES)) +
          geom_col(aes(fill=enriched)) +
          coord_flip() +
          labs(x="Pathway", y="Normalized Enrichment Score",
               title='SCLC-P') +  theme_minimal() + theme(legend.title = element_blank(), axis.text = element_text(size = 15, face = 'bold')) +  scale_fill_manual(values = cl, na.value = cl["NA"])) 
  dev.off()
  
}

### Dotplot ####
df1 = read.csv('20231225_SCLC-P_GEMMGRN_top5.csv')
df1$Group = 'SCLC-P'
df = df1
df %>% colnames()

for (i in 1:length(df$leadingEdge_genes)){
  df$N_leadingEdge_genes[i] = df$leadingEdge_genes[i] %>% str_split(',') %>% unlist %>% length() %>% as.numeric()
}

df$pathway %>% unique
df$pathway = factor(df$pathway,rev(df$pathway))

df$padj = as.numeric(df$padj)
df$NES = as.numeric(df$NES)

df$padj_size_log = -log(df$padj)
df$padj_size_log[df$padj >= 0.05] = 0.01

pdf('/Users/jooyoung/Google Drive/My Drive/PNAS_fin/Suppl/Suppl_updated/SF6_SCLC_P_GSEA_dot_top5.pdf', width = 8,height = 8/1.3)
ggplot(df, aes(x=Group, y=pathway,color=NES, size = padj_size_log)) +geom_point()  +scale_colour_gradient2(low = lighten("lightgray",0.4),mid = lighten("lightgray",0.8),high = 'red3',midpoint = 0, name = 'NES',guide=guide_colorbar(reverse=F) ) +xlab(NULL) + ylab(NULL) + scale_size(range=c(1,10)) + theme_classic() + theme(text = element_text(size = 15))+ guides(size=guide_legend(title="-log(p.adj)")) # , limits = c(-2.5, 2.5),oob = scales::squish
dev.off()

### DEG overlap ####
humantf = read.table('../DEG_gsea/Homo_sapiens_TF.gene_list.txt')
mus_deg = readRDS('~/Documents/samir_macbook/GEMM/samir_mouse/DEG_GEMM_pheno_tune_pou2f3.rds')
genes_mus = mus_deg$HGNC.symbol
sclc_deg = readRDS('sclcp_deg_formus.rds')
genes_sclc = rownames(sclc_deg)

match = intersect(genes_mus, genes_sclc)

mus_deg = mus_deg[mus_deg$HGNC.symbol %in% match,]
rownames(mus_deg) = mus_deg$HGNC.symbol
mus_deg = mus_deg[match,]
sclc_deg = sclc_deg[match,]

df = cbind(sclc_deg,mus_deg)
colnames(df)[c(7,10)] = paste0('GEMM_',colnames(df)[c(7,10)]) 

df$gene = rownames(df)
df$Type = NA
df$Type[df$SCLC_avg_log2FC > 0.25 & df$GEMM_avg_log2FC > 0.25] = 'POU2F3'
df$Type[df$SCLC_avg_log2FC <  -0.25 & df$GEMM_avg_log2FC <  -0.25] = 'Others'

df$Label = NA
df$Label[df$SCLC_avg_log2FC > 0.4 & df$GEMM_avg_log2FC > 0.4] = df$gene[df$SCLC_avg_log2FC > 0.4 & df$GEMM_avg_log2FC > 0.4]
df$Label[df$SCLC_avg_log2FC < -0.4 & df$GEMM_avg_log2FC < -0.4] = df$gene[df$SCLC_avg_log2FC < -0.4 & df$GEMM_avg_log2FC< -0.4]

df$Label_show = NA
df$Label_show[( (!df$Label %>% is.na()) & (!df$tf %>% is.na() ))] = df$Label[( (!df$Label %>% is.na()) & (!df$tf %>% is.na() ))]
df$Label_show[!df$Label_show %>% is.na()] = 'Transcription factor'
df$dotsize = ifelse(df$Type  %in% c('POU2F3','Others'), 4, 1)

tmp = df %>% filter((Type %in% c( "POU2F3",   "Others") ) & Label_show  %>% is.na() )
tmp = tmp %>% filter(!Label %>% is.na())
genes = tmp$HGNC.symbol[(tmp$Type == 'Others') & ( tmp$SCLC_avg_log2FC < -0.9 & tmp$GEMM_avg_log2FC < -2.5)]
genes = c(genes,tmp$HGNC.symbol[(tmp$Type == 'POU2F3') & ( tmp$SCLC_avg_log2FC  > 0.8 & tmp$GEMM_avg_log2FC > 0.8 )])

genes = c(genes,tmp$gene[tmp$Type == 'Others'] %>% sample(10))
genes = c(genes, tmp$gene[tmp$Type == 'POU2F3'] %>% sample(10))
genes = genes %>% unique
df$Label[(! df$Label %in% genes) & (df$Label_show %>% is.na())] = NA

for (i in df$gene){
  if ( ! (df[df$gene == i,'Label_show'] %>% is.na())){
    df[df$gene == i, 'Label'] = paste0("bold('",i,"')")
  }
}

pdf('/Users/jooyoung/Google Drive/My Drive/PNAS_fin/Suppl/SCLC_P_dotplot.pdf', width = 9/1.3,height = 8/1.3)
ggplot(df, aes(x = GEMM_avg_log2FC, y = SCLC_avg_log2FC, label = Label, col = Type)) + geom_point(aes(size=dotsize)) + scale_color_manual(values = c("#80461B",lighten("#5E716A", 0.7)),na.value=lighten("lightgray", 0.2) ) + new_scale_color()  + geom_text_repel(aes(col = Label_show), show.legend = F, parse = TRUE ,box.padding = 0.25, max.overlaps = Inf) + scale_color_manual(values = c('black'))+ theme_classic()  + ggtitle('P group vs Others',subtitle = 'POU2F3: log2FC > 0.4 for both \nOthers: log2FC < -0.4 for both') + xlab('Others vs GEMM-P') + ylab('Others vs SCLC-P') + scale_size_continuous(range = c(0.5, 2.5)) +guides(size = "none")
dev.off()
