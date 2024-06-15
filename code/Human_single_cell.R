library(Seurat)
library(dplyr)
library(ggplot2)
library(ggpubr)
library(ggrepel)
library(cowplot)
library(patchwork)
library(stringr)
library(grid)
library(ComplexHeatmap)
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
library(readxl)
set.seed(1234)

patient_col = c(c(c(darken('#EEDC82'),lighten('#F2A900',0.2),'#EEED09','#6F4E37','#CCA01D',darken('#CFFF04'),'gold','#CC7722')),c("#CDDEEB", "#89CFF0", "#0000FF", "#7FFFD4", "#367588", "#1ca9c9", "#0C2340","#A0B6AC" ,"#6F8078",'#CCCCFF' ), c( '#FF6EC7',darken('#FF2400', 0.1),darken('#FF2400', 0.4)))

scale_color_patient <- function(...){
  ggplot2:::manual_scale(
    'color', 
    values = setNames( patient_col ,tumor$patient %>% levels()), na.value =  setNames(patient_col ,tumor$patient %>% levels())["NA"],
    ...
  )
}

# Sup Fig 3A-C & Data processing & GEO rds data generation ####
lists = list.files(pattern = '^raw.+rds$')

## percent.mt < 25, nCount > 200, doubletFinder ####
for (i in lists){
  print(i)
  a = readRDS(i)
  
  print('Filter')
  # percentmt, nCount 
  a <- subset(a, subset = nFeature_RNA > 200 &  percent.mt < 25)
  
  print('Doublet')
  # DoubletFInder 
  a <- NormalizeData(a)
  a <- FindVariableFeatures(a, selection.method = "vst", nfeatures = 2000)
  a <- ScaleData(a)
  a <- RunPCA(a)
  a <- RunUMAP(a, dims = 1:30)
  
  number = length(colnames(a))
  
  rate = ifelse(number < 1650,0.008, ifelse(number < 3300, 0.016, ifelse(number < 4950, 0.024, ifelse(number < 6600, 0.032, ifelse(number < 8250, 0.04, ifelse(number < 9900, 0.048, ifelse(number < 11550, 0.056, ifelse(number < 13200, 0.064, ifelse(number < 14850, 0.072, 0.080)))))))))
  
  sweep.res <- paramSweep_v3(a, 1:30, sct = F) 
  sweep.stats_kidney <- summarizeSweep(sweep.res, GT = FALSE)
  bcmvn_kidney <- find.pK(sweep.stats_kidney)
  idx = which(bcmvn_kidney$BCmetric == (bcmvn_kidney$BCmetric %>% max()))
  pk = bcmvn_kidney$pK[idx]  %>% as.character() %>% as.numeric()
  
  ## Homotypic Doublet Proportion Estimate ------------
  annotations <- a@meta.data$seurat_clusters
  homotypic.prop <- modelHomotypic(annotations)  #  homotypic doublet proportions are modeled as the sum of squared annotation frequencies. 
  
  # https://kb.10xgenomics.com/hc/en-us/articles/360001378811-What-is-the-maximum-number-of-cells-that-can-be-profiled-
  nExp_poi <- round(rate*nrow(a@meta.data)) # * 앞 에다가 % 적어주기 gse16: 2.4%
  nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))
  a <- doubletFinder_v3(a, pN = 0.25, pK = pk, nExp = nExp_poi, PCs = 1:30)
  a  <- doubletFinder_v3( a , PCs = 1:30, pN = 0.25, pK = pk, nExp = nExp_poi.adj, reuse.pANN = colnames(a@meta.data)[colnames(a@meta.data) %>% startsWith('pANN')] , sct = FALSE)
  
  
  print('Doublet Filterout')
  # Filter Doublet
  column = colnames(a@meta.data) %>% startsWith('DF') %>% which()
  cells = rownames(a@meta.data)[a@meta.data[,column] == 'Singlet']
  a = subset(a, cells = cells)
  
  print('Save')
  saveRDS(a, paste0('Filtered_',i %>% str_remove('raw_')))
  assign(i %>% str_remove('raw_') %>% str_remove('.rds'), a)
  rm(a)
}

## Merge all the objects in the list ####
seuralist = ls()[ls() %>% startsWith('H')]
seuralist = list.files(pattern = '^Fil.+rds$')
seurat_object_list = list()
for (i in 1:length(seuralist)) {
  # seu = get(seuralist[i])
  seu = readRDS(seuralist[i])
  seu <- RenameCells(seu,
                     add.cell.id = seuralist[i] %>% str_remove('Filtered_') %>% str_remove('.rds'))
  seurat_object_list[[i]] <- seu
  
}
rm(seu)
gc()
merged_combined <- Merge_Seurat_List(seurat_object_list)
merged_combined$sample %>% unique
meta = merged_combined@meta.data
meta$barcode = rownames(meta)
meta$sample %>% unique
meta = meta[,! colnames(meta) %>% str_detect('DF|pAN')]
merged_combined@meta.data = meta
merged_combined$sample %>% unique
merged_combined$patient = merged_combined$sample %>% word(sep = '_',start = 1)

## Unify patient name ####
merged_combined$patient = merged_combined$patient %>% str_remove('A|B') 
merged_combined$patient_prev = merged_combined$patient
merged_combined$patient = merged_combined$patient %>% as.character()
merged_combined$patient = ifelse(merged_combined$patient == 'HP95','MSK-HP01',ifelse(merged_combined$patient_prev  == 'HP96','MSK-HP02',ifelse(merged_combined$patient_prev == 'HP97','MSK-HP03',ifelse(merged_combined$patient_prev == 'HP99','MSK-HP04',ifelse(merged_combined$patient_prev == 'HP100','MSK-HP05',ifelse(merged_combined$patient_prev == 'HP101','MSK-HP06',ifelse(merged_combined$patient_prev == 'HMP23','MSK-HP07',ifelse(merged_combined$patient_prev == 'HMP24','MSK-HP08',ifelse(merged_combined$patient_prev == 'HMP05','MSK-HP09',ifelse(merged_combined$patient_prev == 'HMP08','MSK-HP10',ifelse(merged_combined$patient_prev == 'HMP11','MSK-HP11',ifelse(merged_combined$patient_prev == 'HMP13','MSK-HP12',ifelse(merged_combined$patient_prev == 'HMP14', 'MSK-HP13',ifelse(merged_combined$patient_prev == 'HMP19','MSK-HP14',ifelse(merged_combined$patient_prev == 'HMP20','MSK-HP15',ifelse(merged_combined$patient_prev == 'HMP22','MSK-HP16',ifelse(merged_combined$patient_prev == 'HMP25','MSK-HP17',ifelse(merged_combined$patient_prev == 'HMP26','MSK-HP18',ifelse(merged_combined$patient_prev == 'HMP04','MSK-HP19',ifelse(merged_combined$patient_prev == 'HMP16','MSK-HP20','MSK-HP21'))))))))))))))))))))

merged_combined$patient = as.character(merged_combined$patient)
merged_combined$patient = factor(merged_combined$patient, levels = c("MSK-HP01", "MSK-HP02", "MSK-HP03", "MSK-HP04", "MSK-HP05", "MSK-HP06", "MSK-HP07", "MSK-HP08", "MSK-HP09", "MSK-HP10", "MSK-HP11", "MSK-HP12", "MSK-HP13", "MSK-HP14","MSK-HP15","MSK-HP16", "MSK-HP17", "MSK-HP18", "MSK-HP19", "MSK-HP20", "MSK-HP21"))


## Add subtype ####
sub = readxl::read_xlsx('/Users/jooyoung/OneDrive - 고려대학교/PRAD_Samir_Choi/Single_CELL_DATA_05092023.v4.xlsx') %>% as.data.frame()
sub = sub[sub$Group != 'Normal',]
sub=sub[(sub$`Sample ID`) %in% (merged_combined$sample %>% word(sep = '_',start = 1) %>% unique),]

(merged_combined$sample %>% word(sep = '_',start = 1) %>% unique)  %in% sub$`Sample ID` 

colnames(sub)[1] = 'sample'

test = full_join(merged_combined@meta.data, sub, by = 'sample')
dim(test)
colnames(test)
test = test[,c(1:7,9:11)]
colnames(test)[8] = 'subtype_large'
test$subtype = test$subtype_large
test$subtype = test$subtype  %>% gsub(".*-","",.)
test$subtype[test$subtype_large %>% str_detect('DNPC')] = 'CRPC'
test$subtype %>% table
dim(test)
View(test)
test = distinct(test)
test$barcode [!test$barcode %in% rownames(merged_combined@meta.data)]
test = test[! test$barcode %>% is.na(),]

which(test$barcode != rownames(merged_combined@meta.data))
which(test$barcode == rownames(merged_combined@meta.data)) %>% length()

merged_combined@meta.data = test

merged_combined = DietSeurat(merged_combined)
saveRDS(merged_combined,'/Users/jooyoung/OneDrive - 고려대학교/samir/rds_jy/msk.integrated.jy.rds')  
rm(merged_combined)


## Remove uninformative genes and regress cell cycle ####
msk = readRDS('/Users/jooyoung/OneDrive - 고려대학교/samir/rds_jy/msk.integrated.jy.rds')
dim(msk)
msk@meta.data %>% dim
msk@meta.data  = msk@meta.data[,-9]

count = msk@assays$RNA@counts
meta = msk@meta.data
rm(msk)
gc()
idx = ! rownames(count) %>% str_detect('^MT-|^MTMR|^MTND|NEAT1|EGFP|TDTOMATO|TMSB4X|TMSB10|^RPS|^RPL|^MRP|^FAU$|UBA52|MALAT1')
count = count[idx,]
dim(count)
dim(meta)
rownames(meta) = meta$barcode
which(colnames(count) != meta$barcode)

msk = CreateSeuratObject(count = count,meta.data = meta )
msk = NormalizeData(msk)

s.genes <- cc.genes$s.genes ; g2m.genes <- cc.genes$g2m.genes; msk <- CellCycleScoring(msk, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE); msk <- ScaleData(msk, features = rownames(msk),vars.to.regress = c("S.Score", "G2M.Score"))
msk <- FindVariableFeatures(msk, selection.method = "vst", nfeatures = 2000)
msk <- RunPCA(msk, features = VariableFeatures(object = msk))
msk <- FindNeighbors(msk, dims = 1:30)
msk <- FindClusters(msk, resolution = 0.5)
msk <- RunUMAP(msk, dims = 1:30)

msk$subtype[msk$subtype %>% is.na] = 'CSPC'
msk$subtype = factor(msk$subtype, levels = c('CSPC','CRPC','NEPC'))
saveRDS(msk,'/Users/jooyoung/OneDrive - 고려대학교/samir/rds_jy/msk.integrated.remove.cellcycle.jy.rds')

## Coarse annotation  ####
msk = readRDS('/Users/jooyoung/OneDrive - 고려대학교/samir/rds_jy/msk.integrated.remove.cellcycle.jy.rds')
markers = c('PTPRC','COL1A2','TAGLN','ACTA2','PECAM1','EPCAM','YAP1','EZH2','AR','KLK2','KLK3','FOLH1','STEAP1','TACSTD2','DLL3','CHGA','CHGB')
for (i in markers){
  a = FeaturePlot(msk, i,cols = c('lightgray','red3'),label = T , raster=FALSE)
  b = VlnPlot_scCustom(msk,i) + theme(axis.text.x = element_text(size = 13,angle = 45))  + theme(plot.title = element_text(size = 25))
  c = cowplot::plot_grid(a,b, rel_widths = c(1,1.4))
  assign(paste0('g_',i),c)
}
pic = ls()[ls() %>% startsWith("g_")]
names = c()
for(l in pic){names = append(l %>% strsplit('g_') %>%unlist() %>% .[2],names)}
names = sort(names)
plot_list = list()
for (k in 1:length(pic )){
  print(pic[k])
  plot_list[[k]] = get(pic [k])
}
idx = match(markers, names)
plot_grid_args = c(plot_list[idx],list(ncol = 4))
d =   do.call(plot_grid,plot_grid_args)

pdf('Dim10_Markergene_exp.pdf', width = 11*2.4*4, height = 6.5*1.3*ceiling(length(markers)/4))
print(d)
dev.off()

saveRDS(msk,'/Users/jooyoung/OneDrive - 고려대학교/samir/rds_jy/msk.integrated.remove.cellcycle.jy.rds')

## Immune cells ####
msk = readRDS('msk.integrated.remove.cellcycle.jy.rds')
im = subset(msk, seurat_clusters %in% c(3,42,34,25,35,1,28,15,10,31,29,22)); rm(msk)
im = DietSeurat(im)
im = NormalizeData(im)
s.genes <- cc.genes$s.genes ; g2m.genes <- cc.genes$g2m.genes; im <- CellCycleScoring(im, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE); im <- ScaleData(im, features = rownames(im),vars.to.regress = c("S.Score", "G2M.Score"))

im <- FindVariableFeatures(im, selection.method = "vst", nfeatures = 2000)
im <- ScaleData(im)
im <- RunPCA(im)
ElbowPlot(im, ndims = 50)
im <- RunUMAP(im, dims = 1:20)
im <- FindNeighbors(im, dims = 1:20); gc()
im <- FindClusters(im, resolution = 0.2); gc()
saveRDS(im,'msk.integrated.remove.cellcycle.im.rds')

markers = c('PTPRC','CD3D','CD79A','LYZ','CD14')
for (i in markers){
  a = FeaturePlot(im, i,cols = c('lightgray','red4'), pt.size = 0.1, label = T, repel = T ) 
  b = VlnPlot_scCustom(im,i) + theme(axis.text.x = element_text(size = 13,angle = 45))  + theme(plot.title = element_text(size = 25), axis.text.x = element_text(size = 15))
  c = cowplot::plot_grid(a,b, rel_widths = c(1,1.2))
  assign(paste0('g_',i),c)
}
pic = ls()[ls() %>% startsWith('g_')]
names = c()
for(l in pic){names = append(l %>% strsplit('g_') %>%unlist() %>% .[2],names)}
names = sort(names)
plot_list = list()
for (k in 1:length(pic )){
  print(pic[k])
  plot_list[[k]] = get(pic [k])
}
idx = match(markers, names)
plot_grid_args = c(plot_list[idx],list(ncol = 3))
d =   do.call(plot_grid,plot_grid_args)

png(paste0('20231015_immune_marker.png'), width = 70*11*3*2, height = 70*5*2*ceiling(length(markers)/3))
print(d)
dev.off()

## Myeloid ####
mye = subset(im, subset = seurat_clusters %in% c(10,11,16,1,9,5,17,13))
mye = DietSeurat(mye)
mye = NormalizeData(mye)
s.genes <- cc.genes$s.genes ; g2m.genes <- cc.genes$g2m.genes; mye <- CellCycleScoring(mye, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE); mye <- ScaleData(mye, features = rownames(mye),vars.to.regress = c("S.Score", "G2M.Score"))
mye <- FindVariableFeatures(mye, selection.method = "vst", nfeatures = 2000)
mye <- RunPCA(mye)
ElbowPlot(mye, ndims = 50)
mye =DietSeurat(mye)
mye = RunFastMNN(object.list = SplitObject(mye, split.by = "subtype"), features = length(rownames(mye))) 
data <- as.matrix(mye@reductions$mnn@cell.embeddings[,1:30])
clusterings = list()
n = 1
for (k in seq(20,50,5)){
  Rphenograph_out <- Rphenograph(data, k = k)
  pheno = membership(Rphenograph_out[[2]]) %>% as.factor()
  names(pheno) = rownames(mye@reductions$mnn@cell.embeddings)
  clusterings[[n]] = pheno
  n = n+1
}
rand_score = c()
for (n in 1:((length(clusterings)) -1) ){
  print(n)
  score = genieclust::adjusted_rand_score(clusterings[[n]], clusterings[[n+1]])
  rand_score = c(rand_score, score)
}
rand_score
Rphenograph_out <- Rphenograph(data, k = 30) 
pheno = membership(Rphenograph_out[[2]]) 
mye$pheno_cluster = pheno 
mye$barcode <- colnames(mye)

mye$subtype = factor(mye$subtype, levels= c('CSPC','CRPC','NEPC'))
saveRDS(mye,'/Users/jooyoung/Dropbox/msk.integrated.remove.cellcycle.mye.subtype.rds')

## Lymphoid ####
lym = subset(im, subset = seurat_clusters %in% c(10,11,16,1,9,5,17,13), invert = T)
lym = DietSeurat(lym)
lym = NormalizeData(lym)
s.genes <- cc.genes$s.genes ; g2m.genes <- cc.genes$g2m.genes; lym <- CellCycleScoring(lym, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE); lym <- ScaleData(lym, features = rownames(lym),vars.to.regress = c("S.Score", "G2M.Score"))
lym <- FindVariableFeatures(lym, selection.method = "vst", nfeatures = 2000)
lym <- RunPCA(lym)
ElbowPlot(lym, ndims = 50)
lym =DietSeurat(lym)
lym = RunFastMNN(object.list = SplitObject(lym, split.by = "subtype")) 
data <- as.matrix(lym@reductions$mnn@cell.embeddings[,1:30])
clusterings = list()
n = 1
for (k in seq(20,50,5)){
  Rphenograph_out <- Rphenograph(data, k = k)
  pheno = membership(Rphenograph_out[[2]]) %>% as.factor()
  names(pheno) = rownames(lym@reductions$mnn@cell.embeddings)
  clusterings[[n]] = pheno
  n = n+1
}
rand_score = c()
for (n in 1:((length(clusterings)) -1) ){
  print(n)
  score = genieclust::adjusted_rand_score(clusterings[[n]], clusterings[[n+1]])
  rand_score = c(rand_score, score)
}
rand_score
Rphenograph_out <- Rphenograph(data, k = 30) 
pheno = membership(Rphenograph_out[[2]]) 
lym$pheno_cluster = pheno 
lym$barcode <- colnames(lym)

lym$subtype = factor(lym$subtype, levels= c('CSPC','CRPC','NEPC'))
saveRDS(lym,'/Users/jooyoung/Dropbox/msk.integrated.remove.cellcycle.lym.subtype.rds')

## Stromal cell ####
msk = readRDS('msk.integrated.remove.cellcycle.jy.rds')
strom = subset(msk, seurat_clusters %in% c(13,9,16,37,8,38,24,11)); rm(msk)
strom = DietSeurat(strom)
strom = NormalizeData(strom)
s.genes <- cc.genes$s.genes ; g2m.genes <- cc.genes$g2m.genes; strom <- CellCycleScoring(strom, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE); strom <- ScaleData(strom, features = rownames(strom),vars.to.regress = c("S.Score", "G2M.Score"))
strom <- FindVariableFeatures(strom, selection.method = "vst", nfeatures = 2000)
strom <- RunPCA(strom)
ElbowPlot(strom, ndims = 50)
strom=DietSeurat(strom)
strom = RunFastMNN(object.list = SplitObject(strom, split.by = "subtype")) 

data <- as.matrix(strom@reductions$mnn@cell.embeddings[,1:30])
clusterings = list()
n = 1
for (k in seq(20,50,5)){
  Rphenograph_out <- Rphenograph(data, k = k)
  pheno = membership(Rphenograph_out[[2]]) %>% as.factor()
  names(pheno) = rownames(strom@reductions$mnn@cell.embeddings)
  clusterings[[n]] = pheno
  n = n+1
}
rand_score = c()
for (n in 1:((length(clusterings)) -1) ){
  print(n)
  score = genieclust::adjusted_rand_score(clusterings[[n]], clusterings[[n+1]])
  rand_score = c(rand_score, score)
}
rand_score
Rphenograph_out <- Rphenograph(data, k = 30) 
pheno = membership(Rphenograph_out[[2]]) 
strom$pheno_cluster = pheno 
strom$barcode <- colnames(strom)

strom$subtype = factor(strom$subtype, levels= c('CSPC','CRPC','NEPC'))
saveRDS(strom,'/Users/jooyoung/Dropbox/msk.integrated.remove.cellcycle.strom.subtype.rds')


## Tumor cells ####
msk = readRDS('/Users/jooyoung/OneDrive - 고려대학교/samir/rds_jy/msk.integrated.remove.cellcycle.jy.rds')
epi = subset(msk, seurat_clusters %in% c(13,9,3,42,34,25,16,37,8,38,24,11,35,1,28,15,10,31,29,22), invert = T)
epi = DietSeurat(epi)
epi = NormalizeData(epi)
s.genes <- cc.genes$s.genes ; g2m.genes <- cc.genes$g2m.genes; epi <- CellCycleScoring(epi, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE); epi <- ScaleData(epi, features = rownames(epi),vars.to.regress = c("S.Score", "G2M.Score"))
epi <- FindVariableFeatures(epi, selection.method = "vst", nfeatures = 2000)
HVG = VariableFeatures(object = epi)
genes = readxl::read_xlsx('/Users/jooyoung/OneDrive - 고려대학교/samir/rds/Table S18.JAK_STAT_FGFR_Misc_Signatures.xlsx'); genes = genes %>% as.list() %>% unlist(); genes = genes[! genes %>% is.na() ];genes[genes %>% str_detect('orf')];genes[genes == 'Cl orf172'] = 'Clorf172';genes[genes == 'Fl 1R'] = 	'FOLR1';genes[genes == 'FL 1R'] = 	'FOLR1';genes[genes == 'IL11RA1'] = 	'IL11RA';genes[genes == 'TGIF'] = 	'TGIF1';genes[genes == 'GPR56'] = 'ADGRG1';genes[genes == 'GPR110'] = 'ADGRF1';genes[genes == 'INADL'] ='PATJ';genes[genes == "TACSTDI"] ='TACSTD2';genes[genes == "TMEM30B"] ='TMEM3OB';genes[genes == "RBM35A"] = 'ESRP1';genes[genes == "MTAC2D1"] ='TC2N';genes[!genes %in% rownames(epi)];genes = genes[genes %in% rownames(epi)]
genes_hvg = c(genes , HVG) %>% unique()
epi <- RunPCA(epi, features = genes_hvg)
epi <- FindNeighbors(epi, dims = 1:30); gc()
epi <- FindClusters(epi, resolution = 0.5); gc()
epi <- RunUMAP(epi, dims = 1:30)

saveRDS(epi,'/Users/jooyoung/OneDrive - 고려대학교/samir/rds_jy/msk.integrated.remove.cellcycle.epi.rds')

markers = c('EPCAM','YAP1','EZH2','AR','NKX3-1','KLK2','KLK3','FOLH1','STEAP1','STEAP2','TACSTD2','ETV6','ETV1','DLL3','CHGA','CHGB','NEUROD1','ASCL1')
for (i in markers){
  a = FeaturePlot_scCustom(epi,order = F,aspect_ratio = 1, i) +scale_colour_gradientn(colours = rev(brewer.pal(n = 11, 'Spectral'))) # , limits = c(0, 4),oob = scales::squish
  b = VlnPlot_scCustom(epi,i) + theme(axis.text.x = element_text(size = 13,angle = 45))  + theme(plot.title = element_text(size = 25))
  c = cowplot::plot_grid(a,b, rel_widths = c(1,1.4))
  assign(paste0('g_',i),c)
}
pic = ls()[ls() %>% startsWith("g_")]
names = c()
for(l in pic){names = append(l %>% strsplit('g_') %>%unlist() %>% .[2],names)}
names = sort(names)
plot_list = list()
for (k in 1:length(pic )){
  print(pic[k])
  plot_list[[k]] = get(pic [k])
}
idx = match(markers, names)
plot_grid_args = c(plot_list[idx],list(ncol = 3))
d =   do.call(plot_grid,plot_grid_args)

pdf('20231017_Tumormarker_exp.pdf', width = 11*2*3, height = 6.5*1.3*ceiling(length(markers)/3))
print(d)
dev.off()

tumor = subset(epi, seurat_clusters %in% c(0,4,3,15,14,8,23,24,19,10,5), invert = T); rm(msk); gc() 
DimPlot(tumor)
tumor=DietSeurat(tumor)
tumor = RunFastMNN(object.list = SplitObject(tumor, split.by = "subtype"))
tumor <- RunUMAP(tumor, dim = 1:30, reduction = 'mnn')
tumor <- FindNeighbors(tumor, dim = 1:30, reduction = 'mnn'); gc()
tumor <- FindClusters(tumor, resolution = 0.2); gc()
tumor$subtype = factor(tumor$subtype, levels = c('CSPC','CRPC','NEPC'))
saveRDS(tumor,'/Users/jooyoung/OneDrive - 고려대학교/samir/rds_jy/msk.integrated.remove.cellcycle.tumor.fastmnn.subtype.rds')


## Remove Hepatocytes ####
tumor = readRDS('/Users/jooyoung/OneDrive - 고려대학교/samir/rds_jy/msk.integrated.remove.cellcycle.tumor.fastmnn.subtype.rds')
markers = c('CRP', 'VTN', 'A1CF', 'GC','AMBP','ALB','YAP1','EZH2','AR','NKX3-1','KLK2','KLK3','FOLH1','STEAP1','STEAP2','TACSTD2','ETV6','ETV1','DLL3','CHGA','CHGB','NEUROD1','ASCL1')

tumor = subset(tumor, pheno_cluster == 6, invert = T)
tumor = DietSeurat(tumor)
s.genes <- cc.genes$s.genes ; g2m.genes <- cc.genes$g2m.genes; tumor <- CellCycleScoring(tumor, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE); tumor <- ScaleData(tumor, features = rownames(tumor),vars.to.regress = c("S.Score", "G2M.Score"))
gc()
tumor <- FindVariableFeatures(tumor, selection.method = "vst", nfeatures = 2000)
HVG = VariableFeatures(object = tumor)
genes = readxl::read_xlsx('/Users/jooyoung/OneDrive - 고려대학교/samir/rds/Table S18.JAK_STAT_FGFR_Misc_Signatures.xlsx'); genes = genes %>% as.list() %>% unlist(); genes = genes[! genes %>% is.na() ];genes[genes %>% str_detect('orf')];genes[genes == 'Cl orf172'] = 'Clorf172';genes[genes == 'Fl 1R'] = 	'FOLR1';genes[genes == 'FL 1R'] = 	'FOLR1';genes[genes == 'IL11RA1'] = 	'IL11RA';genes[genes == 'TGIF'] = 	'TGIF1';genes[genes == 'GPR56'] = 'ADGRG1';genes[genes == 'GPR110'] = 'ADGRF1';genes[genes == 'INADL'] ='PATJ';genes[genes == "TACSTDI"] ='TACSTD2';genes[genes == "TMEM30B"] ='TMEM3OB';genes[genes == "RBM35A"] = 'ESRP1';genes[genes == "MTAC2D1"] ='TC2N';genes[!genes %in% rownames(tumor)];genes = genes[genes %in% rownames(tumor)]
genes_hvg = c(genes , HVG) %>% unique()
tumor <- RunPCA(tumor, features = genes_hvg)
rm(HVG); rm(genes)
ElbowPlot(tumor, ndims = 50) #20

tumor = DietSeurat(tumor)
tumor = RunFastMNN(object.list = SplitObject(tumor, split.by = "subtype"))
tumor$subtype = factor(tumor$subtype, levels = c('CSPC','CRPC','NEPC'))

saveRDS(msk,'msk.integrated.remove.cellcycle.tumor.cells.rds')


## Phenograph clustering for tumor cells -  ####
tumor = readRDS('msk.integrated.remove.cellcycle.tumor.cells.rds')
data <- as.matrix(tumor@reductions$mnn@cell.embeddings[,1:20])
clusterings = list()
n = 1
for (k in seq(20,50,5)){
  Rphenograph_out <- Rphenograph(data, k = k)
  pheno = membership(Rphenograph_out[[2]]) %>% as.factor()
  names(pheno) = rownames(tumor@reductions$mnn@cell.embeddings)
  clusterings[[n]] = pheno
  n = n+1
}
rand_score = c()
for (n in 1:((length(clusterings)) -1) ){
  print(n)
  score = genieclust::adjusted_rand_score(clusterings[[n]], clusterings[[n+1]])
  rand_score = c(rand_score, score)
}
rand_score

Rphenograph_out <- Rphenograph(data, k = 30) 
pheno = membership(Rphenograph_out[[2]]) %>% as.factor()
tumor$pheno_cluster = pheno
tumor$barcode <- colnames(tumor)
saveRDS(tumor, 'msk.integrated.remove.cellcycle.tumor.cells.rds')

library(Matrix)
df = Matrix(tumor@assays$RNA@data, sparse = TRUE) 
writeMM(df, file='tumor_fastmnn_subtype_nohepa.mtx')
write.table(data.frame('gene'=rownames(df)),file='tumor_gene_names.csv',quote=F,row.names=F,col.names=F)
write.csv(tumor@meta.data, file='tumor_fastmnn_subtype_nohepa_meta.csv', quote=F, row.names=F)
write.csv(tumor@reductions$mnn@cell.embeddings, file='tumor_fastmnn_subtype_nohepa_mnn.csv', quote=F, row.names=F) 

# UMAP with python 
umap = read.csv('tumor_fastmnn_subtype_nohepa_pagaumap.csv')
rownames(umap) = umap$barcode
umap$barcode = NULL
colnames(umap) = c('X1','X2')

which(rownames(umap) != colnames(tumor))
tumor[['umap']] = CreateDimReducObject(embeddings = umap %>% as.matrix(), key = 'UMAP_', assay = 'RNA')

pdf('20231027_tumorpheno_Subtype_fastmnn_subtype_nohepa.pdf',width = 5.7 , height = 5)
DimPlot_scCustom(tumor, reduction = 'umap', group.by = 'subtype' , raster = FALSE, colors_use = c('gold','blue3','red3'))
dev.off()

# Figure 2C ####
shannon_entropy <- function(p) {
  if (p == 0){
    return(0)
  } else {
    return(  -(p * log(p)))
  }
}

## Cancer Entropy ####
cancer = readRDS('msk.integrated.remove.cellcycle.tumor.cells.rds')
cancer$subtype

CSPC = c()
CRPC = c()
NEPC = c()

for (j in c(cancer$subtype %>% unique)){
  print(j)
  seuraa = subset(cancer, subset = subtype == j)
  hcs = c()
  names = c()
  for (i in c(unique(seuraa$pheno_cluster)) %>% as.character()){
    seura = subset(seuraa , subset = pheno_cluster == i )
    if (dim(seura)[2] > 100){
      print(i)
      names = c(names,rep(i,100))
      meta = seura@meta.data
      for (n in 1:100){
        cells = sample(colnames(seura), 100,replace = F )
        df = meta[cells,]
        freq = (table(df$patient)/100) %>% as.data.frame.array()
        hcs = c(hcs, apply(freq, 1, shannon_entropy) %>% sum )
      }
    }
  }
  
  if (j == 'CSPC'){
    CSPC = hcs %>% as.data.frame(); CSPC$cluster = names
  } else if (j == 'CRPC'){
    CRPC = hcs %>% as.data.frame(); CRPC$cluster = names
  } else {
    NEPC = hcs %>% as.data.frame(); NEPC$cluster = names
  }
}


CSPC$subtype = 'CSPC' 
CRPC$subtype = 'CRPC' 
NEPC$subtype = 'NEPC' 

df = rbind(CSPC,rbind(CRPC,NEPC))
df$subtype = factor(df$subtype, c('CSPC','CRPC','NEPC'))
colnames(df)[1] = 'HC'

saveRDS(df,'cancer_fastmnn_subtype_nohepa_df.rds')

## NonCancer Entropy ####
### Myeloid ####
mye = readRDS('../msk.integrated.remove.cellcycle.mye.subtype.rds')
h = 'mye'
hcs = c()
for (j in mye$pheno_cluster %>% unique){
  seura = subset(mye , subset = pheno_cluster == j)
  
  if (dim(seura)[2] > 100){
    print(j)
    meta = seura@meta.data
    
    for (n in 1:100){
      cells = sample(colnames(seura), 100 ,replace = T)
      df = meta[cells,]
      freq = (table(df$patient)/100) %>% as.data.frame.array()
      hcs = c(hcs, apply(freq, 1, shannon_entropy) %>% sum )
    }
  }
}

hcs = hcs %>% as.data.frame()
colnames(hcs) = 'HC'
hcs$subtype = h

df = hcs
saveRDS(df,'mye_df')

### Lymphoid ####
lym = readRDS('../msk.integrated.remove.cellcycle.lym.subtype.rds')

h = 'lym'
hcs = c()
for (j in lym$pheno_cluster %>% unique){
  seura = subset(lym , subset = pheno_cluster == j)
  
  if (dim(seura)[2] > 100){
    print(j)
    meta = seura@meta.data
    
    for (n in 1:100){
      cells = sample(colnames(seura), 100 ,replace = T)
      df = meta[cells,]
      freq = (table(df$patient)/100) %>% as.data.frame.array()
      hcs = c(hcs, apply(freq, 1, shannon_entropy) %>% sum )
    }
  }
}

hcs = hcs %>% as.data.frame()
colnames(hcs) = 'HC'
hcs$subtype = h

df = hcs
saveRDS(df,'lym_df.rds')

### Stromal ####
strom = readRDS('../msk.integrated.remove.cellcycle.strom.subtype.rds')

h = 'stromal'
hcs = c()
for (j in strom$pheno_cluster %>% unique){
  seura = subset(strom , subset = pheno_cluster == j)
  
  if (dim(seura)[2] > 100){
    print(j)
    meta = seura@meta.data
    
    for (n in 1:100){
      cells = sample(colnames(seura), 100 ,replace = T)
      df = meta[cells,]
      freq = (table(df$patient)/100) %>% as.data.frame.array()
      hcs = c(hcs, apply(freq, 1, shannon_entropy) %>% sum )
    }
  }
}

hcs = hcs %>% as.data.frame()
colnames(hcs) = 'HC'
hcs$subtype = h

df = hcs
saveRDS(df,'strom_df.rds')

## Figure 2C, Sup Fig 3H,F ####
noncancer1 = readRDS('strom_df.rds')
noncancer2= readRDS('lym_df.rds')
noncancer3= readRDS('mye_df.rds')
cancer = readRDS('cancer_fastmnn_subtype_nohepa_df.rds')

df = rbind(cancer, noncancer1,noncancer2, noncancer3)
df %>% colnames()
df$subtype

df$subtype = factor(df$subtype, c('CSPC', 'CRPC' ,'NEPC', 'stromal','mye','lym'))

p_label <- function(p_value) {
  if (p_value < 0.0001) {
    return('***')
  } else if (p_value < 0.001){
    return('**')
  }else {
    return(paste0("p = ", p_value))
  }
}

stat.test <- df %>%
  wilcox_test(HC ~ subtype, p.adjust.method = 'bonferroni', alternative = 'greater') %>%
  add_significance() %>%
  add_xy_position(x="condition",step.increase=1) # , alternative = 'greater'

stat.test = stat.test[c(1,2,6),]

stat.test$manual_position <-1.1 * c(1.0,1.13,1.02)
stat.test$label <- stat.test$p.adj.signif

pdf('20231222_shannon_subtype_corrected.pdf',width = 4.5*1.2,height = 5*1.2)
ggplot(df, aes(x = subtype, y = HC)) + geom_violin(aes(fill = subtype))+geom_boxplot( aes(fill = subtype),alpha = 0.3, outlier.shape = NA)   + theme_classic() + NoLegend()  + ylab('Shannon Entropy')+xlab('Subtype') + scale_fill_manual( values = c('darkgoldenrod2', '#2171B5','#CB181D','lightgray','lightgray','lightgray')) + ggsignif::geom_signif(data=stat.test,aes(xmin=group1,xmax=group2,annotations=label,y_position=manual_position),manual=TRUE, inherit.aes=FALSE)
dev.off()


# Sup Fig 3E ####
library(genefu)
library(org.Hs.eg.db)
library(AnnotationDbi)

data(scmod2.robust)
data(pam50.robust)
data(scmgene.robust)
data(sig.ggi)
data(scmod1.robust)
data(sig.genius)

## Get data ####
cancer = readRDS('msk.integrated.remove.cellcycle.tumor.cells.rds')
M = cancer@assays$RNA@data
genes = rownames(cancer) %>% as.data.frame()
rownames(genes) = genes$.
colnames(genes) = 'probe'
genes$NCBI.gene.symbol  = rownames(genes)

## Map geneID ####
entrez_ids <- select(org.Hs.eg.db, keys=genes$NCBI.gene.symbol, keytype="SYMBOL", columns="ENTREZID")
entrez_ids = entrez_ids[! (entrez_ids$SYMBOL %>% duplicated() ),]
rownames(entrez_ids) = entrez_ids$SYMBOL
entrez_ids = entrez_ids[genes$NCBI.gene.symbol,]
genes$EntrezGene.ID = entrez_ids$ENTREZID

## Start learning ####
verbose = F
centroids <- sbt.model$centroids
centroids = centroids[,c("Basal"  ,  "LumA",   "LumB")]
annot.centroids <- sbt.model$centroids.map
method.cor <- sbt.model$method.cor
method.centroids <- sbt.model$method.centroids
std <- sbt.model$std
mq <- sbt.model$rescale.q
mins <- sbt.model$mins

number.cluster <- ncol(centroids)

if(is.null(dimnames(centroids)[[2]])) { name.cluster <- paste("cluster", 1:ncol(centroids), sep=".") } else { name.cluster <- dimnames(centroids)[[2]] }

gt <- nrow(centroids)

centroids.gid <- as.character(annot.centroids[ ,"EntrezGene.ID"])
names(centroids.gid) <- as.character(annot.centroids[ , "probe"])
myx <- !duplicated(centroids.gid) & !is.na(centroids.gid)
centroids.gid <- centroids.gid[myx]
annot.centroids <- annot.centroids[myx, , drop=FALSE]
centroids <- centroids[myx, , drop=FALSE]
gid <- as.character(annot[ ,"EntrezGene.ID"])
names(gid) <- as.character(annot[ ,"probe"])

rr <- geneid.map(geneid1=gid, data1=data, geneid2=centroids.gid, verbose=FALSE)
nn <- match(rr$geneid2, centroids.gid)
nn <- nn[!is.na(nn)]
centroids.gid <- centroids.gid[nn]
annot.centroids <- annot.centroids[nn, ]
centroids <- centroids[nn, , drop=FALSE]
data <- rr$data1

centroids.map <- cbind("probe"=dimnames(data)[[2]], "probe.centroids"=dimnames(centroids)[[1]], "EntrezGene.ID"=as.character(annot[dimnames(data)[[2]], "EntrezGene.ID"]))
dimnames(centroids.map)[[1]] <- dimnames(data)[[2]]
gm <- nrow(centroids)

switch(std,
       "scale"={
         data <- scale(data, center=TRUE, scale=TRUE)
         if(verbose) { message("standardization of the gene expressions") }
       }, 
       "robust"={
         data <- apply(data, 2, function(x) { return((rescale(x, q=mq, na.rm=TRUE) - 0.5) * 2) })
         if(verbose) { message("robust standardization of the gene expressions") }
       }, 
       "none"={ if(verbose) { message("no standardization of the gene expressions") } })

## apply the nearest centroid classifier to classify the samples again
ncor <- t(apply(X=data, MARGIN=1, FUN=function(x, y, method.cor) { 
  rr <- array(NA, dim=ncol(y), dimnames=list(colnames(y)))
  if (sum(complete.cases(x, y)) > 3) {
    rr <- cor(x=x, y=y, method=method.cor, use="complete.obs")
  }
  return (rr)
}, y=centroids, method.cor=method.cor))

nproba <- t(apply(X=ncor, MARGIN=1, FUN=function (x) {
  rr <- array(NA, dim=length(x), dimnames=list(names(x)))
  x[!is.na(x) & x < 0] <- 0
  if (!all(is.na(x))) {
    rr <- x / sum(x, na.rm=TRUE)
  }
  return (rr)
}))

dimnames(ncor) <- dimnames(nproba) <- list(dimnames(data)[[1]], name.cluster)
ncl <- apply(X=ncor, MARGIN=1, FUN=function(x) { return(order(x, decreasing=TRUE, na.last=TRUE)[1]) })
names(ncl) <- dimnames(data)[[1]]

ncln <- name.cluster[ncl]
names(ncln) <- dimnames(data)[[1]]

myx <- sort(unique(ncl))
myx <- myx[!is.na(myx)]
name.cluster2 <- name.cluster[myx]
number.cluster2 <- length(myx)
ps.res <- ncl2 <- NULL
if(do.prediction.strength) {
  ## compute the clustering and cut the dendrogram
  ## hierarchical clustering with correlation-based distance and average linkage
  hcl <- amap::hcluster(x=data, method="correlation", link="average")
  mins.ok <- stop.ok <- FALSE
  nbc <- number.cluster2
  nclust.best <- 1
  while(!mins.ok && !stop.ok) { ## until each cluster contains at least mins samples
    cl <- cutree(tree=hcl, k=nbc)
    tt <- table(cl)
    if(sum(tt >= mins) >= number.cluster2) {
      if(nbc > number.cluster2) { ## put NA for clusters with less than mins samples
        td <- names(tt)[tt < mins]
        cl[is.element(cl, td)] <- NA
        ## rename the clusters
        ucl <- sort(unique(cl))
        ucl <- ucl[!is.na(ucl)]
        cl2 <- cl
        for(i in 1:number.cluster2) { cl2[cl == ucl[i] & !is.na(cl)] <- i }
        cl <- cl2
      }
      nclust.best <- number.cluster2
      mins.ok <- TRUE
    } else {
      if(sum(tt >= mins) > nclust.best) {
        nbc.best <- nbc
        nclust.best <- sum(tt >= mins)
      }
      nbc <- nbc + 1
      if(nbc > (nrow(data) - (number.cluster2 * mins))) {
        warning(sprintf("impossible to find %i main clusters with at least %i individuals!", number.cluster2, mins))
        stop.ok <- TRUE
      }
    }
    if(stop.ok) { ## no convergence for the clustering with mininmum set of individuals
      cl <- cutree(tree=hcl, k=nbc.best)
      tt <- table(cl)
      td <- names(tt)[tt < mins]
      cl[is.element(cl, td)] <- NA
      ## rename the clusters
      ucl <- sort(unique(cl))
      ucl <- ucl[!is.na(ucl)]
      cl2 <- cl
      for(i in 1:nclust.best) { cl2[cl == ucl[i] & !is.na(cl)] <- i }
      cl <- cl2
    }
  }
  ## compute the centroids
  ## take the core samples in each cluster to compute the centroid
  ## not feasible due to low intra correlation within clusters!!!
  ## minimal pairwise cor of approx 0.3
  #cl2 <- cutree(tree=hcl, h=0.7)
  #table(cl, cl2) to detect which core cluster of samples for which cluster.
  cl.centroids <- matrix(NA, nrow=ncol(data), ncol=nclust.best, dimnames=list(dimnames(data)[[2]], paste("cluster", 1:nclust.best, sep=".")))
  for(i in 1:ncol(cl.centroids)) {
    switch(method.centroids, 
           "mean"={ cl.centroids[ ,i] <- apply(X=data[cl == i & !is.na(cl), ,drop=FALSE], MARGIN=2, FUN=mean, na.rm=TRUE, trim=0.025) }, 
           "median"={ cl.centroids[ ,i] <- apply(X=data[cl == i & !is.na(cl), ,drop=FALSE], MARGIN=2, FUN=median, na.rm=TRUE) }, 
           "tukey"={ cl.centroids[ ,i] <- apply(X=data[cl == i & !is.na(cl), ,drop=FALSE], MARGIN=2, FUN=tbrm, na.rm=TRUE, C=9) })
  }
  #apply the nearest centroid classifier to classify the samples again
  ncor2 <- t(apply(X=data, MARGIN=1, FUN=function(x, y, z) { return(cor(x, y, method=z, use="complete.obs")) }, y=cl.centroids, z=method.cor))
  nproba2 <- t(apply(X=ncor2, MARGIN=1, FUN=function(x) { return(abs(x) / sum(abs(x), na.rm=TRUE)) }))
  dimnames(ncor2) <- dimnames(nproba2) <- list(dimnames(data)[[1]], dimnames(cl.centroids)[[2]])
  ncl2 <- apply(X=ncor2, MARGIN=1, FUN=function(x) { return(order(x, decreasing=TRUE)[1]) })
  names(ncl2) <- dimnames(data)[[1]]
  ## rename clusters since we do not expect to get the same id per cluster
  ## this avoids a warning in ps.cluster
  uncl <- sort(unique(ncl))
  uncl <- uncl[!is.na(uncl)]
  nclt <- ncl
  for(mm in 1:length(uncl)) {
    nclt[ncl == uncl[mm]] <- mm
  }
  uncl2 <- sort(unique(ncl2))
  uncl2 <- uncl2[!is.na(uncl2)]
  ncl2t <- ncl2
  for(mm in 1:length(uncl2)) {
    ncl2t[ncl2 == uncl2[mm]] <- mm
  }
  #prediction strength
  ps.res <- ps.cluster(cl.tr=ncl2t, cl.ts=nclt, na.rm=TRUE)
  ## put NA for clusters which are potentially not present in the dataset
  tt <- rep(NA, length(name.cluster))
  names(tt) <- name.cluster
  tt[name.cluster2] <- ps.res$ps.cluster
  ps.res$ps.cluster <- tt
}

df = list("subtype"=ncln, "subtype.proba"=nproba, "cor"=ncor, "prediction.strength"=ps.res, "subtype.train"=ncl2, "profiles"=data, "centroids.map"=centroids.map)
cancer$PAM50_pseudo_nonormher2 = ifelse(cancer$pheno_cluster %in% c(2,4,12,19),'Basal',ifelse(cancer$pheno_cluster %in% c(10,16,17,29,30,31), 'LumB','LumA'))


pdf('20240318_PAM50_pseudo_nonormher.pdf', width = 8, height = 8)
DimPlot_scCustom(cancer, group.by = 'PAM50_pseudo_nonormher2',aspect_ratio = 1)
dev.off()

# Sup Fig 3G ####
library(Seurat)
library(dplyr)
library(vegan)

pc = readRDS('/Users/jooyoung/Dropbox/samir/rds_jy/msk.integrated.remove.cellcycle.lym.mnn.subtype.rds')
pc <- FindNeighbors(pc, reduction = "mnn", dims = 1:30, k.param = 30, return.neighbor=TRUE)

shannon <- setNames(data.frame(matrix(ncol = 2, nrow = 0)), c("cell_id","shannon"))
for (cell in 1:length(colnames(pc))) {
  top30 <- TopNeighbors(pc[["RNA.nn"]],cell=colnames(pc)[cell],n=30)
  vector <- table(pc@meta.data[top30,'patient']) %>% as.vector()
  fraction <- vector/sum(vector)
  H <- diversity(fraction)
  shannon[cell,'cell_id'] <- colnames(pc)[cell]
  shannon[cell,'shannon'] <- H
}
row.names(shannon) <- shannon$cell_id
shannon$cell_id <- NULL
pc = AddMetaData(pc, metadata = shannon)

saveRDS(pc, '/Users/jooyoung/Dropbox/samir/rds_jy/msk.integrated.remove.cellcycle.lym.mnn.subtype.rds')


pc = readRDS('/Users/jooyoung/Dropbox/samir/rds_jy/msk.integrated.remove.cellcycle.mye.mnn.subtype.rds')
pc <- FindNeighbors(pc, reduction = "mnn", dims = 1:30, k.param = 30, return.neighbor=TRUE)
shannon <- setNames(data.frame(matrix(ncol = 2, nrow = 0)), c("cell_id","shannon"))
for (cell in 1:length(colnames(pc))) {
  top30 <- TopNeighbors(pc[["RNA.nn"]],cell=colnames(pc)[cell],n=30)
  vector <- table(pc@meta.data[top30,'patient']) %>% as.vector()
  fraction <- vector/sum(vector)
  H <- diversity(fraction)
  shannon[cell,'cell_id'] <- colnames(pc)[cell]
  shannon[cell,'shannon'] <- H
}
row.names(shannon) <- shannon$cell_id
shannon$cell_id <- NULL
pc = AddMetaData(pc, metadata = shannon)

saveRDS(pc, '/Users/jooyoung/Dropbox/samir/rds_jy/msk.integrated.remove.cellcycle.mye.mnn.subtype.rds')


pc = readRDS('/Users/jooyoung/Dropbox/samir/rds_jy/msk.integrated.remove.cellcycle.strom.mnn.subtype.rds')
pc <- FindNeighbors(pc, reduction = "mnn", dims = 1:30, k.param = 30, return.neighbor=TRUE)
shannon <- setNames(data.frame(matrix(ncol = 2, nrow = 0)), c("cell_id","shannon"))
for (cell in 1:length(colnames(pc))) {
  top30 <- TopNeighbors(pc[["RNA.nn"]],cell=colnames(pc)[cell],n=30)
  vector <- table(pc@meta.data[top30,'patient']) %>% as.vector()
  fraction <- vector/sum(vector)
  H <- diversity(fraction)
  shannon[cell,'cell_id'] <- colnames(pc)[cell]
  shannon[cell,'shannon'] <- H
}
row.names(shannon) <- shannon$cell_id
shannon$cell_id <- NULL
pc = AddMetaData(pc, metadata = shannon)
saveRDS(pc, '/Users/jooyoung/Dropbox/samir/rds_jy/msk.integrated.remove.cellcycle.strom.mnn.subtype.rds')

pc = readRDS('/Users/jooyoung/Dropbox/samir/rds_jy/msk.integrated.remove.cellcycle.tumor.cells.rds')
pc <- FindNeighbors(pc, reduction = "mnn", dims = 1:20, k.param = 30, return.neighbor=TRUE)
shannon <- setNames(data.frame(matrix(ncol = 2, nrow = 0)), c("cell_id","shannon"))
for (cell in 1:length(colnames(pc))) {
  top30 <- TopNeighbors(pc[["RNA.nn"]],cell=colnames(pc)[cell],n=30)
  vector <- table(pc@meta.data[top30,'patient']) %>% as.vector()
  fraction <- vector/sum(vector)
  H <- diversity(fraction)
  shannon[cell,'cell_id'] <- colnames(pc)[cell]
  shannon[cell,'shannon'] <- H
}
row.names(shannon) <- shannon$cell_id
shannon$cell_id <- NULL
pc = AddMetaData(pc, metadata = shannon)
saveRDS(pc, '/Users/jooyoung/Dropbox/samir/rds_jy/msk.integrated.remove.cellcycle.tumor.cells.rds')

