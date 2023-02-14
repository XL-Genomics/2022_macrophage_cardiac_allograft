####--------------------------------------------------------------------------------------------------------------------
####  Pediatric Acute Rejection
####  2022-08-03 by Xiao LI (Texas Heart Institute, US)
####--------------------------------------------------------------------------------------------------------------------


####--------------------------------------------------------------------------------------------------------------------
####  Initiate  ####
####--------------------------------------------------------------------------------------------------------------------
Ver <- '1'
Step <- 'STEP11_Annotation_CM'
Project <- '2022_acute_rejection_dturaga'
Code_dir <- paste0('~/Documents/Bioinformatics/project/', Project, '/code/human_v', Ver, '/')

source(Sys.readlink(paste0(Code_dir, 'src/bioinformatics.R')))
source(Sys.readlink(paste0(Code_dir, 'src/scRNAseq.R')))
source(Sys.readlink(paste0(Code_dir, 'src/geneset_human.R')))
source(paste0(Code_dir, 'human_v', Ver, '.helper_functions.R'))

InitiateProject('Rivendell', Ver, Step, 'human', Project, 'ithil')
####--------------------------------------------------------------------------------------------------------------------


####--------------------------------------------------------------------------------------------------------------------
####  Load data  ####
####--------------------------------------------------------------------------------------------------------------------
srt <- readRDS('integrated/STEP07.scvi_integration.srt.rds')
meta <- readRDS('integrated/STEP09.annotated.srt_meta.rds')

srt <- AddMetaData(srt, metadata = meta)
####--------------------------------------------------------------------------------------------------------------------


####--------------------------------------------------------------------------------------------------------------------
####  Immune cell annotation  ####
####--------------------------------------------------------------------------------------------------------------------
cm.srt <- srt[, srt$Cell_type %in% c('CM')]

DimPlot2(cm.srt, group.by = 'sample')

## Re-embed by scVI
cm.srt <- RunUMAP(cm.srt, reduction = 'scANVI', dims = 1:30,
                  reduction.name = 'sub_scANVI_umap', reduction.key = 'subSCANVIUMAP_')
cm.srt <- FindNeighbors(cm.srt, dims = 1:30, reduction = 'scANVI', force.recalc = T) %>%
        FindClusters(resolution = c(0.1, 0.2, 0.5))

DimPlot2(cm.srt, reduction = 'sub_scANVI_umap', group.by = 'name2', cols = mycol_10)
DimPlot2(cm.srt, reduction = 'sub_scANVI_umap', group.by = 'Doublet_SC')
DimPlot2(cm.srt, reduction = 'sub_scANVI_umap', group.by = 'Cell_state_previous', label = T, cols = mycol_20) +
        DimPlot2(cm.srt, reduction = 'sub_scANVI_umap', group.by = 'CBN_snn_res.0.1', label = T) +
        DimPlot2(cm.srt, reduction = 'sub_scANVI_umap', group.by = 'CBN_snn_res.0.2', label = T) +
        DimPlot2(cm.srt, reduction = 'sub_scANVI_umap', group.by = 'CBN_snn_res.0.5', label = T)

## Annotate
Idents(cm.srt) <- 'CBN_snn_res.0.1'
mk <- FindAllMarkers(cm.srt, only.pos = T, return.thresh = 0.01, logfc.threshold = 0.5)
p <- MarkerHeatmap(cm.srt, mk, n_cells = 500, top = 10)
PlotPDF('01.1.heat.maker_all_cluster', 15, 15)
p
dev.off()

PlotPDF('01.2.dim.all_cluster', 10, 10)
DimPlot2(cm.srt, reduction = 'sub_scANVI_umap', label = T, cols = mycol_10)
dev.off()

PlotPDF('01.3.bar.doublet_in_all_cluster', 5, 5)
CountCellBarPlot(cm.srt, group.var = 'CBN_snn_res.0.1', stack.var = 'Doublet_SC', percentage = T)
dev.off()

PlotPDF('01.4.heat.scANVI_mapping', 12, 8)
p1 <- MappingHeatmap(cm.srt, que_var = 'CBN_snn_res.0.1', ref_var = 'Cell_state_scANVI')
p2 <- MappingHeatmap(cm.srt, que_var = 'CBN_snn_res.0.1', ref_var = 'Cell_state_SingleR')
p1/p2
dev.off()

## Re-annotate
cm.srt$Cell_state <- NA
cm.srt$Cell_state[cm.srt$CBN_snn_res.0.1%in% c(0, 2)] <- 'CM1'
cm.srt$Cell_state[cm.srt$CBN_snn_res.0.1 == 1] <- 'CM2'
cm.srt$Cell_state[cm.srt$CBN_snn_res.0.1 == 3] <- 'CM3'

PlotPDF('02.dim.all_cluster', 10, 10)
DimPlot2(cm.srt, reduction = 'sub_scANVI_umap', label = T, cols = mycol_10, group.by = 'Cell_state')
dev.off()

Idents(cm.srt) <- 'Cell_state'
levels(cm.srt) <- c('CM1', 'CM2', 'CM3')
mk <- FindAllMarkers(cm.srt, only.pos = T, return.thresh = 0.01, logfc.threshold = 0.5)
p <- MarkerHeatmap(cm.srt, mk, n_cells = 500, top = 10)
PlotPDF('03.heat.maker_all_cluster', 15, 15)
p
dev.off()
####--------------------------------------------------------------------------------------------------------------------
saveRDS(cm.srt, 'integrated/STEP11.cm_cells.srt.rds')
####--------------------------------------------------------------------------------------------------------------------

## Testing:
Yap_target1 <- as.vector(read.table('/Volumes/shire/data/other/yap_target_cm_atac_2019_DevCell.csv', header = T)[, 1])
Yap_target1 <- ConvertGeneList(Yap_target1, from = 'mouse', to = 'human')
cm.srt <- AddModuleScore2(cm.srt, features = list(Yap_target1), names = 'Yap_score', assay = 'CBN')
FeaturePlot2(cm.srt, features = 'Yap_score', reduction = 'sub_scANVI_umap')
CountCellBarPlot(cm.srt, group.var = 'sample', stack.var = 'Cell_state', percentage = T)$plot +
        CountCellBarPlot(cm.srt, group.var = 'sample', stack.var = 'Cell_state', percentage = F)$plot
