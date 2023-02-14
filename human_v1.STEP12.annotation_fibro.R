####--------------------------------------------------------------------------------------------------------------------
####  Pediatric Acute Rejection
####  2022-08-03 by Xiao LI (Texas Heart Institute, US)
####--------------------------------------------------------------------------------------------------------------------


####--------------------------------------------------------------------------------------------------------------------
####  Initiate  ####
####--------------------------------------------------------------------------------------------------------------------
Ver <- '1'
Step <- 'STEP12_Annotation_FB'
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
fb.srt <- srt[, srt$Cell_type %in% c('FB')]

DimPlot2(fb.srt, group.by = 'sample')

## Re-embed by scVI
fb.srt <- RunUMAP(fb.srt, reduction = 'scANVI', dims = 1:25,
                  reduction.name = 'sub_scANVI_umap', reduction.key = 'subSCANVIUMAP_')
fb.srt <- FindNeighbors(fb.srt, dims = 1:25, reduction = 'scANVI', force.recalc = T) %>%
        FindClusters(resolution = c(0.1, 0.2, 0.5))

DimPlot2(fb.srt, reduction = 'sub_scANVI_umap', group.by = 'name2', cols = mycol_10)
DimPlot2(fb.srt, reduction = 'sub_scANVI_umap', group.by = 'Doublet_SC')

fb.srt$Cell_state_previous[fb.srt$Cell_state_previous %in%
                                   names(table(fb.srt$Cell_state_previous))[table(fb.srt$Cell_state_previous)<20]] <- NA

DimPlot2(fb.srt, reduction = 'sub_scANVI_umap', group.by = 'Cell_state_previous', label = T, cols  = mycol_20) +
        DimPlot2(fb.srt, reduction = 'sub_scANVI_umap', group.by = 'CBN_snn_res.0.1', label = T)+
        DimPlot2(fb.srt, reduction = 'sub_scANVI_umap', group.by = 'CBN_snn_res.0.2', label = T)+
        DimPlot2(fb.srt, reduction = 'sub_scANVI_umap', group.by = 'CBN_snn_res.0.5', label = T)


FeaturePlot2(fb.srt, features = markers_lvl1_min, reduction = 'sub_scANVI_umap')

## Annotate
Idents(fb.srt) <- 'CBN_snn_res.0.5'
mk <- FindAllMarkers(fb.srt, only.pos = T, return.thresh = 0.01)
p <- MarkerHeatmap(fb.srt, mk, n_cells = 400, top = 20)
PlotPDF('01.1.heat.maker_all_cluster', 15, 15)
p
dev.off()

PlotPDF('01.2.dim.all_cluster', 10, 10)
DimPlot2(fb.srt, reduction = 'sub_scANVI_umap', label = T, cols = mycol_20)
dev.off()

PlotPDF('01.3.bar.doublet_in_all_cluster', 5, 5)
CountCellBarPlot(fb.srt, group.var = 'CBN_snn_res.0.5', stack.var = 'Doublet_SC', percentage = T)
dev.off()

PlotPDF('01.4.heat.scANVI_mapping', 12, 12)
p1 <- MappingHeatmap(fb.srt, que_var = 'CBN_snn_res.0.5', ref_var = 'Cell_state_scANVI', ref.disp.min = 0.1)
p2 <- MappingHeatmap(fb.srt, que_var = 'CBN_snn_res.0.5', ref_var = 'Cell_state_SingleR', ref.disp.min = 0.1)
p1/p2
dev.off()

## Re-annotate
fb.srt$Cell_state <- NA
fb.srt$Cell_state[fb.srt$CBN_snn_res.0.5 %in% c(0, 2, 4)] <- 'FB1'
fb.srt$Cell_state[fb.srt$CBN_snn_res.0.5 %in% c(3)] <- 'FB2'
fb.srt$Cell_state[fb.srt$CBN_snn_res.0.5 %in% c(5)] <- 'FB3'
fb.srt$Cell_state[fb.srt$CBN_snn_res.0.5 %in% c(7, 10)] <- 'FB4'
fb.srt$Cell_state[fb.srt$CBN_snn_res.0.5 %in% c(1, 8, 6, 9, 11)] <- 'Doublet'

PlotPDF('02.dim.all_cluster', 10, 10)
DimPlot2(fb.srt, reduction = 'sub_scANVI_umap', label = T, cols = mycol_10, group.by = 'Cell_state')
dev.off()

Idents(fb.srt) <- 'Cell_state'
levels(fb.srt) <- c('FB1',
                    'FB2',
                    'FB3',
                    'FB4',
                    'Doublet')
mk <- FindAllMarkers(fb.srt, only.pos = T, return.thresh = 0.01, logfc.threshold = 0.5)
p <- MarkerHeatmap(fb.srt, mk, n_cells = 400, top = 20)
PlotPDF('03.heat.maker_all_cluster', 15, 15)
p
dev.off()

####--------------------------------------------------------------------------------------------------------------------
saveRDS(fb.srt, 'integrated/STEP12.fb_cells.srt.rds')
####--------------------------------------------------------------------------------------------------------------------
