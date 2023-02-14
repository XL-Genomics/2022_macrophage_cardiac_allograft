####--------------------------------------------------------------------------------------------------------------------
####  Pediatric Acute Rejection
####  2022-08-03 by Xiao LI (Texas Heart Institute, US)
####--------------------------------------------------------------------------------------------------------------------


####--------------------------------------------------------------------------------------------------------------------
####  Initiate  ####
####--------------------------------------------------------------------------------------------------------------------
Ver <- '1'
Step <- 'STEP14_Annotation_PC'
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
pc.srt <- srt[, srt$Cell_type %in% c('PC')]

DimPlot2(pc.srt, group.by = 'sample')

## Re-embed by scANVI
pc.srt <- RunUMAP(pc.srt, reduction = 'scANVI', dims = 1:30, reduction.name = 'sub_scANVI_umap', reduction.key = 'subSCANVIUMAP_')
pc.srt <- FindNeighbors(pc.srt, dims = 1:30, reduction = 'scANVI', force.recalc = T) %>%
        FindClusters(resolution = c(0.1, 0.2, 0.5))

DimPlot2(pc.srt, reduction = 'sub_scANVI_umap', group.by = 'name2', cols = mycol_10)
DimPlot2(pc.srt, reduction = 'sub_scANVI_umap', group.by = 'Doublet_SC')

pc.srt$Cell_state_previous[pc.srt$Cell_state_previous %in%
                                   names(table(pc.srt$Cell_state_previous))[table(pc.srt$Cell_state_previous)<20]] <- NA

DimPlot2(pc.srt, reduction = 'sub_scANVI_umap', group.by = 'Cell_state_previous', label = T, cols = mycol_20)+
DimPlot2(pc.srt, reduction = 'sub_scANVI_umap', group.by = 'CBN_snn_res.0.1', label = T)+
DimPlot2(pc.srt, reduction = 'sub_scANVI_umap', group.by = 'CBN_snn_res.0.2', label = T)+
DimPlot2(pc.srt, reduction = 'sub_scANVI_umap', group.by = 'CBN_snn_res.0.5', label = T)

## Annotate
Idents(pc.srt) <- 'CBN_snn_res.0.5'
mk <- FindAllMarkers(pc.srt, only.pos = T, return.thresh = 0.01, logfc.threshold = 0.5)
p <- MarkerHeatmap(pc.srt, mk, n_cells = 500, top = 10)
PlotPDF('01.1.heat.maker_all_cluster', 15, 15)
p
dev.off()

PlotPDF('01.2.dim.all_cluster', 10, 10)
DimPlot2(pc.srt, reduction = 'sub_scANVI_umap', label = T, cols = mycol_10)
dev.off()

PlotPDF('01.3.bar.doublet_in_all_cluster', 5, 5)
CountCellBarPlot(pc.srt, group.var = 'CBN_snn_res.0.5', stack.var = 'Doublet_SC', percentage = T)
dev.off()

PlotPDF('01.4.heat.scANVI_mapping', 12, 8)
p1 <- MappingHeatmap(pc.srt, que_var = 'CBN_snn_res.0.5', ref_var = 'Cell_state_scANVI' , ref.disp.min = 0.05)
p2 <- MappingHeatmap(pc.srt, que_var = 'CBN_snn_res.0.5', ref_var = 'Cell_state_SingleR', ref.disp.min = 0.05)
p1/p2
dev.off()

## Re-annotate
pc.srt$Cell_state <- NA
pc.srt$Cell_state[pc.srt$CBN_snn_res.0.5 %in% c(0, 1, 2, 5)] <- 'PC1'
pc.srt$Cell_state[pc.srt$CBN_snn_res.0.5 == 3] <- 'PC2'
pc.srt$Cell_state[pc.srt$CBN_snn_res.0.5 == 7] <- 'PC3'
pc.srt$Cell_state[pc.srt$CBN_snn_res.0.5 == 8] <- 'PC4'
pc.srt$Cell_state[pc.srt$CBN_snn_res.0.5 %in% c(4, 6)] <- 'Doublet'

PlotPDF('02.dim.all_cluster', 10, 10)
DimPlot2(pc.srt, reduction = 'sub_scANVI_umap', label = T, cols = mycol_10, group.by = 'Cell_state')
dev.off()


Idents(pc.srt) <- 'Cell_state'
levels(pc.srt) <- c('PC1',
                    'PC2',
                    'PC3',
                    'PC4',
                    'Doublet')
mk <- FindAllMarkers(pc.srt, only.pos = T, return.thresh = 0.01, logfc.threshold = 0.5)
p <- MarkerHeatmap(pc.srt, mk, n_cells = 500, top = 10)
PlotPDF('03.heat.maker_all_cluster', 15, 15)
p
dev.off()

####--------------------------------------------------------------------------------------------------------------------
saveRDS(pc.srt, 'integrated/STEP14.pc_cells.srt.rds')
####--------------------------------------------------------------------------------------------------------------------
