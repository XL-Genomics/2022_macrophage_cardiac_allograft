####--------------------------------------------------------------------------------------------------------------------
####  Pediatric Acute Rejection
####  2022-08-03 by Xiao LI (Texas Heart Institute, US)
####--------------------------------------------------------------------------------------------------------------------


####--------------------------------------------------------------------------------------------------------------------
####  Initiate  ####
####--------------------------------------------------------------------------------------------------------------------
Ver <- '1'
Step <- 'STEP15_Annotation_SMC'
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
smc.srt <- srt[, srt$Cell_type %in% c('SMC')]

DimPlot2(smc.srt, group.by = 'sample')

## Re-embed by scANVI
smc.srt <- RunUMAP(smc.srt, reduction = 'scANVI', dims = 1:30,
                   reduction.name = 'sub_scANVI_umap', reduction.key = 'subSCANVIUMAP_')
smc.srt <- FindNeighbors(smc.srt, dims = 1:30, reduction = 'scANVI', force.recalc = T) %>%
        FindClusters(resolution = c(0.1, 0.2, 0.5))

DimPlot2(smc.srt, reduction = 'sub_scANVI_umap', group.by = 'name2', cols = mycol_10)
DimPlot2(smc.srt, reduction = 'sub_scANVI_umap', group.by = 'Doublet_SC')

smc.srt$Cell_state_previous[smc.srt$Cell_state_previous %in%
                                    names(table(smc.srt$Cell_state_previous))[table(smc.srt$Cell_state_previous)<20]] <- NA

DimPlot2(smc.srt, reduction = 'sub_scANVI_umap', group.by = 'Cell_state_previous', label = T, cols = mycol_20)+
DimPlot2(smc.srt, reduction = 'sub_scANVI_umap', group.by = 'CBN_snn_res.0.1', label = T)+
DimPlot2(smc.srt, reduction = 'sub_scANVI_umap', group.by = 'CBN_snn_res.0.2', label = T)+
DimPlot2(smc.srt, reduction = 'sub_scANVI_umap', group.by = 'CBN_snn_res.0.5', label = T)

## Annotate
Idents(smc.srt) <- 'CBN_snn_res.0.1'
mk <- FindAllMarkers(smc.srt, only.pos = T, return.thresh = 0.01, logfc.threshold = 0.5)
p <- MarkerHeatmap(smc.srt, mk, n_cells = 500, top = 10)
PlotPDF('01.1.heat.maker_all_cluster', 15, 15)
p
dev.off()

PlotPDF('01.2.dim.all_cluster', 10, 10)
DimPlot2(smc.srt, reduction = 'sub_scANVI_umap', label = T, cols = mycol_10)
dev.off()

PlotPDF('01.3.bar.doublet_in_all_cluster', 5, 5)
CountCellBarPlot(smc.srt, group.var = 'CBN_snn_res.0.1', stack.var = 'Doublet_SC', percentage = T)
dev.off()

PlotPDF('01.4.heat.scANVI_mapping', 12, 8)
p1 <- MappingHeatmap(smc.srt, que_var = 'CBN_snn_res.0.1', ref_var = 'Cell_state_scANVI' , ref.disp.min = 0.05)
p2 <- MappingHeatmap(smc.srt, que_var = 'CBN_snn_res.0.1', ref_var = 'Cell_state_SingleR', ref.disp.min = 0.05)
p1/p2
dev.off()

## Re-annotate
smc.srt$Cell_state <- NA
smc.srt$Cell_state[smc.srt$CBN_snn_res.0.1 == 0] <- 'SMC1'
smc.srt$Cell_state[smc.srt$CBN_snn_res.0.1 == 1] <- 'SMC2'
smc.srt$Cell_state[smc.srt$CBN_snn_res.0.1 == 2] <- 'SMC3'

PlotPDF('02.dim.all_cluster', 10, 10)
DimPlot2(smc.srt, reduction = 'sub_scANVI_umap', label = T, cols = mycol_10, group.by = 'Cell_state')
dev.off()


Idents(smc.srt) <- 'Cell_state'
levels(smc.srt) <- c('SMC1',
                     'SMC2',
                     'SMC3')
mk <- FindAllMarkers(smc.srt, only.pos = T, return.thresh = 0.01, logfc.threshold = 0.5)
p <- MarkerHeatmap(smc.srt, mk, n_cells = 500, top = 10)
PlotPDF('03.heat.maker_all_cluster', 15, 15)
p
dev.off()

####--------------------------------------------------------------------------------------------------------------------
saveRDS(smc.srt, 'integrated/STEP15.smc_cells.srt.rds')
####--------------------------------------------------------------------------------------------------------------------
