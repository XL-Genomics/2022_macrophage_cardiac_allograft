####--------------------------------------------------------------------------------------------------------------------
####  Pediatric Acute Rejection
####  2022-08-03 by Xiao LI (Texas Heart Institute, US)
####--------------------------------------------------------------------------------------------------------------------


####--------------------------------------------------------------------------------------------------------------------
####  Initiate  ####
####--------------------------------------------------------------------------------------------------------------------
Ver <- '1'
Step <- 'STEP13_Annotation_EC'
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
ec.srt <- srt[, srt$Cell_type %in% c('EC')]

DimPlot2(ec.srt, group.by = 'sample')

## Re-embed by scANVI
ec.srt <- RunUMAP(ec.srt, reduction = 'scANVI', dims = 1:30,
                  reduction.name = 'sub_scANVI_umap', reduction.key = 'subSCANVIUMAP_')
ec.srt <- FindNeighbors(ec.srt, dims = 1:30, reduction = 'scANVI', force.recalc = T) %>%
        FindClusters(resolution = c(0.1, 0.2, 0.5, 0.8))

ec.srt$Cell_state_previous[ec.srt$Cell_state_previous %in%
                                   names(table(ec.srt$Cell_state_previous))[table(ec.srt$Cell_state_previous)<20]] <- NA

DimPlot2(ec.srt, reduction = 'sub_scANVI_umap', group.by = 'name2', cols = mycol_10)
DimPlot2(ec.srt, reduction = 'sub_scANVI_umap', group.by = 'Doublet_SC')
DimPlot2(ec.srt, reduction = 'sub_scANVI_umap', group.by = 'Cell_state_previous', label = T, cols = mycol_20)+
        DimPlot2(ec.srt, reduction = 'sub_scANVI_umap', group.by = 'CBN_snn_res.0.2', label = T)+
        DimPlot2(ec.srt, reduction = 'sub_scANVI_umap', group.by = 'CBN_snn_res.0.5', label = T)+
        DimPlot2(ec.srt, reduction = 'sub_scANVI_umap', group.by = 'CBN_snn_res.0.8', label = T)

FeaturePlot2(ec.srt, features = markers_lvl1_min, reduction = 'sub_scANVI_umap')

## Annotate
Idents(ec.srt) <- 'CBN_snn_res.0.2'
mk <- FindAllMarkers(ec.srt, only.pos = T, return.thresh = 0.01, logfc.threshold = 0.5)
p <- MarkerHeatmap(ec.srt, mk, n_cells = 500, top = 10)
PlotPDF('01.1.heat.maker_all_cluster', 15, 15)
p
dev.off()

PlotPDF('01.2.dim.all_cluster', 10, 10)
DimPlot2(ec.srt, reduction = 'sub_scANVI_umap', label = T, cols = mycol_10)
dev.off()

PlotPDF('01.3.bar.doublet_in_all_cluster', 5, 5)
CountCellBarPlot(ec.srt, group.var = 'CBN_snn_res.0.2', stack.var = 'Doublet_SC', percentage = T)
dev.off()

PlotPDF('01.4.heat.scANVI_mapping', 12, 8)
p1 <- MappingHeatmap(ec.srt, que_var = 'CBN_snn_res.0.2', ref_var = 'Cell_state_scANVI' , ref.disp.min = 0.05)
p2 <- MappingHeatmap(ec.srt, que_var = 'CBN_snn_res.0.2', ref_var = 'Cell_state_SingleR', ref.disp.min = 0.05)
p1/p2
dev.off()

## Re-annotate
ec.srt$Cell_state <- NA
ec.srt$Cell_state[ec.srt$CBN_snn_res.0.2 == 0] <- 'EC1_cap'
ec.srt$Cell_state[ec.srt$CBN_snn_res.0.2 == 1] <- 'EC2_veno'
ec.srt$Cell_state[ec.srt$CBN_snn_res.0.5 == 4] <- 'EC3_art'
ec.srt$Cell_state[ec.srt$CBN_snn_res.0.2 == 4] <- 'EC4_endoc'
ec.srt$Cell_state[ec.srt$CBN_snn_res.0.2 == 6] <- 'EC5_lec'
ec.srt$Cell_state[ec.srt$CBN_snn_res.0.2 %in% c(2, 3, 5, 7)] <- 'Doublet'

PlotPDF('02.dim.all_cluster', 10, 10)
DimPlot2(ec.srt, reduction = 'sub_scANVI_umap', label = T, cols = mycol_10, group.by = 'Cell_state')
dev.off()

Idents(ec.srt) <- 'Cell_state'
levels(ec.srt) <- c('EC1_cap', 'EC2_veno', 'EC3_art', 'EC4_endoc', 'EC5_lec', 'Doublet')
mk <- FindAllMarkers(ec.srt, only.pos = T, return.thresh = 0.01, logfc.threshold = 0.5)
p <- MarkerHeatmap(ec.srt, mk, n_cells = 500, top = 10)
PlotPDF('03.heat.maker_all_cluster', 15, 15)
p
dev.off()

####--------------------------------------------------------------------------------------------------------------------
saveRDS(ec.srt, 'integrated/STEP13.ec_cells.srt.rds')
####--------------------------------------------------------------------------------------------------------------------
