####--------------------------------------------------------------------------------------------------------------------
####  Pediatric Acute Rejection
####  2022-08-03 by Xiao LI (Texas Heart Institute, US)
####--------------------------------------------------------------------------------------------------------------------


####--------------------------------------------------------------------------------------------------------------------
####  Initiate  ####
####--------------------------------------------------------------------------------------------------------------------
Ver <- '1'
Step <- 'STEP10_Annotation_Immune'
Project <- '2022_acute_rejection_dturaga'
Code_dir <- paste0('~/Documents/Bioinformatics/project/', Project, '/code/human_v', Ver, '/')

source(Sys.readlink(paste0(Code_dir, 'src/bioinformatics.R')))
source(Sys.readlink(paste0(Code_dir, 'src/scRNAseq.R')))
source(Sys.readlink(paste0(Code_dir, 'src/geneset_human.R')))
source(paste0(Code_dir, 'human_v', Ver, '.helper_functions.R'))

InitiateProject('Rivendell', Ver, Step, 'human', Project, 'ithil')

suppressMessages(library('SingleR'))
suppressMessages(library('celldex'))
suppressMessages(library('BiocParallel'))
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
imm.srt <- srt[, srt$Cell_type %in% c('Mye', 'Lym')]
DimPlot2(imm.srt, group.by = 'name')

## Re-embed by SCT-Harmony
imm.srt <- ProcessSrt_sct(imm.srt, do.umap = F, assay = 'CBN', regress_cc = F)
imm.srt <- ProcessSrt_hmn(imm.srt, do.umap = T, assay = 'SCT', haromize.by = 'sample')
imm.srt <- ProcessSrt_clust(imm.srt, resolution = c(0.1, 0.2, 0.5), reduction = 'harmony')

Table(imm.srt$Cell_state_previous)
imm.srt$tmp <- imm.srt$Cell_state_previous
imm.srt$tmp[imm.srt$tmp %in% names(Table(imm.srt$tmp))[Table(imm.srt$tmp)<30]] <- NA
imm.srt$tmp <- droplevels(imm.srt$tmp)

PlotPDF('01.1.umap.meta_cluster', 20, 12)
DimPlot2(imm.srt, reduction = 'hmn_umap', group.by = 'name', cols = mycol_10) +
        DimPlot2(imm.srt, reduction = 'hmn_umap', group.by = 'Doublet_SC') +
        DimPlot2(imm.srt, reduction = 'hmn_umap', group.by = 'tmp', label = T, cols = mycol_20) +
        DimPlot2(imm.srt, reduction = 'hmn_umap', group.by = 'SCT_snn_res.0.1', label = T, cols = mycol_20) +
        DimPlot2(imm.srt, reduction = 'hmn_umap', group.by = 'SCT_snn_res.0.2', label = T, cols = mycol_20) +
        DimPlot2(imm.srt, reduction = 'hmn_umap', group.by = 'SCT_snn_res.0.5', label = T, cols = mycol_20)
dev.off()

## Annotate
Idents(imm.srt) <- 'SCT_snn_res.0.5'
# mk <- FindAllMarkers(imm.srt, only.pos = T, return.thresh = 0.01)
# p <- MarkerHeatmap(imm.srt, mk)
# PlotPDF('01.2.heat.maker_all_cluster', 15, 15)
# p
# dev.off()

PlotPDF('01.2.dim.all_cluster', 10, 10)
DimPlot2(imm.srt, reduction = 'hmn_umap', label = T, cols = mycol_20)
dev.off()

PlotPDF('01.3.bar.doublet_in_all_cluster', 5, 5)
CountCellBarPlot(imm.srt, group.var = 'SCT_snn_res.0.5', stack.var = 'Doublet_SC', percentage = T)
dev.off()

PlotPDF('01.4.heat.scANVI_mapping', 12, 6)
MappingHeatmap(imm.srt, que_var = 'SCT_snn_res.0.5', ref_var = 'Cell_state_scANVI') |
        MappingHeatmap(imm.srt, que_var = 'SCT_snn_res.0.5', ref_var = 'Cell_type_scANVI')
dev.off()

genes <- c('CD2', 'CD8A',
           'NKG7', 'GZMA',
           'CD163', 'LYVE1',
           'FCGR3A', 'CD14',
           'CD1C', 'CLEC4A',
           'CD83', 'CCR7',
           'CD79A', 'MS4A1',
           'IGHG1', 'JCHAIN',
           'KIT', 'CPA3',
           'PECAM1', 'TNNT2', 'DMD', 'GSN', 'TAGLN', 'RGS5')

PlotPDF('01.5.dot.manual_marker_inspection', 8, 16)
DimPlot2(imm.srt, reduction = 'hmn_umap', group.by = 'SCT_snn_res.0.5', label = T)/
        DotPlot2(imm.srt, features = genes)
dev.off()

imm.srt$Cell_state <- NA
imm.srt$Cell_state[imm.srt$SCT_snn_res.0.5 %in% c(0, 2, 3, 7, 14)] <- 'Mono/MP'
imm.srt$Cell_state[imm.srt$SCT_snn_res.0.5 %in% c(1, 4)] <- 'NK/T'
imm.srt$Cell_state[imm.srt$SCT_snn_res.0.5 %in% c(12)] <- 'BC'
imm.srt$Cell_state[imm.srt$SCT_snn_res.0.5 %in% c(13)] <- 'Mast'
imm.srt$Cell_state[imm.srt$SCT_snn_res.0.5 %in% c(10)] <- 'Mito'
imm.srt$Cell_state[imm.srt$SCT_snn_res.0.5 %in% c(5, 6, 8, 9, 11)] <- 'Doublet'

PlotPDF('01.6.dim.all_cluster', 10, 10)
DimPlot2(imm.srt, reduction = 'hmn_umap', label = T, cols = mycol_10, group.by = 'Cell_state')
dev.off()

Idents(imm.srt) <- 'Cell_state'
mk <- FindAllMarkers(imm.srt, only.pos = T, return.thresh = 0.01)
p <- MarkerHeatmap(imm.srt, mk)
PlotPDF('01.7.heat.maker_all_cluster', 15, 15)
p
dev.off()

imm.srt <- DietSeurat(imm.srt, dimreducs = names(imm.srt@reductions), scale.data = F, assays = names(imm.srt@assays))
####--------------------------------------------------------------------------------------------------------------------
saveRDS(imm.srt, 'integrated/STEP10.immune_cells.srt.rds')
####--------------------------------------------------------------------------------------------------------------------


####--------------------------------------------------------------------------------------------------------------------
####  Myeloid cell annotation  ####
####--------------------------------------------------------------------------------------------------------------------
mye.srt <- imm.srt[, imm.srt$Cell_state %in% c('Mono/MP', 'Mast')]
DimPlot2(mye.srt, group.by = 'Cell_state', reduction = 'hmn_umap')

## Re-embed UMAP
mye.srt <- RunUMAP(mye.srt, dims = 1:16, seed.use = 505, reduction = 'harmony', min.dist = 0.5,
                   reduction.name = 'hmn_umap', reduction.key = 'hmnumap_')

## Make embedding consistant with previous version
x <- mye.srt@reductions$hmn_umap@cell.embeddings[, 1]
y <- mye.srt@reductions$hmn_umap@cell.embeddings[, 2]
mye.srt@reductions$hmn_umap@cell.embeddings[, 1] <- x
mye.srt@reductions$hmn_umap@cell.embeddings[, 2] <- -y

p1 <- DimPlot2(mye.srt, reduction = 'hmn_umap', group.by = 'tmp', cols = mycol_20)
p2 <- DimPlot2(mye.srt, reduction = 'hmn_umap', group.by = 'name', cols = mycol_10)

PlotPDF('02.1.umap.myeloid_cells', 5, 10)
p1/p2
dev.off()

mye.srt <- FindNeighbors(mye.srt, reduction = 'harmony', dims = 1:16) |> FindClusters(res = c(0.1, 0.2, 0.5))
DimPlot2(mye.srt, reduction = 'hmn_umap', group.by = 'SCT_snn_res.0.5', cols = mycol_20, label = T) /
        DimPlot2(mye.srt, reduction = 'hmn_umap', group.by = 'tmp', cols = mycol_20,  label = T)

mye.srt$Cell_state[mye.srt$SCT_snn_res.0.5 %in% c(5)] <- 'MP1'
mye.srt$Cell_state[mye.srt$SCT_snn_res.0.5 %in% c(0, 1, 2, 8, 6)] <- 'MP2'
mye.srt$Cell_state[mye.srt$SCT_snn_res.0.5 %in% c(3, 4, 7, 9)] <- 'Mono'
mye.srt$Cell_state[mye.srt$SCT_snn_res.0.5 %in% c(10)] <- 'Mast'

DimPlot2(mye.srt, reduction = 'hmn_umap', group.by = 'Cell_state', cols = mycol_10,  label = T)
####--------------------------------------------------------------------------------------------------------------------
saveRDS(mye.srt, 'integrated/STEP10.myeloid_cells.srt.rds')
####--------------------------------------------------------------------------------------------------------------------


####--------------------------------------------------------------------------------------------------------------------
####  Lymphoid cell annotation  ####
####--------------------------------------------------------------------------------------------------------------------
##  Celldex ref mapping
que.sce <- as.SingleCellExperiment(DietSeurat(imm.srt))
ref <- BlueprintEncodeData()
pred_cellsubtype <- SingleR(test = que.sce,
                            ref = ref,
                            labels = factor(ref$label.fine),
                            aggr.ref = T,
                            BPPARAM = SnowParam(4))
imm.srt$Cell_state_SingleR_BE <- pred_cellsubtype$pruned.labels
p <- MappingHeatmap(imm.srt, que_var = 'Cell_state', ref_var = 'Cell_state_SingleR_BE', ref.disp.min = 0.1)
PlotPDF('03.1.heat.ref_mapping_to_celldex_be', 5, 5)
p
dev.off()

lym.srt <- imm.srt[, imm.srt$Cell_state %in% c('NK/T')]
lym.srt <- RunUMAP(lym.srt, dims = 1:35, seed.use = 505, reduction = 'harmony', min.dist = 0.4,
                   reduction.name = 'hmn_umap', reduction.key = 'hmnumap_')
lym.srt <- FindNeighbors(lym.srt, reduction = 'harmony', dims = 1:35) |>
        FindClusters(res = c(0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.8, 1))

DimPlot2(lym.srt, reduction = 'hmn_umap', cols = mycol_20, group.by = 'SCT_snn_res.0.8') /
        DimPlot2(lym.srt, reduction = 'hmn_umap', cols = mycol_20, group.by = 'tmp') /
        MappingHeatmap(lym.srt, que_var = 'SCT_snn_res.0.8', ref_var = 'Cell_state_SingleR_BE', ref.disp.min = 0.1)

lym.srt$Cell_state[lym.srt$SCT_snn_res.0.8 %in% c(7, 4, 6)] <- 'CD4 TC'
lym.srt$Cell_state[lym.srt$SCT_snn_res.0.8 %in% c(1, 0, 9)] <- 'CD8 TC'
lym.srt$Cell_state[lym.srt$SCT_snn_res.0.8 %in% c(2, 3)] <- 'NK1'
lym.srt$Cell_state[lym.srt$SCT_snn_res.0.8 %in% c(5)] <- 'NK2'
lym.srt$Cell_state[lym.srt$SCT_snn_res.0.8 %in% c(8)] <- 'NK3'
lym.srt$Cell_state <- factor(lym.srt$Cell_state, levels = sort(U(lym.srt$Cell_state)))

PlotPDF('03.2.dim.all_cluster', 10, 10)
DimPlot2(lym.srt, reduction = 'hmn_umap', label = T, cols = mycol_10, group.by = 'Cell_state')
dev.off()

p <- MappingHeatmap(lym.srt, que_var = 'Cell_state', ref_var = 'Cell_state_SingleR_BE', ref.disp.min = 0.1)
PlotPDF('03.3.heat.ref_mapping_to_celldex_be', 10, 10)
p
dev.off()

Idents(lym.srt) <- 'Cell_state'
mk <- FindAllMarkers(lym.srt, only.pos = T, return.thresh = 0.01, test.use = 'MAST')
p <- MarkerHeatmap(lym.srt, mk)
PlotPDF('03.4.heat.maker_all_cluster', 15, 15)
p
dev.off()
####--------------------------------------------------------------------------------------------------------------------
saveRDS(lym.srt, 'integrated/STEP10.lymphoid_cells.srt.rds')
####--------------------------------------------------------------------------------------------------------------------
