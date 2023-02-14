####--------------------------------------------------------------------------------------------------------------------
####  Pediatric Acute Rejection
####  2022-08-03 by Xiao LI (Texas Heart Institute, US)
####--------------------------------------------------------------------------------------------------------------------


####--------------------------------------------------------------------------------------------------------------------
####  Initiate  ####
####--------------------------------------------------------------------------------------------------------------------
Ver <- '1'
Step <- 'STEP09_Annotation'
Project <- '2022_acute_rejection_dturaga'
Code_dir <- paste0('~/Documents/Bioinformatics/project/', Project, '/code/human_v', Ver, '/')

source(Sys.readlink(paste0(Code_dir, 'src/bioinformatics.R')))
source(Sys.readlink(paste0(Code_dir, 'src/scRNAseq.R')))
source(Sys.readlink(paste0(Code_dir, 'src/geneset_human.R')))
source(paste0(Code_dir, 'human_v', Ver, '.helper_functions.R'))

InitiateProject('Rivendell', Ver, Step, 'human', Project, 'ithil')

Color_cell_type <- c(mycol_40[c(3, 6, 11, 14, 35, 17, 22, 26)], 'grey85')
Color_cell_state <- c(mycol_BuGr[1:4],
                      mycol_30[4:6],
                      mycol_30[7:9], mycol_30[22:24],
                      mycol_40[c(14, 35, 17, 22, 26)],
                      'grey80')
Color_sample_pair <- clr_darken(shift = 0, colorRampPalette(c(brewer.pal(10,"Paired")))(10))
####--------------------------------------------------------------------------------------------------------------------


####--------------------------------------------------------------------------------------------------------------------
####  Load data  ####
####--------------------------------------------------------------------------------------------------------------------
srt <- readRDS('integrated/STEP07.scvi_integration.srt.rds')
srt.bkp <- srt

pred_celltype <- readRDS('analysis/STEP06.singler_pred_celltype.rds')
pred_cellstate <- readRDS('analysis/STEP06.singler_pred_cellsubtype.rds')
marker_score <- readRDS('analysis/STEP06.ref_marker_score.srt_meta.rds')

srt$Cell_type_SingleR <- pred_celltype$pruned.labels
srt$Cell_state_SingleR <- pred_cellstate$pruned.labels
srt <- AddMetaData(srt, metadata = marker_score)

colnames(srt@meta.data)
####--------------------------------------------------------------------------------------------------------------------


####--------------------------------------------------------------------------------------------------------------------
####  Use Louvain res0.8 and scVI umap  ####
####--------------------------------------------------------------------------------------------------------------------
# DimPlot2(srt, reduction = 'scANVI_umap', group.by = 'Cluster_leiden_res0.8', cols = mycol_30, label = T, repel = T)

srt <- DietSeurat(srt, counts = T, data = T, assays = 'CBN', scale.data = F, dimreducs = c('scANVI_umap', 'scANVI'))
srt@reductions <- srt@reductions[c('scANVI', 'scANVI_umap')]

srt <- FindNeighbors(srt, reduction = 'scANVI', dims = 1:30) |> FindClusters(res = c(0.5, 0.8))

DimPlot2(srt, reduction = 'scANVI_umap', group.by = 'CBN_snn_res.0.5', cols = mycol_30, label = T, repel = T) /
DimPlot2(srt, reduction = 'scANVI_umap', group.by = 'CBN_snn_res.0.8', cols = mycol_40, label = T, repel = T)

Idents(srt) <- 'CBN_snn_res.0.5'
####--------------------------------------------------------------------------------------------------------------------


####--------------------------------------------------------------------------------------------------------------------
####  Label doublet clusters  ####
####--------------------------------------------------------------------------------------------------------------------
p1 <- CountCellBarPlot(srt, group.var = 'CBN_snn_res.0.5', stack.var = 'Doublet_SC', percentage = T)$plot
p2 <- VlnPlot2(srt, features = 'Doublet_SC_score', sort = T, group.by = 'CBN_snn_res.0.5')
p3 <- DimPlot2(srt, reduction = 'scANVI_umap', group.by = 'CBN_snn_res.0.5', cols = mycol_30, label = T, repel = T)
p4 <- DimPlot2(srt, reduction = 'scANVI_umap', group.by = 'Doublet_SC', cols = mycol_10, label = F, repel = F)

PlotPDF('01.1.bar.doublet_rich_clusters', 20, 15)
print((p1 + p2)/(p3 + p4))
dev.off()

srt$non_ambiguous <- T
srt$non_ambiguous[srt$CBN_snn_res.0.5 %in% c(8, 15, 16, 18, 19)] <- F

DimPlot2(srt, group.by = 'non_ambiguous')
CountCellBarPlot(srt, group.var = 'sample', stack.var = 'non_ambiguous', percentage = T)$plot
####--------------------------------------------------------------------------------------------------------------------


####--------------------------------------------------------------------------------------------------------------------
####  Annotate Cell Type  ####
####--------------------------------------------------------------------------------------------------------------------
srt$Cell_type_SingleR <- factor(srt$Cell_type_SingleR)
ref_order <- c('Atrial_Cardiomyocyte', 'Ventricular_Cardiomyocyte', 'Fibroblast', 'Mesothelial', 'Endothelial',
              'Smooth_muscle_cells', 'Pericytes', 'Myeloid', 'Lymphoid', 'Neuronal', 'Adipocytes')
# p1 <- MappingHeatmap(srt, que_var = 'CBN_snn_res.0.5', ref_var = 'Cell_type_SingleR',
#                      percentage = T, log10_scale = F, ref.disp.min = 0,
#                      que_order = c(
#                              0, 6, 14,
#                              2, 13, 17,
#                              12, 22,
#                              4, 7, 18,
#                              5, 9, 15, 21,
#                      ),
#                      ref_order = c(ref_order, 'NotAssigned', NA))

prev_v.srt <- readRDS('../human_v0/integrated/STEP18.annotated.srt.rds')
O(Cells(prev_v.srt), Cells(srt))
srt$Cell_type_previous <- factor('New', levels = c(levels(prev_v.srt$Cell_type), 'New'))
srt$Cell_type_previous[Cells(prev_v.srt)] <- prev_v.srt$Cell_type
srt$Cell_state_previous <- factor('New', levels = c(levels(prev_v.srt$Cell_state), 'New'))
srt$Cell_state_previous[Cells(prev_v.srt)] <- prev_v.srt$Cell_state

p1 <- MappingHeatmap(srt, que_var = 'CBN_snn_res.0.5', ref_var = 'Cell_type_previous',
                     percentage = T, log10_scale = F, ref.disp.min = 0,
                     que_order = c(
                             0, 3, 6, 14,
                             2, 13, 16, 17, 19,
                             1, 7, 20, 24, 25,
                             4, 18,
                             12, 22,
                             9,
                             5, 15, 21,
                             11,
                             23,
                             10,
                             8
                     ),
                     ref_order = levels(srt$Cell_type_previous))

p2 <- MappingHeatmap(srt, que_var = 'CBN_snn_res.0.5', ref_var = 'Cell_type_scANVI',
                     percentage = T, log10_scale = F, ref.disp.min = 0,
                     que_order = c(
                             0,
                             3, 6,
                             2, 13, 16, 17, 19,
                             10, 14, 25,
                             1, 7, 8, 20, 24,
                             12, 22,
                             4, 18,
                             5, 15, 21,
                             9,
                             11,
                             23

                     ),
                     ref_order = ref_order
)
p3 <- DotPlot2(srt, features = grep('^Markerset_', colnames(srt@meta.data), value = T))

p <- p1 | p2 | p3
PlotPDF('02.1.dot_heat.celltype_prediction', 16, 10)
print(p)
dev.off()

## Annotation based on scANVI integration with Teichmann data (and previous annotation)
srt$Cell_type <- as.vector(srt$CBN_snn_res.0.5)
srt$Cell_type[srt$CBN_snn_res.0.5 %in% c(0, 3, 6, 14)] <- 'CM'
srt$Cell_type[srt$CBN_snn_res.0.5 %in% c(2, 13, 16, 17, 19)] <- 'FB'
srt$Cell_type[srt$CBN_snn_res.0.5 %in% c(1, 7, 20, 24, 25, 8)] <- 'EC'
srt$Cell_type[srt$CBN_snn_res.0.5 %in% c(12, 22)] <- 'SMC'
srt$Cell_type[srt$CBN_snn_res.0.5 %in% c(4, 18)] <- 'PC'
srt$Cell_type[srt$CBN_snn_res.0.5 %in% c(9)] <- 'Lym'
srt$Cell_type[srt$CBN_snn_res.0.5 %in% c(5, 15, 21)] <- 'Mye'
srt$Cell_type[srt$CBN_snn_res.0.5 %in% c(11)] <- 'Neu'
srt$Cell_type[srt$CBN_snn_res.0.5 %in% c(23)] <- 'Adi'
srt$Cell_type[srt$CBN_snn_res.0.5 %in% c(10)] <- 'Mito'

p1 <- DimPlot2(srt, group.by = 'Cell_type', cols = mycol_20, label = T)
PlotPDF('02.2.umap.annotated', 6, 5)
p1
dev.off()
####--------------------------------------------------------------------------------------------------------------------
df <- srt@meta.data[,!grepl('^Cluster_leiden_|^Pred_cell|^Markerset_|non_ambiguous|LowQual|CBN_snn_res|seurat_clusters',
                            colnames(srt@meta.data))]
saveRDS(df, 'integrated/STEP09.annotated.srt_meta.rds')
####--------------------------------------------------------------------------------------------------------------------

