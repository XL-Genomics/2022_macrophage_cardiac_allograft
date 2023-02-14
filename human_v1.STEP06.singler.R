####--------------------------------------------------------------------------------------------------------------------
####  Pediatric Acute Rejection
####  2022-08-03 by Xiao LI (Texas Heart Institute, US)
####--------------------------------------------------------------------------------------------------------------------


####--------------------------------------------------------------------------------------------------------------------
####  Initiate  ####
####--------------------------------------------------------------------------------------------------------------------
Ver <- '1'
Step <- 'STEP06_SingleR'
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
srt <- readRDS('integrated/STEP05.merged.clean_cbn.srt.rds')
####--------------------------------------------------------------------------------------------------------------------


####--------------------------------------------------------------------------------------------------------------------
####  Mapping with 2020_Nature_STeichmann annotation ####
####--------------------------------------------------------------------------------------------------------------------
## Prepare reference
ref.srt <- readRDS('/Volumes/shire/data/scrnaseq/2020_Nature_STeichmann/matrix_public/raw.srt.rds')

dev.off()
PlotPDF('01.1.heat.2020_Nature_STeichmann_annotation', 13, 4)
pheatmap(log10(as.matrix(Table(ref.srt$cell_type, ref.srt$cell_states))+1), cellwidth = 10, cellheight = 10)
dev.off()

ref.srt <- ref.srt[, ref.srt$cell_type != 'doublets']
ref.srt$Ref_cell_type <- droplevels(ref.srt$cell_type)
ref.srt$Ref_cell_subtype <- droplevels(ref.srt$cell_states)
levels(ref.srt$Ref_cell_subtype)[levels(ref.srt$Ref_cell_subtype)=="nan"] <- "NotAssigned"
Idents(ref.srt) <- 'Ref_cell_subtype'

ref.sce <- as.SingleCellExperiment(DietSeurat(ref.srt))
que.sce <- as.SingleCellExperiment(DietSeurat(srt))
gc()

## SingleR mapping at cell type level
system.time({pred_celltype <- SingleR(test = que.sce,
                                     ref = ref.sce,
                                     labels = factor(ref.sce$Ref_cell_type),
                                     aggr.ref = T,
                                     BPPARAM = SnowParam(4))
}) ## elapsed 1250.677s
####--------------------------------------------------------------------------------------------------------------------
saveRDS(pred_celltype, 'analysis/STEP06.singler_pred_celltype.rds')
####--------------------------------------------------------------------------------------------------------------------

## SingleR mapping at cell subtype level
system.time({pred_cellsubtype <- SingleR(test = que.sce,
                                        ref = ref.sce,
                                        labels = factor(ref.sce$Ref_cell_subtype),
                                        aggr.ref = T,
                                        BPPARAM = SnowParam(4))
}) ## elapsed 1746.813s
####--------------------------------------------------------------------------------------------------------------------
saveRDS(pred_cellsubtype, 'analysis/STEP06.singler_pred_cellsubtype.rds')
####--------------------------------------------------------------------------------------------------------------------

srt$Pred_celltype <- pred_celltype$pruned.labels
srt$Pred_cellsubtype <- pred_cellsubtype$pruned.labels

p <- DimPlot2(srt, group.by = 'Pred_celltype', raster = T, reduction = 'hmn_umap',
              cols = mycol_20, label = T, repel = T)
PlotPDF('02.1.dim.pred_celltype', 12, 14)
print(p)
dev.off()

p <- DimPlot2(srt, group.by = 'Pred_cellsubtype', raster = T, reduction = 'hmn_umap',
              cols = c(mycol_60, 'black'), label = T, repel = T)
PlotPDF('02.2.dim.pred_cellstate', 12, 14)
print(p)
dev.off()


###--------------------------------------------------------------------------------------------------------------------
####  Expression of 2020_Nature_STeichmann markers  ####
####--------------------------------------------------------------------------------------------------------------------
Idents(ref.srt) <- 'Ref_cell_type'
ref.marker <- readRDS('/Volumes/shire/data/scrnaseq/2020_Nature_STeichmann/matrix_public/celltype_markers.rds')
ref.marker.high <-  ref.marker |> filter(pct.1 >= 0.3, p_val_adj <= 1e-8, avg_log2FC >= 10)
genes <- split(ref.marker.high$gene, ref.marker.high$cluster)
srt <- AddModuleScore2(srt, features = genes, name = c(
        'Markerset_Adipo',
        'Markerset_aCM',
        'Markerset_Endo',
        'Markerset_Fibro',
        'Markerset_Lymph',
        'Markerset_Meso',
        'Markerset_Myel',
        'Markerset_Neuro',
        'Markerset_NotAssigned',
        'Markerset_Peric',
        'Markerset_Smooth',
        'Markerset_vCM'
))

p <- FeaturePlot2(srt, features = grep('Markerset_', colnames(srt@meta.data), value = T),
                  raster = T, reduction = 'hmn_umap', ncol = 6, min.cutoff = 'q5', max.cutoff = 'q95')
PlotPDF('02.3.feature.merged_umap_lvl1_marker', 24, 12)
print(p)
dev.off()

df <- srt@meta.data[, grep('Markerset_', colnames(srt@meta.data), value = T)]
head(df)
####--------------------------------------------------------------------------------------------------------------------
saveRDS(df, 'analysis/STEP06.ref_marker_score.srt_meta.rds')
####--------------------------------------------------------------------------------------------------------------------
