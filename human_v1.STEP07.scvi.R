####--------------------------------------------------------------------------------------------------------------------
####  Pediatric Acute Rejection
####  2022-08-03 by Xiao LI (Texas Heart Institute, US)
####--------------------------------------------------------------------------------------------------------------------


####--------------------------------------------------------------------------------------------------------------------
####  Initiate  ####
####--------------------------------------------------------------------------------------------------------------------
Ver <- '1'
Step <- 'STEP07_scVI'
Project <- '2022_acute_rejection_dturaga'
Code_dir <- paste0('~/Documents/Bioinformatics/project/', Project, '/code/human_v', Ver, '/')

source(Sys.readlink(paste0(Code_dir, 'src/bioinformatics.R')))
source(Sys.readlink(paste0(Code_dir, 'src/scRNAseq.R')))
source(Sys.readlink(paste0(Code_dir, 'src/geneset_human.R')))
source(paste0(Code_dir, 'human_v', Ver, '.helper_functions.R'))

InitiateProject('Rivendell', Ver, Step, 'human', Project, 'ithil')


suppressMessages(library('reticulate'))
suppressMessages(library('anndata'))
sc <- import("scanpy")
np <- import('numpy')
sce <- import('scanpy.external')
####--------------------------------------------------------------------------------------------------------------------


####--------------------------------------------------------------------------------------------------------------------
####  Load data  ####
####--------------------------------------------------------------------------------------------------------------------
srt <- readRDS('integrated/STEP05.merged.clean_cbn.srt.rds')
DefaultAssay(srt) <- 'CBN'
####--------------------------------------------------------------------------------------------------------------------


####--------------------------------------------------------------------------------------------------------------------
####  Convert to ann.data  ####
####--------------------------------------------------------------------------------------------------------------------
SaveH5ad(srt, path = 'integrated/', name = 'STEP07.merged.clean_cbn.ann',
         assay = 'CBN', raw_count_only = T, verbose = T)

## Following code is for avoiding "_index" in adata.var bug
adata <- read_h5ad('integrated/STEP07.merged.clean_cbn.ann.h5ad')
adata$raw <- NULL
adata$write_h5ad(filename = 'integrated/STEP07.merged.clean_cbn.ann.h5ad') ## replace the original
adata$X[1:10, 1:20]
####--------------------------------------------------------------------------------------------------------------------



####--------------------------------------------------------------------------------------------------------------------
####  Follow jupyter notebook STEP07 and STEP08  ####
####--------------------------------------------------------------------------------------------------------------------


####--------------------------------------------------------------------------------------------------------------------
####  Re-embed based on scVI latent  ####
####--------------------------------------------------------------------------------------------------------------------
# adata <- read_h5ad('integrated/STEP08.scanvi_integration.h5ad') ## didn't work due to bug

obsm <- read.csv(file = 'integrated/STEP08.scanvi_integration.csvs/obsm.csv', sep = ',', header = T)
obs <- read.csv(file = 'integrated/STEP08.scanvi_integration.csvs/obs.csv', sep = ',', header = T)
rownames(obsm) <- str_remove(obs$X, pattern = '-0$')
O(Cells(srt), rownames(obsm))
identical(Cells(srt), rownames(obsm))

####--------------------------------------------------------------------------------------------------------------------
# scvi_mtx <- adata$layers$get('scvi_normalized')
# dim(scvi_mtx)
# H(scvi_mtx)
# O(colnames(scvi_mtx), rownames(srt))
# O(str_remove(rownames(scvi_mtx), pattern = '-0$'), Cells(srt))
# rownames(scvi_mtx) <- str_remove(rownames(scvi_mtx), pattern = '-0$')
# identical(rownames(scvi_mtx), Cells(srt))
# srt[['SCVI']] <- CreateAssayObject(data = t(scvi_mtx), check.matrix = T)

## Save reduced dimensions
srt@reductions$scVI_umap <- srt@reductions$umap
srt@reductions$scVI_umap@cell.embeddings <- as.matrix(obsm[, paste0('X_scVIumap', 1:2)])
colnames(srt@reductions$scVI_umap@cell.embeddings) <- c('SCVIUMAP_1', 'SCVIUMAP_2')
srt@reductions$scVI_umap@assay.used <- 'CBN'
srt@reductions$scVI_umap@key <- 'SCVIUMAP_'
DimPlot2(srt, reduction = 'scVI_umap', group.by = 'name2')

srt@reductions$scVI <- srt@reductions$umap
srt@reductions$scVI@cell.embeddings <- as.matrix(obsm[, paste0('X_scVI', 1:30)])
rownames(srt@reductions$scVI@cell.embeddings) <- Cells(srt)
colnames(srt@reductions$scVI@cell.embeddings) <- paste0('SCVI_', 1:30)
srt@reductions$scVI@assay.used <- 'CBN'
srt@reductions$scVI@key <- 'SCVI_'

srt@reductions$scANVI_umap <- srt@reductions$umap
srt@reductions$scANVI_umap@cell.embeddings <- as.matrix(obsm[, paste0('X_scANVIumap', 1:2)])
colnames(srt@reductions$scANVI_umap@cell.embeddings) <- c('SCANVIUMAP_1', 'SCANVIUMAP_2')
srt@reductions$scANVI_umap@assay.used <- 'CBN'
srt@reductions$scANVI_umap@key <- 'SCANVIUMAP_'
DimPlot2(srt, reduction = 'scANVI_umap', group.by = 'name2')

srt@reductions$scANVI <- srt@reductions$umap
srt@reductions$scANVI@cell.embeddings <- as.matrix(obsm[, paste0('X_scANVI', 1:30)])
rownames(srt@reductions$scANVI@cell.embeddings) <- Cells(srt)
colnames(srt@reductions$scANVI@cell.embeddings) <- paste0('SCANVI_', 1:30)
srt@reductions$scANVI@assay.used <- 'CBN'
srt@reductions$scANVI@key <- 'SCANVI_'

## Save metadata
srt$Cell_type_scANVI <- obs$scANVI_predict_celltype
srt$Cell_state_scANVI <- obs$scANVI_predict_cellstate
srt$Cluster_leiden_res1.5 <- obs$Cluster_leiden_res1.5
srt$Cluster_leiden_res1.0 <- obs$Cluster_leiden_res1.0
srt$Cluster_leiden_res0.8 <- obs$Cluster_leiden_res0.8
srt$Cluster_leiden_res0.5 <- obs$Cluster_leiden_res0.5
srt$Cluster_leiden_res0.2 <- obs$Cluster_leiden_res0.2

PlotPDF('01.umap.by_sample_by_predicted_celltype', 20, 8)
DimPlot2(srt, group.by = 'name', cols = mycol_20) | DimPlot2(srt, group.by = 'Cell_type_scANVI', cols = mycol_20)
dev.off()

PlotPDF('02.umap.doublet', 10, 8)
DimPlot2(srt, group.by = 'Doublet_SC', cols = mycol_10)
dev.off()

PlotPDF('03.umap.leiden', 22, 16)
DimPlot2(srt, group.by = c(paste0('Cluster_leiden_res', c('0.2', '0.5', '0.8', '1.0', '1.5')), 'Doublet_SC'),
         cols = mycol_40, ncol = 3, label = T, repel = T)
dev.off()

p <- MappingHeatmap(srt, que_var = 'Cluster_leiden_res0.8', ref_var = 'Cell_type_scANVI', ref.disp.min = 0,
                    que_order = c(
                            4,
                            3, 5,
                            2, 17,
                            18, 13,
                            0, 9, 11, 14,
                            7,
                            1,
                            6, 15,
                            8, 16,
                            10,
                            12
                    ),
                    ref_order = c(
                            'Atrial_Cardiomyocyte',
                            'Ventricular_Cardiomyocyte',
                            'Fibroblast',
                            'Mesothelial',
                            'Endothelial',
                            'Smooth_muscle_cells',
                            'Pericytes',
                            'Myeloid',
                            'Lymphoid',
                            'Neuronal',
                            'Adipocytes'
                    )
)
p
PlotPDF('04.heat.leiden_0.8_vs_scANVI_celltype', 5, 8)
p
dev.off()
####--------------------------------------------------------------------------------------------------------------------
saveRDS(srt, 'integrated/STEP07.scvi_integration.srt.rds')
####--------------------------------------------------------------------------------------------------------------------


## Test
srt
srt2 <- RunUMAP(srt, reduction = 'scANVI', reduction.name = 'scANVI_umap2', reduction.key = 'SCANVIUMAP2_', dims = 1:30)
DimPlot2(srt2, reduction = 'scANVI_umap2', group.by = 'Cluster_leiden_res0.8', label = T, cols = mycol_30)+
DimPlot2(srt2, reduction = 'scANVI_umap2', group.by = 'sample', label = T, cols = mycol_10)
