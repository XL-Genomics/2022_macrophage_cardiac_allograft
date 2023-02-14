####--------------------------------------------------------------------------------------------------------------------
####  Pediatric Acute Rejection
####  2022-08-03 by Xiao LI (Texas Heart Institute, US)
####--------------------------------------------------------------------------------------------------------------------


####--------------------------------------------------------------------------------------------------------------------
####  Initiate  ####
####--------------------------------------------------------------------------------------------------------------------
Ver <- '1'
Step <- 'STEP05_Cleanup'
Project <- '2022_acute_rejection_dturaga'
Code_dir <- paste0('~/Documents/Bioinformatics/project/', Project, '/code/human_v', Ver, '/')

source(Sys.readlink(paste0(Code_dir, 'src/bioinformatics.R')))
source(Sys.readlink(paste0(Code_dir, 'src/scRNAseq.R')))
source(Sys.readlink(paste0(Code_dir, 'src/geneset_human.R')))
source(paste0(Code_dir, 'human_v', Ver, '.helper_functions.R'))

InitiateProject('Rivendell', Ver, Step, 'human', Project, 'ithil')
####--------------------------------------------------------------------------------------------------------------------


####--------------------------------------------------------------------------------------------------------------------
####  Global Functions  ####
####--------------------------------------------------------------------------------------------------------------------
CheckMerged <- function(srt, reduction = 'umap', assay = 'RNA'){
        srt@reductions[[1]] <- srt@reductions[[reduction]]
        p1 <- DimPlot2(srt, group.by = 'sample', reduction = reduction, split.by = 'sample',
                       cols = mycol_20, label = F, pt.size = 0.1, ncol = 4) &
                guides(col = guide_legend(ncol = 1)) &
                theme(aspect.ratio = 1)
        p3 <-  DimPlot2(srt, group.by = 'sample', reduction = reduction,
                        cols = mycol_20, label = F, pt.size = 0.1) +
                labs(x = "UMAP 1", y = "UMAP 2",
                     title = paste0(ncol(srt), " Cells x ", nrow(srt), " Genes")) +
                guides(col = guide_legend(ncol = 1))
        p <- VlnPlot2(srt,
                      features = paste0(c('nFeature_', 'nCount_', 'pct_mito_'), assay), group.by = 'sample', cols = mycol_20)
        p4 <- wrap_plots(list((
                p[[1]] + theme(axis.text.x = element_blank())),
                (p[[2]] + theme(axis.text.x = element_blank())),
                p[[3]]), ncol = 1) & theme(aspect.ratio = 0.5)
        p5 <- FeaturePlot2(srt,
                           reduction = reduction,
                           features = intersect(markers_lvl1, rownames(srt)),
                           ncol = ceiling(L(intersect(markers_lvl1, rownames(srt))) / 4))
        return(list(
                p1,
                # p2,
                p3,
                p4,
                p5
        ))
}
####--------------------------------------------------------------------------------------------------------------------


####--------------------------------------------------------------------------------------------------------------------
####  Load data  ####
####--------------------------------------------------------------------------------------------------------------------
merged.srt <- readRDS('integrated/STEP03.merged.flt.srt.rds')
merged.srt@meta.data <- readRDS('integrated/STEP04.merged.dlt.srt_meta.rds')

## Organize metadata
meta_data <- merged.srt@meta.data

## Add new meta data
meta_col <- c(
        'sample',
        'orig.name',
        'study',
        'method',
        'platform',
        'protocol',
        'data_process',
        'tissue',
        'enrichment',
        'preparation',
        'diagnosis',
        'sex_recipient',
        'age_recipient',
        'sex_donor',
        'age_donor',
        'duration',
        'recipient',
        'replicate',
        'name',
        'name2',
        'S.Score',
        'G2M.Score',
        'Phase',
        'nCount_CBN',
        'nFeature_CBN',
        'pct_mito_CBN',
        'nCount_SCT',
        'nFeature_SCT',
        'pct_mito_SCT',
        'Doublet_SC',
        'Doublet_SC_score'
)
meta_data <- meta_data[, meta_col]
merged.srt@meta.data <- meta_data

PlotPNG('01.1.dim.merge_scrublet_doublets', 10, 10)
DimPlot2(merged.srt, reduction = 'hmn_umap', group.by = 'Doublet_SC',
         raster = T, cols = c('grey75', 'red'), pt.size = 0.01)
dev.off()
####--------------------------------------------------------------------------------------------------------------------


####--------------------------------------------------------------------------------------------------------------------
####  Cleanup SCT object  ####
####--------------------------------------------------------------------------------------------------------------------
Table(merged.srt$Doublet_SC)

## Remove genes that have at least 1 UMI in less than 10 cells
x <- rowSums(merged.srt@assays$CBN@counts >= 1)
hist(x, breaks = 10000, xlim = c(0, 100))
sum(x < 10)
genes_keep <- rownames(merged.srt@assays$CBN@counts)[x >= 10]
merged.clean.srt <- DietSeurat(merged.srt,
                               counts = T, data = T, scale.data = T,
                               features = genes_keep,
                               assays = 'SCT',
                               dimreducs = names(merged.srt@reductions))

## Plot for visual inspections
plot.list <- CheckMerged(merged.clean.srt, reduction = 'hmn_umap', assay = 'SCT')
PlotPDF('02.1.merge.clean.quick_check.split_dim', 20, 20)
plot.list[[1]]
dev.off()
PlotPDF('02.2.merge.clean.quick_check.dim', 10, 10)
plot.list[[2]]
dev.off()
PlotPDF('02.3.merge.clean.quick_check.vln', 10, 10)
plot.list[[3]]
dev.off()
PlotPDF('02.4.merge.clean.quick_check.feature', 20, 20)
plot.list[[4]]
dev.off()

####--------------------------------------------------------------------------------------------------------------------
####  Cleanup CBN object  ####
####--------------------------------------------------------------------------------------------------------------------
DefaultAssay(merged.srt) <- 'CBN'
merged.clean2.srt <- DietSeurat(merged.srt,
                               counts = T, data = T, scale.data = F,
                               features = genes_keep,
                               assays = 'CBN',
                               dimreducs = names(merged.srt@reductions))
merged.clean2.srt <- NormalizeData(merged.clean2.srt)
####--------------------------------------------------------------------------------------------------------------------
saveRDS(merged.clean2.srt, 'integrated/STEP05.merged.clean_cbn.srt.rds')
####--------------------------------------------------------------------------------------------------------------------
