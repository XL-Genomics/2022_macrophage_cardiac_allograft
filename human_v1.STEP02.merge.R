####--------------------------------------------------------------------------------------------------------------------
####  Pediatric Acute Rejection
####  2022-08-03 by Xiao LI (Texas Heart Institute, US)
####--------------------------------------------------------------------------------------------------------------------


####--------------------------------------------------------------------------------------------------------------------
####  Initiate  ####
####--------------------------------------------------------------------------------------------------------------------
Ver <- '1'
Step <- 'STEP02_Merge'
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
        p1 <- DimPlot2(srt, group.by = 'study', reduction = reduction, split.by = 'study',
                      cols = mycol_10, label = F, pt.size = 0.1, ncol = 7) &
                guides(col = guide_legend(ncol = 1)) &
                theme(aspect.ratio = 1)
        p3 <-  DimPlot2(srt, group.by = 'study', reduction = reduction,
                       cols = mycol_10, label = F, pt.size = 0.1) +
                labs(x = "UMAP 1", y = "UMAP 2",
                     title = paste0(ncol(srt), " Cells x ", nrow(srt), " Genes")) +
                guides(col = guide_legend(ncol = 1))
        p <- VlnPlot2(srt,
                      features = paste0(c('nFeature_', 'nCount_', 'pct_mito_'), assay), group.by = 'sample')
        p4 <- wrap_plots(list((
                p[[1]] + theme(axis.text.x = element_blank())),
                (p[[2]] + theme(axis.text.x = element_blank())),
                p[[3]]), ncol = 1) & theme(aspect.ratio = 0.5)
        p5 <- FeaturePlot2(srt,
                           reduction = reduction,
                           features = intersect(markers_lvl1_min, rownames(srt)),
                           ncol = ceiling(L(intersect(markers_lvl1_min, rownames(srt))) / 4))
        return(list(
                p1,
                # p2,
                p3,
                p4,
                p5
        ))
}

Preprocess <- function(srt_obj, assay = 'RNA', ...) {
        srt.out <- srt_obj |>
                NormalizeData(verbose = F) |>
                CellCycleScoring(s.features = s.genes, g2m.features = g2m.genes, set.ident = F) |>
                PercentageFeatureSet(pattern = '^MT-', col.name = paste0('pct_mito_', assay), assay = assay)
        srt.out@meta.data[, paste0('pct_mito_', assay)][is.nan(srt.out@meta.data[, paste0('pct_mito_', assay)])] <- 0
        return(srt.out)
}
SCTransform.my <- function(srt_obj, assay = 'RNA', ...){
        srt.out <- srt_obj |>
                SCTransform(assay = assay,
                            method = "glmGamPoi",
                            seed.use = 505,
                            return.only.var.genes = F,
                            vars.to.regress = c(paste0('nFeature_', assay),
                                                paste0('nCount_', assay),
                                                paste0('pct_mito_', assay),
                                                'S.Score',
                                                'G2M.Score'),
                            ...)
        return(srt.out)
}

RunUMAP.my <- function(srt_obj, var.toal = 0.75, haromize.by, assay, ...){
        srt.out <- RunPCA(srt_obj, seed.use = 505, assay = assay)
        dimN <- FindDimNumber(srt_obj = srt.out, var.toal = var.toal, reduction = 'pca')
        srt.out <- RunUMAP(srt.out, dims = 1:dimN, seed.use = 505) |>
                RunHarmony(group.by.vars = haromize.by, assay.use = assay)
        dimN <- FindDimNumber(srt_obj = srt.out, var.toal = var.toal, reduction = 'harmony')
        srt.out <- RunUMAP(srt.out, dims = 1:dimN, seed.use = 505, reduction = 'harmony',
                           reduction.name = 'hmn_umap', reduction.key = 'hmnumap_', ...)
        return(srt.out)
}
####--------------------------------------------------------------------------------------------------------------------


####--------------------------------------------------------------------------------------------------------------------
####  Load sample metadata  ####
####--------------------------------------------------------------------------------------------------------------------
sample_meta.df <- read.csv(paste0(Docu_dir, 'pediatric_sample_meta.csv'))
studies <- U(sample_meta.df$Study)
names(studies) <- U(sample_meta.df$Study_id)
studies_cellbender <- studies[studies %in% sample_meta.df$Study[sample_meta.df$platform %in% c('10X', 'Drop-seq')]]
####--------------------------------------------------------------------------------------------------------------------


####--------------------------------------------------------------------------------------------------------------------
####  Merge raw (via Seurat)  ####
####--------------------------------------------------------------------------------------------------------------------
rds_name <- paste0(names(studies), '.', studies, '.raw.srt.rds')
data.list <- list()
ncells <- 0
for(i in 1:L(rds_name)) {
        data.list[[i]] <- readRDS(paste0('individual/', rds_name[i]))
        message(data.list[[i]]$study[1], ' ', ncol(data.list[[i]]))
        ncells <- ncells + ncol(data.list[[i]])
}
merged.raw.srt <- merge(data.list[[1]], data.list[2:L(data.list)])
c(ncells, ncol(merged.raw.srt))
rm(data.list)
gc()
####--------------------------------------------------------------------------------------------------------------------
saveRDS(merged.raw.srt@meta.data, 'integrated/STEP02.merged.raw.srt_meta.rds')
saveRDS(merged.raw.srt, 'integrated/STEP02.merged.raw.srt.rds')
####--------------------------------------------------------------------------------------------------------------------


####--------------------------------------------------------------------------------------------------------------------
####  Merge cbn datasets (via Seurat)  ####
####--------------------------------------------------------------------------------------------------------------------
rds_name <- paste0(names(studies), '.', studies, '.cbn.srt.rds')
data.list <- list()
ncells <- 0
for(i in 1:L(rds_name)) {
        data.list[[i]] <- readRDS(paste0('individual/', rds_name[i]))
        message(data.list[[i]]$study[1], ' ', ncol(data.list[[i]]))
        ncells <- ncells + ncol(data.list[[i]])
}
merged.cbn.srt <- merge(data.list[[1]], data.list[2:L(data.list)])
c(ncells, ncol(merged.cbn.srt)) ## 140536 140536
rm(data.list)
gc()

## Intersect datasets
consensus_cell <- intersect(Cells(merged.raw.srt), Cells(merged.cbn.srt))
merged.cbn.srt <- merged.cbn.srt[, consensus_cell]
consensus_gene <- intersect(rownames(merged.raw.srt), rownames(merged.cbn.srt))
merged.cbn.srt <- DietSeurat(merged.cbn.srt, features = consensus_gene)

rm(merged.raw.srt)
gc()

## Preprocess
merged.cbn.srt$pct_mito <- NULL
merged.cbn.srt <- RenameAssays(merged.cbn.srt, 'RNA' = 'CBN')
assay <- 'CBN'
merged.cbn.srt <- Preprocess(merged.cbn.srt,  assay)

## Merged SCT workflow
merged.cbn.srt <- SCTransform.my(merged.cbn.srt, assay = assay, conserve.memory = T)
merged.cbn.srt <- RunUMAP.my(merged.cbn.srt, assay = 'SCT', var.toal = 0.75, haromize.by = 'sample')
merged.cbn.srt

merged.cbn.srt$nCount_RNA <- NULL
merged.cbn.srt$nFeature_RNA <- NULL
merged.cbn.srt <- PercentageFeatureSet(merged.cbn.srt, pattern = '^MT-', col.name = paste0('pct_mito_SCT'), assay = 'SCT')
gc()


## Plot for visual inspections
plot.list <- CheckMerged(merged.cbn.srt, reduction = 'umap', assay = 'SCT')
PlotPDF('01.1.merge_cellbender_consensus.quick_check.split_dim', 20, 20)
plot.list[[1]]
dev.off()
PlotPDF('01.2.merge_cellbender_consensus.quick_check.dim', 10, 10)
plot.list[[2]]
dev.off()
PlotPDF('01.3.merge_cellbender_consensus.quick_check.vln', 10, 10)
plot.list[[3]]
dev.off()
PlotPDF('01.4.merge_cellbender_consensus.quick_check.feature', 20, 20)
plot.list[[4]]
dev.off()

## Plot for visual inspections
plot.list <- CheckMerged(merged.cbn.srt, reduction = 'hmn_umap', assay = 'SCT')
PlotPDF('02.1.merge_cellbender_consensus.quick_check.split_dim', 20, 20)
plot.list[[1]]
dev.off()
PlotPDF('02.2.merge_cellbender_consensus.quick_check.dim', 10, 10)
plot.list[[2]]
dev.off()
PlotPDF('02.3.merge_cellbender_consensus.quick_check.feature', 20, 20)
plot.list[[4]]
dev.off()
####--------------------------------------------------------------------------------------------------------------------


####--------------------------------------------------------------------------------------------------------------------
####  Organize meta data  ####
####--------------------------------------------------------------------------------------------------------------------
study_order <- c(
        '2019_Unpub7_JMartin',
        '2021_Circulation_EPorrello',
        '2022_Unpub2_DTuraga',
        '2022_Unpub3_DTuraga'
)
sample_order <- str_sort(U(merged.cbn.srt$sample))
method_order <- c(
        'Nucleus'
)
platform_order <- c(
        '10X'
)
protocol_order <- c(
        'GEM_3p_v3'
)
tissue_order <- c(
        'Left_ventricle',
        'Right_ventricle'
)
enrichment_order <- c(
        'All_cells'
)
sex_recipient_order <- c(
        'F',
        'M'
)
sex_donor_order <- c(
        'F',
        'M',
        'Unknown',
        'None'
)
age_recipient_order <- sort(U(merged.cbn.srt$age_recipient))
age_donor_order <- c(
        7,
        20,
        'Unknown',
        'None'
)
diagnosis_order <- c(
        'Healthy',
        'PGF',
        'CR'
)
recipient_order <- c(
        '3_198',
        '13_235',
        'Sample12',
        'Sample13',
        'P136',
        'P182',
        'P170'
)
duration_order <- sort(U(merged.cbn.srt$duration))

## Edit processed info
merged.cbn.srt$data_process <- 'CellBender'
merged.cbn.srt$preparation <- NA
merged.cbn.srt$preparation[merged.cbn.srt$study == '2019_Unpub7_JMartin'] <- 'Density_gradient_nucleus_isolation'
merged.cbn.srt$preparation[merged.cbn.srt$study == '2021_Circulation_EPorrello'] <- 'Unknown'
merged.cbn.srt$preparation[merged.cbn.srt$study == '2022_Unpub2_DTuraga'] <- 'FACS_nucleus_isolation'
merged.cbn.srt$preparation[merged.cbn.srt$study == '2022_Unpub3_DTuraga'] <- 'FACS_nucleus_isolation'

## Set merged data ordering
merged.cbn.srt$study <- factor(merged.cbn.srt$study, levels = study_order)
merged.cbn.srt$sample <- factor(merged.cbn.srt$sample, levels = sample_order)
merged.cbn.srt$method <- factor(merged.cbn.srt$method, levels = method_order)
merged.cbn.srt$platform <- factor(merged.cbn.srt$platform, levels = platform_order)
merged.cbn.srt$protocol <- factor(merged.cbn.srt$protocol, levels = protocol_order)
merged.cbn.srt$tissue <- factor(merged.cbn.srt$tissue, levels = tissue_order)
merged.cbn.srt$enrichment <- factor(merged.cbn.srt$enrichment, levels = enrichment_order)
merged.cbn.srt$sex_recipient <- factor(merged.cbn.srt$sex_recipient, levels = sex_recipient_order)
merged.cbn.srt$age_recipient <- factor(merged.cbn.srt$age_recipient, levels = age_recipient_order)
merged.cbn.srt$sex_donor <- factor(merged.cbn.srt$sex_donor, levels = sex_donor_order)
merged.cbn.srt$age_donor <- factor(merged.cbn.srt$age_donor, levels = age_donor_order)
merged.cbn.srt$diagnosis <- factor(merged.cbn.srt$diagnosis, levels = diagnosis_order)
merged.cbn.srt$duration <- factor(merged.cbn.srt$duration, levels = duration_order)
merged.cbn.srt$orig.ident <- NULL
merged.cbn.srt$name <- revalue(merged.cbn.srt$sample, replace = c(
        P1_S001 = 'Control 1',
        P1_S002 = 'Control 1',
        P1_S003 = 'Control 1',
        P1_S004 = 'Control 1',
        P1_S005 = 'Control 2',
        P1_S006 = 'Control 2',
        P2_S001 = 'Control 3',
        P2_S002 = 'Control 4',
        P3_S001 = '5d PT',
        P3_S002 = '15m PT',
        P4_S001 = '12y PT'
))
merged.cbn.srt$name <- factor(merged.cbn.srt$name, levels = c(
        'Control 1',
        'Control 2',
        'Control 3',
        'Control 4',
        '5d PT',
        '15m PT',
        '12y PT'
))

merged.cbn.srt$name2 <- revalue(merged.cbn.srt$name, replace = c(
        'Control 1' = 'Control',
        'Control 2' = 'Control',
        'Control 3' = 'Control',
        'Control 4' = 'Control'
))

## Checking basic stat of the merged dataset
stat <- c('study', 'name', 'name2', 'sex_recipient')
col_choice <- list(mycol_10,
                   mycol_20,
                   mycol_10,
                   c(mycol_10[c(2, 1)])
                   )
p.list <- list()
for(i in 1:4){
        data <- data.frame(group = factor(names(table(merged.cbn.srt@meta.data[, stat[i]])),
                                          levels = levels(merged.cbn.srt@meta.data[, stat[i]])),
                           value = table(merged.cbn.srt@meta.data[, stat[i]]))
        p.list[[i]] <- ggplot(data, aes(x = "", y = value.Freq, fill = group)) +
                geom_bar(stat = "identity", width = 1) +
                coord_polar("y", start = 0, direction = -1) +
                labs(title = stat[i]) +
                scale_fill_manual(values = col_choice[[i]]) +
                WhiteBackground() +
                NoAxes() +
                theme(aspect.ratio = 1, legend.title = element_blank())
}

PlotPDF('04.merged_consensus.meta_stat', 10, 10, onefile = T)
p.list
dev.off()

merged.cbn.srt <- DietSeurat(merged.cbn.srt, assays = c('SCT', 'CBN'), dimreducs = names(merged.cbn.srt@reductions))
####--------------------------------------------------------------------------------------------------------------------
saveRDS(merged.cbn.srt@meta.data, 'integrated/STEP02.merged.cbn.srt_meta.rds')
saveRDS(merged.cbn.srt, 'integrated/STEP02.merged.cbn.srt.rds')
####--------------------------------------------------------------------------------------------------------------------
