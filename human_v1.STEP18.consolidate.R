####--------------------------------------------------------------------------------------------------------------------
####  Pediatric Acute Rejection
####  2022-08-03 by Xiao LI (Texas Heart Institute, US)
####--------------------------------------------------------------------------------------------------------------------


####--------------------------------------------------------------------------------------------------------------------
####  Initiate  ####
####--------------------------------------------------------------------------------------------------------------------
Ver <- '1'
Step <- 'STEP18_Consolidate_Annotation'
Project <- '2022_acute_rejection_dturaga'
Code_dir <- paste0('~/Documents/Bioinformatics/project/', Project, '/code/human_v', Ver, '/')

source(Sys.readlink(paste0(Code_dir, 'src/bioinformatics.R')))
source(Sys.readlink(paste0(Code_dir, 'src/scRNAseq.R')))
source(Sys.readlink(paste0(Code_dir, 'src/geneset_human.R')))
source(paste0(Code_dir, 'human_v', Ver, '.helper_functions.R'))

InitiateProject('Rivendell', Ver, Step, 'human', Project, 'ithil')

Color_cell_type <- c(mycol_40[c(3, 6, 11, 14, 35, 17, 22, 26, 30)], 'grey50', 'grey80', 'grey90')
####--------------------------------------------------------------------------------------------------------------------


####--------------------------------------------------------------------------------------------------------------------
####  Load data  ####
####--------------------------------------------------------------------------------------------------------------------
srt <- readRDS('integrated/STEP07.scvi_integration.srt.rds')
meta <- readRDS('integrated/STEP09.annotated.srt_meta.rds')

srt <- AddMetaData(srt, metadata = meta)

srt$Cell_state <- NA
srt$Cell_type <- as.vector(srt$Cell_type)

Table(srt$Cell_state, srt$Cell_type)

srt@reductions$sub_scANVI_umap <- srt@reductions$scVI_umap
srt@reductions$sub_scANVI_umap@cell.embeddings[,1] <- NA
srt@reductions$sub_scANVI_umap@cell.embeddings[,2] <- NA
srt@reductions$sub_scANVI_umap@key <- 'subSCANVIUMAP_'
colnames(srt@reductions$sub_scANVI_umap@cell.embeddings) <- c('subSCANVIUMAP_1', 'subSCANVIUMAP_2')
####--------------------------------------------------------------------------------------------------------------------


####--------------------------------------------------------------------------------------------------------------------
####  Consolidate immune cells  ####
####--------------------------------------------------------------------------------------------------------------------
imm.srt <- readRDS('integrated/STEP10.immune_cells.srt.rds')

srt$Cell_state[Cells(imm.srt)] <- imm.srt$Cell_state
Table(srt$Cell_state, srt$Cell_type)

mye.srt <- readRDS('integrated/STEP10.myeloid_cells.srt.rds')
srt$Cell_state[Cells(mye.srt)] <- mye.srt$Cell_state
srt$Cell_type[srt$Cell_state %in% c('MP1', 'MP2', 'Mono', 'Mast', 'Mito')] <- 'Mye'
Table(srt$Cell_state, srt$Cell_type)
srt@reductions$sub_scANVI_umap@cell.embeddings[Cells(mye.srt),] <- mye.srt@reductions$hmn_umap@cell.embeddings

lym.srt <- readRDS('integrated/STEP10.lymphoid_cells.srt.rds')
srt$Cell_state[Cells(lym.srt)] <- as.vector(lym.srt$Cell_state)
srt$Cell_type[srt$Cell_state %in% c('CD4 TC', 'CD8 TC', 'NK1', 'NK2', 'NK3', 'BC')] <- 'Lym'
Table(srt$Cell_state, srt$Cell_type)
srt@reductions$sub_scANVI_umap@cell.embeddings[Cells(lym.srt),] <- lym.srt@reductions$hmn_umap@cell.embeddings

DimPlot2(srt, reduction = 'sub_scANVI_umap', group.by = 'Cell_state', split.by = 'Cell_type', ncol = 3)
####--------------------------------------------------------------------------------------------------------------------


####--------------------------------------------------------------------------------------------------------------------
####  Consolidate CMs  ####
####--------------------------------------------------------------------------------------------------------------------
cm.srt <- readRDS('integrated/STEP11.cm_cells.srt.rds')

srt$Cell_state[Cells(cm.srt)] <- cm.srt$Cell_state

Table(srt$Cell_state, srt$Cell_type)

srt@reductions$sub_scANVI_umap@cell.embeddings[Cells(cm.srt),] <- cm.srt@reductions$sub_scANVI_umap@cell.embeddings
DimPlot2(srt, reduction = 'sub_scANVI_umap', group.by = 'Cell_state', split.by = 'Cell_type', ncol = 3)
####--------------------------------------------------------------------------------------------------------------------


####--------------------------------------------------------------------------------------------------------------------
####  Consolidate FBs  ####
####--------------------------------------------------------------------------------------------------------------------
fb.srt <- readRDS('integrated/STEP12.fb_cells.srt.rds')

srt$Cell_state[Cells(fb.srt)] <- fb.srt$Cell_state

Table(srt$Cell_state, srt$Cell_type)

srt@reductions$sub_scANVI_umap@cell.embeddings[Cells(fb.srt),] <- fb.srt@reductions$sub_scANVI_umap@cell.embeddings
DimPlot2(srt, reduction = 'sub_scANVI_umap', group.by = 'Cell_state', split.by = 'Cell_type', ncol = 3)
####--------------------------------------------------------------------------------------------------------------------


####--------------------------------------------------------------------------------------------------------------------
####  Consolidate ECs  ####
####--------------------------------------------------------------------------------------------------------------------
ec.srt <- readRDS('integrated/STEP13.ec_cells.srt.rds')

srt$Cell_state[Cells(ec.srt)] <- ec.srt$Cell_state

Table(srt$Cell_state, srt$Cell_type)

srt@reductions$sub_scANVI_umap@cell.embeddings[Cells(ec.srt),] <- ec.srt@reductions$sub_scANVI_umap@cell.embeddings
DimPlot2(srt, reduction = 'sub_scANVI_umap', group.by = 'Cell_state', split.by = 'Cell_type', ncol = 3)
####--------------------------------------------------------------------------------------------------------------------


####--------------------------------------------------------------------------------------------------------------------
####  Consolidate PCs  ####
####--------------------------------------------------------------------------------------------------------------------
pc.srt <- readRDS('integrated/STEP14.pc_cells.srt.rds')

srt$Cell_state[Cells(pc.srt)] <- pc.srt$Cell_state

Table(srt$Cell_state, srt$Cell_type)

srt@reductions$sub_scANVI_umap@cell.embeddings[Cells(pc.srt),] <- pc.srt@reductions$sub_scANVI_umap@cell.embeddings

DimPlot2(srt, reduction = 'sub_scANVI_umap', group.by = 'Cell_state', split.by = 'Cell_type', ncol = 3)
####--------------------------------------------------------------------------------------------------------------------


####--------------------------------------------------------------------------------------------------------------------
####  Consolidate SMCs  ####
####--------------------------------------------------------------------------------------------------------------------
smc.srt <- readRDS('integrated/STEP15.smc_cells.srt.rds')

srt$Cell_state[Cells(smc.srt)] <- smc.srt$Cell_state

Table(srt$Cell_state, srt$Cell_type)

srt@reductions$sub_scANVI_umap@cell.embeddings[Cells(smc.srt),] <- smc.srt@reductions$sub_scANVI_umap@cell.embeddings
DimPlot2(srt, reduction = 'sub_scANVI_umap', group.by = 'Cell_state', split.by = 'Cell_type', ncol = 3)
####--------------------------------------------------------------------------------------------------------------------


####--------------------------------------------------------------------------------------------------------------------
####  Organize order  ####
####--------------------------------------------------------------------------------------------------------------------
Table(srt$Cell_state, srt$Cell_type)
srt$Cell_state[srt$Cell_type %in% c("Adi" ,'Neu', 'Mito')] <-
        srt$Cell_type[srt$Cell_type %in% c("Adi" ,'Neu', 'Mito')]

srt$Cell_type <- factor(srt$Cell_type,
                        levels = c('CM', 'FB', 'EC', 'PC', 'SMC', 'Lym', 'Mye', 'Neu', 'Adi', 'Mito', 'Doublet'))

srt$Cell_state <- factor(srt$Cell_state, levels = c(
        'CM1', 'CM2', 'CM3',
        'FB1', 'FB2', 'FB3', 'FB4',
        'EC1_cap', 'EC2_veno', 'EC3_art', 'EC4_endoc', 'EC5_lec',
        'PC1', 'PC2', 'PC3', 'PC4',
        'SMC1', 'SMC2', 'SMC3',
        'CD4 TC', 'CD8 TC', 'NK1', 'NK2', 'NK3', 'BC',
        'Mono', 'MP1', 'MP2', 'Mast',
        'Neu', 'Adi',
        'Mito', 'Doublet'
))
Table(srt$Cell_state, srt$Cell_type)

dev.off()
PlotPDF('01.1.heat.celltype_vs_cellstate', 5, 5)
pheatmap(log(Table(srt$Cell_state, srt$Cell_type)+1), cluster_rows = F, cluster_cols = F)
dev.off()

p1 <- DimPlot2(srt, group.by = 'Cell_type', cols = Color_cell_type)
p2 <- DimPlot2(srt, group.by = 'Cell_state', cols = c(mycol_30, 'grey20', 'grey80', 'grey90'), label = T)
PlotPDF('01.2.umap.annotated', 12, 5)
p1+p2
dev.off()

srt$Cell_type_non_ambig <- T
srt$Cell_type_non_ambig[srt$Cell_type %in% c('Doublet', 'Mito')] <- F
srt$Cell_type_non_ambig[srt$Cell_state %in% c('Doublet', 'Mito')] <- F
Table(srt$Cell_type_non_ambig, srt$Cell_type)
Table(srt$Cell_type_non_ambig, srt$Cell_state)
####--------------------------------------------------------------------------------------------------------------------


####--------------------------------------------------------------------------------------------------------------------
####  Re-embed without ambiguous cells
####--------------------------------------------------------------------------------------------------------------------
clean.srt <- srt[, srt$Cell_type_non_ambig]
clean.srt <- RunUMAP(clean.srt, reduction = 'scVI', dims = 1:30, n.neighbors = 50, min.dist = 0.5,
                     reduction.name = 'clean_umap', reduction.key = 'cUMAP_')
DimPlot2(clean.srt, reduction = 'clean_umap', group.by = 'Cell_type',
         cols = Color_cell_type, label = T, raster = T) /
DimPlot2(clean.srt, reduction = 'clean_umap', group.by = 'Cell_state',
         cols = c(mycol_30, 'grey'), label = T, raster = T)

## Clean umap without doublets or mito-high cells
srt@reductions$clean_umap <- srt@reductions$scVI_umap
srt@reductions$clean_umap@cell.embeddings[,1] <- NA
srt@reductions$clean_umap@cell.embeddings[,2] <- NA
srt@reductions$clean_umap@cell.embeddings[Cells(clean.srt), 1:2] <- clean.srt@reductions$clean_umap@cell.embeddings[,1:2]
srt@reductions$clean_umap@key <- 'cUMAP_'
colnames(srt@reductions$clean_umap@cell.embeddings) <- c('cUMAP_1', 'cUMAP_2')

## Full umap with doublets and mito-high cells
srt@reductions$full_umap <- srt@reductions$scANVI_umap
srt@reductions$full_umap@key <- 'fUMAP_'
colnames(srt@reductions$full_umap@cell.embeddings) <- c('fUMAP_1', 'fUMAP_2')

## Sub-umap for CM, FB, NK/T, myeloid, EC, PC and SMC
srt@reductions$sub_umap <- srt@reductions$sub_scANVI_umap
srt@reductions$sub_umap@key <- 'subUMAP_'
colnames(srt@reductions$sub_umap@cell.embeddings) <- c('subUMAP_1', 'subUMAP_2')

all_reductions <- srt@reductions

srt2 <- DietSeurat(srt, counts = T, data = T, scale.data = F,
                   assays = 'CBN',
                   dimreducs = c('sub_umap', 'full_umap', 'clean_umap'))


p1 <- DimPlot2(srt2, reduction = 'full_umap', group.by = 'Cell_type', cols = Color_cell_type, label = T)
p2 <- DimPlot2(srt2, reduction = 'full_umap', group.by = 'Cell_state',
               cols = c(mycol_30, 'grey20', 'grey80', 'grey90'), label = T)
p3 <- DimPlot2(srt2, reduction = 'clean_umap', group.by = 'Cell_type', cols = Color_cell_type, label = T)
p4 <- DimPlot2(srt2, reduction = 'clean_umap', group.by = 'Cell_state',
               cols = c(mycol_30, 'grey20', 'grey80', 'grey90'), label = T)
p <- wrap_plots(p1, p2, p3, p4, ncol = 2)
PlotPDF('02.1.umap.annotation_clean_umap', 12, 10)
p
dev.off()

i <- c('CM', 'FB', 'EC', 'PC', 'SMC', 'Lym', 'Mye')
plist <- list()
for (j in 1:L(i)){
        message(i[j])
        plist[[j]] <-  DimPlot2(srt2[, srt$Cell_type == i[j]],
                                reduction = 'sub_umap', group.by = 'Cell_state', cols = mycol_10, label = T)
}
p <- wrap_plots(plist, ncol = 3)
PlotPDF('02.2.umapannotation_subumap', 15, 15)
p
dev.off()
####--------------------------------------------------------------------------------------------------------------------


####--------------------------------------------------------------------------------------------------------------------
####  Markers  ####
####--------------------------------------------------------------------------------------------------------------------
clean.srt <- srt2[, srt2$Cell_type_non_ambig]
clean.srt$Cell_type <- droplevels(clean.srt$Cell_type)
Table(clean.srt$Cell_type)

Idents(clean.srt) <- 'Cell_type'
Cell_type_marker <- FindAllMarkers(clean.srt,
                                   assay = 'CBN',
                                   test.use = 'wilcox',
                                   only.pos = T,
                                   max.cells.per.ident = 1000,
                                   logfc.threshold = 1,
                                   random.seed = 505)

Cell_type_marker.dedup <- Cell_type_marker |> group_by(cluster) |> top_n(25, wt = avg_log2FC)
Cell_type_marker.dedup <- Cell_type_marker.dedup[!duplicated(Cell_type_marker.dedup$gene),]
Table(Cell_type_marker.dedup$cluster)

p <- MarkerHeatmap(clean.srt, marker.df = Cell_type_marker.dedup, top = 20)
PlotPDF('03.1.heat.cell_type_marker', 15, 20)
print(p)
dev.off()

clean.srt$Cell_state <- droplevels(clean.srt$Cell_state)
Table(clean.srt$Cell_state)
Idents(clean.srt) <- 'Cell_state'
Cell_state_marker <- FindAllMarkers(clean.srt,
                                    assay = 'CBN',
                                    test.use = 'wilcox',
                                    only.pos = T,
                                    max.cells.per.ident = 1000,
                                    logfc.threshold = 0.25,
                                    random.seed = 505)

Cell_state_marker.dedup <- Cell_state_marker |> group_by(cluster) |> top_n(24,  wt = avg_log2FC)
Cell_state_marker.dedup <- Cell_state_marker.dedup[!duplicated(Cell_state_marker.dedup$gene),]
sort(Table(Cell_state_marker.dedup$cluster))

p <- MarkerHeatmap(clean.srt, marker.df = Cell_state_marker.dedup, top = 5)
PlotPDF('03.2.heat.cell_state_marker', 15, 25)
print(p)
dev.off()

srt2@misc$marker <- list()
srt2@misc$marker$Cell_type_marker <- Cell_type_marker
srt2@misc$marker$Cell_state_marker <- Cell_state_marker
####--------------------------------------------------------------------------------------------------------------------


####--------------------------------------------------------------------------------------------------------------------
saveRDS(all_reductions, 'integrated/STEP18.all_reductions.srt_dimreducs.rds')
####--------------------------------------------------------------------------------------------------------------------


####--------------------------------------------------------------------------------------------------------------------
####  Donor recipient origin  ####
####--------------------------------------------------------------------------------------------------------------------
P136_freemuxlet <- read.table(gzfile(
        paste0('/Volumes/shire/data/scrnaseq/2022_Unpub2_DTuraga/matrix/popscle/',
               '2022_Unpub2_DTuraga_P136_s1/2022_Unpub2_DTuraga_P136_s1.freemux.clust1.samples.gz')),
        header = T)
rownames(P136_freemuxlet) <- paste0('2022_Unpub2_DTuraga:P3_S001:', P136_freemuxlet$BARCODE)
P170_freemuxlet <- read.table(gzfile(
        paste0('/Volumes/shire/data/scrnaseq/2022_Unpub2_DTuraga/matrix/popscle/',
               '2022_Unpub2_DTuraga_P170_s1/2022_Unpub2_DTuraga_P170_s1.freemux.clust1.samples.gz')),
        header = T)
rownames(P170_freemuxlet) <- paste0('2022_Unpub2_DTuraga:P3_S002:', P170_freemuxlet$BARCODE)
P182_freemuxlet <- read.table(gzfile(
        paste0('/Volumes/shire/data/scrnaseq/2022_Unpub3_DTuraga/matrix/popscle/',
               '2022_Unpub3_DTuraga_P182/2022_Unpub3_DTuraga_P182.freemux.clust1.samples.gz')),
        header = T)
rownames(P182_freemuxlet) <- paste0('2022_Unpub3_DTuraga:P4_S001:', P182_freemuxlet$BARCODE)
freemuxlet <- rbind(P136_freemuxlet, P170_freemuxlet, P182_freemuxlet)
O(Cells(srt2), rownames(freemuxlet))

freemuxlet <- freemuxlet[Cells(srt2), ]

srt2$Fmux_doublet <- 'NA'
srt2$Fmux_doublet <- factor(freemuxlet$DROPLET.TYPE, levels = c('SNG', 'DBL', 'NA'))
srt2$Fmux_doublet[is.na(srt2$Fmux_doublet)] <- 'NA'

srt2$Fmux_genotype <- 'NA'
srt2$Fmux_genotype <- revalue(freemuxlet$BEST.GUESS, c('0,0' = 'Donor',
                                                       '1,1' = 'Recipient',
                                                       '1,0' = 'Ambiguous'))
srt2$Fmux_genotype <- factor(srt2$Fmux_genotype, levels = c('Donor', 'Recipient', 'Ambiguous', 'NA'))
srt2$Fmux_genotype[is.na(srt2$Fmux_genotype)] <- 'NA'

srt2$Fmux_genotype_clean <- srt2$Fmux_genotype
srt2$Fmux_genotype_clean[(! srt2$Cell_type_non_ambig) & srt2$Fmux_doublet != 'NA'] <- 'Ambiguous'

Table(srt2$Fmux_genotype_clean, srt2$Fmux_doublet)


#### Correcting a critical error in meta data:
srt2$name <- revalue(srt2$name, replace = c(
        '15m PT' = '12y PT',
        '12y PT' = '15m PT'
))
srt2$name <- factor(srt2$name, levels = c(paste("Control", 1:4), "5d PT", "15m PT", "12y PT"))
Table(srt2$name, srt2$recipient)

srt2$name2 <- revalue(srt2$name2, replace = c(
        '15m PT' = '12y PT',
        '12y PT' = '15m PT'
))
srt2$name2 <- factor(srt2$name2, levels = c("Control", "5d PT", "15m PT", "12y PT"))
Table(srt2$name2, srt2$recipient)
####--------------------------------------------------------------------------------------------------------------------
saveRDS(srt2, 'integrated/STEP18.annotated.srt.rds')
####--------------------------------------------------------------------------------------------------------------------
