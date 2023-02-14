####--------------------------------------------------------------------------------------------------------------------
####  Pediatric Acute Rejection
####  2022-08-03 by Xiao LI (Texas Heart Institute, US)
####--------------------------------------------------------------------------------------------------------------------


####--------------------------------------------------------------------------------------------------------------------
####  Initiate  ####
####--------------------------------------------------------------------------------------------------------------------
Ver <- '1'
Step <- 'STEP22_Analysis_Myeloid_SCENIC'
Project <- '2022_acute_rejection_dturaga'
Code_dir <- paste0('~/Documents/Bioinformatics/project/', Project, '/code/human_v', Ver, '/')

source(Sys.readlink(paste0(Code_dir, 'src/bioinformatics.R')))
source(Sys.readlink(paste0(Code_dir, 'src/scRNAseq.R')))
source(Sys.readlink(paste0(Code_dir, 'src/geneset_human.R')))
source(paste0(Code_dir, 'human_v', Ver, '.helper_functions.R'))

InitiateProject('Rivendell', Ver, Step, 'human', Project, 'ithil')

mycol_genotype <- c('grey85', 'cyan3', 'deeppink2' ,'lightskyblue4')
mycol_group <- c('deeppink2', 'cyan3', 'lightskyblue4')

library('SCENIC')
library('AUCell')
library('RcisTarget')
library('KernSmooth')
library('BiocParallel')
library('data.table')
library('grid')
library('ComplexHeatmap')
####--------------------------------------------------------------------------------------------------------------------


####--------------------------------------------------------------------------------------------------------------------
full.srt <- readRDS('integrated/STEP18.annotated.srt.rds')
full.srt$Origin <- factor(full.srt$Fmux_genotype_clean, levels = c('Donor', 'Recipient', 'Ambiguous', 'Control'))
full.srt$Origin[is.na(full.srt$Origin)] <- 'Control'
mye.srt <- full.srt[, full.srt$Cell_type == 'Mye' &
                            full.srt$Origin != 'Ambiguous' &
                            full.srt$Cell_state != 'Mast' &
                            full.srt$Cell_type_non_ambig]
mye.srt$Cell_state <- droplevels(mye.srt$Cell_state)
mye.srt$Origin <- droplevels(mye.srt$Origin)
mye.srt$Cell_state_orig <- as.vector(mye.srt$Cell_state)
mye.srt$Cell_state_orig[mye.srt$Origin == 'Donor' & mye.srt$Cell_state == 'MP2'] <- 'Donor MP2'
mye.srt$Cell_state_orig[mye.srt$Origin == 'Recipient' & mye.srt$Cell_state == 'MP2'] <- 'Recipient MP2'
mye.srt$Cell_state_orig[mye.srt$Origin == 'Control' & mye.srt$Cell_state == 'MP2'] <- 'Control MP2'
mye.srt$Cell_state_orig <- factor(mye.srt$Cell_state_orig, levels = c(
        'Control MP2',
        'Donor MP2',
        'Recipient MP2',
        'Mono',
        'MP1'))
Table(mye.srt$Cell_state_orig)

mye_p136.srt <- mye.srt[, mye.srt$name2 == '5d PT']
mye_p136.srt$Cell_state_orig <- droplevels(mye_p136.srt$Cell_state_orig)
####--------------------------------------------------------------------------------------------------------------------


####--------------------------------------------------------------------------------------------------------------------
####  Control + 5d + 15m + 12y, Mono + MP1 + MP2  ####
####--------------------------------------------------------------------------------------------------------------------
Idents(mye.srt) <- 'Cell_state'

scenicOptions <- initializeScenic(org = "hgnc",
                                  dbDir = "~/Documents/Bioinformatics/r/db/cisTarget/",
                                  dbs = c("hg38__refseq-r80__10kb_up_and_down_tss.mc9nr.feather"),
                                  nCores = 16)

exprMat <- as.matrix(GetAssayData(mye.srt, assay = 'CBN', slot = 'count'))
genesKept <- geneFiltering(exprMat, scenicOptions)
exprMat_filtered <- exprMat[genesKept, ]
runCorrelation(exprMat_filtered, scenicOptions)
exprMat_filtered_log <- log2(exprMat_filtered + 1)
runGenie3(exprMat_filtered_log, scenicOptions) ## Extremely time consuming

exprMat_log <- log2(exprMat+1)
scenicOptions@settings$dbs <- c('hg38__refseq-r80__10kb_up_and_down_tss.mc9nr.feather')
runSCENIC_1_coexNetwork2modules(scenicOptions)
runSCENIC_2_createRegulons(scenicOptions, coexMethod=c("top50perTarget"))
scenicOptions@settings$nCores <- 1
runSCENIC_3_scoreCells(scenicOptions, exprMat_log)

cellInfo <- data.frame(Group = Idents(mye.srt))
regulonAUC <- loadInt(scenicOptions, "aucell_regulonAUC")
regulonAUC <- regulonAUC[onlyNonDuplicatedExtended(rownames(regulonAUC)),
                         Cells(mye.srt)]


regulonActivity_byCellType <- sapply(split(rownames(cellInfo), cellInfo$Group),
                                     function(cells) rowMeans(getAUC(regulonAUC)[, cells]))
regulonActivity_byCellType_Scaled <- t(scale(t(regulonActivity_byCellType), center = T, scale=T))
data <- regulonActivity_byCellType_Scaled
data[data < 0] <- 0

dev.off()
PlotPDF('01.heatmap.scenic_mye_all_cellstate', 6, 16)
pheatmap::pheatmap(data,
                   fontsize_row = 10,
                   color = mycol_RnBw,
                   breaks = seq(0, 1.5, length.out = 100),
                   treeheight_row = 20, treeheight_col = 20, border_color = NA, cluster_cols = F,
                   clustering_distance_rows = "euclidean", clustering_method = "ward.D2")
dev.off()


regulonActivity_byCellType_Scaled.nonExtend <-
        regulonActivity_byCellType_Scaled[! grepl('extended', rownames(regulonActivity_byCellType_Scaled)), ]
data <- regulonActivity_byCellType_Scaled.nonExtend
data[data<0] <- 0

dev.off()
PlotPDF('02.heatmap.scenic_mye_all_cellstate_nonextended', 6, 16)
pheatmap::pheatmap(data,
                   fontsize_row = 10,
                   color = mycol_RnBw,
                   breaks = seq(0, 1.5, length.out = 100),
                   treeheight_row = 20, treeheight_col = 20, border_color = NA, cluster_cols = F,
                   clustering_distance_rows = "euclidean", clustering_method = "ward.D2")
dev.off()
####--------------------------------------------------------------------------------------------------------------------
saveRDS(scenicOptions, 'analysis/STEP22.scenic_all_mye.scenicoptions.rds')
saveRDS(regulonAUC, 'analysis/STEP22.scenic_all_mye.regulonauc.rds')
####--------------------------------------------------------------------------------------------------------------------


####--------------------------------------------------------------------------------------------------------------------
####  5dPT ResMP vs MoMP2 vs Mono vs MP1  ####
####--------------------------------------------------------------------------------------------------------------------
Idents(mye_p136.srt) <- 'Cell_state_orig'

scenicOptions <- initializeScenic(org = "hgnc",
                                  dbDir = "~/Documents/Bioinformatics/r/db/cisTarget/",
                                  dbs = c("hg38__refseq-r80__10kb_up_and_down_tss.mc9nr.feather"),
                                  nCores = 24)

exprMat <- as.matrix(GetAssayData(mye_p136.srt, assay = 'CBN', slot = 'count'))
genesKept <- geneFiltering(exprMat, scenicOptions)
exprMat_filtered <- exprMat[genesKept, ]
runCorrelation(exprMat_filtered, scenicOptions)
exprMat_filtered_log <- log2(exprMat_filtered + 1)
runGenie3(exprMat_filtered_log, scenicOptions) ## Extremely time consuming

exprMat_log <- log2(exprMat+1)
scenicOptions@settings$dbs <- c('hg38__refseq-r80__10kb_up_and_down_tss.mc9nr.feather')
runSCENIC_1_coexNetwork2modules(scenicOptions)
runSCENIC_2_createRegulons(scenicOptions, coexMethod=c("top50perTarget"))
scenicOptions@settings$nCores <- 1
runSCENIC_3_scoreCells(scenicOptions, exprMat_log)

cellInfo <- data.frame(Group = Idents(mye_p136.srt))
regulonAUC <- loadInt(scenicOptions, "aucell_regulonAUC")
regulonAUC <- regulonAUC[onlyNonDuplicatedExtended(rownames(regulonAUC)), ]
regulonActivity_byCellType <- sapply(split(rownames(cellInfo), cellInfo$Group),
                                     function(cells) rowMeans(getAUC(regulonAUC)[, cells]))
regulonActivity_byCellType_Scaled <- t(scale(t(regulonActivity_byCellType), center = T, scale=T))
data <- regulonActivity_byCellType_Scaled
data[data < 0] <- 0

dev.off()
PlotPDF('03.heatmap.scenic_mye_all_cell_origin', 6, 16)
pheatmap(data,
         fontsize_row = 10,
         color = mycol_RnBw,
         breaks = seq(0, 1.5, length.out = 100),
         treeheight_row = 20, treeheight_col = 20, border_color = NA, cluster_cols = F,
         clustering_distance_rows = "euclidean", clustering_method = "ward.D2")
dev.off()


regulonActivity_byCellType_Scaled.nonExtend <-
        regulonActivity_byCellType_Scaled[! grepl('extended', rownames(regulonActivity_byCellType_Scaled)), ]
data <- regulonActivity_byCellType_Scaled.nonExtend
data[data < 0] <- 0

dev.off()
PlotPDF('04.heatmap.scenic_mye_all_cell_origin_nonextended', 6, 16)
pheatmap::pheatmap(data,
                   fontsize_row = 10,
                   color = mycol_RnBw,
                   breaks = seq(0, 1.5, length.out = 100),
                   treeheight_row = 20, treeheight_col = 20, border_color = NA, cluster_cols = F,
                   clustering_distance_rows = "euclidean", clustering_method = "ward.D2")
dev.off()
####--------------------------------------------------------------------------------------------------------------------
saveRDS(scenicOptions, 'analysis/STEP22.scenic_p136_mye_with_prediction.scenicoptions.rds')
saveRDS(regulonAUC, 'analysis/STEP22.scenic_p136_mye_with_prediction.regulonauc.rds')
####--------------------------------------------------------------------------------------------------------------------
