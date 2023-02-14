####--------------------------------------------------------------------------------------------------------------------
####  Pediatric Acute Rejection
####  2022-08-03 by Xiao LI (Texas Heart Institute, US)
####--------------------------------------------------------------------------------------------------------------------


####--------------------------------------------------------------------------------------------------------------------
####  Initiate  ####
####--------------------------------------------------------------------------------------------------------------------
Ver <- '1'
Step <- 'STEP21_Analysis_immune_velocity'
Project <- '2022_acute_rejection_dturaga'
Code_dir <- paste0('~/Documents/Bioinformatics/project/', Project, '/code/human_v', Ver, '/')

source(Sys.readlink(paste0(Code_dir, 'src/bioinformatics.R')))
source(Sys.readlink(paste0(Code_dir, 'src/scRNAseq.R')))
source(Sys.readlink(paste0(Code_dir, 'src/geneset_human.R')))
source(paste0(Code_dir, 'human_v', Ver, '.helper_functions.R'))

InitiateProject('Rivendell', Ver, Step, 'human', Project, 'ithil')
library('velocyto.R')
####--------------------------------------------------------------------------------------------------------------------


####--------------------------------------------------------------------------------------------------------------------
####  Prepare loom file from 10x cellranger run on 'Bobbyd' server  ####
####--------------------------------------------------------------------------------------------------------------------
## Shell command:
## conda activate velocyto
## cd ~/work/rejection/matrix/
## velocyto run10x \
##      -@ 24 \
##      --samtools-memory 10240 \
##      -m ~/work/genome/velocyto/hg38_rmsk.gtf \
##      ./2022_Unpub2_DTuraga_P136_s1 \
##      ~/work/genome/index_cellranger/refdata-gex-GRCh38-2020-A/genes/genes.gtf
####--------------------------------------------------------------------------------------------------------------------


####--------------------------------------------------------------------------------------------------------------------
####  P136 sample  ####
####--------------------------------------------------------------------------------------------------------------------
##  Convert Loom to Seurat
loom.data <- ReadVelocity(paste0('/Volumes/shire/data/scrnaseq/2022_Unpub2_DTuraga/',
                                 'matrix/velocyto/2022_Unpub2_DTuraga_P136_s1.loom'))
srtLoom <- as.Seurat(loom.data)

##  Load annotated Seurat
srt <- readRDS('integrated/STEP18.annotated.srt.rds')
srt <- srt[, srt$name == '5d PT' & srt$Fmux_genotype_clean != 'Ambiguous']
DimPlot(srt, reduction = 'clean_umap')

##  Match cell names
names <- Cells(srtLoom)
names <- str_replace(names, '2022_Unpub2_DTuraga_P136_s1', '2022_Unpub2_DTuraga:P3_S001')
names <- str_replace(names, 'x$', '-1')
srtLoom <- RenameCells(srtLoom, new.names = names)
O(Cells(srtLoom), Cells(srt))
srtLoom <- srtLoom[, Cells(srt)]
srtLoom <- AddMetaData(srtLoom, metadata = srt@meta.data)

##  Process Loom converted Seurat and Run Velocity
dimNum <- 10
srtLoom <- SCTransform(srtLoom, assay = "spliced", seed.use = 505) %>%
        RunPCA(verbose = F, seed.use = 505) %>%
        RunUMAP(dims = 1:dimNum, seed.use = 505)
DimPlot(srtLoom, group.by = 'Cell_type')

srtLoom <- RunVelocity(srtLoom, deltaT = 1, kCells = 25, fit.quantile = 0.02, ncores = 16)

##  Plot myeloid cells
colors <- mycol_10[3:6]
mye.srt <- srt[, srt$Cell_type == 'Mye' & srt$Cell_type_non_ambig]
DimPlot2(mye.srt, group.by = 'Cell_state', cols = colors, reduction = 'sub_umap')

srtLoom_imm <- srtLoom[, Cells(mye.srt)]
srtLoom_imm$Cell_state <- droplevels(srtLoom_imm$Cell_state)

names(colors) <- levels(srtLoom_imm$Cell_state)
cell.colors <- colors[srtLoom_imm$Cell_state[colnames(srtLoom_imm)]]
names(cell.colors) <- colnames(srtLoom_imm)

mye.vlc <- show.velocity.on.embedding.cor(emb = Embeddings(mye.srt, reduction = "sub_umap"),
                                          vel = Tool(srtLoom_imm, slot = "RunVelocity"),
                                          scale = 'sqrt',
                                          n.cores = 16)

PlotPDF('01.velocity.p136_myeloid_sub_umap_projection', 8, 8)
show.velocity.on.embedding.cor(emb = Embeddings(mye.srt, reduction = "sub_umap"),
                               vel = Tool(srtLoom_imm, slot = "RunVelocity"),
                               cc = mye.vlc$cc,
                               cell.colors = ac(cell.colors, alpha = 0.5),
                               cell.border.alpha = 0.04,
                               show.grid.flow = T,
                               cex = 1,
                               n = 1000,
                               grid.n = 40,
                               arrow.scale = 2,
                               arrow.lwd = 1.5,
                               min.grid.cell.mass = 6,
                               do.par = F,
                               asp = 0.75)
dev.off()
####--------------------------------------------------------------------------------------------------------------------
saveRDS(mye.vlc, 'analysis/STEP21.p136_myeloid_velocity_embedding.vlc.rds')
saveRDS(srtLoom_imm, 'analysis/STEP21.p136_myeloid_velocity.srt.rds')
saveRDS(srtLoom, 'analysis/STEP21.p136_allcell_velocity.srt.rds')
####--------------------------------------------------------------------------------------------------------------------



####--------------------------------------------------------------------------------------------------------------------
####  p182 sample  ####
####--------------------------------------------------------------------------------------------------------------------
##  Convert Loom to Seurat
loom.data <- ReadVelocity(paste0('/Volumes/shire/data/scrnaseq/2022_Unpub3_DTuraga/',
                                 'matrix/velocyto/2022_Unpub3_DTuraga_P182.loom'))
srtLoom <- as.Seurat(loom.data)

##  Load annotated Seurat
srt <- readRDS('integrated/STEP18.annotated.srt.rds')
srt <- srt[, srt$name == '15m PT']
DimPlot(srt, reduction = 'clean_umap')

##  Match cell names
names <- Cells(srtLoom)
names <- str_replace(names, '2022_Unpub3_DTuraga_P182', '2022_Unpub3_DTuraga:P4_S001')
names <- str_replace(names, 'x$', '-1')
srtLoom <- RenameCells(srtLoom, new.names = names)
O(Cells(srtLoom), Cells(srt))
srtLoom <- srtLoom[, Cells(srt)]
srtLoom <- AddMetaData(srtLoom, metadata = srt@meta.data)

##  Process Loom converted Seurat and Run Velocity
dimNum <- 10
srtLoom <- SCTransform(srtLoom, assay = "spliced", seed.use = 505) %>%
        RunPCA(verbose = F, seed.use = 505) %>%
        RunUMAP(dims = 1:dimNum, seed.use = 505)
DimPlot(srtLoom, group.by = 'Cell_type')

srtLoom <- RunVelocity(srtLoom, deltaT = 1, kCells = 25, fit.quantile = 0.02, ncores = 8)

##  Plot myeloid cells
colors <- mycol_10[3:6]
mye.srt <- srt[, srt$Cell_type == 'Mye' & srt$Cell_type_non_ambig]
DimPlot2(mye.srt, group.by = 'Cell_state', cols = colors, reduction = 'sub_umap')

srtLoom_imm <- srtLoom[, Cells(mye.srt)]
srtLoom_imm$Cell_state <- droplevels(srtLoom_imm$Cell_state)

names(colors) <- levels(srtLoom_imm$Cell_state)
cell.colors <- colors[srtLoom_imm$Cell_state[colnames(srtLoom_imm)]]
names(cell.colors) <- colnames(srtLoom_imm)

mye.vlc <- show.velocity.on.embedding.cor(emb = Embeddings(mye.srt, reduction = "sub_umap"),
                                          vel = Tool(srtLoom_imm, slot = "RunVelocity"),
                                          scale = 'sqrt',
                                          n.cores = 16)

PlotPDF('02.velocity.p182_myeloid_sub_umap_projection', 8, 8)
show.velocity.on.embedding.cor(emb = Embeddings(mye.srt, reduction = "sub_umap"),
                               vel = Tool(srtLoom_imm, slot = "RunVelocity"),
                               cc = mye.vlc$cc,
                               cell.colors = ac(cell.colors, alpha = 0.5),
                               cell.border.alpha = 0.04,
                               show.grid.flow = T,
                               cex = 1,
                               n = 1000,
                               grid.n = 30,
                               arrow.scale = 2,
                               arrow.lwd = 1.5,
                               min.grid.cell.mass = 4,
                               do.par = F,
                               asp = 0.75)
dev.off()
####--------------------------------------------------------------------------------------------------------------------
saveRDS(mye.vlc, 'analysis/STEP21.p182_myeloid_velocity_embedding.vlc.rds')
saveRDS(srtLoom_imm, 'analysis/STEP21.p182_myeloid_velocity.srt.rds')
saveRDS(srtLoom, 'analysis/STEP21.p182_allcell_velocity.srt.rds')
####--------------------------------------------------------------------------------------------------------------------


####--------------------------------------------------------------------------------------------------------------------
####  p170 sample  ####
####--------------------------------------------------------------------------------------------------------------------
##  Convert Loom to Seurat
loom.data <- ReadVelocity(paste0('/Volumes/shire/data/scrnaseq/2022_Unpub2_DTuraga/',
                                 'matrix/velocyto/2022_Unpub2_DTuraga_P170_s1.loom'))
srtLoom <- as.Seurat(loom.data)

##  Load annotated Seurat
srt <- readRDS('integrated/STEP18.annotated.srt.rds')
srt <- srt[, srt$name == '12y PT']
DimPlot(srt, reduction = 'clean_umap')

##  Match cell names
names <- Cells(srtLoom)
names <- str_replace(names, '2022_Unpub2_DTuraga_P170_s1', '2022_Unpub2_DTuraga:P3_S002')
names <- str_replace(names, 'x$', '-1')
srtLoom <- RenameCells(srtLoom, new.names = names)
O(Cells(srtLoom), Cells(srt))
srtLoom <- srtLoom[, Cells(srt)]
srtLoom <- AddMetaData(srtLoom, metadata = srt@meta.data)

##  Process Loom converted Seurat and Run Velocity
dimNum <- 10
srtLoom <- SCTransform(srtLoom, assay = "spliced", seed.use = 505) %>%
        RunPCA(verbose = F, seed.use = 505) %>%
        RunUMAP(dims = 1:dimNum, seed.use = 505)
DimPlot(srtLoom, group.by = 'Cell_type')

srtLoom <- RunVelocity(srtLoom, deltaT = 1, kCells = 25, fit.quantile = 0.02, ncores = 8)

##  Plot myeloid cells
colors <- mycol_10[3:6]
mye.srt <- srt[, srt$Cell_type == 'Mye' & srt$Cell_type_non_ambig]
DimPlot2(mye.srt, group.by = 'Cell_state', cols = colors, reduction = 'sub_umap')

srtLoom_imm <- srtLoom[, Cells(mye.srt)]
srtLoom_imm$Cell_state <- droplevels(srtLoom_imm$Cell_state)

names(colors) <- levels(srtLoom_imm$Cell_state)
cell.colors <- colors[srtLoom_imm$Cell_state[colnames(srtLoom_imm)]]
names(cell.colors) <- colnames(srtLoom_imm)

mye.vlc <- show.velocity.on.embedding.cor(emb = Embeddings(mye.srt, reduction = "sub_umap"),
                                          vel = Tool(srtLoom_imm, slot = "RunVelocity"),
                                          scale = 'sqrt',
                                          n.cores = 16)

PlotPDF('03.velocity.p170_myeloid_sub_umap_projection', 8, 8)
show.velocity.on.embedding.cor(emb = Embeddings(mye.srt, reduction = "sub_umap"),
                               vel = Tool(srtLoom_imm, slot = "RunVelocity"),
                               cc = mye.vlc$cc,
                               cell.colors = ac(cell.colors, alpha = 0.5),
                               cell.border.alpha = 0.04,
                               show.grid.flow = T,
                               cex = 1,
                               n = 1000,
                               grid.n = 30,
                               arrow.scale = 2,
                               arrow.lwd = 1.5,
                               min.grid.cell.mass = 2,
                               do.par = F,
                               asp = 0.75)
dev.off()
####--------------------------------------------------------------------------------------------------------------------
saveRDS(mye.vlc, 'analysis/STEP21.p170_myeloid_velocity_embedding.vlc.rds')
saveRDS(srtLoom_imm, 'analysis/STEP21.p170_myeloid_velocity.srt.rds')
saveRDS(srtLoom, 'analysis/STEP21.p170_allcell_velocity.srt.rds')
####--------------------------------------------------------------------------------------------------------------------


####--------------------------------------------------------------------------------------------------------------------
####  Plot P136 velocity on diffusion map  ####
####--------------------------------------------------------------------------------------------------------------------
srt <- readRDS('integrated/STEP18.annotated.srt.rds')
mye.vlc <- readRDS('analysis/STEP21.p136_myeloid_velocity_embedding.vlc.rds')
srtLoom_imm <- readRDS('analysis/STEP21.p136_myeloid_velocity.srt.rds')
mye.dif <- readRDS('analysis/STEP20.myeloid_diffusion_result.rds')

colors <- mycol_10[3:6]
mye.srt <- srt[, srt$Cell_type == 'Mye' &
                       srt$Cell_state != 'Mast'&
                       srt$Cell_type_non_ambig &
                       srt$name == '5d PT' &
                       srt$Fmux_genotype_clean != 'Ambiguous']

mye.srt@reductions$diffusion <- mye.srt@reductions$sub_umap
mye.srt@reductions$diffusion@cell.embeddings <- mye.dif[Cells(mye.srt), 1:2]*100
colnames(mye.srt@reductions$diffusion@cell.embeddings) <- c('DC_1', 'DC_2')
mye.srt@reductions$diffusion@assay.used <- 'CBN'
mye.srt@reductions$diffusion@key <- 'DC_'

DimPlot2(mye.srt, group.by = 'Cell_state', cols = colors, reduction = 'sub_umap')/
        DimPlot2(mye.srt, group.by = 'Cell_state', cols = colors, reduction = 'diffusion')

names(colors) <- levels(srtLoom_imm$Cell_state)
cell.colors <- colors[srtLoom_imm$Cell_state[colnames(srtLoom_imm)]]
names(cell.colors) <- colnames(srtLoom_imm)
mye.vlc <- show.velocity.on.embedding.cor(emb = Embeddings(mye.srt, reduction = "diffusion"),
                                          vel = Tool(srtLoom_imm, slot = "RunVelocity"),
                                          scale = 'sqrt',
                                          n.cores = 16)

PlotPDF('04.velocity.p136_myeloid_diffusion_projection', 8, 8)
show.velocity.on.embedding.cor(emb = Embeddings(mye.srt, reduction = "diffusion"),
                               vel = Tool(srtLoom_imm, slot = "RunVelocity"),
                               cc = mye.vlc$cc,
                               cell.colors = ac(cell.colors, alpha = 0.5),
                               cell.border.alpha = 0.04,
                               show.grid.flow = T,
                               cex = 1,
                               n = 1000,
                               grid.n = 25,
                               arrow.scale = 2,
                               arrow.lwd = 1.5,
                               min.grid.cell.mass = 6,
                               do.par = F,
                               asp = 0.75)
dev.off()
####--------------------------------------------------------------------------------------------------------------------

