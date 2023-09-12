####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
####  Pediatric Acute Rejection
####  2022-08-03 by Xiao LI (Texas Heart Institute, US)
####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
####  Initiate  ####
####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~----
Ver <- '1'
Step <- 'STEP30_scVelo'
Project <- '2022_acute_rejection_dturaga'
Code_dir <- paste0('~/Documents/Bioinformatics/project/', Project, '/code/human_v', Ver, '/')

source(Sys.readlink(paste0(Code_dir, 'src/bioinformatics.R')))
source(Sys.readlink(paste0(Code_dir, 'src/scRNAseq.R')))
source(Sys.readlink(paste0(Code_dir, 'src/geneset_human.R')))
source(paste0(Code_dir, 'human_v', Ver, '.helper_functions.R'))

InitiateProject('Rivendell', Ver, Step, 'human', Project, 'bree')

library('velocyto.R')
library('anndata')

colors <- mycol_10


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

####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
####  Load Data  ####
####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~----
sample_meta.df <- read.csv(paste0(Docu_dir, 'pediatric_sample_meta.csv'))
rownames(sample_meta.df) <- sample_meta.df$Name_on_disk

full.srt <- readRDS('integrated/STEP18.annotated.srt.rds')
dims <- readRDS('integrated/STEP18.all_reductions.srt_dimreducs.rds')
identical(rownames(dims$scVI), colnames(full.srt))
full.srt@reductions$scVI <- dims$scVI

U(full.srt@meta.data[, c('name', 'recipient')])
Idents(full.srt) <- 'Cell_type'
full.srt <- full.srt[, full.srt$Cell_type_non_ambig]
full.srt$Cell_type_pub <- droplevels(revalue(full.srt$Cell_type, replace = c(
        'CM' = 'Cardiomyocyte',
        'FB' = 'Fibroblast',
        'EC' = 'Endothelial cell',
        'PC' = 'Pericyte',
        'SMC' = 'Smooth muscle cell',
        'Lym' = 'Lympoid cell',
        'Mye' = 'Myeloid cell',
        'Neu' = 'Neuronal cell',
        'Adi' = 'Adipocyte')))
full.srt$Origin <- factor(full.srt$Fmux_genotype_clean, levels = c('Donor', 'Recipient', 'Ambiguous', 'Control cells'))
full.srt$Origin[is.na(full.srt$Origin)] <- 'Control cells'

mye.srt <- full.srt[, full.srt$Cell_type == 'Mye' &
                            full.srt$Origin != 'Ambiguous']
mye.srt$Cell_state <- droplevels(mye.srt$Cell_state)
mye.srt$Origin <- droplevels(mye.srt$Origin)
Idents(mye.srt) <- 'Cell_state'
mye.srt$Cell_state_orig <- as.vector(mye.srt$Cell_state)
mye.srt$Cell_state_orig[mye.srt$Origin == 'Donor' & mye.srt$Cell_state == 'MP2'] <- 'Donor MP2'
mye.srt$Cell_state_orig[mye.srt$Origin == 'Recipient' & mye.srt$Cell_state == 'MP2'] <- 'Recipient MP2'
mye.srt$Cell_state_orig[mye.srt$Origin == 'Control cells' & mye.srt$Cell_state == 'MP2'] <- 'Control MP2'
mye.srt$Cell_state_orig <- factor(mye.srt$Cell_state_orig, levels = c(
        'Control MP2',
        'Donor MP2',
        'Recipient MP2',
        'Mono',
        'MP1',
        'Mast'
))
## reduce white space in umap
x <- mye.srt@reductions$sub_umap@cell.embeddings[, 1]
x[x>7.5] <- x[x>7.5]-3
mye.srt@reductions$sub_umap@cell.embeddings[, 1] <- x
####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
####  Modify Seurat for scVelo  ####
####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~----
full.srt@reductions$mac_umap <- full.srt@reductions$sub_umap
full.srt@reductions$mac_umap@cell.embeddings[,1] <- NA
full.srt@reductions$mac_umap@cell.embeddings[,2] <- NA
full.srt@reductions$mac_umap@cell.embeddings[Cells(mye.srt), ] <- mye.srt@reductions$sub_umap@cell.embeddings
colnames(full.srt@reductions$mac_umap@cell.embeddings) <- c('macUMAP_1', 'macUMAP_2')

full.srt <- RenameAssays(full.srt, 'CBN' = 'RNA')

## rename cells to match loom files
library_order <- c(
        '2019_Unpub7_JMartin_ctrl_s4',
        '2019_Unpub7_JMartin_ctrl_s5',
        '2019_Unpub7_JMartin_ctrl_s11',
        '2019_Unpub7_JMartin_ctrl_s12',
        '2019_Unpub7_JMartin_ctrl_s6',
        '2019_Unpub7_JMartin_ctrl_s7',
        '2021_Circulation_EPorrello_Y2',
        '2021_Circulation_EPorrello_Y3',
        '2022_Unpub2_DTuraga_P136_s1',
        '2022_Unpub2_DTuraga_P170_s1',
        '2022_Unpub3_DTuraga_P182'
        )
rownames(sample_meta.df) <- sample_meta.df$Name_on_disk
sample_id_seq <- sample_meta.df[library_order, 'Sample_id']
barcode <- split(full.srt$orig.name,  as.vector(full.srt$sample))
barcode <- barcode[sample_id_seq]
id <- split(Cells(full.srt),  as.vector(full.srt$sample))
id <- id[sample_id_seq]
for(i in 1:L(barcode)){
        barcode[[i]] <- str_split(barcode[[i]], pattern = '-', simplify = T)[,1]
        barcode[[i]] <- paste0(barcode[[i]], '_', i)
        names(barcode[[i]]) <- as.vector(id[[i]])
}
new_name <- unlist(barcode)
new_name <- new_name[paste0(full.srt$sample, '.', Cells(full.srt))]
head(new_name)

full.srt <- RenameCells(full.srt, new.names = new_name)
####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
####  Export Seurat data for scVelo  ####
####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~----
## Save metadata table:
full.srt$barcode <- Cells(full.srt)
full.srt$UMAP_1 <- full.srt@reductions$clean_umap@cell.embeddings[,1]
full.srt$UMAP_2 <- full.srt@reductions$clean_umap@cell.embeddings[,2]
full.srt$SUBUMAP_1 <- full.srt@reductions$sub_umap@cell.embeddings[,1]
full.srt$SUBUMAP_2 <- full.srt@reductions$sub_umap@cell.embeddings[,2]
full.srt$MACUMAP_1 <- full.srt@reductions$mac_umap@cell.embeddings[,1]
full.srt$MACUMAP_2 <- full.srt@reductions$mac_umap@cell.embeddings[,2]
full.srt$group1 <- full.srt$name
full.srt$group2 <- full.srt$name2
meta <- full.srt@meta.data[, c('sample', 'group1', 'group2', 'Cell_type', 'Cell_state', 'Origin',
                               'barcode', 'UMAP_1', 'UMAP_2', 'SUBUMAP_1', 'SUBUMAP_2', 'MACUMAP_1', 'MACUMAP_2')]
WriteCSV(meta, title = 'full.srt_meta')

# write dimesnionality reduction matrix
scVI <- full.srt@reductions$scVI@cell.embeddings
WriteCSV(as.data.frame(scVI), title = 'full_scVI.srt_dim')

# write gene names
WriteCSV(as.data.frame(rownames(full.srt)), title = 'full.gene_names', col.names = F)

## Write expression counts matrix ## HUGE FILE!
counts_matrix <- GetAssayData(full.srt, assay = 'RNA', slot = 'counts')
writeMM(counts_matrix, file = paste0(Meta_dir, 'full.counts_mtx'))

## Mannually gzip all 4 txt files and move to data drive rdata/analysis
####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
####  Follow Python script: human_v1.STEP30.scvelo.ipynb  ####
####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~----


####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
####  Extract latent time  ####
####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~----
library('reticulate')
library('anndata')
use_condaenv("cellrank", required = TRUE)
scv <- import("scvelo")
mac_adata <- read_h5ad('analysis/PART24.main_mac_scvelo.h5ad')
mac_adata
mac_ltime <- mac_adata$obs$latent_time
names(mac_ltime) <- mac_adata$obs_names

fb_adata <- read_h5ad('analysis/PART24.fb_scvelo.h5ad')
fb_adata
fb_ltime <- fb_adata$obs$latent_time
names(fb_ltime) <- fb_adata$obs_names

ec_adata <- read_h5ad('analysis/PART24.ec_scvelo.h5ad')
ec_adata
ec_ltime <- ec_adata$obs$latent_time
names(ec_ltime) <- ec_adata$obs_names
####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
saveRDS(mac_ltime, 'analysis/PART24.main_mac_scvelo.latent_time_vector.rds')
saveRDS(fb_ltime, 'analysis/PART24.fb_scvelo.latent_time_vector.rds')
saveRDS(ec_ltime, 'analysis/PART24.ec_scvelo.latent_time_vector.rds')
####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
####  Fit time_correlated genes  ####
####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~----
fit_genes <- data.frame(LH=mac_adata$var['fit_likelihood'],
                        Gene=rownames(mac_adata$var['fit_likelihood']))
fit_genes <- fit_genes[order(fit_genes$fit_likelihood, decreasing = T),]
fit_genes$Top <- F
fit_genes$Top[fit_genes$fit_likelihood > 0.1] <- T
Table(fit_genes$Top)

fit_genes <- data.frame(LH=ec_adata$var['fit_likelihood'],
                        Gene=rownames(ec_adata$var['fit_likelihood']))
fit_genes <- fit_genes[order(fit_genes$fit_likelihood, decreasing = T),]
fit_genes$Top <- F
fit_genes$Top[fit_genes$fit_likelihood > 0.1] <- T
Table(fit_genes$Top)
####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
saveRDS(fit_genes, 'analysis/PART24.main_mac_scvelo.latent_time_fit_gene_df.rds')
####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
