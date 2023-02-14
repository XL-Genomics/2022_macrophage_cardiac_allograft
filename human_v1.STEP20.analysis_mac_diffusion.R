####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
####  Pediatric Acute Rejection
####  2022-08-03 by Xiao LI (Texas Heart Institute, US)
####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
####  Initiate  ####
####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~----
Ver <- '1'
Step <- 'STEP20_Analysis_Mac_Diffusion'
Project <- '2022_acute_rejection_dturaga'
Code_dir <- paste0('~/Documents/Bioinformatics/project/', Project, '/code/human_v', Ver, '/')

source(Sys.readlink(paste0(Code_dir, 'src/bioinformatics.R')))
source(Sys.readlink(paste0(Code_dir, 'src/scRNAseq.R')))
source(Sys.readlink(paste0(Code_dir, 'src/geneset_human.R')))
source(paste0(Code_dir, 'human_v', Ver, '.helper_functions.R'))

InitiateProject('Rivendell', Ver, Step, 'human', Project, 'ithil')

library('monocle3')
library('destiny')
suppressMessages(library('reticulate'))
suppressMessages(library('anndata'))
sc <- import("scanpy")
np <- import('numpy')
sce <- import('scanpy.external')
####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
####  Load data  ####
####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~----
full.srt <- readRDS('integrated/STEP18.annotated.srt.rds')
mye.srt <- full.srt[, full.srt$Cell_type == 'Mye' & full.srt$Cell_state != 'Mast' & full.srt$Cell_type_non_ambig]
mye.srt$Cell_type <- droplevels(mye.srt$Cell_type)
mye.srt$Cell_state <- droplevels(mye.srt$Cell_state)

DimPlot2(mye.srt, reduction = 'sub_umap', group.by = 'Cell_state')
DimPlot2(mye.srt, reduction = 'sub_umap', group.by = 'name2')

# saveRDS(mye.srt, 'tmp/STEP20.mye.srt.rds')
mye.srt <- readRDS('tmp/STEP20.mye.srt.rds')
mye.srt.bkp <- mye.srt
mye.srt <- mye.srt.bkp

mye.srt_orig <- readRDS('integrated/STEP10.myeloid_cells.srt.rds')
####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
####  Myeloid diffusion  ####
####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~----
mye.srt <- ProcessSrt_sct(mye.srt, do.umap = F, assay = 'CBN')
## Compute diffusion map on genes
dif <- scater::calculateDiffusionMap(mye.srt@assays$SCT@data)

mye.srt@reductions$diffusion <- mye.srt@reductions$sub_umap
mye.srt@reductions$diffusion@cell.embeddings <- dif[,1:2]*100
mye.srt@reductions$diffusion@assay.used <- 'SCT'
mye.srt@reductions$diffusion@key <- 'DC'

DimPlot2(mye.srt, reduction = 'diffusion', group.by = 'Cell_state')/
        DimPlot2(mye.srt, reduction = 'diffusion', group.by = 'Fmux_genotype_clean')

cells <- DownsampleByMeta(mye.srt, meta_var = 'name2', down_to_min_group = T, random = T, n = NULL)
DimPlot2(mye.srt[, unlist(cells)],
         reduction = 'diffusion', group.by = 'Cell_state', split.by = 'name2', ncol = 2) /
        DimPlot2(mye.srt[, unlist(cells)],
                 reduction = 'diffusion', group.by = 'Fmux_genotype_clean', split.by = 'name2', ncol = 2)

## Compute diffusion pseudotime
mye.sce <- as.SingleCellExperiment(mye.srt)
mye.dm <- DiffusionMap(mye.sce, n_eigs = 2)
mye.dpt <- DPT(mye.dm, tips = 3)

p <- plot(mye.dpt, root = 3)
mye.srt$DPT <- p$data$Colour

DimPlot2(mye.srt, reduction = 'diffusion', group.by = 'Cell_state')/
        FeaturePlot2(mye.srt, reduction = 'diffusion', features = 'DPT')

cells <- DownsampleByMeta(mye.srt, meta_var = 'name2', down_to_min_group = T, random = T, n = NULL)
DimPlot2(mye.srt[, unlist(cells)],
         reduction = 'diffusion', group.by = 'Cell_state', split.by = 'name2', ncol = 4) /
        FeaturePlot2(mye.srt[, unlist(cells)], reduction = 'diffusion', features = 'DPT', split.by = 'name2')
####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
saveRDS(dif, 'analysis/STEP20.myeloid_diffusion_result.rds')
saveRDS(mye.dpt, 'analysis/STEP20.myeloid_dpt_result.rds')
####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
####  Myeloid diffusion, control+P136  ####
####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~----
mye.srt <- mye.srt[, mye.srt$name2 %in% c('Control', '5d PT')]
mye.srt$name2 <- droplevels(mye.srt$name2)
mye.srt <- ProcessSrt_sct(mye.srt, do.umap = F, assay = 'CBN')

## Compute diffusion map on genes
dif <- scater::calculateDiffusionMap(mye.srt@assays$SCT@data)
mye.srt@reductions$diffusion <- mye.srt@reductions$sub_umap
mye.srt@reductions$diffusion@cell.embeddings <- dif[,1:2]*100
mye.srt@reductions$diffusion@assay.used <- 'SCT'
mye.srt@reductions$diffusion@key <- 'DC'

dif <- mye.srt@reductions$diffusion@cell.embeddings
mye.srt$Diffusion <- ((dif[,1] + abs(min(dif[,1])))^2 + (dif[,2] + abs(min(dif[,2])))^2)^0.5
mye.srt$Diffusion_rank <- rank(mye.srt$Diffusion)

DimPlot2(mye.srt, reduction = 'diffusion', group.by = 'Cell_state')/
        FeaturePlot2(mye.srt, reduction = 'diffusion', features = 'Diffusion')


cells <- DownsampleByMeta(mye.srt, meta_var = 'name2', down_to_min_group = T, random = T, n = NULL)
DimPlot2(mye.srt[, unlist(cells)],
         reduction = 'diffusion', group.by = 'Cell_state', split.by = 'name2', ncol = 2) /
        DimPlot2(mye.srt[, unlist(cells)],
                 reduction = 'diffusion', group.by = 'Fmux_genotype_clean', split.by = 'name2', ncol = 2)
####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
saveRDS(dif, 'analysis/STEP20.myeloid_ctrl_p136_diffusion_result.rds')
####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
####  Myeloid diffusion, P136  ####
####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~----
p136.sce <- as.SingleCellExperiment(mye.srt[, mye.srt$name == '5d PT' & mye.srt$Fmux_genotype_clean != 'Ambiguous'])
p136.dm <- DiffusionMap(p136.sce, n_eigs = 2)
p136.dpt <- DPT(p136.dm, tips = 1)

p <- plot(p136.dpt, root = 3)
mye.srt$DPT <- NA
mye.srt$DPT[mye.srt$name == '5d PT' & mye.srt$Fmux_genotype_clean != 'Ambiguous'] <- rank(p$data$Colour)

DimPlot2(mye.srt[, mye.srt$name == '5d PT' & mye.srt$Fmux_genotype_clean != 'Ambiguous'],
         reduction = 'sub_umap', group.by = 'Cell_state')/
        FeaturePlot2(mye.srt[, mye.srt$name == '5d PT' & mye.srt$Fmux_genotype_clean != 'Ambiguous'],
                     reduction = 'sub_umap', features = 'DPT')
pgf_mo_mp2_dpt <- mye.srt$DPT[mye.srt$name == '5d PT' &
                                      mye.srt$Fmux_genotype_clean != 'Ambiguous' &
                                      mye.srt$Cell_state %in% c('Mono', 'MP2')]
p <- plot(p136.dpt, root = 1)
mye.srt$DPT <- NA
mye.srt$DPT[mye.srt$name == '5d PT' & mye.srt$Fmux_genotype_clean != 'Ambiguous'] <- rank(p$data$Colour)

DimPlot2(mye.srt[, mye.srt$name == '5d PT' & mye.srt$Fmux_genotype_clean != 'Ambiguous'],
         reduction = 'sub_umap', group.by = 'Cell_state')/
        FeaturePlot2(mye.srt[, mye.srt$name == '5d PT' & mye.srt$Fmux_genotype_clean != 'Ambiguous'],
                     reduction = 'sub_umap', features = 'DPT')
pgf_mo_mp1_dpt <- mye.srt$DPT[mye.srt$name == '5d PT' &
                                      mye.srt$Fmux_genotype_clean != 'Ambiguous' &
                                      mye.srt$Cell_state %in% c('Mono', 'MP1')]
sub_mye.srt <- mye.srt[, mye.srt$name == '5d PT' &
                               mye.srt$Fmux_genotype_clean != 'Ambiguous' &
                               mye.srt$Cell_state %in% c('Mono', 'MP1')]
ave_exp <- AverageExpression(sub_mye.srt, group.by = 'Cell_state', assays = 'CBN')$CBN
genes_keep <- rownames(ave_exp)[rowMins(ave_exp) >= 0.1]
L(genes_keep) ## 7710 highly expressed genes
ave_exp <- as.data.frame(ave_exp[genes_keep, ])
ave_exp$PST_PCC <- NA
ave_exp$PST_SCC <- NA
expr_mat <- sub_mye.srt@assays$CBN@data[genes_keep, ]
for(i in 1:nrow(ave_exp)){
        if(i %% 100 == 0){message(i)}
        ave_exp$PST_SCC[i] <- cor(x = expr_mat[rownames(ave_exp)[i], ], y = sub_mye.srt$DPT,
                                  method = 'spearman')
        ave_exp$PST_PCC[i] <- cor(x = expr_mat[rownames(ave_exp)[i], ], y = sub_mye.srt$DPT,
                                  method = 'pearson')
}
ave_exp$PST_PCC_zscore <- as.vector(scale(ave_exp$PST_PCC))
ave_exp$PST_SCC_zscore <- as.vector(scale(ave_exp$PST_SCC))

pos_gene <- rownames(ave_exp)[ave_exp$PST_PCC_zscore >= 2]
pos_gene <- pos_gene[order(ave_exp[pos_gene, 'PST_PCC_zscore'], decreasing = T)]
L(pos_gene)
neg_gene <- rownames(ave_exp)[ave_exp$PST_PCC_zscore <= -2]
neg_gene <- neg_gene[order(ave_exp[neg_gene, 'PST_PCC_zscore'], decreasing = F)]
L(neg_gene)

sub_mye.srt <- RunALRA(sub_mye.srt, genes.use = rownames(sub_mye.srt), assay = 'CBN')
sub_mye.srt <- ScaleData(sub_mye.srt, features = c(pos_gene, neg_gene), assay = 'alra')
expr_mat <- sub_mye.srt@assays$alra@scale.data[c(pos_gene, neg_gene),
                                               colnames(sub_mye.srt)[order(sub_mye.srt$DPT)]]
FlattenExpr <- function(expr_scale, genes){
        expr_flat <- expr_scale[genes,]
        expr_flat <- as.data.frame(matrix(expr_flat, ncol = 1))
        expr_flat <- cbind(expr_flat,
                           rep(1:ncol(expr_scale), each = length(genes)),
                           rep(genes, times = ncol(expr_scale)))
        colnames(expr_flat) <- c("Expression", "Cells", "Genes")
        expr_flat <- expr_flat[order(expr_flat$Genes), ]
        expr_flat <- cbind(expr_flat,
                           colMeans(expr_scale[genes, ]),
                           colMedians(expr_scale[genes, ]),
                           colQuantiles(expr_scale[genes, ], probs = 0.75),
                           colQuantiles(expr_scale[genes, ], probs = 0.25))
        colnames(expr_flat)[4:7] <- c("Mean Expr", "Median Expr",
                                      "Upper Quartile", "Lower Quartile")
        return(expr_flat)
}
module <- list(Pos = pos_gene, Neg = neg_gene)
####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
saveRDS(module, 'analysis/STEP20.p136_mono_mp1_dpt_coupled_genes.rds')
####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
expr_mat <- sub_mye.srt@assays$alra@scale.data[c(pos_gene, neg_gene),
                                               colnames(sub_mye.srt)[order(sub_mye.srt$DPT)]]
for(i in 1:nrow(expr_mat)){expr_mat[i,] <- smooth.spline(expr_mat[i, ], spar = 0.9)$y}
df_all <- mapply(function(x, y){
        expr_mat <- t(apply(expr_mat[x,], 1, Range01)) ## normalize to [0-1]
        df <- FlattenExpr(expr_mat[, seq(1, ncol(sub_mye.srt), 10)], x)
        df$DPT <- Range01(df$Cells)
        df$mod_id <- names(module)[y]
        return(df)},
        x = module, y = seq_along(module), SIMPLIFY = F)
df_all <- bind_rows(df_all)
name <- paste0("Module ", 1:2, ' (', lapply(module, function(x) length(x)), ' genes)')
p <- ggplot(data = df_all) +
        geom_ribbon(aes(x = DPT, ymax = `Upper Quartile`, ymin = `Lower Quartile`, fill = mod_id),
                    alpha = 0.1, show.legend = F) +
        geom_line(aes(x = DPT, y = `Median Expr`, color = mod_id),
                  show.legend = T, size = 1, alpha = 1) +
        labs(y = 'Normalized expression') +
        scale_color_manual(values = mycol_10, labels = name) +
        scale_fill_manual(values = mycol_10, labels = name)  +
        scale_y_continuous(breaks = seq(0, 1, 0.5), limits = c(0, 1)) +
        scale_x_continuous(limits = c(0.15, 1)) +
        theme_classic() +
        theme(aspect.ratio = 1, axis.title.x = element_blank(), axis.text.x = element_blank(),
              text = element_text(color = "black"),
              line = element_line(color = "black"),
              axis.line = element_line(color = "black"),
              legend.key.size = unit(10, 'pt'),
              legend.text = element_text(size = 10),
              legend.box.margin = margin(rep(1,1)), legend.title = element_blank(),
              legend.direction = 'vertical', legend.position = 'right',
              #panel.grid.major.x = element_line(color = 'black', linetype = 'dotted')
        )
p
####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
saveRDS(p136.dpt, 'analysis/STEP20.myeloid_p136_dpt_result.rds')
saveRDS(pgf_mo_mp1_dpt, 'analysis/STEP20.p136_mono_mp1_dpt_rank.rds')
saveRDS(pgf_mo_mp2_dpt, 'analysis/STEP20.p136_mono_mp2_dpt_rank.rds')
####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
####  Myeloid monocle pseudotime  ####
####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~----
library('scran')

DefaultAssay(mye.srt) <- 'CBN'
mye.srt2 <- DietSeurat(mye.srt, assays = 'CBN')
mye.srt2 <- RenameAssays(mye.srt2, 'CBN' = 'RNA')
mye.srt_list <- SplitObject(mye.srt2, split.by = 'name')
mye.mnn <- RunFastMNN(mye.srt_list, reduction.name = 'fMNN', reduction.key = 'fmnn_', reconstructed.assay = 'MNN')
mye.mnn <- ScaleData(mye.mnn, assay = 'MNN', vars.to.regress = c('S.Score', 'G2M.Score')) |> RunPCA(assay = 'MNN')
mye.mnn <- RunHarmony(mye.mnn, group.by.vars = 'sample')
mye.mnn <- RunUMAP(mye.mnn, reduction = 'pca', dims = 1:10)
DimPlot2(mye.mnn, group.by = 'Cell_state')

mye.cds <- new_cell_data_set(expression_data = as(as.matrix(GetAssayData(mye.mnn,
                                                                         assay = "MNN",
                                                                         slot = "scale.data")), 'sparseMatrix'),
                             cell_metadata = data.frame(mye.mnn@meta.data),
                             gene_metadata = data.frame(gene_short_name = rownames(mye.mnn@assays$MNN),
                                                        row.names = row.names(mye.mnn@assays$MNN)))
mye.cds@metadata$variable_genes <- rownames(mye.mnn@assays$MNN)  ## Adding Seurat feature genes info

### Define the root cells
root_cell <- "2022_Unpub3_DTuraga:P4_S001:TTTGACTAGACCATGG-1"
mye.cds <- preprocess_cds(mye.cds, num_dim = 8, method = "PCA", norm_method = "size_only",
                          pseudo_count = 0, use_genes = mye.cds@metadata$variable_genes)
mye.cds <- reduce_dimension(mye.cds, preprocess_method = "PCA", reduction_method = "UMAP",
                            umap.min_dist = 0.3, umap.n_neighbors = 30)

p1 <- plot_cells(mye.cds, color_cells_by = "Cell_state", cell_size = 0.5,
                 show_trajectory_graph = F, group_label_size = 5)
p2 <- plot_cells(mye.cds, color_cells_by = "name2", cell_size = 0.5,
                 show_trajectory_graph = F, group_label_size = 5)
p1/p2
mye.cds <- cluster_cells(mye.cds, resolution = 0.001, random_seed = 5128)

### Map pseudotime
mye.cds <- learn_graph(mye.cds, use_partition = F, close_loop = T,
                       learn_graph_control = list(minimal_branch_len = 20))
mye.cds <- order_cells(mye.cds, reduction_method = "UMAP")

p3 <- plot_cells(mye.cds, color_cells_by = "pseudotime", trajectory_graph_segment_size = 2)
p4 <- plot_cells(mye.cds, color_cells_by = "Cell_state", trajectory_graph_segment_size = 2,
                 label_cell_groups = T, group_label_size = 5)
p3/p4
####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
####  P136 Myeloid diffusion map + Slingshot trajectory  ####
####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~----
mye.srt <- mye.srt.bkp
mye.srt <- mye.srt[, mye.srt$name2 %in% c('Control', '5d PT') & mye.srt$Fmux_genotype_clean != 'Ambiguous']
mye.srt$name2 <- droplevels(mye.srt$name2)
dif <- readRDS('analysis/STEP20.myeloid_ctrl_p136_diffusion_result.rds')
mye.srt@reductions$diffusion <- mye.srt@reductions$sub_umap
mye.srt@reductions$diffusion@cell.embeddings <- dif[,1:2]*100
mye.srt@reductions$diffusion@assay.used <- 'SCT'
mye.srt@reductions$diffusion@key <- 'DC'

p136.srt <- mye.srt[, mye.srt$name2 == '5d PT']
DimPlot2(p136.srt, reduction = 'diffusion', group.by = 'Cell_state') /
        DimPlot2(p136.srt, reduction = 'diffusion', group.by = 'Fmux_genotype_clean')

#### Trajectory inference
library('slingshot')
# x <- mye.srt@reductions$diffusion@cell.embeddings
# fit0 <- principal_curve(x, smoother = 'lowess', maxit = 2)
# plot(x)
# points(fit0$s)
p136.sce <- as.SingleCellExperiment(p136.srt)
p136.sce <- slingshot(p136.sce,
                      clusterLabels = 'Cell_state',
                      reducedDim = 'DIFFUSION', thresh = 0.0001,
                      start.clus = 'Mono'
                      )
colors <- colorRampPalette(brewer.pal(11,'Spectral')[-6])(100)
plotcol <- colors[cut(p136.sce$slingPseudotime_1, breaks=100)]
PlotPDF('1.1.p136_mono_mp1_slingshot_trajectory', 6, 6)
plot(reducedDims(p136.sce)$DIFFUSION, col = plotcol, pch=16, asp = 1)
lines(SlingshotDataSet(p136.sce), lwd=2, col='black')
dev.off()

plotcol <- colors[cut(p136.sce$slingPseudotime_2, breaks=100)]
PlotPDF('1.2.p136_mono_mp2_slingshot_trajectory', 6, 6)
plot(reducedDims(p136.sce)$DIFFUSION, col = plotcol, pch=16, asp = 1)
lines(SlingshotDataSet(p136.sce), lwd=2, col='black')
dev.off()

p136.srt$PST1 <- p136.sce$slingPseudotime_2
p136.srt$PST2 <- p136.sce$slingPseudotime_1

## Mono-MP1
p136_momp1.srt <- p136.srt[,!is.na(p136.srt$PST1)]
ave_exp <- AverageExpression(p136_momp1.srt, group.by = 'Cell_state', assays = 'CBN')$CBN
genes_keep <- rownames(ave_exp)[rowMins(ave_exp) >= 0.1]
L(genes_keep) ## 7350 highly expressed genes
ave_exp <- as.data.frame(ave_exp[genes_keep, ])
ave_exp$PST_PCC <- NA
expr_mat <- p136_momp1.srt@assays$CBN@data[genes_keep, Cells(p136_momp1.srt)]
for(i in 1:nrow(ave_exp)){
        if(i %% 100 == 0){message(i)}
        ave_exp$PST_PCC[i] <- cor(x = expr_mat[rownames(ave_exp)[i], ],
                                  y = p136_momp1.srt$PST1,
                                  method = 'pearson')
}
ave_exp$PST_PCC_zscore <- as.vector(scale(ave_exp$PST_PCC))

pos_gene <- rownames(ave_exp)[ave_exp$PST_PCC_zscore >= 2]
pos_gene <- pos_gene[order(ave_exp[pos_gene, 'PST_PCC_zscore'], decreasing = T)]
L(pos_gene) ## 312
c('NFKB1', 'RELB', 'STAT1', 'IRF1') %in% pos_gene

neg_gene <- rownames(ave_exp)[ave_exp$PST_PCC_zscore <= -2]
neg_gene <- neg_gene[order(ave_exp[neg_gene, 'PST_PCC_zscore'], decreasing = F)]
L(neg_gene) ## 52

gl <- list(MoMP1 = list(Pos = pos_gene, Neg = neg_gene))

## Mono-MP2
p136_momp2.srt <- p136.srt[,!is.na(p136.srt$PST2)]
ave_exp <- AverageExpression(p136_momp2.srt, group.by = 'Cell_state', assays = 'CBN')$CBN
genes_keep <- rownames(ave_exp)[rowMins(ave_exp) >= 0.1]
L(genes_keep) ## 6509 highly expressed genes
ave_exp <- as.data.frame(ave_exp[genes_keep, ])
ave_exp$PST_PCC <- NA
expr_mat <- p136_momp2.srt@assays$CBN@data[genes_keep, Cells(p136_momp2.srt)]
for(i in 1:nrow(ave_exp)){
        if(i %% 100 == 0){message(i)}
        ave_exp$PST_PCC[i] <- cor(x = expr_mat[rownames(ave_exp)[i], ],
                                  y = p136_momp2.srt$PST2,
                                  method = 'pearson')
}
ave_exp$PST_PCC_zscore <- as.vector(scale(ave_exp$PST_PCC))

pos_gene <- rownames(ave_exp)[ave_exp$PST_PCC_zscore >= 2]
pos_gene <- pos_gene[order(ave_exp[pos_gene, 'PST_PCC_zscore'], decreasing = T)]
L(pos_gene) ## 239

neg_gene <- rownames(ave_exp)[ave_exp$PST_PCC_zscore <= -2]
neg_gene <- neg_gene[order(ave_exp[neg_gene, 'PST_PCC_zscore'], decreasing = F)]
L(neg_gene) ## 82

gl <- list(MoMP1 = list(Pos = pos_gene, Neg = neg_gene))
gl$MoMP2 <- list(Pos = pos_gene, Neg = neg_gene)
####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
saveRDS(p136.sce, 'analysis/STEP20.p136_slingshot.sce.rds')
saveRDS(gl, 'analysis/STEP20.p136_mono_mp1_mp2_sling_pst_coupled_genes.list.rds')
####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
####  Myeloid diffusion using SCANPY  ####
####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~----
mye.srt_orig <- mye.srt_orig[, Cells(mye.srt)]
DefaultAssay(mye.srt_orig) <- 'CBN'
for(i in 1:L(mye.srt_orig@reductions)){mye.srt_orig@reductions[[i]]@assay.used <- 'CBN'}

identical(Cells(mye.srt_orig), Cells(mye.srt))
mye.srt_orig@meta.data <- mye.srt@meta.data
SaveH5ad(mye.srt_orig, path = 'analysis/', name = 'STEP20.mye.ann',
         assay = 'CBN', raw_count_only = T, verbose = T)
## Following code is for avoiding "_index" in adata.var bug
adata <- read_h5ad('analysis/STEP20.mye.ann.h5ad')
adata$raw <- NULL
adata$write_h5ad(filename = 'analysis/STEP20.mye.ann.h5ad') ## replace the original
adata$X[1:10, 1:20]
####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

