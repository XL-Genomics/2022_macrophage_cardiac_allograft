####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
####  Pediatric Acute Rejection
####  2022-08-03 by Xiao LI (Texas Heart Institute, US)
####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~----


####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
####  Initiate  ####
####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
Ver <- '1'
Step <- 'STEP91_Plot_Circulation_Revision'
Project <- '2022_acute_rejection_dturaga'
Code_dir <- paste0('~/Documents/Bioinformatics/project/', Project, '/code/human_v', Ver, '/')

source(Sys.readlink(paste0(Code_dir, 'src/bioinformatics.R')))
source(Sys.readlink(paste0(Code_dir, 'src/scRNAseq.R')))
source(Sys.readlink(paste0(Code_dir, 'src/geneset_human.R')))
source(paste0(Code_dir, 'human_v', Ver, '.helper_functions.R'))

InitiateProject('Rivendell', Ver, Step, 'human', Project, 'bree')

suppressMessages(library('ggtern'))
suppressMessages(library('ggdendro'))
suppressMessages(library('ggupset'))
suppressMessages(library('velocyto.R'))
suppressMessages(library('CellChat'))
suppressMessages(library('slingshot'))

Color_cell_type <- c(mycol_40[c(3, 6, 11, 14, 35, 17, 22, 26, 30)], 'grey60', 'grey80')
Color_cell_state <- c(mycol_40[1:4],
                      mycol_40[5:8],
                      mycol_30[16:17], mycol_30[13:14],  mycol_30[7:8],
                      mycol_30[10:11],
                      mycol_30[19:21],
                      mycol_40[29:32],
                      mycol_10[9],
                      'grey85')
Color_cell_state_all <- c(Color_cell_state, paste0('grey', seq(20, 80, 10)))
Color_cell_state_mye <- mycol_10[c(5, 3, 4, 7)]
Color_cell_state_orig_mye <- c('deeppink2', 'cyan3', mycol_10[5])
Color_controls <- mycol_60[c(6, 4, 2, 1)]
Color_transplant <- mycol_14[c(1, 10, 9)]
Color_donor <- c(Color_controls, Color_transplant)
Color_condition <- c('grey85', Color_transplant)
Color_genotype <- c('grey85', 'cyan3', 'deeppink2' ,'lightskyblue4')
Color_genotype2 <- c('deeppink2', 'cyan3', 'black', 'grey85')
Color_genotype3 <- c('deeppink2', 'cyan3', 'lightskyblue4')
Color_genotype4 <- c("#427FC2", "#F7941D", 'black', 'grey85')
Color_mp2_group <- c(Color_genotype2[1:2], Color_genotype4[1:2])
####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~----


####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
####  Edit previously saved objects  ####
####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#### All cell types, all samples -- full.srt  ####
full.srt <- readRDS('integrated/STEP19.annotated_alra.srt.rds')
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

#### All Myeloid  -- mye.srt  ####
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

#### P136 Myeloid Velocyto  -- velo_all.srt  ####
velo_pgf.vlc <- readRDS('analysis/STEP21.p136_myeloid_velocity_embedding.vlc.rds')
velo_pgf.srt <- readRDS('analysis/STEP21.p136_myeloid_velocity.srt.rds')
velo_all.srt <- mye.srt
velo_all.srt@reductions$sub_umap <- full.srt[, Cells(velo_all.srt)]@reductions$sub_umap
velo_all.srt$Cell_state <- full.srt$Cell_state[Cells(velo_all.srt)]
velo_all.srt$Cell_state <- droplevels(velo_all.srt$Cell_state)
Idents(velo_all.srt) <- 'Cell_state'
velo_all.srt$Cell_state_orig <- 'Control'
velo_all.srt$Cell_state_orig[velo_all.srt$name2 != 'Control'] <-
        as.vector(mye.srt$Cell_state_orig[Cells(velo_all.srt)[velo_all.srt$name2 != 'Control']])
velo_all.srt$Cell_state_orig[is.na(velo_all.srt$Cell_state_orig)] <- 'Control'
velo_all.srt$Cell_state_orig <- factor(velo_all.srt$Cell_state_orig,
                                       levels = c('Control', 'Donor MP2', 'Recipient MP2', 'Mono', 'MP1'))

#### P136 Myeloid (w/o Mast cell) Diffusion --  momp_p136.srt ####
momp_p136.srt <- mye.srt[, mye.srt$name == '5d PT' & mye.srt$Cell_state != 'Mast']
velo_pgf_momp.srt <- velo_pgf.srt[, Cells(momp_p136.srt)]
dif <- readRDS('analysis/STEP20.myeloid_ctrl_p136_diffusion_result.rds')
momp_p136.srt@reductions$diffusion <- momp_p136.srt@reductions$sub_umap
momp_p136.srt@reductions$diffusion@cell.embeddings[Cells(momp_p136.srt), ] <- dif[Cells(momp_p136.srt), 1:2]
momp_p136.srt@reductions$diffusion@assay.used <- 'CBN'
momp_p136.srt@reductions$diffusion@key <- 'DC_'
colnames(momp_p136.srt@reductions$diffusion@cell.embeddings) <- c('DC_1', 'DC_2')
Idents(momp_p136.srt) <- 'Cell_state_orig'

#### P136 Myeloid (w/o Mast cell) SlingShot --  momp_p136.sce ####
momp_p136.sce <- readRDS('analysis/STEP20.p136_slingshot.sce.rds')
momp_p136.sce$Cell_state_orig <- droplevels(momp_p136.srt$Cell_state_orig)
## add slingshot pseudotime to momp_p136.srt
momp_p136.srt$PST_slingshot <- momp_p136.sce$slingPseudotime_2
momp_p136.srt$PST_slingshot[is.na(momp_p136.srt$PST_slingshot)] <-
        momp_p136.sce$slingPseudotime_1[is.na(momp_p136.srt$PST_slingshot)]
momp_p136.srt$PST_slingshot_rank <- rank(momp_p136.srt$PST_slingshot)

#### Allograft MP2  -- mp2.srt  ####
mp2.srt <- mye.srt[, mye.srt$Cell_state == 'MP2' & mye.srt$name2 != 'Control' &
                           ! paste(mye.srt$Cell_state_orig, mye.srt$name2) %in% c('Donor MP2 15m PT', 'Donor MP2 12y PT')]
Table(mp2.srt$Cell_state_orig, mp2.srt$name2)
mp2.srt$Group <- factor(paste(revalue(mp2.srt$Cell_state_orig,
                                      replace = c('Donor MP2' = 'ResMP', 'Recipient MP2' = 'MoMP2')),
                              mp2.srt$name2),
                        levels = c('ResMP 5d PT', 'MoMP2 5d PT', 'MoMP2 15m PT', 'MoMP2 12y PT'))
Table(mp2.srt$Group)
Idents(mp2.srt) <- 'Group'

#### Allograft Late MoMP2 + ResMP  -- mp2_sub.srt  ####
mp2_sub.srt <- mp2.srt[, mp2.srt$Group != 'MoMP2 5d PT']
mp2_sub.srt$Group2 <- revalue(mp2_sub.srt$Group, replace = c(
        'ResMP 5d PT' = 'ResMP',
        'MoMP2 15m PT' = 'Late MoMP2',
        'MoMP2 12y PT' = 'Late MoMP2'
))
mp2_sub.srt$Group2 <- factor(mp2_sub.srt$Group2, levels = c('ResMP', 'Late MoMP2'))
Idents(mp2_sub.srt) <- 'Group2'
mp2_sub.srt <- FindVariableFeatures(mp2_sub.srt, nfeatures = 5000)

#### MP2 CellChat  -- mp2_cc_list  ####
mp2_cc_list <- list(
        '5d PT'  = readRDS('analysis/STEP24.p136_alra_cellchat.cch.rds'),
        '15m PT' = readRDS('analysis/STEP24.p182_alra_cellchat.cch.rds'),
        '12y PT' = readRDS('analysis/STEP24.p170_alra_cellchat.cch.rds'),
        'LR' = readRDS('analysis/STEP24.final_cellchat_lr_pairs.list.rds')
)

#### MP with DCM  -- mp_dcm.srt  ####
mp_dcm.srt <- readRDS('analysis/STEP25.allo_dcm_ctrl_mp_compare.srt.rds')

#### NK with DCM  -- nk_dcm.srt  ####
nk_dcm.srt <- readRDS('analysis/STEP25.allo_dcm_ctrl_nk_compare.srt.rds')

#### P136 CellChat  -- p136_cc_list  ####
p136_cc_list <- list(
        'All'  = readRDS('analysis/STEP26.p136_alra_group_by_all_cell_states.cellchat.rds'),
        'Simple' = readRDS('analysis/STEP26.p136_alra_group_by_all_cell_states_simplified.cellchat.rds')
)

#### All Lymphoid -- lym.srt  ####
lym.srt <- full.srt[, full.srt$Cell_type == 'Lym' &
                            full.srt$Origin != 'Ambiguous' &
                            full.srt$Cell_state != 'BC']
lym.srt$Cell_state <- droplevels(lym.srt$Cell_state)
lym.srt$Origin <- droplevels(lym.srt$Origin)
Idents(lym.srt) <- 'Cell_state'

#### All cell types, all samples with doublets -- full_dbl.srt  ####
full_dbl.srt <- readRDS('integrated/STEP18.annotated_with_ambig.srt.rds')
full_raw.srt <- readRDS('integrated/STEP02.merged.cbn.srt.rds')
# full_raw.srt <- NormalizeData(full_raw.srt) |>
#         FindVariableFeatures() |>
#         ScaleData() |>
#         RunPCA() |>
#         RunUMAP(dims = 1:50)
# lq_cell <- readRDS('analysis/STEP03.cell_filtered.low_quality_cells.rds')
# full_raw.srt$LowQ <- F
# full_raw.srt@meta.data[unlist(lq_cell), 'LowQ'] <- T
# full_raw.srt$Ambiguous <- NA
# full_raw.srt@meta.data[Cells(full_dbl.srt)[!full_dbl.srt$Cell_type_non_ambig], 'Ambiguous'] <- T
# full_raw.srt@meta.data[Cells(full_dbl.srt)[full_dbl.srt$Cell_type_non_ambig], 'Ambiguous'] <- F
# saveRDS(full_raw.srt@meta.data, 'analysis/STEP91.full_raw.srt_meta.rds')
# saveRDS(full_raw.srt@reductions, 'analysis/STEP91.full_raw.srt_dr.rds')
full_raw.srt@meta.data <- readRDS('analysis/STEP91.full_raw.srt_meta.rds')
full_raw.srt@reductions <- readRDS('analysis/STEP91.full_raw.srt_dr.rds')
####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~----


####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
####  Supp Table 3  Unassigned Cells - R1Q3   ####
####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
mtx1 <- Table(full.srt$Cell_type[full.srt$name2 == '5d PT'], full.srt$Origin[full.srt$name2 == '5d PT'])
mtx2 <- Table(full.srt$Cell_type[full.srt$name2 == '15m PT'], full.srt$Origin[full.srt$name2 == '15m PT'])
mtx3 <- Table(full.srt$Cell_type[full.srt$name2 == '12y PT'], full.srt$Origin[full.srt$name2 == '12y PT'])
data <- as.data.frame(as.matrix(rbind(mtx1, mtx2, mtx3)[, 1:3]))
data$Sample <- rep(c('5d PT', '15m PT', '12y PT'), each = 11)
data$Cell_type <- rownames(data)[1:11]
data <- data[!data$Cell_type %in% c('Mito', 'Doublet'), ]
WriteCSV(data, 'STEP91.1.cell_count_per_origin_cell_type_and_sample')
####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~----


####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
####  Supp Fig 1  QC metrics - R1Q1   ####
####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
####  Plot 1  Violin Per cell UMI nGene and mito% QC  ####
lib_meta.df <- read.csv(paste0(Docu_dir, 'pediatric_sample_meta.csv'))
tmp.df <- readRDS('integrated/STEP02.merged.cbn.srt_meta.rds')
cutoffs <- readRDS('analysis/PART03.cell_filtered.low_quality_cutoffs.df.rds')
tmp.df <- inner_join(x = tmp.df, y = lib_meta.df, by = c('sample' = 'Sample_id'))
tmp.df <- inner_join(x = tmp.df, y = cutoffs, by = c('sample' = 'name'))
tmp.df$group3 <- paste(tmp.df$name, tmp.df$Replicate)
tmp.df$group3 <- str_replace(tmp.df$group3, pattern = '_0', replacement = ' ')
tmp.df$group3 <- str_replace(tmp.df$group3, pattern = '_', replacement = ' ')
tmp.df$group3[tmp.df$sample == 'P1_S003'] <- 'Control 1 Rep3'
tmp.df$group3[tmp.df$sample == 'P1_S004'] <- 'Control 1 Rep4'
tmp.df$group3 <- factor(tmp.df$group3, labels = str_sort(U(tmp.df$group3), numeric = T))
data2 <- tmp.df[!duplicated(tmp.df$sample),
                c('nCount_max', 'nCount_min', 'nFeature_min', 'nFeature_max', 'sample', 'group3')]

p1.1 <- ggplot(tmp.df) +
        geom_violin(aes(x = group3, y = log10(nCount_CBN), fill = group3)) +
        labs(x = '', y = 'Log10 UMI count per nucleus') +
        geom_point(data = data2,
                   aes(x = group3, y = log10(nCount_max)), shape = 3, size = 2) +
        geom_point(data = data2,
                   aes(x = group3, y = log10(nCount_min)), shape = 3, size = 2)
p1.2 <- ggplot(tmp.df) +
        geom_violin(aes(x = group3, y = nFeature_CBN/1e3, fill = group3)) +
        scale_y_continuous(limits = c(0, 15)) +
        labs(x = '', y = 'Gene count per nucleus') +
        geom_point(data = data2,
                   aes(x = group3, y = nFeature_max/1e3), shape = 3, size = 2) +
        geom_point(data = data2,
                   aes(x = group3, y = nFeature_min/1e3), shape = 3, size = 2)
p1.3 <- ggplot(tmp.df) +
        geom_violin(aes(x = group3, y = pct_mito_CBN, fill = group3), scale = 'width') +
        scale_y_continuous(limits = c(0, 10)) +
        labs(x = '', y = 'Mitochondrial UMI % per nucleus') +
        geom_point(data = data2,
                   aes(x = group3, y = 5), shape = 3, size = 2)

p1 <- p1.1/p1.2/p1.3 &
        scale_fill_manual(values = c(Color_transplant, Color_controls[c(1, 1, 1, 1, 2, 2, 3, 4)])) &
        theme_classic() +
        theme(aspect.ratio = 0.33) &
        RotatedAxis()
p1
PlotPDF('1.1.vln.ncount_nfeature_pctmito', 10, 10)
p1
dev.off()

####  Plot 2  UMAP Ambiguous cells before filter  ####
pct1 <- Table(full_raw.srt$LowQ)[['TRUE']]/ncol(full_raw.srt)
p2.1 <- DimPlot2(full_raw.srt, group.by = 'LowQ', cols = c('grey90', 'black'),
                 reduction = 'umap', raster = F) +
        labs(x = 'UMAP 1', y = 'UMAP 2', title = paste0('LowQ Cells (', round(pct1*100, 1),'%)'))
pct2 <- Table(full_raw.srt[, !full_raw.srt$LowQ]$Ambiguous)[['TRUE']]/ncol(full_raw.srt[, !full_raw.srt$LowQ])
p2.2 <- DimPlot2(full_raw.srt[, !full_raw.srt$LowQ], group.by = 'Ambiguous', cols = c('grey90', 'black'),
               reduction = 'umap', raster = F) +
        labs(x = 'UMAP 1', y = 'UMAP 2', title = paste0('Ambiguous Cells (', round(pct2*100, 1),'%)'))
p2.1 + p2.2
PlotPDF('1.2.umap.global_cell_type_with_ambiguous', 20, 10)
p2.1 + p2.2
dev.off()

####  Plot 3  UMAP Before and after integration  ####
p3 <- DimPlot2(full_raw.srt[, Cells(full.srt)],
               group.by = 'name2', cols = Color_condition, raster = F, pt.size = 0.1, order = T, reduction = 'umap') +
        labs(title = 'Condition', x = 'UMAP1', y = 'UMAP2')
p3
PlotPDF('1.3.umap.pc_umap', 10, 10)
p3
dev.off()

####  Plot 4  Alignment Score  ####
post_int <- AlignmentScore(seurat = full.srt, batch = 'sample', reduction = 'clean_umap',
                           k = NULL, downsample = 100, iter = 10)
pre_int <- AlignmentScore(seurat = tmp.srt, batch = 'sample', reduction = 'umap',
                          k = NULL, downsample = 100, iter = 10)
df <- data.frame(Score = c(pre_int, post_int),
                 Group = factor(rep(c('Pre', 'Post'), each = 10), levels = c('Pre', 'Post')))
p4 <- ggplot(df) +
        geom_boxplot(aes(x = Group, y = Score), outlier.shape = NA) +
        scale_y_continuous(limits = c(0, 1)) +
        theme_Publication()
p4
PlotPDF('1.4.box.alignment_score', 5, 5)
p4
dev.off()
####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~----


####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
####  Supp Fig 2  LYVE1 and CD163 - R1Q3   ####
####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
####  Plot 1 Violin LYVE1 CD163 CD163L1 CYBB expression  ####
tmp.srt <- mye.srt[, mye.srt$Cell_state_orig %in% c('Donor MP2', 'Recipient MP2', 'Mono')]
tmp.srt <- DropMetaLevels(tmp.srt)
p1 <- VlnPlot2(tmp.srt, features = c('MRC1', 'MERTK', 'CD163', 'LYVE1', 'CYBB', 'CD163L1'),
               cols = c(Color_genotype[c(3, 2)], Color_cell_state_mye[1]) , same.y.lims = 6,
               group.by = 'Cell_state_orig', pt.size = -1, ncol = 3) + RestoreLegend()
p1
PlotPDF('2.1.vln.lyve1_cd163_exp', 8, 8)
p1
dev.off()

gl <- c('MRC1', 'MERTK', 'CD163', 'LYVE1', 'CYBB', 'CD163L1')
pl <- list()
for(i in 1:6){
        pl[[i]] <- BoxPlot(tmp.srt, feature = gl[i], cols = c(Color_genotype[c(1, 3, 2)], Color_cell_state_mye[2]),
                           group.by = 'Cell_state_orig')}
p2 <- wrap_plots(pl, ncol = 3)   
p2
PlotPDF('2.1.vln.lyve1_cd163_exp', 8, 8)
p1
dev.off()

####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~----


####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
####  RNAScope Probes   - R1Q4 & R1Q6   ####
####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
deg <- readRDS('analysis/STEP23.mp2_resmp_late_momp2_deg.mk.rds')
deg$gene <- rownames(deg)
genes <- rownames(deg)[deg$avg_log2FC < 0]
tmp <- deg
tmp$p_val_adj[tmp$p_val_adj < 1e-20] <- 1e-20

g <- c('MRC1', 'ITGAM', 'PTPRC',
       'LYVE1', 'CYBB', 'MARCO',
       'IGSF21', 'CD163L1', 'PDGFB')
tmp.srt <- full.srt[, unlist(DownsampleByMeta(full.srt, meta_var = 'name2', down_to_min_group = T))]
p <- FeaturePlot2_Dark(tmp.srt, features = g,
                       split.by = 'name2', reduction = 'clean_umap', raster = T, pt.size = 3)
PlotPDF('3.1.feat.candidate_probe_gene', 10, 20)
p
dev.off()

tmp2.srt <- lym.srt[, lym.srt$Cell_state %in% c('NK1', 'NK2', 'NK3')]
levels(tmp2.srt) <- c('NK1', 'NK2', 'NK3')
mk <- FindAllMarkers(tmp2.srt, assay = 'CBN', test.use = 'MAST',
                     only.pos = T, logfc.threshold = 0.25, return.thresh = 0.001)
mk_dedup <- mk |> group_by(cluster) |> top_n(n = 30, wt = avg_log2FC)
mk_dedup <- mk_dedup[! mk_dedup$gene %in% mk_dedup$gene[duplicated(mk_dedup$gene)], ]

g <- c('CRTAM', 'CCL4L2', 'CCL4', 'XCL1', 'XCL2', 'IFNG')
tmp.srt <- full.srt[, unlist(DownsampleByMeta(full.srt, meta_var = 'name2', down_to_min_group = T))]
p <- FeaturePlot2_Dark(tmp.srt, features = g,
                       split.by = 'name2', reduction = 'clean_umap', raster = T, pt.size = 3)
p
PlotPDF('3.2.heat.nk_candidate_probe_gene', 15, 10)
p
dev.off()

p1 <- FeaturePlot2_Dark(full.srt, features = 'CRTAM', reduction = 'clean_umap', raster = F,
                        pt.size = 0.2, min.cutoff = 3.5)
p2 <- FeaturePlot2_Dark(lym.srt, features = 'CRTAM', reduction = 'sub_umap', raster = F,
                        pt.size = 1, min.cutoff = 3.5)
p1 | p2
PlotPDF('3.3.feat.crtam_nk3_marker', 10, 5)
p1 | p2
dev.off()

g <- c('NR4A3', 'JARID2', 'RNF19A', 'CD82')
FeaturePlot2_Dark(lym.srt, features = g, reduction = 'sub_umap', raster = F,
                  pt.size = 0.2, min.cutoff = 'q5')



g <- c('CYBB', 'VAV1', 'CD163L1',  'PDGFB', 'CD163')
tmp.srt <- full.srt[, unlist(DownsampleByMeta(full.srt, meta_var = 'name2', down_to_min_group = T))]
p <- FeaturePlot2_Dark(tmp.srt, features = c(g, 'CRTAM'),
                       split.by = 'name2', reduction = 'clean_umap', raster = T, pt.size = 3)
p
Idents(tmp.srt) <- 'Cell_type'
VlnPlot2(tmp.srt, features = g, group.by = 'name2', idents = 'Mye', ncol = 2)
DotPlot2(tmp.srt, features = g, group.by = 'Cell_type')
Idents(tmp.srt) <- 'Cell_state'
VlnPlot2(tmp.srt, features = 'CRTAM', group.by = 'name2', idents = c('NK1', 'NK2', 'NK3'), ncol = 1, pt.size = 0.1)

DefaultAssay(momp_p136.srt) <- 'alra'
tmp.srt <- momp_p136.srt
tmp.srt$cd163_cybb <- 'Neg'
tmp.srt$cd163_cybb[tmp.srt@assays$CBN@data['CYBB', ] > 0.5] <- 'CYBB'
tmp.srt$cd163_cybb[tmp.srt@assays$CBN@data['CD163', ] > 4] <- 'CD163'
tmp.srt$cd163_cybb[tmp.srt@assays$CBN@data['CD163', ] > 4 & tmp.srt@assays$CBN@data['CYBB', ] > 0.5] <- 'pos'
DimPlot2(tmp.srt, reduction = 'sub_umap', group.by = 'cd163_cybb', 
         cols = c('yellow', 'lightblue', 'grey90', 'black'), pt.size = 2)

FeatureScatter(momp_p136.srt, feature1 = 'CYBB',feature2 =  'CD163', group.by = 'Cell_state_orig')

DefaultAssay(mye.srt) <- 'CBN'
FeaturePlot2_Dark(mye.srt, features = c('CYBB', 'CD163L1', 'CD163'), split.by = 'name2',
                 reduction = 'sub_umap', raster = T, pt.size = 3, min.cutoff = 1)

DotPlot2(mye.srt, features = c('CYBB', 'CD163L1', 'CD163'), group.by = 'name2')


DotPlot2(full.srt[, full.srt$name2 == '5d PT'], features = c('CD163L1', 'CD163')) +
        DotPlot2(momp_p136.srt, features = c('CD163L1', 'CD163'), 
                 group.by = 'Cell_state_orig', col.min = 0, col.max = 2)
####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~----


####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
####  RNAScope Quantification   - R1Q4 & R1Q6   ####
####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
library(ggpubr)
#### CD163L1+CD163+  ####
data <- data.frame('Sample' = factor(rep(c('5dPT', '15mPT', '12yPT'), each = 5), levels = c('5dPT', '15mPT', '12yPT')),
                   'Total_DAPI' = c(188, 163, 173, 122, 145,
                                    309, 118, 279, 326, 345,
                                    144, 119, 178, 115, 105),
                   'CD163L1_DAPI' = c(0, 0, 0, 0, 0,
                                      11, 14, 11, 13, 13,
                                      7, 4, 12, 10, 6))
data$Percentage <- data$CD163L1_DAPI/data$Total_DAPI*100
p1.1 <- ggplot(data, aes(x = Sample, y = Percentage)) +
        geom_boxplot(outlier.shape = NA) +
        geom_quasirandom(size = 1.5, color = 'black') +
        stat_compare_means(comparisons = list(c("5dPT", "15mPT"), c("5dPT", "12yPT"), c("15mPT", "12yPT")),
                           method = "t.test",  paired = F,  label = "p.format", hide.ns = T)
data2 <- data |> group_by(Sample) |> mutate(`Mean Percentage` = sum(CD163L1_DAPI)*100/sum(Total_DAPI))
data2 <- U(data2[, c(1, 5)])
p1.2 <- ggplot(data2) +
        geom_bar(aes(x = Sample, y = `Mean Percentage`), stat = 'identity')
p <- p1.1 + p1.2 &
        scale_y_continuous(limits = c(0, 15)) &
        theme_Publication() &
        theme(aspect.ratio = 3)
p
PlotPDF('3.4.rec_mp2_cd163l2', 8, 5)
p
dev.off()

#### CYBB+CD163+  ####
data <- data.frame('Sample' = factor(rep(c('5dPT', '15mPT', '12yPT'), each = 5), levels = c('5dPT', '15mPT', '12yPT')),
                   'Total_DAPI' = c(162, 106, 126, 189, 163,
                                    336, 290, 313, 353, 318,
                                    123, 88, 90, 212, 180),
                   'CYBB_DAPI' = c(14, 8, 10, 23, 13,
                                   0, 0, 0, 0, 0,
                                   0, 2, 2, 7, 4))
data$Percentage <- data$CYBB_DAPI/data$Total_DAPI*100
p2.1 <- ggplot(data, aes(x = Sample, y = Percentage)) +
        geom_boxplot(outlier.shape = NA) +
        geom_quasirandom(size = 1.5, color = 'black') +
        stat_compare_means(comparisons = list(c("5dPT", "15mPT"), c("5dPT", "12yPT"), c("15mPT", "12yPT")),
                           method = "t.test",  paired = F,  label = "p.format", hide.ns = T)
data2 <- data |> group_by(Sample) |> mutate(`Mean Percentage` = sum(CYBB_DAPI)*100/sum(Total_DAPI))
data2 <- U(data2[, c(1, 5)])
p2.2 <- ggplot(data2) +
        geom_bar(aes(x = Sample, y = `Mean Percentage`), stat = 'identity')
p <- p2.1 + p2.2 &
        scale_y_continuous(limits = c(0, 15)) &
        theme_Publication() &
        theme(aspect.ratio = 3)
p
PlotPDF('3.5.don_mp2_cybb', 8, 5)
p
dev.off()

p1 <- DotPlot2(mye.srt[, mye.srt$name2 != 'Control'], features = c('CYBB', 'CD163L1', 'CD163'), 
               group.by = 'name2',  scale = 'radius',
              col.min = -2) +
        scale_color_gradient(low = 'grey85', high = mycol_BuGr[1], limits = c(-2, 1.5))
p1
p2 <- VlnPlot2(mye.srt[, mye.srt$name2 != 'Control'], features = c('CYBB', 'CD163L1', 'CD163'), group.by = 'name2',
               cols = Color_transplant) /
        VlnPlot2(mye.srt[, mye.srt$Cell_state_orig %in% c('Donor MP2', 'Recipient MP2')], 
                 features = c('CYBB', 'CD163L1', 'CD163'), group.by = 'Cell_state_orig',
                 cols = Color_cell_state_orig_mye)
p2
PlotPDF('3.6.dot_vln.cybb_cd163l1', 10, 4)
p1+p2
dev.off()

####  PDGFB+CD163+  ####
data <- data.frame('Sample' = factor(rep(c('Cont', '5dPT', '15mPT', '12yPT'), each = 5),
                                     levels = c('Cont', '5dPT', '15mPT', '12yPT')),
                   'Total_DAPI' = c(447, 367, 450, 410, 400,
                                    142, 189, 165, 171, 164,
                                    344, 353, 316, 328, 302,
                                    112, 120, 178, 121, 171),
                   'PDGFB_DAPI' = c(0, 3, 2, 4, 0,
                                    0, 0, 0, 0, 0,
                                    3, 4, 2, 0, 2,
                                    7, 9, 17, 8, 8))
data$Percentage <- data$PDGFB_DAPI/data$Total_DAPI*100
p3.1 <- ggplot(data, aes(x = Sample, y = Percentage)) +
        geom_boxplot(outlier.shape = NA) +
        geom_quasirandom(size = 1.5, color = 'black') +
        stat_compare_means(comparisons = list(c('Cont', '5dPT'), c('Cont', '15mPT'), c('Cont', '12yPT'), 
                                              c("5dPT", "15mPT"), c("5dPT", "12yPT"), 
                                              c("15mPT", "12yPT")),
                           method = "t.test",  paired = F,  label = "p.format", hide.ns = T)
data2 <- data |> group_by(Sample) |> mutate(`Mean Percentage` = sum(PDGFB_DAPI)*100/sum(Total_DAPI))
data2 <- U(data2[, c(1, 5)])
p3.2 <- ggplot(data2) +
        geom_bar(aes(x = Sample, y = `Mean Percentage`), stat = 'identity')
p <- p3.1 + p3.2 &
        scale_y_continuous(limits = c(0, 15)) &
        theme_Publication() &
        theme(aspect.ratio = 1)
p
PlotPDF('3.7.fibrosis_pdgfb_cd163', 8, 5)
p
dev.off()


#### CRTAM+cTNT  #### 
data <- data.frame('Sample' = factor(rep(c('Cont', '15mPT', '12yPT', 'DCM-HF'), each = 5),
                                     levels = c('Cont', '15mPT', '12yPT', 'DCM-HF')),
                   'Total_DAPI' = c(382, 302, 284, 340, 388,
                                    342, 317, 279, 245, 329,
                                    151, 133, 143, 161, 163,
                                    165, 134, 153, 133, 161),
                   'CRTAM_DAPI' = c(0, 0, 0, 0, 0,
                                    4, 7, 11, 3, 3,
                                    0, 0, 0, 1, 0,
                                    1, 3, 4, 1, 5))
data$Percentage <- data$CRTAM_DAPI/data$Total_DAPI*100
p4.1 <- ggplot(data, aes(x = Sample, y = Percentage)) +
        geom_boxplot(outlier.shape = NA) +
        geom_quasirandom(size = 1.5, color = 'black') +
        stat_compare_means(comparisons = list(c('Cont', '15mPT'), c('Cont', '12yPT'), c('Cont', 'DCM-HF'), 
                                              c("15mPT", "12yPT"), c("15mPT", "DCM-HF"), 
                                              c("12yPT", "DCM-HF")),
                           method = "t.test",  paired = F,  label = "p.format", hide.ns = T)
data2 <- data |> group_by(Sample) |> mutate(`Mean Percentage` = sum(CRTAM_DAPI)/sum(Total_DAPI))
data2 <- U(data2[, c(1, 5)])
p4.2 <- ggplot(data2) +
        geom_bar(aes(x = Sample, y = `Mean Percentage`), stat = 'identity')
p <- p4.1 + p4.2 &
        theme_Publication() &
        theme(aspect.ratio = 1)
p
PlotPDF('3.8.nk3_crtam', 8, 5)
p
dev.off()


####  PDGFB+CD163+  ####
data <- data.frame('Sample' = factor(rep(c('Cont', '5dPT', '15mPT', '12yPT'), each = 3),
                                     levels = c('Cont', '5dPT', '15mPT', '12yPT')),
                   'Total_DAPI' = c(1213, 1396, 1260,
                                    751, 922, 732, 
                                    1047, 891, 1160, 
                                    807, 712, 838),
                   'Total_CD163' = c(80, 77, 74, 
                                     242, 243, 281, 
                                     209, 229, 186, 
                                     116, 113, 120),
                   'PDGFB_CD163' = c(0, 0, 0,
                                     0, 0, 0,
                                     2, 2, 1,
                                     6, 7, 6))
data$Percentage <- data$PDGFB_CD163/data$Total_CD163*100
p3.1 <- ggplot(data, aes(x = Sample, y = Percentage)) +
        geom_boxplot(outlier.shape = NA) +
        geom_quasirandom(size = 1.5, color = 'black') +
        stat_compare_means(comparisons = list(c('Cont', '5dPT'), c('Cont', '15mPT'), c('Cont', '12yPT'), 
                                              c("5dPT", "15mPT"), c("5dPT", "12yPT"), 
                                              c("15mPT", "12yPT")),
                           method = "wilcox",  paired = F,  label = "p.format", hide.ns = F)
data2 <- data |> group_by(Sample) |> mutate(`Mean Percentage` = sum(PDGFB_CD163)*100/sum(Total_CD163))
data2 <- U(data2[, c(1, 6)])
p3.2 <- ggplot(data2) +
        geom_bar(aes(x = Sample, y = `Mean Percentage`), stat = 'identity')
p <- p3.1 + p3.2 &
        scale_y_continuous(limits = c(0, 8)) &
        theme_Publication() &
        theme(aspect.ratio = 1)
p
PlotPDF('3.9.fibrosis_pdgfb_cd163_v2', 8, 5)
p
dev.off()
####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~----


####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
####  Fig 4 & Supp Fig 4  MP2 Phagocytosis and Fibrotic Signaling in DCM - R2Q5   ####
####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
dcm.srt <- readRDS('external/2022_NatCardRes_KLavine.dcm_donor_scrna.srt.rds')
grp1 <- read.table('external/GOBP_PHAGOCYTOSIS.v7.5.1.grp',
                   comment.char = '#', header = T)[,1]
TRMP_gene <- c('CD163L1', 'MRC1', 'MAF', 'SIGLEC1', 'CCL8', 'CCL14', 'LILRB5', 'LYVE1', 'IL2RA', 'PDGFC', 'WLS', 'DAB2',
               'NRP1', 'SCN9A', 'FGF13', 'GDF15', 'IGF1', 'FMOD', 'SLIT3', 'EGFL7', 'ECM1', 'SDC3')
MoMP_gene <- c('CCL17', 'CCR5', 'CCR2', 'CXCL9', 'CXCL3', 'CXCL2', 'CXCL10', 'CSF2RA', 'CCL5', 'CCL20', 'TREM1', 'TET2',
               'NFKBIE', 'IL1R1', 'IL1R2', 'NFKB1', 'REL', 'MAP2K3', 'IL1A', 'NLRP3', 'TRAF1', 'MAPK6', 'SRC', 'IL1B',
               'NOD2', 'JAK3', 'IRAK2', 'MAP3K8', 'RELB', 'SOCS3', 'MYD88', 'MMP9', 'IL20', 'AREG', 'PTX3', 'IL23A',
               'IL10', 'OSM', 'TIMP1', 'EREG', 'IL27') ## These genes are from the K. Lavine's Nature Medicine Paper

####  Plot 1 UMAP DCM MP2   ####
dcm_mp.srt <- dcm.srt[, dcm.srt@active.ident %in% c('Monocytes', 'Macrophages')]
dcm_mp.srt <- AddModuleScore2(dcm_mp.srt, features = list(TRMP_gene, MoMP_gene), names = c('TRMP', 'MoMP'), nbin = 20)
p1 <- FeaturePlot2(dcm_mp.srt, features = 'TRMP', min.cutoff = -2, max.cutoff = 2, pt.size = 0.5,
                   split.by = 'Condition', ncol = 1) +
        NoLegend()

dcm_mp.srt$RMP_bin <- 'Other'
dcm_mp.srt$RMP_bin[dcm_mp.srt$Condition == 'DCM' & dcm_mp.srt$TRMP > 0.5] <- 'DCM MP2'
dcm_mp.srt$RMP_bin[dcm_mp.srt$Condition == 'Donor' & dcm_mp.srt$TRMP > 0.5] <- 'Donor MP2'

p2 <- DimPlot2(dcm_mp.srt[, dcm_mp.srt$Condition == 'Donor'], group.by = 'RMP_bin',
               cols = c('grey10', 'grey90'), pt.size = 0.5,
               split.by = 'Condition', ncol = 1)
p3 <- DimPlot2(dcm_mp.srt[, dcm_mp.srt$Condition == 'DCM'], group.by = 'RMP_bin',
               cols = c('seagreen3', 'grey90'), pt.size = 0.5,
               split.by = 'Condition', ncol = 1)
p <- wrap_plots(p1[[2]], p1[[1]], p2, p3, ncol = 2)

PlotPDF('4.1.umap_feat.dcm_mp2', 8, 8)
p
dev.off()

####  Plot 2 Feature DCM MP2   ####
rmp.srt <- dcm_mp.srt[, dcm_mp.srt$RMP_bin != 'Other']
DefaultAssay(rmp.srt) <- 'RNA'
rmp.srt <- DietSeurat(rmp.srt, assays = 'RNA')
rmp.srt$Group <- rmp.srt$Condition
rmp.srt <- RenameAssays(rmp.srt, 'RNA' = 'CBN')
rmp.srt <- NormalizeData(rmp.srt)
tmp.srt <- merge(mp2_sub.srt, rmp.srt)
tmp.srt <- AddModuleScore2(tmp.srt, features = list(grp1), names = 'Phago')
p3 <- BoxPlot(tmp.srt, feature = 'Phago', group.by = 'Group')
data <- tmp.srt@meta.data[, c('Phago', 'Group')]
t.test(as.vector(data[data$Group == 'DCM', 1]),
       as.vector(data[data$Group == 'Donor', 1]))$p.value
p3

PlotPDF('4.2.box.mp2_phago_geneset_expr', 6, 4)
p3
dev.off()

####  Plot 3 Box Pro-Fibrosis Gene  ####
DefaultAssay(dcm_mp.srt) <- 'RNA'
dcm_mp.srt <- NormalizeData(dcm_mp.srt)
dcm_mp.srt <- RunALRA(dcm_mp.srt, genes.use = rownames(dcm_mp.srt), assay = 'RNA')
rmp.srt <- dcm_mp.srt[, dcm_mp.srt$RMP_bin != 'Other']
rmp.srt <- DietSeurat(rmp.srt, assays = c('RNA', 'alra'))
rmp.srt$Group <- rmp.srt$Condition
rmp.srt <- RenameAssays(rmp.srt, 'RNA' = 'CBN')
tmp.srt <- merge(mp2_sub.srt, rmp.srt)

Idents(tmp.srt) <- 'Group'
levels(tmp.srt) <- c( 'ResMP 5d PT','MoMP2 12y PT','MoMP2 15m PT','Donor','DCM')
p4 <- DotPlot2(tmp.srt, features = c('PDGFB', 'WNT5B', 'TGFB2', 'IGF1', 'TGFA'), assay = 'alra',
               col.max = 0.5)
p4
PlotPDF('4.3.dot.dcm_mp2_fibrotic', 6, 4)
p4
dev.off()

####  Plot 4 Box Fibrosis quantification  ####
data <- read_xlsx(paste0('~/Documents/Bioinformatics/project/2022_acute_rejection_dturaga/',
                         'doc/human_v1/Fibrosis_Area_Quantification.xlsx'), sheet = 2)
data <- data*100
data <- melt(data)
data <- data[! is.na(data$value) & data$variable %in% c('3B26D', 'P136', 'P182', 'P170', 'P146', 'P154'), ]
data$Sample <- factor(revalue(data$variable, c(
        '3B26D' = 'Control',
        'P136' = '5dPT',
        'P182' = '15mPT',
        'P170' = '12yPT',
        'P146' = 'DCM HF 1',
        'P154' = 'DCM HF 2'
)), levels = c('Control', '5dPT', '15mPT', '12yPT', 'DCM HF 1', 'DCM HF 2'))
p5.1 <- ggplot(data) +
        geom_boxplot(aes(x = Sample, y = value), outlier.shape = NA) +
        ggbeeswarm::geom_quasirandom(aes(x = Sample, y = value, color = Sample)) +
        scale_color_manual(values = c('grey50', Color_transplant, 'grey20', 'grey20')) +
        scale_y_continuous(limits = c(0, 100)) +
        theme_classic() +
        theme_Publication() +
        theme(aspect.ratio = 2) +
        NoLegend()
p5.1

data <- data.frame(Sample = factor(c('Control', '5dPT', '15mPT', '12yPT', 'DCM HF 1', 'DCM HF 2'),
                                   levels = c('Control', '5dPT', '15mPT', '12yPT', 'DCM HF 1', 'DCM HF 2')),
                   Pct = c(5.62058, 7.5048, 19.5001, 19.695529, 11.8514, 11.1726))
p5.2 <- ggplot(data, aes(color = Sample, y = Pct, x = Sample)) +
        geom_point(size = 5) +
        scale_color_manual(values = c('grey50', Color_transplant, 'grey20', 'grey20')) +
        scale_y_continuous(limits = c(0, 20), breaks = seq(0, 20, 5)) +
        labs(y = 'Percent Area') +
        theme_classic() +
        theme(aspect.ratio = 2) +
        NoLegend()
p5.1 | p5.2
PlotPDF('4.4.dot.fibrotic_area_quantification', 8, 4)
p5.1 | p5.2
dev.off()
####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~----



####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
####  Fig 4  Donor MP2 vs Recipient MP2 - R1Q3   ####
####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
tmp.srt <- mp2.srt
Idents(tmp.srt) <- 'Origin'
DimPlot2(tmp.srt, reduction = 'sub_umap', group.by = c('Origin', 'Group'), cols = mycol_10)

####  Plot 1 umap Donor MP2 vs Recipient MP2  ####
p1 <- DimPlotSplit(mp2.srt, reduction = 'sub_umap',  split_by = 'Origin', split_order = levels(mp2.srt$Origin),
                   ncol = 2, cols.highlight = Color_mp2_group)
p1[[3]] <- NULL
PlotPDF('5.1.umap.donor_mp2_vs_recipient_mp2', 8, 4)
p1
dev.off()


deg <- FindMarkers(tmp.srt, assay = 'CBN', test.use = 'bimod',
                   ident.1 = 'Recipient', ident.2 = 'Donor', logfc.threshold = 0.25)
deg$gene <- rownames(deg)
deg <- deg[deg$p_val_adj <= 0.01, ]
deg <- deg[, c('gene', 'avg_log2FC', 'p_val_adj')]
colnames(deg) <- c('Gene_symbol', 'Mean_log2_fold_change', 'Adjusted_p_value')
deg$Gene_symbol <- paste0("'", deg$Gene_symbol, "'")
deg <- deg[order(deg$Mean_log2_fold_change, decreasing = T),]
deg$Differentially_expressed <- 'Recipient MP2 High'
deg$Differentially_expressed[deg$Mean_log2_fold_change < 0] <- 'Donor MP2 High'
saveRDS(deg, 'analysis/STEP91.mp2_donor_vs_recipient_deg.mk.rds')


####  Table 1 Donor MP2 vs Recipient MP2 DEG  ####
deg <- readRDS('analysis/STEP91.mp2_donor_vs_recipient_deg.mk.rds')
WriteCSV(deg, 'STEP91.5.1.table.donor_mp2_vs_recipient_mp2_deg')



####  Plot 2 volcano Donor MP2 vs Recipient MP2 DEG  ####
deg <- readRDS('analysis/STEP91.mp2_donor_vs_recipient_deg.mk.rds')
genes <- rownames(deg)[(deg$p_val_adj <= 0.01 & abs(deg$avg_log2FC) > 1) | deg$p_val_adj <= 1e-5]
tmp <- deg
tmp$p_val_adj[tmp$p_val_adj < 1e-20] <- 1e-20
p2 <- MarkerVolcano(tmp, label = T, label_genes = genes)
p2
PlotPDF('5.2.vol.donor_mp2_vs_recipient_mp2', 10, 10)
p2
dev.off()



tmp.srt <- AddModuleScore2(tmp.srt, features = genelist, names = names(genelist), assay = 'CBN', return_z = T)
BoxPlot(tmp.srt, feature = 'Donor', group.by = 'Group') + BoxPlot(tmp.srt, feature = 'Recipient', group.by = 'Group')


genelist <- split(str_remove_all(deg$Gene_symbol, "'"), deg$Differentially_expressed)
names(genelist) <- c('Donor', 'Recipient')
deg_orig <- readRDS('analysis/STEP23.mp2_resmp_late_momp2_deg.mk.rds')
deg_orig_hi <- deg_orig[deg_orig$p_val_adj < 0.05 & abs(deg_orig$avg_log2FC) >= 0.25 , ]
deg_orig_hi$gene <- rownames(deg_orig_hi)
deg_orig_hi$cluster <- 'Late MoMP2'
deg_orig_hi$cluster[deg_orig_hi$avg_log2FC < 0] <- 'ResMP'
genelist_orig <- split(deg_orig_hi$gene, deg_orig_hi$cluster)
O(genelist$Donor, genelist_orig$ResMP)
O(genelist$Recipient, genelist_orig$`Late MoMP2`)

p3.1 <- Plot2WayVenn(all_genes = rownames(tmp.srt),
                     x = genelist$Donor, genelist_orig$ResMP, x_name = 'New', y_name = 'Old')
p3.2 <- Plot2WayVenn(all_genes = rownames(tmp.srt),
                     x = genelist$Recipient, genelist_orig$`Late MoMP2`, x_name = 'New', y_name = 'Old')
p3 <- wrap_plots(p3.1, p3.2)
p3
PlotPDF('5.3.venn.comparison_with_previous', 10, 5)
p3
dev.off()


####  Plot 2 bar Donor MP2 vs Recipient MP2 deg GO enrich ####
enrich <- ModuleEnrichment(module_list = genelist, human_or_mouse = 'human')
saveRDS(enrich, 'analysis/STEP91.donor_mp2_recipient_mp2_deg.enrich.rds')
enrich <- readRDS('analysis/STEP91.donor_mp2_recipient_mp2_deg.enrich.rds')
donor_go <- enrich$GO$GO_Donor
donor_go <- donor_go[donor_go$p.adjust < 0.1, ]
recipient_go <- enrich$GO$GO_Recipient
recipient_go <- recipient_go[recipient_go$p.adjust < 0.1, ]

go_choice_1 <- c(
        'phagocytosis',
        'Fc-gamma receptor signaling pathway',
        'negative regulation of leukocyte activation',
        'vesicle organization',
        'toll-like receptor signaling pathway'
)
gene_choice_1 <- list()
for(i in 1:L(go_choice_1)){
        pool <- unlist(str_split(donor_go$geneID[donor_go$Description == go_choice_1[i]], '/'))
        pool_sort <- pool[order(deg[pool, 'Mean_log2_fold_change'], decreasing = F)]
        gene_choice_1[[i]] <- head(pool_sort, 5)
}
df <- donor_go[donor_go$Description %in% go_choice_1, ]
df$Description <- factor(df$Description, levels = go_choice_1)
df <- df[order(df$Description),]
df$gene <- NA
for(i in 1:L(go_choice_1)){df$gene[i] <- paste(gene_choice_1[[i]], collapse = '  ')}
p4.1 <- ggplot(df, aes(x = -log10(p.adjust), y = Description))+
        geom_bar(stat = 'identity', position = "dodge", fill = mycol_10[1], alpha = 0.4) +
        geom_text(aes(label = gene, x = max(-log10(p.adjust))/100, y = Description),
                  color = 'black', size = 3, hjust = 0, face = 'italic') +
        scale_y_discrete(limit = rev) +
        labs(x = '-Log10 adjusted p value', y = 'Enriched GO terms', title = 'Donor MP2') +
        theme_classic() +
        theme(aspect.ratio = 1)

go_choice_2 <- c(
        'actin filament polymerization' ,
        'proteoglycan metabolic process',
        'chondroitin sulfate biosynthetic process',
        'regulation of tolerance induction'
)
gene_choice_2 <- list()
for(i in 1:L(go_choice_2)){
        pool <- unlist(str_split(recipient_go$geneID[recipient_go$Description == go_choice_2[i]], '/'))
        pool_sort <- pool[order(deg[pool, 'Mean_log2_fold_change'], decreasing = T)]
        gene_choice_2[[i]] <- head(pool_sort, 5)
}
df <- recipient_go[recipient_go$Description %in% go_choice_2, ]
df$Description <- factor(df$Description, levels = go_choice_2)
df <- df[order(df$Description),]
df$gene <- NA
for(i in 1:L(go_choice_2)){df$gene[i] <- paste(gene_choice_2[[i]], collapse = '  ')}
p4.2 <- ggplot(df, aes(x = -log10(p.adjust), y = Description))+
        geom_bar(stat = 'identity', position = "dodge", fill = mycol_10[1], alpha = 0.4) +
        geom_text(aes(label = gene, x = max(-log10(p.adjust))/100, y = Description),
                  color = 'black', size = 3, hjust = 0, face = 'italic') +
        scale_y_discrete(limit = rev) +
        labs(x = '-Log10 adjusted p value', y = 'Enriched GO terms', title = 'Recipient MP2') +
        theme_classic() +
        theme(aspect.ratio = 0.8)
p4.1/p4.2
PlotPDF('5.4.bar.resmp_vs_late_momp2_top_go_enrich', 10, 10)
p4.1/p4.2
dev.off()


####  Plot 4 box phago glyco functional geneset  ####
# grp1 <- read.table('external/GOBP_PHAGOCYTOSIS.v7.5.1.grp', comment.char = '#', header = T)[,1]
# grp1 <- intersect(grp1, VariableFeatures(mp2_sub.srt)) ## 61
grp2 <- read.table('external/KEGG_FC_GAMMA_R_MEDIATED_PHAGOCYTOSIS.v7.5.1.grp', comment.char = '#', header = T)[,1]
grp2 <- intersect(grp2, VariableFeatures(mp2_sub.srt)) ## 27
grp3 <- read.table('external/KEGG_GLYCOSAMINOGLYCAN_BIOSYNTHESIS_CHONDROITIN_SULFATE.v7.5.1.grp',
                   comment.char = '#', header = T)[,1]
grp3 <- intersect(grp3, VariableFeatures(mp2_sub.srt)) ## 8

tmp.srt <- AddModuleScore2(mp2.srt, features = list(grp2, grp3),
                           names = c('KEGG_phago', 'KEGG_glyco'),
                           return_z = T)
p5.1 <- BoxPlot(tmp.srt, group.by = 'Group', feature = 'KEGG_phago', cols = Color_mp2_group) + NoLegend()
p5.2 <- BoxPlot(tmp.srt, group.by = 'Group', feature = 'KEGG_glyco', cols = Color_mp2_group)
p5 <- p5.1 + p5.2
p5
data <- tmp.srt@meta.data[, c('KEGG_phago', 'KEGG_glyco', 'Group')]
pval <- list()
for(i in 1:2){
        pval[[i]] <- c(t.test(as.vector(data[data$Group == 'ResMP 5d PT', i]),
                              as.vector(data[data$Group == 'MoMP2 15m PT', i]))$p.value,
                       t.test(as.vector(data[data$Group == 'ResMP 5d PT', i]),
                              as.vector(data[data$Group == 'MoMP2 12y PT', i]))$p.value)
}
all(unlist(pval) < 0.001)
PlotPDF('5.5.box.mp2_phago_geneset_expr', 8, 4)
p5
dev.off()
####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~----


####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
####  Fig 6  NK3 cells not in HF - R1Q6   ####
####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
que.srt <- lym.srt[, lym.srt$Cell_state %in% c('NK1', 'NK2', 'NK3')]
que.srt$Study <- 'Allograft'
que.srt$Group <- revalue(que.srt$name2, replace = c('Control' = 'HTx_Control', 
                                                    '5d PT' = 'HTx_PT',
                                                    '15m PT'  = 'HTx_PT',
                                                    '12y PT' = 'HTx_PT'))
que.srt$Donor <- paste0(que.srt$Group, '_', que.srt$recipient)
que.srt@reductions$final_umap <- que.srt@reductions$sub_umap
que.srt <- RenameAssays(que.srt, 'CBN' = 'RNA')

####  Load NK Data
# ref1.srt <- readRDS('/Volumes/shire/data/scrnaseq/2022_Nature_JMartin/matrix_public/')
#### DCM (Keonig et al) data ####
ref1.srt <- nk_dcm.srt
ref1.srt <- ref1.srt[, ref1.srt$Study == 'Lavine']
ref1.srt$Group <- revalue(ref1.srt$name2, replace = c('Donor' = 'DCM_Control'))
ref1.srt$Donor <- paste0(ref1.srt$Group, '_', ref1.srt$name)
ref1.srt$Study <- 'DCM'

#### MI data ####
ref2.srt <- readRDS('/Volumes/shire/data/scrnaseq/2022_Nature_RKramann/matrix_public/Lymphoid_snRNA_snATAC.Rds')
ref2.srt.bkp <- ref2.srt
ref2.srt <- DietSeurat(ref2.srt, assays = 'RNA', dimreducs = names(ref2.srt))
ref2.srt$Group <- revalue(ref2.srt$patient_group, replace = c('group_1' = 'MI_Myogenic',
                                                              'group_2' = 'MI_Ischemic',
                                                              'group_3' = 'MI_Fibrotic',
                                                              'control' = 'MI_Control'))
ref2.srt$Donor <- paste0(ref2.srt$Group, '_', ref2.srt$donor)
ref2.srt$Study <- 'MI'
ref2.srt$donor <- ref2.srt$patient
ref2.srt@reductions$final_umap <- ref2.srt@reductions$umap_harmony_v2
ref2.srt$Cell_state <- ref2.srt$annotation
DimPlot2(ref2.srt, reduction = 'umap_harmony_v2', group.by = c('Group', 'annotation'), cols = mycol_10)

#### Cardiomyopathy data  ####
ref3.srt <- readRDS('/Volumes/shire/data/scrnaseq/2022_Science_CSeidman/matrix_public/lymphoid.rds')
ref3.srt.bkp <- ref3.srt
genes <- ConvertGeneID(rownames(ref3.srt), species = 'human')
mtx <- GetAssayData(ref3.srt)
mtx <- mtx[genes$GENEID, ]
identical(rownames(mtx), genes$GENEID)
rownames(mtx) <- genes$SYMBOL
mtx <- mtx[!duplicated(rownames(mtx)), ]
tmp.srt <-  CreateSeuratObject(counts = mtx, meta.data = ref3.srt@meta.data, min.cells = 1, min.features = 1)
tmp.srt@reductions <- ref3.srt@reductions
tmp.srt -> ref3.srt
ref3.srt$Study <- 'Myopathy'
ref3.srt$donor <- ref3.srt$donor_id
ref3.srt@reductions$final_umap <- ref3.srt@reductions$umap
ref3.srt$Cell_state <- ref3.srt$cell_states
ref3.srt$Group <- revalue(ref3.srt$disease, replace = c('non-compaction cardiomyopathy' = 'CM_NCCM',
                                                        'dilated cardiomyopathy' = 'CM_DCM',
                                                        'arrhythmogenic right ventricular cardiomyopathy' = 'CM_ACM',
                                                        'normal' = 'CM_Control'))
ref3.srt$Donor <- paste0(ref3.srt$Group, '_', ref3.srt$donor)
Idents(ref3.srt) <- 'cell_states'
DimPlot2(ref3.srt, reduction = 'umap', label = T, cols = mycol_30)

#### HTx endocardial biopsy data  ####
# srt1 <- CreateSeuratObject(ReadCB_h5('~/Downloads/matrix_public/GSM6176135_RV_Control_1/KA_1_output_filtered.h5'))
# srt1$Sample <- 'Control_1_RV'
# srt2 <- CreateSeuratObject(ReadCB_h5('~/Downloads/matrix_public/GSM6176136_RV_Control_2/KA_2_output_filtered.h5'))
# srt2$Sample <- 'Control_2_RV'
# srt3 <- CreateSeuratObject(ReadCB_h5('~/Downloads/matrix_public/GSM6176137_RV_CAV_1/KA_3_output_filtered.h5'))
# srt3$Sample <- 'CAV_1_RV'
# srt4 <- CreateSeuratObject(ReadCB_h5('~/Downloads/matrix_public/GSM6176138_LV_CAV_2/KA_4_output_filtered.h5'))
# srt4$Sample <- 'CAV_2_LV'
# srt5 <- CreateSeuratObject(ReadCB_h5('~/Downloads/matrix_public/GSM6176139_RV_Control_3/KA_5_output_filtered.h5'))
# srt5$Sample <- 'Control_3_RV'
# srt6 <- CreateSeuratObject(ReadCB_h5('~/Downloads/matrix_public/GSM6176140_RV_CAV_2/KA_6_output_filtered.h5'))
# srt6$Sample <- 'CAV_2_RV'
# srt7 <- CreateSeuratObject(ReadCB_h5('~/Downloads/matrix_public/GSM6176141_RV_CAV_3/KA_7_output_filtered.h5'))
# srt7$Sample <- 'CAV_3_RV'
# srt8 <- CreateSeuratObject(ReadCB_h5('~/Downloads/matrix_public/GSM6176142_LV_CAV_3/KA_8_output_filtered.h5'))
# srt8$Sample <- 'CAV_3_LV'
# srt9 <- CreateSeuratObject(ReadCB_h5('~/Downloads/matrix_public/GSM6176143_RV_CAV_4/KA_9_output_filtered.h5'))
# srt9$Sample <- 'CAV_4_RV'
# srt10 <- CreateSeuratObject(ReadCB_h5('~/Downloads/matrix_public/GSM6176144_LV_CAV_1/KA_10_output_filtered.h5'))
# srt10$Sample <- 'CAV_1_LV'
ref4.srt <- merge(srt1, list(srt2, srt3, srt4, srt5, srt6, srt7, srt8, srt9, srt10))
ref4.srt <- ProcessSrt_std(ref4.srt, do.umap = F) 
VlnPlot2(ref4.srt, features = c('nCount_RNA', 'nFeature_RNA', 'pct_mito_RNA'), group.by = 'Sample')
ref4.srt <- ref4.srt[, ref4.srt$nCount_RNA < 20000 &
                             ref4.srt$nCount_RNA > 200 &
                             ref4.srt$nFeature_RNA < 5000 &
                             ref4.srt$nFeature_RNA > 100 &
                             ref4.srt$pct_mito_RNA < 1]
ref4.srt <- ProcessSrt_std(ref4.srt, do.umap = F) |> ProcessSrt_hmn(haromize.by = 'Sample', do.umap = T)
ref4.srt <- FindNeighbors(ref4.srt, reduction = 'harmony') |> FindClusters(resolution = 0.5)
FeaturePlot2(ref4.srt, features = markers_lvl2_immune) + DimPlot2(ref4.srt, label = T, cols = mycol_20)
ref4.srt <- DietSeurat(ref4.srt, dimreducs = names(ref4.srt@reductions))
saveRDS(ref4.srt, 'individual/2023_CircHeartFail_JMoslehi.full.srt')
ref4_sub.srt <- ref4.srt[, ref4.srt$RNA_snn_res.0.5 == 1] ## NK + TC
ref4_sub.srt <- ProcessSrt_std(ref4_sub.srt, do.umap = F) |> ProcessSrt_hmn(haromize.by = 'Sample', do.umap = T)
ref4_sub.srt <- FindNeighbors(ref4_sub.srt, reduction = 'harmony') |> FindClusters(resolution = 0.5)
FeaturePlot2(ref4_sub.srt, features = c('GZMB', 'KLRF1', 'PRF1', 'GNLY', 'KLRC4', 'FCRL3', 'KLRC2')) + 
                     DimPlot2(ref4_sub.srt, label = T, cols = mycol_10)
ref4_sub.srt <- DietSeurat(ref4_sub.srt, dimreducs = names(ref4_sub.srt@reductions))
saveRDS(ref4_sub.srt, 'individual/2023_CircHeartFail_JMoslehi.TC_NK_subset.srt')
ref4_nk.srt <- ref4_sub.srt[, ref4_sub.srt$RNA_snn_res.0.5 == 2]
ref4_nk.srt$Study <- 'EMB'
ref4_nk.srt$Cell_state <- 'NK'
ref4_nk.srt$Group <- paste0('EMB_', str_split(ref4_nk.srt$Sample, pattern = '_', simplify = T)[, 1])
ref4_nk.srt$Donor <- paste0('EMB_', str_remove_all(ref4_nk.srt$Sample, pattern = '_LV|_RV'))
ref4_nk.srt <- DietSeurat(ref4_nk.srt, dimreducs = names(ref4_nk.srt@reductions))
saveRDS(ref4_nk.srt, 'individual/2023_CircHeartFail_JMoslehi.NK_subset.srt')

#### Pediatric CHD HF data  ####
load('~/Downloads/ImmuneCells_Pediatric_snRNA_GEO')
ref5.srt <- ImmuN4
ref5_sub.srt <- ref5.srt[, ref5.srt@active.ident == 'Tcells']
ref5_sub.srt <- ProcessSrt_std(ref5_sub.srt, do.umap = F) |> ProcessSrt_hmn(haromize.by = 'orig.ident', do.umap = T)
ref5_sub.srt <- FindNeighbors(ref5_sub.srt, reduction = 'harmony') |> FindClusters(resolution = 0.5)
FeaturePlot2(ref5_sub.srt, features = c('FCGR3A', 'NCAM1', 'GZMB', 'KLRF1', 'PRF1', 'GNLY', 'KLRC4', 'FCRL3', 'KLRC2')) + 
        DimPlot2(ref5_sub.srt, label = T, cols = mycol_10, reduction = 'hmn_umap')
ref5_nk.srt <- ref5_sub.srt[, ref5_sub.srt$RNA_snn_res.0.5 == 4]
Table(ref5_nk.srt$diagnosis)
ref5_nk.srt$Study <- 'CHD'
ref5_nk.srt$Cell_state <- 'NK'
ref5_nk.srt$Group <- 'CHD_HF'
ref5_nk.srt$Group[ref5_nk.srt$diagnosis == 'Donor'] <- 'CHD_Control'
ref5_nk.srt$Donor <- paste0(ref5_nk.srt$Group, '_', ref5_nk.srt$labID)
ref5_nk.srt <- ref5_nk.srt[, ref5_nk.srt$labID %in% c('P8',
                                                      'P26',
                                                      'P28',
                                                      'P33',
                                                      'P40',
                                                      'P64')]

#### Pediatric HLHS HF data  ####
ref6.srt <- readRDS('../../../2022_hlhs_dturaga/rdata/human_v2/integrated/PART19.annotated.srt.rds')
ref6.srt <- RenameAssays(ref6.srt, 'CBN' = 'RNA')
ref6.srt <- ref6.srt[, ref6.srt$Cell_state == 'NK' & ref6.srt$palliation %in% c('II', 'III')]
ref6.srt$Study <- 'HLHS'
ref6.srt$Group <- 'HLHS_HF'
ref6.srt$Donor <- paste0(ref6.srt$Group, '_', ref6.srt$donor)


#### Define NK3 Signature ####
nk_markers <- FindAllMarkers(lym.srt, assay = 'CBN', only.pos = T)
nk_markers <- nk_markers[nk_markers$p_val_adj < 0.001, ]
nk_markers_hi <- nk_markers[nk_markers$avg_log2FC > 1, ]
gl <- split(nk_markers_hi$gene, nk_markers_hi$cluster)
nk3_mk <- gl$NK3[! gl$NK3 %in% unlist(gl[1:4])]
nk3_mk <- RemoveRiboMito(nk3_mk, human_or_mouse = 'human') 
L(nk3_mk) ## 146


#### Score NK3 signatures ####
combine.srt <- merge(x = que.srt,
                     y = list(ref1.srt,
                              ref2.srt[, ref2.srt$annotation %in% c('NK', 'NK_T')], 
                              ref3.srt[, ref3.srt$cell_states %in% c('NK_CD16hi', 'NK_CD16hiIFNGhi', 'NK_CD56hi')],
                              ref4_nk.srt,
                              ref5_nk.srt,
                              ref6.srt), 
                     merge.dr = 'final_umap')
combine.srt$Group <- factor(combine.srt$Group, levels = c(
        'HTx_Control',
        'HTx_PT',
        'CHD_HF',
        'HLHS_HF',
        'MI_Ischemic',
        'MI_Fibrotic',
        'MI_Myogenic',
        'CM_Control',
        'CM_DCM',
        'CM_ACM',
        'CM_NCCM',
        'DCM_Control',
        'DCM',
        'EMB_Control',
        'EMB_CAV'
))
combine.srt <- AddModuleScore2(combine.srt, features = list(nk3_mk), names = 'NK3_sig', return_z = T)
combine.srt$NK3_sig[combine.srt$NK3_sig < 0] <- 0
data <- combine.srt@meta.data[, c('NK3_sig', 'Group', 'Donor')] |> group_by(Donor) |> mutate(Mean_NK3 = mean(NK3_sig))
data <- U(data[, c('Group', 'Donor', 'Mean_NK3')])
# data <- data[! data$Group %in% c('CM_Control'), ]
p6 <- ggplot(data, aes(x = Group, y = Mean_NK3)) +
        geom_boxplot(outlier.shape = NA) +
        geom_point() +
        scale_y_continuous(limits = c(0, 1)) +
        theme_classic() +
        RotatedAxis()
p6
PlotPDF('6.1.box.nk3_signature_expression_across_diseases', 8, 4)
p6
dev.off()

gene1 <- read.table('external/HALLMARK_ALLOGRAFT_REJECTION.v7.5.1.grp', comment.char = '#', header = T)[,1]
tmp.srt <- nk_dcm.srt[, nk_dcm.srt$Study == 'Lavine']
tmp.srt$name2 <- factor(tmp.srt$name2, levels = c('Donor', 'DCM'))
tmp.srt <- AddModuleScore2(tmp.srt, features = list(gene), names = 'ALLO', return_z = T)
p16.1 <- BoxPlot(tmp.srt, feature = 'ALLO', group.by = 'name2') +
        ggtitle(label = t.test(tmp.srt$ALLO[tmp.srt$name2 == 'DCM'],
                               tmp.srt$ALLO[tmp.srt$name2 == 'Donor'])$p.value)
tmp.srt <- lym.srt
tmp.srt$group <- NA
tmp.srt$group[tmp.srt$name != 'Control' & tmp.srt$Cell_state == 'NK3'] <- 'NK3'
tmp.srt$group[tmp.srt$name2 == 'Control' & tmp.srt$Cell_state %in% c('NK1', 'NK2', 'NK3')] <- 'Control'
tmp.srt <- tmp.srt[, !is.na(tmp.srt$group)]
tmp.srt <- AddModuleScore2(tmp.srt, features = list(gene), names = 'ALLO', return_z = T)
p16.2 <- BoxPlot(tmp.srt, feature = 'ALLO', group.by = 'group') +
        ggtitle(label = t.test(tmp.srt$ALLO[tmp.srt$group == 'NK3'],
                               tmp.srt$ALLO[tmp.srt$group == 'Control'])$p.value)
tmp.srt <- ref4_nk.srt
tmp.srt$Group <- factor(tmp.srt$Group, levels = c('EMB_Control', 'EMB_CAV'))
tmp.srt <- AddModuleScore2(tmp.srt, features = list(gene), names = 'ALLO', return_z = T)
p16.3 <- BoxPlot(tmp.srt, feature = 'ALLO', group.by = 'Group') +
        ggtitle(label = t.test(tmp.srt$ALLO[tmp.srt$Group == 'EMB_CAV'],
                               tmp.srt$ALLO[tmp.srt$Group == 'EMB_Control'])$p.value)
p16 <- p16.1 + p16.2 + p16.3 & theme(aspect.ratio = 4) &  scale_y_continuous(limits = c(-3, 3), breaks = seq(-3, 3))
p16
PlotPDF('6.2.box.nk_rejection_score_vs_dcm_vs_emb', 8, 4)
p16
dev.off()

####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~----


####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
####  Testing   ####
####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
DefaultAssay(dcm.srt) <- 'RNA'
dcm.srt <- NormalizeData(dcm.srt)
x <- merge(full.srt, dcm.srt)
x$name2[is.na(x$name2)] <- x$Condition[is.na(x$name2)]
DotPlot2(x, features = 'VAV1', group.by = 'name2')
VlnPlot2(x, features = 'VAV1', group.by = 'name2', pt.size=0.1)


lym.srt <- AddModuleScore2(lym.srt, features = list(NK_cytotoxicity, 
                                                          NK_HLA_depend_act_receptors, 
                                                          NK_HLA_depend_inh_receptors,
                                                          NK_HLA_independ_act_receptors,
                                                          NK_HLA_independ_inh_receptors),
                              names = c('NK_cytotoxicity',
                                        'NK_HLA_depend_act_receptors',
                                        'NK_HLA_depend_inh_receptors',
                                        'NK_HLA_independ_act_receptors',
                                        'NK_HLA_independ_inh_receptors'), return_z = F)
VlnPlot2(lym.srt, features = c('NK_cytotoxicity',
                                  'NK_HLA_depend_act_receptors',
                                  'NK_HLA_depend_inh_receptors',
                                  'NK_HLA_independ_act_receptors',
                                  'NK_HLA_independ_inh_receptors'), group.by = 'Cell_state')




