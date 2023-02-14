####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
####  Pediatric Acute Rejection
####  2022-08-03 by Xiao LI (Texas Heart Institute, US)
####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
####  Initiate  ####
####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~----
Ver <- '1'
Step <- 'STEP25_Analysis_DCM_Data'
Project <- '2022_acute_rejection_dturaga'
Code_dir <- paste0('~/Documents/Bioinformatics/project/', Project, '/code/human_v', Ver, '/')

source(Sys.readlink(paste0(Code_dir, 'src/bioinformatics.R')))
source(Sys.readlink(paste0(Code_dir, 'src/scRNAseq.R')))
source(Sys.readlink(paste0(Code_dir, 'src/geneset_human.R')))
source(paste0(Code_dir, 'human_v', Ver, '.helper_functions.R'))

InitiateProject('Rivendell', Ver, Step, 'human', Project, 'ithil')
####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
####  Load DCM Data -- Do not re-run  ####
####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~----
# load('/Volumes/shire/data/scrnaseq/2022_NatCardRes_KLavine/matrix_public/GSE183852_DCM_Cells.Robj')
# dcm.srt <- HDCM
#
# DimPlot2(dcm.srt)
# Table(dcm.srt$Names, dcm.srt$Condition)
####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# saveRDS(dcm.srt, 'external/2022_NatCardRes_KLavine.dcm_donor_scrna.srt.rds')
dcm.srt <- readRDS('external/2022_NatCardRes_KLavine.dcm_donor_scrna.srt.rds')
####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
####  Load Allograft Data  ####
####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~----
full.srt <- readRDS('integrated/STEP18.annotated.srt.rds')
full_imm.srt <- full.srt[, full.srt$Cell_type %in% c('Lym', 'Mye') &
                                 full.srt$Cell_type_non_ambig &
                                 full.srt$Fmux_genotype_clean != 'Ambiguous']
####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
####  Global comparative PCA  ####
####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~----
DefaultAssay(full.srt) <- 'CBN'
tmp.srt <- RenameAssays(full.srt, 'CBN' = 'RNA')
tmp.srt <- NormalizeData(tmp.srt)
tmp.srt@meta.data <- tmp.srt@meta.data[, c('name', 'name2', 'age_recipient', 'age_donor', 'sex_recipient', 'sex_donor')]
tmp.srt <- FindVariableFeatures(tmp.srt)

Table(tmp.srt$name)
Table(tmp.srt$name2)

psb <- AggregateExpression(tmp.srt, assays = 'RNA', slot = 'count', group.by = 'name')$RNA
meta <- tmp.srt@meta.data[!duplicated(tmp.srt$name),]
rownames(meta) <- meta$name
meta <- meta[colnames(psb), ]
dds <- DESeq2::DESeqDataSetFromMatrix(psb+1,
                                      colData = meta,
                                      design = ~ name)
rld <- DESeq2::rlog(dds, blind = T)
p <- DESeq2::plotPCA(rld, intgroup = "name", ntop = 500) +
        scale_color_manual(values = c(mycol_40[13:16], mycol_10[1:3])) +
        theme_classic() +
        theme(aspect.ratio = 1)
p
PlotPDF('0.1.pc.all_cell_global_pc_by_donor', 5, 5)
p
dev.off()
####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
saveRDS(p$data, 'analysis/STEP25.allograft_ctrls_all_cell_pseudobulk_pca.df.rds')
####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
####  Immune comparative analysis  ####
####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~----
DefaultAssay(full_imm.srt) <- 'CBN'
tmp.srt <- RenameAssays(full_imm.srt, 'CBN' = 'RNA')
tmp.srt <- NormalizeData(tmp.srt)
tmp.srt$Study <- 'Current'
tmp.srt@meta.data <- tmp.srt@meta.data[, c('Study', 'name', 'name2', 'age_recipient', 'Cell_state')]
tmp.srt <- FindVariableFeatures(tmp.srt)

dcm_imm.srt <- dcm.srt[, dcm.srt$Names %in%  c('NK_Cells', 'T_Cells', 'Monocytes', 'Macrophages')]
DefaultAssay(dcm_imm.srt) <- 'RNA'
dcm_imm.srt <- NormalizeData(dcm_imm.srt)
dcm_imm.srt$Study <- 'Lavine'
dcm_imm.srt@meta.data <- dcm_imm.srt@meta.data[, c('Study', 'orig.ident', 'Condition', 'Age_Group_Tertile', 'Names')]
colnames(dcm_imm.srt@meta.data) <- c('Study', 'name', 'name2', 'age_recipient', 'Cell_state')
dcm_imm.srt$Cell_state <- revalue(dcm_imm.srt$Cell_state, c('NK_Cells' = 'NK',
                                                            'T_Cells' = 'TC',
                                                            'Monocytes' = 'Mono',
                                                            'Macrophages' = 'Mac'))
dcm_imm.srt$name <- revalue(dcm_imm.srt$name, c('HDCM5' = 'Donor1',
                                                'HDCM7' = 'Donor2'))
dcm_imm.srt <- FindVariableFeatures(dcm_imm.srt)

imm_merge <- merge(DietSeurat(tmp.srt, assays = 'RNA'),
                   list(DietSeurat(dcm_imm.srt, assays = 'RNA')))
gene_keep <- rownames(imm_merge)[rowMaxs(imm_merge@assays$RNA@counts)>1]
imm_merge <- DietSeurat(imm_merge, features = gene_keep)
imm_merge$group <- paste0(imm_merge$Study, '_', imm_merge$name)
Table(imm_merge$group)
Table(imm_merge$name)
Table(imm_merge$name2)

psb <- AggregateExpression(imm_merge, assays = 'RNA', slot = 'count', group.by = 'name')$RNA
meta <- imm_merge@meta.data[!duplicated(imm_merge$name),]
rownames(meta) <- meta$name
meta <- meta[colnames(psb), ]
dds <- DESeq2::DESeqDataSetFromMatrix(psb+1,
                                      colData = meta,
                                      design = ~ name)
rld <- DESeq2::rlog(dds, blind = T)
p <- DESeq2::plotPCA(rld, intgroup = "name", ntop = 500) +
        scale_color_manual(values = c(mycol_10[1:3], mycol_40[13:16], mycol_40[25:26], mycol_40[29:33])) +
        theme_classic() +
        theme(aspect.ratio = 1)
p
PlotPDF('1.1.pc.all_immune_global_pc_by_donor', 5, 5)
p
dev.off()

imm_merge$name3 <- revalue(imm_merge$name2, c(Donor = 'Control'))
psb <- AggregateExpression(imm_merge, assays = 'RNA', slot = 'count', group.by = 'name3')$RNA
meta <- imm_merge@meta.data[!duplicated(imm_merge$name3),]
rownames(meta) <- meta$name3
meta <- meta[colnames(psb), ]
dds <- DESeq2::DESeqDataSetFromMatrix(psb+1,
                                      colData = meta,
                                      design = ~ name3)
rld <- DESeq2::rlog(dds, blind = T)
p <- DESeq2::plotPCA(rld, intgroup = "name3", ntop = 500) +
        scale_color_manual(values = mycol_10) +
        theme_classic() +
        theme(aspect.ratio = 1)
p
PlotPDF('1.2.pc.all_immune_global_pc_by_group', 5, 5)
p
dev.off()
####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
####  MP1 analysis  ####
####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~----
mac_merge <- imm_merge[, imm_merge$Cell_state %in% c('Mac', 'MP1', 'MP2') &
                               !imm_merge$name %in% c('12y PT', '15m PT')]
Table(mac_merge$name, mac_merge$Cell_state)
mac_merge <- mac_merge %>%
        NormalizeData() %>%
        CellCycleScoring(s.features = s.genes,
                         g2m.features = g2m.genes,
                         set.ident = F) %>%
        PercentageFeatureSet(pattern = '^MT-', col.name = 'pct_mito_RNA', assay = 'RNA')

deg <- full.srt@misc$marker$Cell_state_marker
deg_hi <- deg[deg$avg_log2FC>0.5 & deg$p_val_adj<0.001 & deg$cluster %in% c('MP1', 'Mac2', 'Mono'), 'gene']
L(deg_hi) ## 2489

mac_merge <- ScaleData(mac_merge, features = deg_hi, vars.to.regress = c('nFeature_RNA',
                                                                         'pct_mito_RNA',
                                                                         'nCount_RNA',
                                                                         'S.Score',
                                                                         'G2M.Score'))
mac_merge <- RunPCA(mac_merge, features = deg_hi)
mac_merge$name4 <- revalue(mac_merge$name, c('Donor1' = 'Control 5',
                                             'Donor2' = 'Control 6',
                                             'HDCM1' = 'DCM 1',
                                             'HDCM3' = 'DCM 2',
                                             'HDCM4' = 'DCM 3',
                                             'HDCM6' = 'DCM 4',
                                             'HDCM8' = 'DCM 5',
                                             '5d PT' = '5d PT MP1'
))
mac_merge$name4 <- factor(mac_merge$name4,
                          levels = c(paste('Control', 1:6), paste('DCM', 1:5), paste0('5d PT MP', 1:2)))
mac_merge$name4[mac_merge$name == '5d PT'] <- paste('5d PT', mac_merge$Cell_state[mac_merge$name == '5d PT'])
Table(mac_merge$name4)

DimPlot2(mac_merge, group.by = 'name4', cols = mycol_20, reduction = 'pca')

mac_merge <- RunHarmony(mac_merge, group.by.vars = c('Study'))
DimPlot2(mac_merge, group.by = 'name4', reduction = 'harmony', cols = mycol_20)

## Plot PCA
df <- as.data.frame(mac_merge@reductions$pca@cell.embeddings[, 1:2])
df$Group <- mac_merge$name4
df$Group2 <- as.vector(df$Group)
df$Group2[df$Group %in% paste('Control', 1:6)] <- 'Control MP'
df$Group2[df$Group %in% paste('DCM', 1:5)] <- 'DCM MP'

df2 <- df |>
        group_by(Group) |>
        summarize(PC_1 = mean(PC_1), PC_2 = mean(PC_2))
df2$Group2 <- as.vector(df2$Group)
df2$Group2[df2$Group %in% paste('Control', 1:6)] <- 'Control MP'
df2$Group2[df2$Group %in% paste('DCM', 1:5)] <- 'DCM MP'

p <- ggplot() +
        geom_point(data = df, mapping = aes(x = PC_1, y = PC_2, color = Group2), alpha = 0.05, size = 0.5) +
        geom_point(data = df2, mapping = aes(x = PC_1, y = PC_2, color = Group2), size = 3) +
        scale_color_manual(values = mycol_10) +
        scale_x_continuous(limits = c(-30, 30)) +
        theme_classic() +
        theme(aspect.ratio = 1)
p
PlotPDF('2.1.pca.mac_group', 5, 5)
p
dev.off()

## Plot Harmonized PCA
df <- as.data.frame(mac_merge@reductions$harmony@cell.embeddings[, 1:2])
df$Group <- mac_merge$name4
df$Group2 <- as.vector(df$Group)
df$Group2[df$Group %in% paste('Control', 1:6)] <- 'Control MP'
df$Group2[df$Group %in% paste('DCM', 1:5)] <- 'DCM MP'

df2 <- df |>
        group_by(Group) |>
        summarize(PC_1 = mean(harmony_1), PC_2 = mean(harmony_2))
df2$Group2 <- as.vector(df2$Group)
df2$Group2[df2$Group %in% paste('Control', 1:6)] <- 'Control MP'
df2$Group2[df2$Group %in% paste('DCM', 1:5)] <- 'DCM MP'

p <- ggplot() +
        geom_point(data = df, mapping = aes(x = harmony_1, y = harmony_2, color = Group2), alpha = 0.05, size = 0.5) +
        geom_point(data = df2, mapping = aes(x = PC_1, y = PC_2, color = Group2), size = 3) +
        scale_color_manual(values = mycol_10) +
        scale_x_continuous(limits = c(-30, 30)) +
        theme_classic() +
        theme(aspect.ratio = 1)
p
PlotPDF('2.2.harmony_pca.mac_group', 5, 5)
p
dev.off()

mac_merge <- DietSeurat(mac_merge, scale.data = F, dimreducs = names(mac_merge@reductions))
####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
saveRDS(mac_merge, 'analysis/STEP25.allo_dcm_ctrl_mp_compare.srt.rds')
####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
####  NK comparative analysis  ####
####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~----
current_nk.srt <- full_imm.srt[, full_imm.srt$Cell_state %in% c('NK1', 'NK2', 'NK3') &
                                       full_imm.srt$name2 %in% c('Control', '5d PT')]
DefaultAssay(current_nk.srt) <- 'CBN'
current_nk.srt <- RenameAssays(current_nk.srt, 'CBN' = 'RNA')
current_nk.srt <- NormalizeData(current_nk.srt)
current_nk.srt$Study <- 'Current'
current_nk.srt@meta.data <- current_nk.srt@meta.data[, c('Study', 'name2', 'name')]

dcm_nk.srt <- dcm.srt[, dcm.srt$Names == 'NK_Cells']
DefaultAssay(dcm_nk.srt) <- 'RNA'
dcm_nk.srt <- NormalizeData(dcm_nk.srt)
dcm_nk.srt$Study <- 'Lavine'
dcm_nk.srt@meta.data <- dcm_nk.srt@meta.data[, c('Study', 'Condition', 'orig.ident')]
colnames(dcm_nk.srt@meta.data) <- c('Study', 'name2', 'name')

nk_merge <- merge(DietSeurat(current_nk.srt, assays = 'RNA'), DietSeurat(dcm_nk.srt, assays = 'RNA'))
gene_keep <- rownames(nk_merge)[rowMaxs(nk_merge@assays$RNA@counts)>1]
nk_merge <- DietSeurat(nk_merge, features = gene_keep)

nk_merge$Cell_state <- 'NK'
nk_merge$Cell_state[Cells(current_nk.srt)] <- as.vector(full_imm.srt$Cell_state[Cells(current_nk.srt)])
nk_merge <- nk_merge %>%
        NormalizeData() %>%
        CellCycleScoring(s.features = s.genes,
                         g2m.features = g2m.genes,
                         set.ident = F) %>%
        PercentageFeatureSet(pattern = '^MT-', col.name = 'pct_mito_RNA', assay = 'RNA')
Table(nk_merge$Cell_state, nk_merge$name)

deg <- full.srt@misc$marker$Cell_state_marker
deg_hi <- deg[deg$avg_log2FC>0.5 & deg$p_val_adj<0.001 & deg$cluster %in% c('NK1', 'NK2', 'NK3'), 'gene']
LU(deg_hi) ## 742
gene_use <- deg_hi

nk_merge <- ScaleData(nk_merge, features = gene_use, vars.to.regress = c('nFeature_RNA',
                                                                         'pct_mito_RNA',
                                                                         'nCount_RNA',
                                                                         'S.Score',
                                                                         'G2M.Score'))
nk_merge <- RunPCA(nk_merge, features = gene_use)
nk_merge$name4 <- revalue(nk_merge$name, c('HDCM5' = 'Control 5',
                                           'HDCM7' = 'Control 6',
                                           'HDCM1' = 'DCM 1',
                                           'HDCM3' = 'DCM 2',
                                           'HDCM4' = 'DCM 3',
                                           'HDCM6' = 'DCM 4',
                                           'HDCM8' = 'DCM 5',
                                           '5d PT' = '5d PT NK1'
))
nk_merge$name4 <- factor(nk_merge$name4,
                          levels = c(paste('Control', 1:6), paste('DCM', 1:5), paste0('5d PT NK', 1:3)))
nk_merge$name4[nk_merge$name == '5d PT'] <- paste('5d PT', nk_merge$Cell_state[nk_merge$name == '5d PT'])
Table(nk_merge$name4)

nk_merge <- RunHarmony(nk_merge, group.by.vars = c('Study'))
DimPlot2(nk_merge, group.by = 'name4', cols = mycol_20, reduction = 'pca') /
        DimPlot2(nk_merge, group.by = 'name4', reduction = 'harmony', cols = mycol_20)

## Plot PCA
df <- as.data.frame(nk_merge@reductions$pca@cell.embeddings[, 1:2])
df$Group <- nk_merge$name4
df$Group2 <- as.vector(df$Group)
df$Group2[df$Group %in% paste('Control', 1:6)] <- 'Control NK'
df$Group2[df$Group %in% paste('DCM', 1:5)] <- 'DCM NK'

df2 <- df |>
        group_by(Group) |>
        summarize(PC_1 = mean(PC_1), PC_2 = mean(PC_2))
df2$Group2 <- as.vector(df2$Group)
df2$Group2[df2$Group %in% paste('Control', 1:6)] <- 'Control NK'
df2$Group2[df2$Group %in% paste('DCM', 1:5)] <- 'DCM NK'

varExplainedPct <- (ExplainedVariance(nk_merge)) # Get % variance explained
p <- ggplot() +
        geom_point(data = df, mapping = aes(x = PC_1, y = PC_2, color = Group2), alpha = 0.2) +
        geom_point(data = df2, mapping = aes(x = PC_1, y = PC_2, color = Group2), size = 3) +
        scale_color_manual(values = mycol_10) +
        #scale_x_continuous(limits = c(-30, 30)) +
        labs(x = paste0('PC1: ', varExplainedPct[1], '% variance'),
             y = paste0('PC2: ', varExplainedPct[2], '% variance')) +
        theme_classic() +
        theme(aspect.ratio = 1)
p
PlotPDF('3.1.pca.nk_sc_by_group', 5, 5)
p
dev.off()

## Plot harmony PCA
df <- as.data.frame(nk_merge@reductions$harmony@cell.embeddings[, 1:2])
df$Group <- nk_merge$name4
df$Group2 <- as.vector(df$Group)
df$Group2[df$Group %in% paste('Control', 1:6)] <- 'Control NK'
df$Group2[df$Group %in% paste('DCM', 1:5)] <- 'DCM NK'

df2 <- df |>
        group_by(Group) |>
        summarize(PC_1 = mean(harmony_1), PC_2 = mean(harmony_2))
df2$Group2 <- as.vector(df2$Group)
df2$Group2[df2$Group %in% paste('Control', 1:6)] <- 'Control NK'
df2$Group2[df2$Group %in% paste('DCM', 1:5)] <- 'DCM NK'

varExplainedPct <- ExplainedVariance(nk_merge, reduction = 'harmony') # Get % variance explained
p <- ggplot() +
        geom_point(data = df, mapping = aes(x = harmony_1, y = harmony_2, color = Group2), alpha = 0.2) +
        geom_point(data = df2, mapping = aes(x = PC_1, y = PC_2, color = Group2), size = 3) +
        scale_color_manual(values = mycol_10) +
        #scale_x_continuous(limits = c(-30, 30)) +
        labs(x = paste0('PC1: ', varExplainedPct[1], '% variance'),
             y = paste0('PC2: ', varExplainedPct[2], '% variance')) +
        theme_classic() +
        theme(aspect.ratio = 1)
p
PlotPDF('3.2.harmony_pca.nk_sc_by_group', 5, 5)
p
dev.off()

nk_merge <- DietSeurat(nk_merge, scale.data = F, dimreducs = names(nk_merge@reductions))
####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
saveRDS(nk_merge, 'analysis/STEP25.allo_dcm_ctrl_nk_compare.srt.rds')
####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
####  Allograft NK representative genes  ####
####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~----
nk_no_donor <- nk_merge[, nk_merge$name2 %in% c('5d PT', 'DCM')]
DefaultAssay(nk_no_donor) <- 'RNA'
Idents(nk_no_donor) <- 'name2'
nk_no_donor <- NormalizeData(nk_no_donor)
nk_deg <- FindMarkers(nk_no_donor, ident.1 = "5d PT", ident.2 = 'DCM', test.use = 'MAST', assay = 'RNA')
nk_deg$gene <- rownames(nk_deg)
p <- MarkerVolcano(mk.df = nk_deg)
p
PlotPDF('4.1.vol.nk_deg_5dpt_vs_dcm', 10, 10)
p
dev.off()

nk_har_up <- rownames(nk_deg)[nk_deg$avg_log2FC >= 0.25 & nk_deg$p_val_adj<0.05]
nk_har_dn <- rownames(nk_deg)[nk_deg$avg_log2FC <= -0.25 & nk_deg$p_val_adj<0.05]

nk_mk <- full.srt@misc$marker$Cell_state_marker[grepl('^NK', full.srt@misc$marker$Cell_state_marker$cluster), ]
nk_mk <- nk_mk[nk_mk$avg_log2F >= 0.25 & nk_mk$p_val_adj < 0.05,]
nk3_mk <- nk_mk[nk_mk$cluster == 'NK3', 'gene']

grp1 <- read.table('external/KEGG_NATURAL_KILLER_CELL_MEDIATED_CYTOTOXICITY.v7.5.1.grp',
                   comment.char = '#', header = T)[,1]
grp2 <- read.table('external/KEGG_ALLOGRAFT_REJECTION.v7.5.1.grp',
                   comment.char = '#', header = T)[,1]
grp3 <- read.table('external/WP_ALLOGRAFT_REJECTION.v7.5.1.grp',
                   comment.char = '#', header = T)[,1]
grp4 <- read.table('external/KEGG_GRAFT_VERSUS_HOST_DISEASE.v7.5.1.grp',
                   comment.char = '#', header = T)[,1]
grp5 <- read.table('external/HALLMARK_ALLOGRAFT_REJECTION.v7.5.1.grp',
                   comment.char = '#', header = T)[,1]

nk3_gene <- intersect(nk_har_up, nk3_mk)
genes1 <- intersect(nk3_gene, grp1)
genes2 <- intersect(nk3_gene, grp2)
genes3 <- intersect(nk3_gene, grp3)
genes4 <- intersect(nk3_gene, grp4)
genes5 <- intersect(nk3_gene, grp5)

Idents(nk_merge) <- 'name4'
p1 <- DotPlot2(nk_merge, features = genes1, idents = c(paste0('DCM', 1:5), paste0('5d PT NK', 1:3)), col.min = 0)
p2 <- DotPlot2(nk_merge, features = genes4, idents = c(paste0('DCM', 1:5), paste0('5d PT NK', 1:3)), col.min = 0)
p3 <- DotPlot2(nk_merge, features = genes5, idents = c(paste0('DCM', 1:5), paste0('5d PT NK', 1:3)), col.min = 0)

PlotPDF('4.2.dot.nk_5dpt_dcm_deg_in_pathways', 20, 10)
wrap_plots(p1, p2, p3, ncol = 3)
dev.off()
####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
saveRDS(list(genes1, genes4, genes5), 'analysis/STEP25.nk_5dpt_vs_dcm_deg_select.rds')
####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
