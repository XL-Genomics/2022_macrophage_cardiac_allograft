####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
####  Pediatric Acute Rejection
####  2022-08-03 by Xiao LI (Texas Heart Institute, US)
####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~----


####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
####  Initiate  ####
####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~----
Ver <- '1'
Step <- 'STEP90_Plot'
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
####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
####  Edit previously saved objects  ####
####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~----
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
mp2_sub.srt$Group2 <- factor(mp2_sub.srt$Group, levels = c('ResMP', 'Late MoMP2'))
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
####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~----


####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
####  Fig 1  Cardiac Allografts Harbor Recipient-derived Immune Cells   ####
####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~----
####  Plot 1 dot sample gender and duration ####
df <- readRDS('analysis/STEP25.allograft_ctrls_all_cell_pseudobulk_pca.df.rds')
df$Age <- c(11, 11, 10, 14, 20, 8, 17) ## Age of the heart (donor age + duration for allograft)
df$Sex <- c('M', 'F', 'F', 'M', 'M', 'F', 'M') ## Sex of the heart (donor sex for allograft)
df$Count <- as.vector(Table(full.srt$name))
p1.1 <- ggplot(df) +
        geom_point(aes(x = PC1, y = PC2, color = group)) +
        scale_color_manual(values = Color_donor) +
        labs(color = '', x = '', y = '')
p1.2 <- ggplot(df) +
        geom_point(aes(x = PC1, y = PC2, color = Age)) +
        scale_color_distiller(palette = 'Spectral', breaks = seq(8, 20, 4)) +
        labs(x = '', y = '')
p1.3 <- ggplot(df) +
        geom_point(aes(x = PC1, y = PC2, color = Sex)) +
        scale_color_manual(values = mycol_10[2:1]) +
        labs(x = 'PC1 39% variance', y = 'PC2 25% variance')
p1.4 <- ggplot(df) +
        geom_point(aes(x = PC1, y = PC2, size = Count)) +
        scale_size(breaks = c(1e3, 1e4, 4e4, 7e4), limits = c(1e2, 7e4)) +
        labs(size = 'Nuclei count', x = '', y = '')
p1 <- wrap_plots(p1.1, p1.2, p1.3, p1.4, ncol = 2) &
        theme(aspect.ratio = 1,
              axis.line = element_line(size = 0.5, color = 'black'),
              axis.text = element_blank(),
              axis.ticks = element_blank(),
              legend.background = element_blank(),
              panel.background = element_blank())
p1
PlotPDF('01.1.pca.sample_metadata', 6, 6)
p1
dev.off()

####  Plot 2 umap color by condition ####
p2 <- DimPlot2(full.srt, group.by = 'name2', cols = Color_condition, raster = F, pt.size = 0.1, order = T) +
        labs(title = 'Condition', x = 'UMAP1', y = 'UMAP2')
p2

PlotPDF('01.2.umap.color_by_condition', 8, 8)
p2
dev.off()

####  Plot 3 umap color by cell type ####
p3 <- DimPlot2(full.srt, group.by = 'Cell_type_pub', label = F, raster = F, pt.size = 0.1,
               cols = Color_cell_type) +
        labs(title = 'Cell type', x = 'UMAP1', y = 'UMAP2')
p3

PlotPDF('01.3.umap.color_by_celltype', 8, 8)
p3
dev.off()

####  Plot 4 dot ref cell type score ####
ref.marker <- readRDS('/Volumes/shire/data/scrnaseq/2020_Nature_STeichmann/matrix_public/celltype_markers.rds')
ref.marker.high <-  ref.marker |> filter(pct.1 >= 0.3, p_val_adj <= 1e-8, avg_log2FC >= 10)
genes <- split(ref.marker.high$gene, ref.marker.high$cluster)
tmp.srt <- AddModuleScore2(full.srt, features = genes, name = c(
        'Adipocyte score',
        'Atrial CM score',
        'Endothelial score',
        'Fibroblast score',
        'Lymphoid score',
        'Mesothelial score',
        'Myeloid score',
        'Neuronal score',
        'NotAssigned score',
        'Pericyte score',
        'Smooth muscle score',
        'Ventricular CM score'
))
p4 <- DotPlot2(tmp.srt, group.by = 'Cell_type_pub', cols = mycol_BuGr[c(100, 1)], col.min = 0,
              features = c(
                      'Ventricular CM score',
                      'Fibroblast score',
                      'Endothelial score',
                      'Pericyte score',
                      'Smooth muscle score',
                      'Lymphoid score',
                      'Myeloid score',
                      'Neuronal score',
                      'Adipocyte score'
              )) +
        labs(x = 'Reference cell types', y = 'Annotation', fill = 'Mean Expr.', size = 'Cell expr. %')
p4
PlotPDF('01.4.dot.cell_type_annotation_to_ref', 7, 5)
p4
dev.off()

####  Plot 5 box cell type fold change vs each control ####
tmp.srt <- full.srt[, full.srt$Cell_type != 'CM']
tmp.srt$tmp <- factor('Non-immune', levels = c('Non-immune', 'Lymphoid', 'Myeloid'))
tmp.srt$tmp[tmp.srt$Cell_type == 'Lym'] <- 'Lymphoid'
tmp.srt$tmp[tmp.srt$Cell_type == 'Mye'] <- 'Myeloid'
df <- CountCellBarPlot(tmp.srt, group.var = 'tmp', stack.var = 'name', percentage = T)$data
df$Fraction <- df$Count/rep(Table(tmp.srt$name), each = L(levels(df$GroupVar)))
mtx1 <- matrix(NA, 4, LU(df$GroupVar))
rownames(mtx1) <- paste0('Control ', 1:4)
colnames(mtx1) <- U(df$GroupVar)
mtx2 <- mtx1
mtx3 <- mtx1
for(i in 1:4){
        for(j in 1:3){
                celltype <- levels(df$GroupVar)[j]
                mtx1[i, j] <- df$Fraction[df$StackVar == '5d PT' & df$GroupVar == celltype]/
                        df$Fraction[df$StackVar == paste('Control', i) & df$GroupVar == celltype]
                mtx2[i, j] <- df$Fraction[df$StackVar == '15m PT' & df$GroupVar == celltype]/
                        df$Fraction[df$StackVar == paste('Control', i) & df$GroupVar == celltype]
                mtx3[i, j] <- df$Fraction[df$StackVar == '12y PT' & df$GroupVar == celltype]/
                        df$Fraction[df$StackVar == paste('Control', i) & df$GroupVar == celltype]
        }
}
mtx_log <- log2(rbind(mtx1, mtx2, mtx3))
df2 <- melt(mtx_log)
colnames(df2) <- c('Patient', 'Cell.Type', 'Log2FC')
df2$Group <- rep(c('5d PT', '15m PT', '12y PT'), each = 4)
df2$Group <- factor(df2$Group, levels = c('5d PT', '15m PT', '12y PT'))
p5 <- ggplot(df2, aes(x = Group, y = Log2FC, group = Group)) +
        geom_boxplot(outlier.shape = NA, width = 0.5) +
        geom_point(aes(color = Patient), position = position_jitterdodge(jitter.width = 0.5)) +
        geom_hline(yintercept = c(0), color = 'red4', linetype = 'dashed') +
        scale_color_manual(values = Color_donor) +
        theme_classic() +
        scale_y_continuous(minor_breaks = seq(-3, 6.5, 0.5), breaks = seq(-3, 6.5, 1), limits = c(-3, 6.5)) +
        theme(aspect.ratio = 3, panel.grid.major.y = element_line()) +
        RotatedAxis() +
        labs(x = '', y = 'Log2 fold change (Allograft vs. Controls)',
             title = 'Cell composition change', color = '') +
        facet_wrap(~Cell.Type)
p5
PlotPDF('01.5.box.celltype_foldchange_har_vs_donors', 5, 5)
p5
dev.off()

####  Plot 6 umap color by graft vs host  ####
p6 <- DimPlot(full.srt, shuffle = T, reduction = 'clean_umap', group.by = 'Origin',
              order = c('Donor', 'Recipient', 'Ambiguous', 'Control cells'),
              cols = rev(Color_genotype2), raster = F, pt.size = 0.1) &
        theme(aspect.ratio = 1, axis.text = element_blank(), axis.ticks = element_blank()) &
        labs(title = 'Donor and Recipient cell origins', subtitle = '(infered by Demuxlet)', x = 'UMAP1', y = 'UMAP2')
p6
PlotPDF('01.6.umap.p136_all_cell_demuxlet_color_by_graft_host', 10, 10)
p6
dev.off()

####  Plot 7 violin plot showing Y chromosome - P182  ####
tmp.srt <- full.srt[, full.srt$name == '15m PT' & full.srt$Origin != 'Ambiguous']
tmp.srt <- AddModuleScore2(tmp.srt, features = list(chrY_genes$chrY_only_genes), assay = 'alra',
                           names = 'Score_chrY_only_genes', return_z = F)
p7 <- VlnPlot2(tmp.srt, features = 'Score_chrY_only_genes', group.by = 'Origin',
               pt.size = 0.1, adjust = 1.5, cols = Color_genotype2) +
        scale_y_continuous(limits = c(-0.3, 1.5))
p7
PlotPDF('01.7.vln.p136_all_cell_demuxlet_color_by_graft_host', 4, 4)
p7
dev.off()

####  Plot 8 bar plot showing consensus percentage - P136  ####
tmp.srt <- full.srt[, full.srt$name == '5d PT']
P136_souporcell <- read.table(
        paste0('/Volumes/shire/data/scrnaseq/2022_Unpub2_DTuraga/matrix/souporcell/',
               '2022_Unpub2_DTuraga_P136_s1/output/clusters.tsv'),
        header = T)
rownames(P136_souporcell) <- paste0('2022_Unpub2_DTuraga:P3_S001:', P136_souporcell$barcode)
tmp.srt$Origin_2 <- revalue(P136_souporcell[Cells(tmp.srt), 'assignment'],
                                      c('1' = 'Donor',
                                        '0' = 'Recipient',
                                        '1/0' = 'Ambiguous',
                                        '0/1' = 'Ambiguous'))
tmp.srt$Origin_2 <- factor(tmp.srt$Origin_2, levels = c('Donor', 'Recipient', 'Ambiguous'))
Table(tmp.srt$Origin_2, tmp.srt$Origin)
data <- data.frame(table(droplevels(tmp.srt$Origin_2), droplevels(tmp.srt$Origin)))
colnames(data) <- c('Souporcell', 'Freemuxlet', 'Count')
data$Aggrement_with_Souporcell <- factor(c('Agree',
                                           rep('Disagree', 3),
                                           'Agree',
                                           rep('Disagree', 3),
                                           'Agree'),
                                         levels = c('Agree', 'Disagree'))
data$Freemuxlet <- factor(data$Freemuxlet, levels = c('Donor', 'Recipient', 'Ambiguous'))
frac <- paste0(
        round(sum(data$Count[data$Aggrement_with_Souporcell == 'Agree'])*100/sum(data$Count), digits = 1), '%, ',
        round(sum(data$Count[data$Aggrement_with_Souporcell != 'Agree'])*100/sum(data$Count), digits = 1), '%')
p8 <- ggplot(data) +
        geom_bar(aes(x = Aggrement_with_Souporcell, y = Count, fill = Freemuxlet), stat = 'identity') +
        scale_fill_manual(values = Color_genotype2) +
        labs(subtitle =  frac) +
        theme_classic() +
        theme(aspect.ratio = 1) +
        RotatedAxis()
p8
PlotPDF('01.8.bar.souporcell_vs_freemuxlet', 4, 4)
p8
dev.off()

####  Plot 9 umap plot showing example alleles - P136   ####
tmp.srt <- full.srt[, full.srt$name == '5d PT']
barcodes <- read.table(gzfile(paste0('/Volumes/shire/data/scrnaseq/2022_Unpub2_DTuraga/matrix/',
                                     '2022_Unpub2_DTuraga_P136_s1/outs/filtered_feature_bc_matrix/barcodes.tsv.gz')),
                       header = F)
vcf <- read.table(paste0('/Volumes/shire/data/scrnaseq/2022_Unpub2_DTuraga/matrix/',
                         'vartrix/2022_Unpub2_DTuraga_P136_s1/outvariants_consensus.txt'))
snv_matrix <- readMM(paste0('/Volumes/shire/data/scrnaseq/2022_Unpub2_DTuraga/matrix/',
                            'vartrix/2022_Unpub2_DTuraga_P136_s1/output_consensus.mat'))
rownames(snv_matrix) <- vcf$V1
colnames(snv_matrix) <- barcodes$V1
snv_matrix_p136 <- snv_matrix[, tmp.srt$orig.name]
## NEAT1 chr11_65429459
tmp.srt$SNV_chr11_65429459 <- factor(snv_matrix_p136['chr11_65429459', ])
tmp.srt$SNV_chr11_65429459 <- revalue(tmp.srt$SNV_chr11_65429459, replace = c(
        '0' = 'No coverage', '1' = 'Ref/Ref', '2' = 'Alt/Alt', '3' = 'Ref/Alt'))
## MALAT1 chr11_65504360
tmp.srt$SNV_chr11_65504360 <- factor(snv_matrix_p136['chr11_65504360', ])
tmp.srt$SNV_chr11_65504360 <- revalue(tmp.srt$SNV_chr11_65504360, replace = c(
        '0' = 'No coverage', '1' = 'Ref/Ref', '2' = 'Alt/Alt', '3' = 'Ref/Alt'))
p9.1 <- DimPlot2(tmp.srt, reduction = 'clean_umap', group.by = 'SNV_chr11_65429459', cols = Color_genotype, order = T) +
        labs(title = 'NEAT1 chr11:65429459 T/A')
p9.2 <- DimPlot2(tmp.srt, reduction = 'clean_umap', group.by = 'SNV_chr11_65504360', cols = Color_genotype, order = T) +
        labs(title = 'MALAT1 chr11:65504360 C/T')
p9.1+p9.2
PlotPDF('01.9.umap.example_alleles', 10, 5)
p9.1+p9.2
dev.off()

####  Plot 10 bar plot cell origin composition   ####
tmp.srt <- full.srt[, full.srt$name2 != 'Control' & full.srt$Origin != 'Ambiguous']
tmp.srt$tmp <- factor('Non-immune', levels = c('Non-immune', 'Lymphoid', 'Myeloid'))
tmp.srt$tmp[tmp.srt$Cell_type == 'Lym'] <- 'Lymphoid'
tmp.srt$tmp[tmp.srt$Cell_type == 'Mye'] <- 'Myeloid'
tmp.srt$tmp2 <- paste(tmp.srt$name, tmp.srt$Origin, sep = ':')
df <- CountCellBarPlot(tmp.srt, group.var = 'tmp', stack.var = 'tmp2', percentage = T)$data
df$Group <- factor(str_split(df$StackVar, pattern = ':', simplify = T)[,1], levels = paste(c('5d', '15m', '12y'), 'PT'))
df$StackVar <- str_split(df$StackVar, pattern = ':', simplify = T)[,2]
df <- df[df$StackVar == 'Donor', ]
df$Fraction <- NA
for(i in 1:nrow(df)){
        df$Fraction[i] <- df$Count[i]/Table(tmp.srt$tmp, tmp.srt$name)[as.vector(df$GroupVar[i]),
                                                                       as.vector(df$Group[i])]
}
p10 <- ggplot(df, aes(fill = StackVar, y = Fraction, x = Group)) +
        geom_bar(position = "stack", stat = "identity") +
        labs(title = 'Freemuxlet genotyping', fill = '', x = 'Cell type', y = 'Cell fraction') +
        theme_minimal() +
        theme(
                panel.grid.major.x =  element_blank(),
                panel.grid.minor.x = element_blank()) +
        RotatedAxis() +
        facet_wrap(~GroupVar)
p10
PlotPDF('01.10.bar.genotype_composition', 5, 5)
p10
dev.off()

####  Plot 11 heatmap cell type markes   ####
mk <- full.srt@misc$marker$Cell_type_marker
mk_dedup <- mk |> group_by(cluster) |> top_n(n = 50, wt = avg_log2FC)
mk_dedup <- mk_dedup[! mk_dedup$gene %in% mk_dedup$gene[duplicated(mk_dedup$gene)], ]
p11 <- MarkerHeatmap(full.srt, marker.df = mk_dedup, disp.min = 0, top = 40, group.cols = Color_cell_type, raster = T)
p11
PlotPDF('01.11.heat.full_cell_type_marker', 15, 15)
p11
dev.off()

####  Table 2 heatmap cell type markes  ####
mk <- full.srt@misc$marker$Cell_type_marker
mk_dedup <- mk |> group_by(cluster) |> top_n(n = 50, wt = avg_log2FC)
mk_dedup <- mk_dedup[! mk_dedup$gene %in% mk_dedup$gene[duplicated(mk_dedup$gene)], ]
mk_dedup <- mk_dedup[, c('gene', 'cluster', 'avg_log2FC', 'p_val_adj')]
colnames(mk_dedup) <- c('Gene_symbol', 'Cell_type', 'Mean_log2_fold_change', 'Adjusted_p_value')
mk_dedup$Gene_symbol <- paste0("'", mk_dedup$Gene_symbol, "'")
WriteCSV(mk_dedup, 'STEP90.01.1.table.full_cell_type_marker')
####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~----


####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
####  Fig 2  Graft MP is tissue-resident MP ####
####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~----
####  Plot 1 umap myeloid color by cell state ####
p1 <- DimPlot2(mye.srt, reduction = 'sub_umap', raster = F, cols = Color_cell_state_mye) +
        labs(x = 'UMAP1', y = "UMAP2", title = 'Myeloid cell subtypes')
p1
PlotPDF('02.1.umap.myeloid_color_by_cell_subtype', 5, 4)
p1
dev.off()

####  Plot 2 heat myeloid cell subtype markers  ####
mk <- FindAllMarkers(mye.srt, logfc.threshold = 0.5, only.pos = T, return.thresh = 0.001, test.use = 'MAST')
mk <- mk[! mk$gene %in% c('NAMPT','RIPOR2', "SIPA1L1", "ATP1B3", 'LRRK2'),]
p2 <- MarkerHeatmap(mye.srt, marker.df = mk, group.cols = Color_cell_state_mye,
                    n_cells = 300, disp.min = 0, top = 50, raster = T)
p2
PlotPDF('02.2.heat.immune_cell_subtype_markers', 10, 10)
print(p2)
dev.off()

####  Table 3 heat myeloid cell subtype markers  ####
mk <- FindAllMarkers(mye.srt, logfc.threshold = 0.5, only.pos = T, return.thresh = 0.001, test.use = 'MAST')
mk <- mk[! mk$gene %in% c('NAMPT','RIPOR2', "SIPA1L1", "ATP1B3", 'LRRK2'),]
mk <- mk[, c('gene', 'cluster', 'avg_log2FC', 'p_val_adj')]
colnames(mk) <- c('Gene_symbol', 'Cell_type', 'Mean_log2_fold_change', 'Adjusted_p_value')
mk$Gene_symbol <- paste0("'", mk$Gene_symbol, "'")
WriteCSV(mk, 'STEP90.02.1.table.myeloid_cell_subtype_marker')

####  Plot 3 box myeloid expressing Bajpai/Dick mac markers  ####
TRMP_gene <- c('CD163L1', 'MRC1', 'MAF', 'SIGLEC1', 'CCL8', 'CCL14', 'LILRB5', 'LYVE1', 'IL2RA', 'PDGFC', 'WLS', 'DAB2',
               'NRP1', 'SCN9A', 'FGF13', 'GDF15', 'IGF1', 'FMOD', 'SLIT3', 'EGFL7', 'ECM1', 'SDC3')
MoMP_gene <- c('CCL17', 'CCR5', 'CCR2', 'CXCL9', 'CXCL3', 'CXCL2', 'CXCL10', 'CSF2RA', 'CCL5', 'CCL20', 'TREM1', 'TET2',
               'NFKBIE', 'IL1R1', 'IL1R2', 'NFKB1', 'REL', 'MAP2K3', 'IL1A', 'NLRP3', 'TRAF1', 'MAPK6', 'SRC', 'IL1B',
               'NOD2', 'JAK3', 'IRAK2', 'MAP3K8', 'RELB', 'SOCS3', 'MYD88', 'MMP9', 'IL20', 'AREG', 'PTX3', 'IL23A',
               'IL10', 'OSM', 'TIMP1', 'EREG', 'IL27') ## These genes are from the K. Lavine's Nature Medicine Paper
genes <- read_excel('external/sciimmunol.abf7777_tables_s1_to_s9/sciimmunol.abf7777_table_s7.xlsx',
                    sheet = 3, col_names = T, skip = 1) ## These genes are from the E. Slava's Science Immun Paper
genes <- genes[genes$avg_logFC > 0.75 & genes$p_val_adj < 0.001, ]
gl <- split(genes$gene, genes$cluster)

tmp.srt <- AddModuleScore2(mye.srt, features = c(list(TRMP_gene, MoMP_gene), gl),
                           names = c('TRMP', 'MoMP', names(gl)), return_z = T)
p3.1 <- BoxPlot(tmp.srt, feature = c('MoMP'), group.by = 'Cell_state', cols = Color_cell_state_mye) + NoLegend()
p3.2 <- BoxPlot(tmp.srt, feature = c('TRMP'), group.by = 'Cell_state', cols = Color_cell_state_mye) + NoLegend()
p3.3 <- BoxPlot(tmp.srt, feature = c(names(gl)[2]), group.by = 'Cell_state', cols = Color_cell_state_mye)
p3.1 + p3.2 + p3.3
PlotPDF('02.3.box.mp2_marker_exp', 10, 3)
p3.1 + p3.2 + p3.3
dev.off()

####  Plot 4 feat density TLF+CCR2+MKI67  ####
p4 <- FeaturePlot3(mye.srt,
                   features = c('TIMD4', 'LYVE1', 'FOLR2', 'CCR2', 'FCGR3A', 'MKI67', 'TOP2A'),
                   adjust = 3, pt.size = 1, ncol = 7)
p4
PlotPDF('02.4.feat.resmp_marker', 22, 3)
p4
dev.off()

####  Plot 5 umap P136 myeloid color by cell origin ####
p5 <- DimPlot2(mye.srt, reduction = 'sub_umap', raster = F, split.by = 'name2',
               cols = Color_genotype2[-c(3)], group.by = 'Origin') +
        labs(x = 'UMAP1', y = "UMAP2", title = 'Myeloid cell subtypes')
p5
PlotPDF('02.5.umap.myeloid_color_by_cell_origin', 10, 4)
p5
dev.off()

####  Plot 6 bar P136 mye cell origin composition  ####
tmp.srt <- full.srt[, full.srt$name2 != 'Control' & full.srt$Cell_type == 'Mye' & full.srt$Origin != 'Ambiguous']
tmp.srt$Origin <- droplevels(tmp.srt$Origin)
p6.1 <- CountCellBarPlot(tmp.srt[, tmp.srt$name == '5d PT'], group.var = 'Cell_state', stack.var = 'Origin',
                         stack.color = Color_genotype2)$plot
p6.2 <- CountCellBarPlot(tmp.srt[, tmp.srt$name == '15m PT'], group.var = 'Cell_state', stack.var = 'Origin',
                         stack.color = Color_genotype2)$plot
p6.3 <- CountCellBarPlot(tmp.srt[, tmp.srt$name == '12y PT'], group.var = 'Cell_state', stack.var = 'Origin',
                         stack.color = Color_genotype2)$plot
p6 <- p6.1 + p6.2 + p6.3 &
        labs(title = 'Freemuxlet genotyping', fill = '', xlab = '', ylab = 'Cell count')
p6
PlotPDF('02.6.bar.p136_mye_genotype_composition', 9, 3)
p6
dev.off()

####  Plot 7 box P136 myeloid expressing Bajpai/Dick mac markers split by origin  ####
TRMP_gene <- c('CD163L1', 'MRC1', 'MAF', 'SIGLEC1', 'CCL8', 'CCL14', 'LILRB5', 'LYVE1', 'IL2RA', 'PDGFC', 'WLS', 'DAB2',
               'NRP1', 'SCN9A', 'FGF13', 'GDF15', 'IGF1', 'FMOD', 'SLIT3', 'EGFL7', 'ECM1', 'SDC3')
genes <- read_excel('external/sciimmunol.abf7777_tables_s1_to_s9/sciimmunol.abf7777_table_s7.xlsx',
                    sheet = 3, col_names = T, skip = 1) ## These genes are from the E. Slava's Science Immun Paper
genes <- genes[genes$avg_logFC > 0.75 & genes$p_val_adj < 0.001, ]
gl <- split(genes$gene, genes$cluster)
tmp.srt <- AddModuleScore2(momp_p136.srt, features = c(list(TRMP_gene), gl),
                           names = c('TRMP', names(gl)), return_z = T)
tmp.srt <- tmp.srt[, tmp.srt$Cell_state %in% c('MP2', 'Mono')]
p7.1 <- BoxPlot(tmp.srt, feature = 'TRMP', group.by = 'Cell_state_orig', cols = Color_cell_state_orig_mye)
p7.2 <- BoxPlot(tmp.srt, feature = names(gl)[2], group.by = 'Cell_state_orig', cols = Color_cell_state_orig_mye)
p7.1 + p7.2
PlotPDF('02.7.box.mp2_marker_exp_in_p136_group_by_origin', 7, 3)
p7.1 + p7.2
dev.off()

####  Plot 8 box P136 MP2 MHC  ####
tmp.srt <- AddModuleScore2(mye.srt, features = mhc_genes,
                           return_z = T, assay = 'alra', names = c('Score_MHC1', 'Score_MHC2'))
tmp.srt <- tmp.srt[, tmp.srt$Cell_state != 'Mast']
p8.1 <- BoxPlot(tmp.srt, feature = 'Score_MHC1', group.by = 'Cell_state_orig',
                cols = c(Color_genotype[c(1,3,2)], Color_cell_state_mye[c(1,2)])) + NoLegend()
p8.2 <- BoxPlot(tmp.srt, feature = 'Score_MHC2', group.by = 'Cell_state_orig',
                cols = c(Color_genotype[c(1,3,2)], Color_cell_state_mye[c(1,2)]))
p8 <- p8.1 + p8.2 & scale_y_continuous(limits = c(-2, 3))
p8
PlotPDF('02.8.box.mp2_mhc_expression_by_origin', 7, 3)
p8
dev.off()

####  Plot 9 violin MRC1 MERTK expression  ####
tmp.srt <- mye.srt[, mye.srt$Cell_state_orig %in% c('Control MP2', 'Donor MP2', 'Recipient MP2', 'MP1')]
p9 <- VlnPlot2(tmp.srt, features = c('MRC1', 'MERTK'),
                cols = Color_genotype, same.y.lims = 6,
                group.by = 'Cell_state_orig', pt.size = -1, ncol = 2) + RestoreLegend()
p9
PlotPDF('02.9.vln.mrc1_mertk_exp', 8, 4)
p9
dev.off()

####  Plot 10 bar Mono MP composition  ####
tmp.srt <- full.srt[, full.srt$Cell_type == 'Mye' & full.srt$Origin != 'Ambiguous' & full.srt$Cell_state != 'Mast']
p10 <- CountCellBarPlot(tmp.srt, group.var = 'name2', stack.var = 'Cell_state',
                        percentage = T, stack.color = Color_cell_state_mye)$plot
p10
PlotPDF('02.10.bar.mono_mp_composition', 4, 4)
p10
dev.off()
####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~----


####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
####  Fig 3  Tissue-resident MP Replenished by Monocyte-derived MP2 in transplant ####
####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~----
####  Plot 1 velocity myeloid umap  ####
colors <- c('grey80', 'deeppink2', 'cyan3', Color_cell_state_mye[c(1, 2)])
names(colors) <- levels(velo_all.srt$Cell_state_orig)
cell.colors <- colors[velo_all.srt$Cell_state_orig[colnames(velo_all.srt)]]
names(cell.colors) <- colnames(velo_all.srt)
PlotPDF('03.1.umap.mac_mono_velocity', 8, 8)
show.velocity.on.embedding.cor(emb = Embeddings(velo_all.srt, reduction = "sub_umap"),
                               vel = Tool(velo_pgf.srt, slot = "RunVelocity"),
                               cc = velo_pgf.vlc$cc,
                               cell.colors = ac(cell.colors, alpha = 0.5),
                               cell.border.alpha = 0.04,
                               show.grid.flow = T,
                               cex = 1,
                               n = 1000,
                               grid.n = 15,
                               arrow.scale = 6,
                               arrow.lwd = 2,
                               fixed.arrow.length = F,
                               min.grid.cell.mass = 4,
                               do.par = F,
                               asp = 1)
dev.off()

####  Plot 2 velocity p136 myeloid diffusion   ####
PlotPDF('03.2.diff.pgf_mono_mac_slingshot_trajectory', 6, 4)
plot(reducedDims(momp_p136.sce)$DIFFUSION,
     col = c('deeppink2', 'cyan3', Color_cell_state_mye[c(1, 2)])[momp_p136.sce$Cell_state_orig],
     asp = 0.4, pch = 19)
lines(slingshot::SlingshotDataSet(momp_p136.sce), lwd=2,  col = 'black')
dev.off()

####  Plot 3 pseudotime   ####
tmp.srt <- momp_p136.srt
p3 <- FeaturePlot2(tmp.srt, features = 'PST_slingshot_rank',
                   reduction = 'diffusion', pt.size = 0.5, order = F) +
        scale_color_viridis_c(option = 'magma') +
        theme(aspect.ratio = 0.4)
p3
PlotPDF('03.3.umap.pgf_mono_mac_slingshot_pseudotime', 5, 2)
p3
dev.off()

####  Plot 4 feat PCC with PBMC monocyte  ####
mye.srt <- FindVariableFeatures(mye.srt)
data("pbmc3k.final")
genes <- intersect(rownames(pbmc3k.final), mye.srt@assays$CBN@var.features)
ave.exp <- AverageExpression(pbmc3k.final, assays = 'RNA')$RNA
ref <- rowMeans(ave.exp[genes, c('CD14+ Mono', 'FCGR3A+ Mono')])
exp <- (as.matrix(mye.srt@assays$CBN@data))[genes,]
pcc <- rep(NA, ncol(mye.srt))
for(i in 1:L(pcc)){pcc[i] <- cor(exp[, i], ref)}
tmp.srt <- mye.srt
tmp.srt$pcc <- pcc
p4 <- FeaturePlot3(tmp.srt, features = 'pcc', reduction = 'sub_umap', pt.size = 0.5, adjust = 1) +
        scale_color_distiller(palette = 'RdYlBu', values = c(0, 0.7, 1))
p4
PlotPDF('03.4.umap.pbmc_mono_similarity', 4, 4)
p4
dev.off()


####  Plot 5 heat mono-mp2 up genes  ####
module <- readRDS('analysis/STEP20.p136_mono_mp1_mp2_sling_pst_coupled_genes.list.rds')
module <- module$MoMP2$Pos
goi <- c('F13A1',
         'RBPJ',
         'NRP1',
         'PDGFC',
         'MAF',
         'MRC1',
         'CD163',
         'RUNX1'
         )
tmp.srt <- mye.srt[, mye.srt$Cell_state %in% c('Mono', 'MP2') & mye.srt$name2 != 'Control']
tmp.srt$Group <- factor('Mono', levels = c("ResMP 5d PT", "MoMP2 15m PT", "MoMP2 12y PT", "MoMP2 5d PT", 'Mono'))
tmp.srt$Group[Cells(mp2.srt)] <- mp2.srt$Group
tmp.srt <- AddModuleScore2(tmp.srt, features = list(module), names = 'UpGene_Score', return_z = T)
df <- AverageExpression(tmp.srt, assays = 'CBN', group.by = 'Group', features = goi)$CBN
for(i in 1:nrow(df)){df[i,] <- scale(df[i,])}

gene_order <- rownames(df)[order(df[, 'ResMP 5d PT'], decreasing = T)]
df <- melt(df[gene_order, ])
colnames(df) <- c('Gene', 'Group', 'Expr')
df$Gene <- factor(df$Gene, levels = gene_order)
df$Expr[df$Expr > 1] <- 1
df$Expr[df$Expr < -1] <- -1
p5 <- ggplot(df) +
        geom_tile(aes(x = Group, y = Gene, fill = Expr)) +
        scale_fill_distiller(palette = 'RdYlBu') +
        scale_y_discrete(limits = rev) +
        labs(x = 'MoMP2 trajectory-correlated genes', y = '', fill = 'Expr') +
        theme_classic() +
        theme(aspect.ratio = 13/4) +
        RotatedAxis()
p5
PlotPDF('03.5.heat.momp2_trajectory_correlated_genes2', 3, 4)
p5
dev.off()


####  Table 4 heat mono-mp2 up genes  ####
module <- readRDS('analysis/STEP20.p136_mono_mp1_mp2_sling_pst_coupled_genes.list.rds')
module <- module$MoMP2$Pos
data <- data.frame('Gene_symbol' = module)
data$Gene_symbol <- paste0("'", data$Gene_symbol, "'")
WriteCSV(data, 'STEP90.03.1.table.momp2_trajectory_positively_correlated_genes')


####  Plot 6 pc allograft MoMP2 similarity to ResMP2  ####
pc_dim <- readRDS('analysis/STEP23.mp2_harmonized_pc.srt_redu.rds')
df <- as.data.frame(pc_dim@cell.embeddings[, c(1, 3)])
df$Group[Cells(mp2.srt)] <- as.vector(mp2.srt$Group)
df$Group <- factor(df$Group, levels = levels(mp2.srt))
df2 <- df |>
        group_by(Group) |>
        summarize(PC_1 = mean(harmony_1), PC_2 = mean(harmony_3))
p6 <- ggplot() +
        geom_point(data = df, mapping = aes(x = harmony_1, y = harmony_3, color = Group), alpha = 0.2) +
        geom_point(data = df2, mapping = aes(x = PC_1, y = PC_2, color = Group), size = 3) +
        scale_color_manual(values = c(Color_genotype2[1:2], Color_genotype4[1:2])) +
        #scale_x_continuous(limits = c(-30, 30)) +
        labs(color = '', x = 'PC1', y = 'PC2') +
        theme_classic() +
        theme(aspect.ratio = 1, axis.text = element_blank(), axis.ticks = element_blank())
p6
PlotPDF('03.6.pca.mp2_pca_origin', 8, 8)
p6
dev.off()

####  Plot 7 tree allograft MoMP2 similarity to ResMP2  ####
pc_dim <- readRDS('analysis/STEP23.mp2_harmonized_pc.srt_redu.rds')
df <- as.data.frame(pc_dim@cell.embeddings[, 1:2])
df$Group[Cells(mp2.srt)] <- as.vector(mp2.srt$Group)
df$Group <- factor(df$Group, levels = levels(mp2.srt))
df2 <- df |> group_by(Group) |> summarize(PC_1 = median(harmony_1), PC_2 = median(harmony_2))
mtx <- as.matrix(df2[, 2:3])
rownames(mtx) <- df2$Group
hc <- hclust(dist(mtx, method = 'euclidean'), method = 'complete')
p7 <- ggdendrogram(hc, rotate = T, size = 2, leaf_labels = T)
p7
PlotPDF('03.7.dendro.mp2_pca_origin', 5, 5)
p7
dev.off()

####  Plot 8 baton MP2 composition allograft   ####
tmp.srt <- full.srt[, full.srt$Cell_type == 'Mye' & full.srt$name2 != 'Control' & full.srt$Cell_type_non_ambig]
tmp.srt$tmp <- paste(tmp.srt$Origin, tmp.srt$Cell_state)
tmp.srt$name <- droplevels(tmp.srt$name)
df <- CountCellBarPlot(tmp.srt, group.var = 'name', stack.var = 'tmp')$data
total <- Table(tmp.srt$name)
df$Pct <- df$Count/total[df$GroupVar]
df2 <- df[df$StackVar %in% c('Donor MP2'),]
p8.1 <- ggplot(df2, aes(x = GroupVar, y = Pct, group = '')) +
        geom_line(color = 'black') +
        geom_point() +
        scale_y_continuous(limits = c(0, 0.2), minor_breaks = seq(0, 0.2, 0.05), breaks = seq(0, 0.2, 0.05)) +
        theme(aspect.ratio = 1.5,
              panel.background = element_blank(),
              panel.grid.major.y = element_line(colour = 'grey'),
              panel.grid.minor.y = element_line(colour = 'grey')) +
        labs(x = '', y = 'Fraction of all Myeloid cells', title = 'Donor MP2 abundance')
df3 <- df[df$StackVar %in% c('Recipient MP2'),]
p8.2 <- ggplot(df3, aes(x = GroupVar, y = Pct, group = '')) +
        geom_line(color = 'black') +
        geom_point() +
        scale_y_continuous(limits = c(0, 1), minor_breaks = seq(0, 1, 0.2), breaks = seq(0, 1, 0.2)) +
        theme(aspect.ratio = 1.5,
              panel.background = element_blank(),
              panel.grid.major.y = element_line(colour = 'grey'),
              panel.grid.minor.y = element_line(colour = 'grey')) +
        labs(x = '', y = 'Fraction of all Myeloid cells', title = 'Recipient MP2 abundance')
p8.1 | p8.2
PlotPDF('03.8.baton.mp2_ratio_in_allograft', 8, 5)
p8.1 | p8.2
dev.off()

####  Plot 9 scatter MP2 proliferation   ####
tmp.srt <- full.srt[, full.srt$Cell_state == 'MP2']
dsp <- DownsampleByMeta(tmp.srt, meta_var = 'name2', down_to_min_group = T, random = T)
tmp.srt <- tmp.srt[, unlist(dsp)]
Idents(tmp.srt) <- 'name2'
p9 <- FeatureScatter(tmp.srt, feature1 = "S.Score", feature2 = "G2M.Score",
                     group.by = 'name2', cols = c('grey20', Color_transplant)) +
        scale_x_continuous(breaks = seq(-0.2, 1, 0.4), limits = c(-0.2, 0.8)) +
        scale_y_continuous(breaks = seq(-0.2, 1.2, 0.4), limits = c(-0.2, 1.2)) +
        theme(panel.grid.major = element_line(colour = 'grey80'),
              aspect.ratio = 1) +
        facet_wrap(~colors)
p9
x <- p9$data
x$Bin <- x$S.Score > 0.4 | x$G2M.Score > 0.4
x <- Table(x$Bin, x$colors)
x <- x[2,]*100/(x[1,]+x[2,])
x ## percentage of proliferative cells
PlotPDF('03.9.scatter.mp2_cell_cycle_score', 6, 6)
p9
dev.off()

####  Plot 10 box myeloid expressing Bajpai/Dick mac markers  ####
TRMP_gene <- c('CD163L1', 'MRC1', 'MAF', 'SIGLEC1', 'CCL8', 'CCL14', 'LILRB5', 'LYVE1', 'IL2RA', 'PDGFC', 'WLS', 'DAB2',
               'NRP1', 'SCN9A', 'FGF13', 'GDF15', 'IGF1', 'FMOD', 'SLIT3', 'EGFL7', 'ECM1', 'SDC3')
genes <- read_excel('external/sciimmunol.abf7777_tables_s1_to_s9/sciimmunol.abf7777_table_s7.xlsx',
                    sheet = 3, col_names = T, skip = 1) ## These genes are from the E. Slava's Science Immun Paper
genes <- genes[genes$avg_logFC > 0.75 & genes$p_val_adj < 0.001, ]
gl <- split(genes$gene, genes$cluster)
module <- readRDS('analysis/STEP20.p136_mono_mp1_mp2_sling_pst_coupled_genes.list.rds')
module <- module$MoMP2$Pos
tmp.srt <- AddModuleScore2(mp2.srt, features = list(U(c(TRMP_gene, gl$MF1)),
                                                    module),
                           names = c('ResMP_Score', 'Module_Score'), return_z = F)
levels(tmp.srt) <- c("ResMP 5d PT", "MoMP2 15m PT", "MoMP2 12y PT", "MoMP2 5d PT")

p10.1 <- BoxPlot(tmp.srt, feature = 'ResMP_Score', group.by = NULL,
               cols = Color_mp2_group[c(1,3,4,2)])
p_mtx <- matrix(NA, 4, 4, dimnames = list(levels(tmp.srt), levels(tmp.srt)))
data <- split(p10.1$data$ResMP_Score, p10$data$ident)
for(i in 1:4){for(j in 1:4){p_mtx[i,j] <- t.test(data[[i]], data[[j]])$p.value}}
p_mtx < 1e-3
p10.2 <- BoxPlot(tmp.srt, feature = 'Module_Score', group.by = NULL,
                 cols = Color_mp2_group[c(1,3,4,2)])
p_mtx <- matrix(NA, 4, 4, dimnames = list(levels(tmp.srt), levels(tmp.srt)))
data <- split(p10.2$data$Module_Score, p10$data$ident)
for(i in 1:4){for(j in 1:4){p_mtx[i,j] <- t.test(data[[i]], data[[j]])$p.value}}
p_mtx < 1e-3
p10 <- p10.1/p10.2
p10
PlotPDF('03.10.box.m2_marker_on_mp2_groups', 5, 10)
p10
dev.off()

####  Plot 11  trend mono-mp2 up genes  ####
module <- readRDS('analysis/STEP20.p136_mono_mp1_mp2_sling_pst_coupled_genes.list.rds')
module <- module$MoMP2
tmp.srt <- momp_p136.srt[, momp_p136.srt$Cell_state %in% c('Mono', 'MP2')]
tmp.srt <- ScaleData(tmp.srt, features = c(module$Pos, module$Neg), assay = 'alra')
expr_mat <- tmp.srt@assays$alra@scale.data[c(module$Pos, module$Neg),
                                           colnames(tmp.srt)[order(tmp.srt$PST_slingshot_rank)]]
for(i in 1:nrow(expr_mat)){expr_mat[i, ] <- smooth.spline(expr_mat[i, ], spar = 1.2)$y}

pos_expr_mat <- expr_mat[module$Pos,]
c('F13A1', 'RBPJ', 'NRP1', 'PDGFC', 'MAF', 'VAV3') %in% module$Pos

par(pty="s")
PlotPDF('03.11.trend.p136_mono_mac2_diffusion_pseudotime_module', 4, 4)
plot(1:ncol(pos_expr_mat), pos_expr_mat[1, ], type = 'l', col = alpha('black', 0.05))
for(i in 2:nrow(pos_expr_mat)){
        points(1:ncol(pos_expr_mat), pos_expr_mat[i, ], type = 'l', col = alpha('black', 0.05))}
dev.off()
####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~----


####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
####  Fig 4  MoMP2 Differ in Phagocytotic and Firbogenic Functions  ####
####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~----
####  Plot 1 umap ResMP and MoMP2  ####
p1 <- DimPlotSplit(mp2.srt, reduction = 'sub_umap',  split_by = 'Group', split_order = levels(mp2.srt$Group),
                   ncol = 2, cols.highlight = Color_mp2_group)
p1[[2]] <- NULL
PlotPDF('04.1.umap.resmp_vs_late_momp2', 4, 4)
p1
dev.off()

p1.1 <- DimPlot2(mp2.srt, reduction = 'sub_umap',  group.by = 'Group',
                 cols = c(Color_mp2_group[1], 'grey80', 'grey80', 'grey80'))
p1.2 <- DimPlot2(mp2.srt, reduction = 'sub_umap',  group.by = 'Group',
                 cols = c('grey80', 'grey80', Color_mp2_group[3], Color_mp2_group[4]))
p1.1 + p1.2
PlotPDF('04.1.umap.resmp_vs_late_momp2_v2', 8, 3)
p1.1 + p1.2
dev.off()

####  Plot 2 bar ResMP vs late MoMP2 deg GO enrich ####
enrich <- readRDS('analysis/STEP23.mp2_resmp_late_momp2_deg.enrich.rds')
res_go <- enrich$GO$GO_ResMP[enrich$GO$GO_ResMP$p.adjust <= 0.05, ]
momp_go <- enrich$GO$`GO_Late MoMP2`[enrich$GO$`GO_Late MoMP2`$p.adjust <= 0.05, ]

go_choice_1 <- c('phagocytosis',
                 'negative regulation of leukocyte activation',
                 'Fc-gamma receptor signaling pathway',
                 'toll-like receptor 4 signaling pathway',
                 'vesicle organization')
gene_choice_1 <- list()
for(i in 1:5){
        pool <- unlist(str_split(res_go$geneID[res_go$Description == go_choice_1[i]], '/'))
        pool_sort <- pool[order(deg[pool, 'avg_log2FC'], decreasing = F)]
        gene_choice_1[[i]] <- head(pool_sort, 5)
}
df <- res_go[res_go$Description %in% go_choice_1, ]
df$Description <- factor(df$Description, levels = go_choice_1)
df <- df[order(df$Description),]
df$gene <- NA
for(i in 1: 5){df$gene[i] <- paste(gene_choice_1[[i]], collapse = '  ')}
p2.1 <- ggplot(df, aes(x = -log10(p.adjust), y = Description))+
        geom_bar(stat = 'identity', position = "dodge", fill = mycol_10[1], alpha = 0.4) +
        geom_text(aes(label = gene, x = 0.5, y = Description), color = 'black', size = 3, hjust = 0, face = 'italic') +
        scale_y_discrete(limit = rev) +
        labs(x = '-Log10 adjusted p value', y = 'Enriched GO terms', title = 'ResMP2') +
        theme_classic()
go_choice_2 <- c(
        'T cell tolerance induction' ,
        'chondroitin sulfate biosynthetic process',
        'response to fatty acid' ,
        'proteoglycan metabolic process' ,
        'membrane assembly'
)
gene_choice_2 <- list()
for(i in 1:5){
        pool <- unlist(str_split(momp_go$geneID[momp_go$Description == go_choice_2[i]], '/'))
        pool_sort <- pool[order(deg[pool, 'avg_log2FC'], decreasing = T)]
        gene_choice_2[[i]] <- head(pool_sort, 5)
}
df <- momp_go[momp_go$Description %in% go_choice_2, ]
df$Description <- factor(df$Description, levels = go_choice_2)
df <- df[order(df$Description),]
df$gene <- NA
for(i in 1: 5){df$gene[i] <- paste(gene_choice_2[[i]], collapse = '  ')}
p2.2 <- ggplot(df, aes(x = -log10(p.adjust), y = Description))+
        geom_bar(stat = 'identity', position = "dodge", fill = mycol_10[1], alpha = 0.4) +
        geom_text(aes(label = gene, x = 0.5, y = Description), color = 'black', size = 3, hjust = 0, face = 'italic') +
        scale_y_discrete(limit = rev) +
        labs(x = '-Log10 adjusted p value', y = 'Enriched GO terms', title = 'Late MoMP2') +
        theme_classic()
p2.1/p2.2
PlotPDF('04.2.bar.resmp_vs_late_momp2_top_go_enrich', 10, 10)
p2.1/p2.2
dev.off()

####  Plot 3 volcano ResMP vs late MoMP2 DEG  ####
deg <- readRDS('analysis/STEP23.mp2_resmp_late_momp2_deg.mk.rds')
deg$gene <- rownames(deg)
genes <- rownames(deg)[abs(deg$avg_log2FC) > 1.5 & deg$p_val_adj <= 0.01]
tmp <- deg
tmp$p_val_adj[tmp$p_val_adj < 1e-20] <- 1e-20
p3 <- MarkerVolcano(tmp, label = T, label_genes = genes)
p3
PlotPDF('04.3.vol.resmp_vs_late_momp2_deg', 10, 10)
p3
dev.off()

####  Table 5 volcano ResMP vs late MoMP2 DEG  ####
deg <- readRDS('analysis/STEP23.mp2_resmp_late_momp2_deg.mk.rds')
deg$gene <- rownames(deg)
deg <- deg[deg$p_val_adj <= 0.01, ]
deg <- deg[, c('gene', 'avg_log2FC', 'p_val_adj')]
colnames(deg) <- c('Gene_symbol', 'Mean_log2_fold_change', 'Adjusted_p_value')
deg$Gene_symbol <- paste0("'", deg$Gene_symbol, "'")
deg <- deg[order(deg$Mean_log2_fold_change, decreasing = T),]
deg$Differentially_expressed <- 'ResMP High'
deg$Differentially_expressed[deg$Mean_log2_fold_change < 0] <- 'Late MoMP2 High'
WriteCSV(deg, 'STEP90.04.1.table.resmp_vs_late_momp2_deg')

####  Plot 4 box phago glyco functional geneset  ####
# grp1 <- read.table('external/GOBP_PHAGOCYTOSIS.v7.5.1.grp', comment.char = '#', header = T)[,1]
# grp1 <- intersect(grp1, VariableFeatures(mp2_sub.srt)) ## 61
grp2 <- read.table('external/KEGG_FC_GAMMA_R_MEDIATED_PHAGOCYTOSIS.v7.5.1.grp', comment.char = '#', header = T)[,1]
grp2 <- intersect(grp2, VariableFeatures(mp2_sub.srt)) ## 27
grp3 <- read.table('external/KEGG_GLYCOSAMINOGLYCAN_BIOSYNTHESIS_CHONDROITIN_SULFATE.v7.5.1.grp',
                   comment.char = '#', header = T)[,1]
grp3 <- intersect(grp3, VariableFeatures(mp2_sub.srt)) ## 8

tmp.srt <- AddModuleScore2(mp2_sub.srt, features = list(grp2, grp3),
                           names = c('KEGG_phago', 'KEGG_glyco'),
                           return_z = T)
p4.1 <- BoxPlot(tmp.srt, group.by = 'Group', feature = 'KEGG_phago', cols = Color_mp2_group[-c(2)]) + NoLegend()
p4.2 <- BoxPlot(tmp.srt, group.by = 'Group', feature = 'KEGG_glyco', cols = Color_mp2_group[-c(2)])
p4 <- p4.1 + p4.2
p4
data <- tmp.srt@meta.data[, c('KEGG_phago', 'KEGG_glyco', 'Group')]
pval <- list()
for(i in 1:2){
        pval[[i]] <- c(t.test(as.vector(data[data$Group == 'ResMP 5d PT', i]),
                              as.vector(data[data$Group == 'MoMP2 15m PT', i]))$p.value,
                       t.test(as.vector(data[data$Group == 'ResMP 5d PT', i]),
                              as.vector(data[data$Group == 'MoMP2 12y PT', i]))$p.value)
}
all(unlist(pval) < 0.001)
PlotPDF('04.4.box.mp2_phago_geneset_expr', 8, 4)
p4
dev.off()

####  Plot 5 vortex cellchat individual LR pair  ####
ResLR <- list(
        'FGF13_FGFR1' = 'FGF',
        'HBEGF_EGFR_ERBB2' = 'EGF',
        'VEGFB_VEGFR1' = 'VEGF'
)
MoMP2LR <- list(
        'PDGFB_PDGFRA' = 'PDGF',
        'PDGFB_PDGFRB' = 'PDGF',
        'WNT5B_FZD6' = 'ncWNT',
        'TGFB2_TGFBR1_TGFBR2' = 'TGFb',
        'IGF1_IGF1R' = 'IGF',
        'TGFA_EGFR' = 'EGF'
)
PlotPDF('04.5.vertex.cellchat_resmp_late_momp2_to_all_cardiac', 4, 4)
for(i in 1:L(ResLR)){
        print(netVisual_individual(mp2_cc_list$`5d PT`,
                                   sources.use = 'MP2',
                                   targets.use = c('CM', 'EC', 'FB', 'PC', 'SMC'),
                                   pairLR.use = names(ResLR)[i],
                                   signaling = ResLR[[i]],
                                   color.use = Color_cell_type,
                                   layout = "circle"))
}
for(i in 1:L(MoMP2LR)){
        print(netVisual_individual(mp2_cc_list$`15m PT`,
                                   sources.use = 'MP2',
                                   targets.use = c('CM', 'EC', 'FB', 'PC', 'SMC'),
                                   pairLR.use = names(MoMP2LR)[i],
                                   signaling = MoMP2LR[[i]],
                                   color.use = Color_cell_type,
                                   layout = "circle"))
}
for(i in 1:L(MoMP2LR)){
        print(netVisual_individual(mp2_cc_list$`12y PT`,
                                   sources.use = 'MP2',
                                   targets.use = c('CM', 'EC', 'FB', 'PC', 'SMC'),
                                   pairLR.use = names(MoMP2LR)[i],
                                   signaling = MoMP2LR[[i]],
                                   color.use = Color_cell_type,
                                   layout = "circle"))
}
dev.off()

####  Plot 6 chord cellchat individual LR pair  ####
PlotPDF('04.6.chord.cellchat_late_momp2_vs_resmp_to_all_cardiac_common_signal', 6, 4)
netVisual_chord_gene(mp2_cc_list$`5d PT`, sources.use = 'MP2', targets.use = c('CM', 'EC', 'FB', 'PC', 'SMC'),
                     slot.name = 'netP', pairLR.use = mp2_cc_list$LR$res_pairLR, color.use = Color_cell_type,
                     scale = T, lab.cex = 0.8, small.gap = 3.5, reduce = 0.005,
                     title.name = paste0("Signaling from 5d PT ResMP to cardiac cells"))
netVisual_chord_gene(mp2_cc_list$`15m PT`, sources.use = 'MP2', targets.use = c('CM', 'EC', 'FB', 'PC', 'SMC'),
                     slot.name = 'netP', pairLR.use = mp2_cc_list$LR$momp2_pairLR, color.use = Color_cell_type,
                     scale = T, lab.cex = 0.8, small.gap = 3.5, reduce = 0.005,
                     title.name = paste0("Signaling from 15m PT MoMP2 to cardiac cells"))
netVisual_chord_gene(mp2_cc_list$`12y PT`, sources.use = 'MP2', targets.use = c('CM', 'EC', 'FB', 'PC', 'SMC'),
                     slot.name = 'netP', pairLR.use = mp2_cc_list$LR$momp2_pairLR, color.use = Color_cell_type,
                     scale = T, lab.cex = 0.8, small.gap = 3.5, reduce = 0.005,
                     title.name = paste0("Signaling from 12y PT MoMP2 to cardiac cells"))
dev.off()

####  Plot 7 box fibrosis score  ####
tmp.srt <- full.srt[, full.srt$Cell_type == 'FB']
gl <- read_excel('external/Kuppe_et_al_Nature_MI_Fibro_marker_SuppData.xlsx', sheet = 2)
gl <- gl$gene[gl$avg_log2FC>0.75]
## 20 genes used for scoring:
## 'POSTN', 'FN1', 'TNC', 'COL1A1', 'COL1A2', 'COL3A1', 'DEC1', 'RUNX1', 'ADAM12', 'KIF26B', 'THBS4', 'FAP', 'COL5A1',
## 'SERPINE1', 'SLC20A1', 'KALRN', 'PRICKLE1', 'CDH11', 'FGF14', 'THBS2'
tmp.srt <- AddModuleScore2(tmp.srt, features = list(gl), names = 'Fibrosis', return_z = T)
p7 <- BoxPlot(tmp.srt, feature = 'Fibrosis', group.by = 'name2', cols = Color_condition)
p7
p_mtx <- matrix(NA, 4, 4, dimnames = list(levels(tmp.srt$name2), levels(tmp.srt$name2)))
data <- split(p7$data$Fibrosis, p7$data$ident)
for(i in 1:4){for(j in 1:4){p_mtx[i,j] <- t.test(data[[i]], data[[j]])$p.value}}
p_mtx < 0.05
PlotPDF('04.7.box.fibrosis_score', 4, 4)
p7
dev.off()

####  Plot 8 upset MoMP2 vs ResMP DEG overlap  ####
deg_all <- readRDS('analysis/STEP23.mp2_3_way_compare_with_resmp.srt_mk.rds')
deg_all <- deg_all[deg_all$cluster != 'MoMP2 5d PT', ]
deg_pos <- deg_all[deg_all$p_val_adj < 0.05 & deg_all$avg_log2FC > 0.5, ]
deg_pos_df <- data.frame(Genes = U(deg_pos$gene), Cluster = NA)
for(i in 1:nrow(deg_pos_df)){
        deg_pos_df$Cluster[i] <- list(as.vector(deg_pos$cluster[deg_pos$gene == deg_pos_df$Genes[i]]))
}
deg_neg <- deg_all[deg_all$p_val_adj < 0.05 & deg_all$avg_log2FC < -0.5, ]
deg_neg_df <- data.frame(Genes = U(deg_neg$gene), Cluster = NA)
for(i in 1:nrow(deg_neg_df)){
        deg_neg_df$Cluster[i] <- list(as.vector(deg_neg$cluster[deg_neg$gene == deg_neg_df$Genes[i]]))
}
deg_df <- rbind(deg_pos_df, deg_neg_df)
p8 <- ggplot(deg_df) +
        geom_bar(aes(x=Cluster)) +
        scale_x_upset(order_by = 'degree') +
        theme_classic() +
        theme(aspect.ratio = 1)
p8
PlotPDF('04.8.upset.momp2_vs_resmp_deg', 4, 4)
p8
dev.off()
####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~----


####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
####  Fig 5  Recipient Monocytes Derived MP1 are Associated with Rejection   ####
####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~----
####  Plot 1  umap immune color by control and allograft time point ####
tmp.srt <- mye.srt[, (mye.srt$Origin == 'Recipient' | mye.srt$name2 == 'Control') & mye.srt$Cell_state != 'Mast']
sub <- DownsampleByMeta(tmp.srt, meta_var = 'name2', down_to_min_group = T, random = T)
tmp.srt <- tmp.srt[, unlist(sub)]
p1 <- DimPlot2(tmp.srt, reduction = 'sub_umap', split.by = 'name2', cols = Color_cell_state_mye, ncol = 2) &
        labs(x = 'UMAP1', y = "UMAP2") &
        scale_x_continuous(limits = c(-6.6, 7.8)) &
        scale_y_continuous(limits = c(-6.8, 6.8))
p1
PlotPDF('05.1.umap.immune_color_by_condition', 4, 4)
p1
dev.off()

####  Plot 2  box immune cell state fold change allograft vs each controls ####
df <- CountCellBarPlot(mye.srt, group.var = 'Cell_state', stack.var = 'name', percentage = T)$data
df$Fraction <- df$Count/rep(Table(full.srt$name), each = LU(df$GroupVar))
mtx1 <- matrix(NA, 4, LU(df$GroupVar), dimnames = list(paste0('Control ', 1:4), U(df$GroupVar)))
mtx2 <- mtx1
mtx3 <- mtx1
for(i in 1:4){
        for(j in 1:4){
                celltype <- levels(df$GroupVar)[j]
                mtx1[i, j] <- df$Fraction[df$StackVar == '5d PT' & df$GroupVar == celltype]/
                        df$Fraction[df$StackVar == paste('Control', i) & df$GroupVar == celltype]
                mtx2[i, j] <- df$Fraction[df$StackVar == '15m PT' & df$GroupVar == celltype]/
                        df$Fraction[df$StackVar == paste('Control', i) & df$GroupVar == celltype]
                mtx3[i, j] <- df$Fraction[df$StackVar == '12y PT' & df$GroupVar == celltype]/
                        df$Fraction[df$StackVar == paste('Control', i) & df$GroupVar == celltype]
        }
}
df2 <- rbind(melt(log2(mtx1)),  melt(log2(mtx2)), melt(log2(mtx3)))
df2$Group <- factor(rep(c('5d PT', '15m PT', '12y PT'), each = 16), levels = c('5d PT', '15m PT', '12y PT'))
colnames(df2) <- c('Patient', 'Cell.Type', 'Log2FC', 'Group')
p2 <- ggplot(df2, aes(x = Cell.Type, y = Log2FC, fill = Group)) +
        geom_boxplot(outlier.shape = NA, width = 0.5) +
        geom_hline(yintercept = c(0), color = 'red4', linetype = 'dashed') +
        scale_fill_manual(values = Color_condition[2:4]) +
        theme_classic() +
        scale_y_continuous(breaks = seq(-4, 8, 2), limits = c(-4, 8)) +
        theme(aspect.ratio = 1, panel.grid.major.y = element_line()) +
        RotatedAxis() +
        labs(x = '', y = 'Log2 fold change (PGF vs. 4 Controls)',
             title = 'Immune cell composition change', color = 'PGF vs.')
p2
PlotPDF('05.2.box.immune_celltype_foldchange_allograft_vs_controls', 5, 4)
print(p2)
dev.off()

####  Plot 3  trend mono-mp1 up genes  ####
module <- readRDS('analysis/STEP20.p136_mono_mp1_mp2_sling_pst_coupled_genes.list.rds')
module <- module$MoMP1
tmp.srt <- momp_p136.srt[, momp_p136.srt$Cell_state %in% c('Mono', 'MP1')]
tmp.srt <- ScaleData(tmp.srt, features = c(module$Pos, module$Neg), assay = 'alra')
expr_mat <- tmp.srt@assays$alra@scale.data[c(module$Pos, module$Neg),
                                           colnames(tmp.srt)[order(tmp.srt$PST_slingshot_rank)]]
for(i in 1:nrow(expr_mat)){expr_mat[i, ] <- smooth.spline(expr_mat[i, ], spar = 1)$y}
df_all <- mapply(function(x, y){
        expr_mat <- t(apply(expr_mat[x,], 1, Range01)) ## normalize to [0-1]
        df <- FlattenExpr(expr_mat[, seq(1, ncol(tmp.srt), 10)], x)
        df$DPT <- Range01(df$Cells)
        df$mod_id <- names(module)[y]
        return(df)},
        x = module, y = seq_along(module), SIMPLIFY = F)
df_all <- bind_rows(df_all)
name <- paste0("Module ", 1:2, ' (', lapply(module, function(x) length(x)), ' genes)')

p3 <- ggplot(data = df_all[df_all$mod_id == 'Pos', ]) +
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
p3
PlotPDF('05.3.trend.p136_mono_mac1_diffusion_pseudotime_module', 4, 4)
p3
dev.off()

c('NFKB1', 'RELB', 'STAT1', 'IRF1') %in% module$Pos

pos_expr_mat <- expr_mat[module$Pos,]
par(pty="s")
PlotPDF('05.3.trend.p136_mono_mac1_diffusion_pseudotime_module_v2', 4, 4)
plot(1:ncol(pos_expr_mat), pos_expr_mat[1, ], type = 'l', col = alpha('black', 0.05))
for(i in 2:nrow(pos_expr_mat)){
        points(1:ncol(pos_expr_mat), pos_expr_mat[i, ], type = 'l', col = alpha('black', 0.05))}
dev.off()

####  Table 6 heat mono-mp1 up genes  ####
module <- readRDS('analysis/STEP20.p136_mono_mp1_mp2_sling_pst_coupled_genes.list.rds')
module <- module$MoMP1$Pos
data <- data.frame('Gene_symbol' = module)
data$Gene_symbol <- paste0("'", data$Gene_symbol, "'")
WriteCSV(data, 'STEP90.05.1.table.momp1_trajectory_positively_correlated_genes')


####  Plot 4  heat scenic MP1 + Mono  ####
## See script 22

####  Plot 5  feat density MP1 cytokine markers  ####
p5 <- FeaturePlot3(mye.srt[, mye.srt$name2 != 'Control'],
                   features = c('TNF', 'STAT1', 'CD80',  'CD83',
                                'IL1A', 'IL1B', 'IL15', 'IL18',
                                'CCL3', 'CCL4', 'CXCL2', 'CXCL3'),
                   ncol = 4, pt.size = 1.8,
                   reduction = 'sub_umap', adjust = 2) &
        NoLegend()
p5[[12]] <- p5[[12]] + RestoreLegend()
p5 <- p5 & theme(axis.line = element_blank())
p5
PlotPDF('05.5.feature.cytotoxic_markers', 9, 13)
print(p5)
dev.off()

####  Plot 6  bar MP1 up genes GO & MSigDB enrichment  ####
deg <- full.srt@misc$marker$Cell_state_marker
deg_hi <- deg[deg$avg_log2FC>1 & deg$p_val_adj<0.001 & deg$cluster == 'MP1', 'gene']
enrich <- ModuleEnrichment(list(deg_hi), human_or_mouse = 'human')
go <- enrich$GO$GO_
go <- go[go$p.adjust <= 0.05, ]
go_choice <- list(
        'cytokine-mediated signaling pathway' = 'CCL4  CCL8  CXCL2  CSF2RA  CSF1R',
        'T cell proliferation' = 'DOCK2  IL15  IL1B  IL18  CD80',
        'response to interferon-gamma' = 'IRF1  IRF8  IFNGR2  CGAS  IFNAR2',
        'response to tumor necrosis factor' = 'BIRC3  JAK2  STAT1  TNFAIP3  NFKB1',
        'natural killer cell mediated cytotoxicity' = 'IL15  FGR  PIK3CD  IL18  SLAMF7'
)
df <- go[go$Description %in% names(go_choice), ]
df$Description <- factor(df$Description, levels = names(go_choice))
df <- df[order(df$Description),]
df$gene <- go_choice
p6.1 <- ggplot(df, aes(x = -log10(p.adjust), y = Description))+
        geom_bar(stat = 'identity', position = "dodge", fill = mycol_10[1], alpha = 0.4) +
        geom_text(aes(label = gene, x = 0.5, y = Description), color = 'black', size = 3, hjust = 0, face = 'italic') +
        scale_y_discrete(limit = rev) +
        labs(x = '-Log10 adjusted p value', y = 'Enriched GO terms') +
        theme_classic()
MSigDB_enrich <- read.table('external/MSigDB_Hallmark_2020_table.txt', sep = '\t', header = T) ## from EnrichR
MSigDB_choice <- list(
        'Interferon Gamma Response' = 'IFI30  TNFSF10  IL15  IRF1  IFNAR2',
        'TNF-alpha Signaling via NF-kB' = 'TNFAIP8  TNFAIP6  NFKB1  RELB  REL',
        'Inflammatory Response' = 'CXCL9  CXCL8  TNFSF10  IL15  IL1B',
        'IL-6/JAK/STAT3 Signaling' = 'CSF3R  CSF2RB  IL6ST  CD38  TLR2',
        'Allograft Rejection' = 'CD86  B2M  PRKCB  IL15  IL18')
df <- MSigDB_enrich[MSigDB_enrich$Term %in% names(MSigDB_choice), ]
df$Description <- factor(df$Term, levels = names(MSigDB_choice))
df <- df[order(df$Description),]
df$gene <- MSigDB_choice
p6.2 <- ggplot(df, aes(x = -log10(Adjusted.P.value), y = Description))+
        geom_bar(stat = 'identity', position = "dodge", fill = mycol_10[1], alpha = 0.4) +
        geom_text(aes(label = gene, x = 0.5, y = Description), color = 'black', size = 3, hjust = 0, face = 'italic') +
        scale_y_discrete(limit = rev) +
        labs(x = '-Log10 adjusted p value', y = 'Enriched MSigDB terms') +
        theme_classic()
wrap_plots(p6.1, p6.2, ncol = 1)
PlotPDF('05.6.bar.mac1_upgene_go_msigdb_enrichment', 5.5, 5)
wrap_plots(p6.1, p6.2, ncol = 1)
dev.off()

####  Plot 7  pca MP subset compare with DCM   ####
df <- as.data.frame(mp_dcm.srt@reductions$pca@cell.embeddings[, 1:2])
df$Group <- mp_dcm.srt$name4
df$Group2 <- factor(df$Group, levels = c('Control MP (n = 6)', '5d PT MP1', '5d PT MP2', 'DCM MP (n = 5)'))
df$Group2[df$Group %in% paste('Control', 1:6)] <- 'Control MP (n = 6)'
df$Group2[df$Group %in% paste('DCM', 1:5)] <- 'DCM MP (n = 5)'

df2 <- df |>
        group_by(Group) |>
        summarize(PC_1 = mean(PC_1), PC_2 = mean(PC_2))
df2$Group2 <- factor(df2$Group, levels = c('Control MP (n = 6)', '5d PT MP1', '5d PT MP2', 'DCM MP (n = 5)'))
df2$Group2[df2$Group %in% paste('Control', 1:6)] <- 'Control MP (n = 6)'
df2$Group2[df2$Group %in% paste('DCM', 1:5)] <- 'DCM MP (n = 5)'
p7 <- ggplot() +
        geom_point(data = df, mapping = aes(x = PC_1, y = PC_2, color = Group2), alpha = 0.1, size = 1) +
        geom_point(data = df2, mapping = aes(x = PC_1, y = PC_2, color = Group2), size = 3) +
        scale_color_manual(values = mycol_10[c(1, 3, 4, 2)]) +
        scale_x_continuous(limits = c(-30, 20)) +
        scale_y_continuous(limits = c(-30, 20)) +
        labs(color = '') +
        theme_classic() +
        theme(aspect.ratio = 1, axis.text = element_blank(), axis.ticks = element_blank())
p7
PlotPDF('05.7.pca.mp_vs_dcm', 5, 5)
p7
dev.off()

####  Plot 8  heat Mac PGF vs DCM selected genes  ####
MSigDB_enrich <- read.table('external/MSigDB_Hallmark_2020_table.txt', sep = '\t', header = T) ## from EnrichR
genelist <- str_split(MSigDB_enrich$Genes[MSigDB_enrich$Term == 'Allograft Rejection'], ';')[[1]]
rejection_marker <- read_excel('external/Piening_table1.xlsx')
genes <- intersect(rejection_marker$Gene_symbol[rejection_marker$logFC>1.5], rownames(mp_dcm.srt))
mp_dcm.srt$tmp2 <- as.vector(mp_dcm.srt$name2)
mp_dcm.srt$tmp2[mp_dcm.srt$tmp2 == '5d PT'] <- as.vector(mp_dcm.srt$Cell_state)[mp_dcm.srt$tmp2 == '5d PT']
mp_dcm.srt$tmp2[mp_dcm.srt$tmp2 == 'Donor'] <- 'Control'
mp_dcm.srt$tmp2 <- factor(mp_dcm.srt$tmp2, levels = c('DCM', 'MP1', 'MP2', 'Control'))
df <- AverageExpression(mp_dcm.srt, features = genelist, group.by = 'tmp2')$RNA
for(i in 1:nrow(df)){df[i,] <- scale(df[i,])}
gene_order <- rownames(df)[order(df[, 'MP1'], decreasing = T)]
gene_order <- gene_order[1:10]
df <- melt(df[gene_order, ])
colnames(df) <- c('Gene', 'Group', 'Expr')
df$Gene <- factor(df$Gene, levels = gene_order)
df$Expr[df$Expr > 1.5] <- 1.5
df$Expr[df$Expr < -0.5] <- -0.5
p8.1 <- ggplot(df) +
        geom_tile(aes(x = Gene, y = Group, fill = Expr)) +
        scale_fill_distiller(palette = 'RdYlBu', limits = c(-0.5, 1.5)) +
        scale_y_discrete(limits = rev) +
        labs(x = 'Allograft rejection-related genes', y = '', fill = 'Expr') +
        theme_classic() +
        theme(aspect.ratio = 5/10) +
        RotatedAxis()
df <- AverageExpression(mp_dcm.srt, features = genes, group.by = 'tmp2')$RNA
for(i in 1:nrow(df)){df[i,] <- scale(df[i,])}
gene_order <- rownames(df)[order(df[, 'MP1'], decreasing = T)]
gene_order <- gene_order[1:10]
df <- melt(df[gene_order, ])
colnames(df) <- c('Gene', 'Group', 'Expr')
df$Gene <- factor(df$Gene, levels = gene_order)
df$Expr[df$Expr > 1.5] <- 1.5
df$Expr[df$Expr < -0.5] <- -0.5
p8.2 <- ggplot(df) +
        geom_tile(aes(x = Gene, y = Group, fill = Expr)) +
        scale_fill_distiller(palette = 'RdYlBu', limits = c(-0.5, 1.5)) +
        scale_y_discrete(limits = rev) +
        labs(x = 'Rejection signaure genes Piening et al.', y = '', fill = 'Expr') +
        theme_classic() +
        theme(aspect.ratio = 5/10) +
        RotatedAxis()
p8.1 | p8.2
PlotPDF('05.8.heat.expression_mp1_high_genes_vs_dcm', 10, 5)
print(p8.1 | p8.2)
dev.off()

####  Plot 9  vortex MP1 signal to TC and Mono  ####
PlotPDF('05.9.chord.pgf_mp1_tc_mono', 8, 8)
netVisual_individual(p136_cc_list$Simple,
                     sources.use = c('MP1'),
                     targets.use = c('CM', 'EC', 'FB', 'PC', 'SMC', 'Neu', 'TC', 'NK', 'MP1', 'MP2', 'Mono'),
                     group = c('CM', 'EC', 'FB', 'PC', 'SMC', 'Neu','TC', 'NK', 'MP1', 'MP2', 'Mono'),
                     cell.order = c(1:6, 8:9, 7, 10:12),
                     pairLR.use = c('CCL4_CCR5', 'CCL3_CCR5', 'IL1A_IL1R2', 'IL1B_IL1R2'),
                     signaling = c('CCL', 'IL1'),
                     #color.use = mycol_20,
                     layout = "circle")
dev.off()

####  Plot 10  vln MP1 signal to TC and Mono  ####
p10 <- CellChat::plotGeneExpression(p136_cc_list$Simple,
                                    signaling = c('CCL', 'IL1'),
                                    features = c('CCL3', 'CCL4', 'CCR5', 'IL1A', 'IL1B', 'IL1R2'))
p10
PlotPDF('05.10.vln.pgf_mp1_tc_mono', 8, 8)
p10
dev.off()
####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~----


####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
####  Fig 6 NK3 is a Hypertoxic NK Subtype Activated by MP1  ####
####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~----
####  Plot 1 umap of lymphoid cells color by cell subtype  ####
p1 <- DimPlot2(lym.srt, reduction = 'sub_umap', group.by = 'Cell_state',
               raster = F, cols = c(mycol_10[1:5])) +
        labs(x = 'UMAP1', y = "UMAP2", title = 'Lymphoid cell subtypes')
p1
PlotPDF('06.1.umap.lym_color_by_cell_subtype', 5, 5)
p1
dev.off()

####  Plot 2  box immune cell state fold change allograft vs each controls ####
df <- Table(lym.srt$Cell_state, lym.srt$name) + 0.1
df <- melt(df)
df$Fraction <- df$value/rep(Table(full.srt$name), each = LU(df$Var1))
mtx1 <- matrix(NA, 4, LU(df$Var1), dimnames = list(paste0('Control ', 1:4), U(df$Var1)))
mtx2 <- mtx1
mtx3 <- mtx1
for(i in 1:4){
        for(j in 1:5){
                celltype <- levels(df$Var1)[j]
                mtx1[i, j] <- df$Fraction[df$Var2 == '5d PT' & df$Var1 == celltype]/
                        df$Fraction[df$Var2 == paste('Control', i) & df$Var1 == celltype]
                mtx2[i, j] <- df$Fraction[df$Var2 == '15m PT' & df$Var1 == celltype]/
                        df$Fraction[df$Var2 == paste('Control', i) & df$Var1 == celltype]
                mtx3[i, j] <- df$Fraction[df$Var2 == '12y PT' & df$Var1 == celltype]/
                        df$Fraction[df$Var2 == paste('Control', i) & df$Var1 == celltype]
        }
}
df2 <- rbind(melt(log2(mtx1)),  melt(log2(mtx2)), melt(log2(mtx3)))
df2$Group <- factor(rep(c('5d PT', '15m PT', '12y PT'), each = 4*5), levels = c('5d PT', '15m PT', '12y PT'))
colnames(df2) <- c('Patient', 'Cell.Type', 'Log2FC', 'Group')
p2 <- ggplot(df2, aes(x = Group, y = Log2FC, fill = Group)) +
        geom_boxplot(outlier.shape = NA, width = 0.5) +
        geom_hline(yintercept = c(0), color = 'red4', linetype = 'dashed') +
        scale_fill_manual(values = Color_condition[2:4]) &
        theme_classic() &
        scale_y_continuous(breaks = seq(0, 12, 2), limits = c(-0.5, 12)) &
        theme(aspect.ratio = 3, panel.grid.major.y = element_line()) &
        RotatedAxis() &
        labs(x = '', y = 'Log2 fold change (Allograft vs. 4 Controls)', fill = '',
             title = 'Lymphoid cell composition change') &
        facet_wrap(~Cell.Type, drop = T, nrow = 1)
p2
PlotPDF('06.2.box.lym_foldchange_allograft_vs_controls', 8, 4)
print(p2)
dev.off()

####  Plot 3 dot lym cell state marker  ####
p3.1 <- DotPlot2(lym.srt, group.by = 'Cell_state',
                 features = c('CD2', 'CD5', 'BCL11B'),
                 assay = 'CBN', cols = c('grey85', mycol_BuGr[1]),
                 idents = c('CD4 TC', 'CD8 TC', 'NK1', 'NK2', 'NK3'),
                 col.min = 0) +
        theme(aspect.ratio = 8/3) +
        labs(title = 'pan-TC\nmarker') +
        NoLegend()
p3.2 <- DotPlot2(lym.srt, group.by = 'Cell_state',
                 features = c('CD4', 'CD40LG', 'IL7R',
                              'CD8B', 'TOX', 'THEMIS'),
                 assay = 'CBN', cols = c('grey85', mycol_BuGr[1]),
                 idents = c('CD4 TC', 'CD8 TC', 'NK1', 'NK2', 'NK3'),
                 col.min = 0) +
        theme(axis.text.y = element_blank(), aspect.ratio = 8/6) +
        labs(title = 'CD4 and CD8 \nTC marker') +
        NoLegend()
p3.3 <- DotPlot2(lym.srt, group.by = 'Cell_state',
                 features = c('KLRC2', 'FCRL3', 'KLRC4',
                              'GNLY', 'PRF1', 'KLRF1',
                              'GZMB', 'IFNG', 'REL'),
                 assay = 'CBN', cols = c('grey85', mycol_BuGr[1]),
                 idents = c('CD4 TC', 'CD8 TC', 'NK1', 'NK2', 'NK3'),
                 col.min = 0)+
        theme(axis.text.y = element_blank(), aspect.ratio = 8/9) +
        labs(title = 'NK subset marker')
p3.1 + p3.2 + p3.3
PlotPDF('06.3.dot.lym_cell_state_marker_select', 11.5, 5)
print(p3.1 + p3.2 + p3.3)
dev.off()

####  Plot 4  box PGF NK KEGG score ####
grp1 <- read.table('external/KEGG_NATURAL_KILLER_CELL_MEDIATED_CYTOTOXICITY.v7.5.1.grp',
                   comment.char = '#', header = T)[,1]
grp2 <- read.table('external/KEGG_ALLOGRAFT_REJECTION.v7.5.1.grp',
                   comment.char = '#', header = T)[,1]
tmp.srt <- AddModuleScore2(lym.srt, assay = 'alra',
                           features = list(grp1, grp2), names = c('Score_NK_KEGG_tox', 'Score_NK_KEGG_allo'))
levels(tmp.srt) <- levels(tmp.srt$Cell_state)
tmp.srt <- tmp.srt[, tmp.srt$name == '5d PT']
p4.1 <- BoxPlot(tmp.srt, feature = 'Score_NK_KEGG_tox', cols = mycol_10, group.by = 'Cell_state') +
        labs(title = 'KEGG: NK-mediated cytotoxicity') +
        scale_y_continuous(limits = c(0, 0.4))
p4.2 <- BoxPlot(tmp.srt, feature = 'Score_NK_KEGG_allo', cols = mycol_10, group.by = 'Cell_state') +
        labs(title = 'KEGG: Allograft rejection')  +
        scale_y_continuous(limits = c(-0.2, 0.3))
p4 <- p4.1 / p4.2 &
        NoLegend() &
        theme(aspect.ratio = 0.75, axis.text.x = element_blank())
p4
PlotPDF('06.4.box.nk_kegg_tox_pw_score', 5, 8)
p4
dev.off()

####  Plot 5  ternary PGF NK cytotoxicity gene expression ####
library('ggtern')
grp <- read.table('external/KEGG_NATURAL_KILLER_CELL_MEDIATED_CYTOTOXICITY.v7.5.1.grp', comment.char = '#', header = T)
df <- AverageExpression(lym.srt,
                        features = grp$KEGG_NATURAL_KILLER_CELL_MEDIATED_CYTOTOXICITY)
df <- df$alra[, c('NK1', 'NK2', 'NK3')]
df <- as.data.frame(df)
df$gene <- rownames(df)
expr <- lym.srt@assays$CBN@data
df$AveExpr <- rowMeans(expr[df$gene, ])
p5 <- ggtern(data = df, aes(NK1, NK3, NK2, size = AveExpr)) +
        geom_mask() +
        geom_point(fill = "red", shape = 21) +
        geom_text(data = df, aes(label = gene),
                  color = 'black', size = 1) +
        theme_rgbg()
p5
PlotPDF('06.5.ternary.nk_cytotoxic_genes', 5, 5)
print(p5)
dev.off()

####  Plot 6  heat nk subtype marker expression ####
tmp.srt <- lym.srt[, lym.srt$Cell_state %in% c('NK1', 'NK2', 'NK3')]
levels(tmp.srt) <- c('NK1', 'NK2', 'NK3')
mk <- FindAllMarkers(tmp.srt, assay = 'CBN', test.use = 'MAST',
                     only.pos = T, logfc.threshold = 0.25, return.thresh = 0.001)
mk_dedup <- mk |> group_by(cluster) |> top_n(n = 30, wt = avg_log2FC)
mk_dedup <- mk_dedup[! mk_dedup$gene %in% mk_dedup$gene[duplicated(mk_dedup$gene)], ]
p6 <- MarkerHeatmap(tmp.srt, marker.df = mk_dedup, disp.min = 0, top = 20, group.cols = mycol_10[3:5], raster = T)
p6
PlotPDF('06.6.heat.nk_subtype_marker', 15, 10)
p6
dev.off()

####  Table 7 heat nk cell subtype markers  ####
tmp.srt <- lym.srt[, lym.srt$Cell_state %in% c('NK1', 'NK2', 'NK3')]
levels(tmp.srt) <- c('NK1', 'NK2', 'NK3')
mk <- FindAllMarkers(tmp.srt, assay = 'CBN', test.use = 'MAST',
                     only.pos = T, logfc.threshold = 0.25, return.thresh = 0.001)
mk_dedup <- mk |> group_by(cluster) |> top_n(n = 30, wt = avg_log2FC)
mk_dedup <- mk_dedup[! mk_dedup$gene %in% mk_dedup$gene[duplicated(mk_dedup$gene)], ]
mk_dedup <- mk_dedup[, c('gene', 'cluster', 'avg_log2FC', 'p_val_adj')]
colnames(mk_dedup) <- c('Gene_symbol', 'Cell_type', 'Mean_log2_fold_change', 'Adjusted_p_value')
mk_dedup$Gene_symbol <- paste0("'", mk_dedup$Gene_symbol, "'")
WriteCSV(mk_dedup, 'STEP90.06.1.table.nk_cell_subtype_marker')

####  Plot 7  pca NK integrated with DCM ####
df <- as.data.frame(nk_dcm.srt@reductions$pca@cell.embeddings[, 1:2])
df$Group <- nk_dcm.srt$name4
df$Group2 <- as.vector(df$Group)
df$Group2[df$Group %in% paste('Control', 1:6)] <- 'Control NK'
df$Group2[df$Group %in% paste('DCM', 1:5)] <- 'DCM NK'
df2 <- df |>
        group_by(Group) |>
        summarize(PC_1 = mean(PC_1), PC_2 = mean(PC_2))
df2$Group2 <- as.vector(df2$Group)
df2$Group2[df2$Group %in% paste('Control', 1:6)] <- 'Control NK'
df2$Group2[df2$Group %in% paste('DCM', 1:5)] <- 'DCM NK'
p7 <- ggplot() +
        geom_point(data = df, mapping = aes(x = PC_1, y = PC_2, color = Group2), alpha = 0.2) +
        geom_point(data = df2, mapping = aes(x = PC_1, y = PC_2, color = Group2), size = 3) +
        scale_color_manual(values = mycol_10[c(2,1,3:5)]) +
        #scale_x_continuous(limits = c(-30, 30)) +
        theme_classic() +
        theme(aspect.ratio = 1, axis.text = element_blank(), axis.ticks = element_blank())
p7
PlotPDF('06.7.pca.nk_sc_by_condition_group', 5, 5)
p7
dev.off()

####  Plot 8  heat NK P136 vs DCM selected genes  ####
genelist <- readRDS('analysis/STEP25.nk_5dpt_vs_dcm_deg_select.rds')
genelist <- U(unlist(genelist[1:3]))
genes <- c('IKBKB', 'VAV3', 'STAT4', 'ST8SIA4', 'ITK', 'GBP2', 'FYN', 'IL4R', 'CRTAM', 'NFAT5')
nk_dcm.srt$tmp2 <- as.vector(nk_dcm.srt$name2)
nk_dcm.srt$tmp2[nk_dcm.srt$tmp2 == '5d PT'] <- as.vector(nk_dcm.srt$name4)[nk_dcm.srt$tmp2 == '5d PT']
df <- AverageExpression(nk_dcm.srt, features = genes, group.by = 'tmp2')$RNA
for(i in 1:nrow(df)){df[i,] <- scale(df[i,])}
gene_order <- rownames(df)[order(df[, '5d PT NK3'])]
df <- melt(df)
colnames(df) <- c('Gene', 'Group', 'Expr')
df$Gene <- factor(df$Gene, levels = gene_order)
p8 <- ggplot(df) +
        geom_tile(aes(x = Gene, y = Group, fill = Expr)) +
        scale_fill_distiller(palette = 'RdYlBu') +
        scale_y_discrete(limits = rev) +
        labs(x = 'Allograft rejection-related genes', y = '', fill = 'Expr') +
        theme_classic() +
        theme(aspect.ratio = 5/10) +
        RotatedAxis()
p8
PlotPDF('06.8.heat.expression_nk3_high_genes_vs_dcm', 8, 5)
print(p8)
dev.off()

####  Plot 9  vortex NK3 to cardiac signaling  ####
PlotPDF('06.9.vortex.p136_nk3_signaling', 8, 8)
netVisual_individual(p136_cc_list$All,
                     pairLR.use = c('FASL_FAS'),
                     signaling = c('FASLG'),
                     top = 0.035,
                     layout = "circle")
netVisual_individual(p136_cc_list$All,
                     pairLR.use = c('CSF2_CSF2RA_CSF2RB', 'IL18_IL18R1_IL18RAP'),
                     signaling = c('IL4', 'IL1'), top = 0.01,
                     layout = "circle")
dev.off()

####  Plot 10  vln MP1 signal to TC and Mono  ####
p10 <- plotGeneExpression(p136_cc_list$All,
                          signaling = c('IL4', 'IL1', 'FASLG'),
                          features = c('FASLG', 'FAS', 'CSF2', 'CSF2RB', 'IL18', 'IL18R1'))
p10
PlotPDF('06.10.vln.p136_nk3_signaling', 8, 8)
p10
dev.off()

####  Plot 11  umap resident lymphoid cells  ####
tmp.srt <- lym.srt[, lym.srt$name == '5d PT']
p11 <- DimPlot(tmp.srt, shuffle = T, reduction = 'sub_umap', group.by = 'Origin',
              order = c('Donor', 'Recipient'),
              cols = Color_genotype2[2:1], raster = F, pt.size = 0.5) &
        theme(aspect.ratio = 1, axis.text = element_blank(), axis.ticks = element_blank()) &
        labs(title = 'Donor and Recipient cell origins', subtitle = '(infered by Demuxlet)', x = 'UMAP1', y = 'UMAP2')
p11
PlotPDF('06.11.umap.p136_lymphoid_origin', 4, 4)
p11
dev.off()

####  Plot 12  vln resident NK1 markers  ####
ResNK_mk <- c('ITGAD', 'CLNK', 'SYK', 'HPGD', 'KLRC1', 'KLRC2', 'PARP8', 'GNLY', 'TXNIP', 'RBPJ')
tmp.srt <- lym.srt[, lym.srt$name == '5d PT' & lym.srt$Cell_state == 'NK1']
Idents(tmp.srt) <- 'Origin'
p12 <- VlnPlot2(tmp.srt, feature = ResNK_mk, pt.size = 1, same.y.lims = T, y.max = 5,
                cols = Color_genotype2[1:2], ncol = 5)
p12
PlotPDF('06.12.vln.p136_resnk_signature', 8, 4)
p12
dev.off()

####  Plot 13  rank resident NK1 markers score  ####
ResNK_mk <- c('ITGAD', 'CLNK', 'SYK', 'HPGD', 'KLRC1', 'KLRC2', 'PARP8', 'GNLY', 'TXNIP', 'RBPJ')
tmp.srt <- AddModuleScore2(lym.srt, features = list(ResNK_mk), names = 'ResNK', return_z = T)
tmp.srt <- tmp.srt[, tmp.srt$name2 == '5d PT' & tmp.srt$Cell_state %in% c('NK1', 'NK2', 'NK3')]
Idents(tmp.srt) <- 'Origin'
df <- tmp.srt@meta.data[, c('ResNK', 'Origin')]
df <- df[order(df$ResNK), ]
df$Rank <- 1:nrow(df)

p13 <- ggplot(df) +
        geom_point(aes(x = Rank,  y = ResNK, color = Origin), size = 0.1) +
        scale_color_manual(values = Color_genotype2[1:2]) +
        theme_classic() +
        theme(aspect.ratio = 0.5)
p13
PlotPDF('06.13.trend.p136_resnk_signature', 8, 4)
p13
dev.off()

####  Plot 14  feat resident NK1 markers score  ####
ResNK_mk <- c('ITGAD', 'CLNK', 'SYK', 'HPGD', 'KLRC1', 'KLRC2', 'PARP8', 'GNLY', 'TXNIP', 'RBPJ')
tmp.srt <- AddModuleScore2(lym.srt, features = list(ResNK_mk), names = 'ResNK', return_z = T)
tmp.srt <- tmp.srt[, tmp.srt$name2 == '5d PT']
p14 <- FeaturePlot2(tmp.srt, features = 'ResNK', pt.size = 1.5, reduction = 'sub_umap',
                    min.cutoff = -2, max.cutoff = 6)
p14
PlotPDF('06.14.feat.p136_resnk_signature', 4, 4)
p14
dev.off()

####  Plot 15  umap T/NK color by control and allograft time point ####
tmp.srt <- lym.srt
tmp.srt$tmp <- factor('Control', levels = c('Control', 'Allograft'))
tmp.srt$tmp[tmp.srt$name2 != 'Control'] <- 'Allograft'
p15 <- DimPlot2(tmp.srt, reduction = 'sub_umap', cols = mycol_10, ncol = 2, split.by = 'tmp') +
        labs(x = 'UMAP1', y = "UMAP2")
p15
PlotPDF('06.15.umap.control_t_nk_color_by_cell_state', 8, 4)
p15
dev.off()

####  Plot 16  box NK rejection signature vs DCM ####
gene1 <- read.table('external/HALLMARK_ALLOGRAFT_REJECTION.v7.5.1.grp', comment.char = '#', header = T)[,1]
gene2 <- read.table('external/KEGG_ALLOGRAFT_REJECTION.v7.5.1.grp', comment.char = '#', header = T)[,1]
gene <- union(gene1, gene2)
tmp.srt <- nk_dcm.srt[, nk_dcm.srt$Study == 'Lavine']
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

p16 <- p16.1 + p16.2 & theme(aspect.ratio = 4) & scale_y_continuous(limits = c(-3, 3))
p16
PlotPDF('06.16.box.nk_rejection_score_vs_dcm', 6, 4)
p16
dev.off()
####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~----


####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#### Test ####
####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~----
