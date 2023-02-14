####--------------------------------------------------------------------------------------------------------------------
####  Pediatric Acute Rejection
####  2022-08-03 by Xiao LI (Texas Heart Institute, US)
####--------------------------------------------------------------------------------------------------------------------


####--------------------------------------------------------------------------------------------------------------------
####  Initiate  ####
####--------------------------------------------------------------------------------------------------------------------
Ver <- '1'
Step <- 'STEP27_Analysis_NK'
Project <- '2022_acute_rejection_dturaga'
Code_dir <- paste0('~/Documents/Bioinformatics/project/', Project, '/code/human_v', Ver, '/')

source(Sys.readlink(paste0(Code_dir, 'src/bioinformatics.R')))
source(Sys.readlink(paste0(Code_dir, 'src/scRNAseq.R')))
source(Sys.readlink(paste0(Code_dir, 'src/geneset_human.R')))
source(paste0(Code_dir, 'human_v', Ver, '.helper_functions.R'))

InitiateProject('Rivendell', Ver, Step, 'human', Project, 'ithil')
####--------------------------------------------------------------------------------------------------------------------


####--------------------------------------------------------------------------------------------------------------------
####  Load data  ####
####--------------------------------------------------------------------------------------------------------------------
full.srt <- readRDS('integrated/STEP19.annotated_alra.srt.rds')

lym.srt <- full.srt[, full.srt$Cell_type %in% c('Lym') & full.srt$Cell_type_non_ambig & full.srt$Cell_state != 'BC']
lym.srt$Cell_type <- droplevels(lym.srt$Cell_type)
lym.srt$Cell_state <- droplevels(lym.srt$Cell_state)
Idents(lym.srt) <- 'Cell_state'
DimPlot2(lym.srt, reduction = 'sub_umap', group.by = 'Cell_state')
####--------------------------------------------------------------------------------------------------------------------

## Plot cell_state markers
mk <- FindAllMarkers(lym.srt, assay = 'CBN', test.use = 'MAST',
                     only.pos = T, logfc.threshold = 0.25, return.thresh = 0.001)
mk_dedup <- mk |> group_by(cluster) |> top_n(n = 20, wt = avg_log2FC)
mk_dedup <- mk_dedup[! mk_dedup$gene %in% mk_dedup$gene[duplicated(mk_dedup$gene)], ]
PlotPDF('1.1.heat.lym_cell_state3_marker', 10, 10)
MarkerHeatmap(lym.srt, marker.df = mk_dedup, disp.min = 0, top = 20, group.cols = mycol_10)
dev.off()

lym.srt$tmp <- factor('TC', levels = c('TC', 'NK'))
lym.srt$tmp[lym.srt$Cell_state %in% c('NK1', 'NK2', 'NK3')] <-  'NK'
Idents(lym.srt) <- 'tmp'
mk_con <- FindMarkers(lym.srt, assay = 'CBN', test.use = 'MAST',
                      ident.1 = 'TC', only.pos = T)
Idents(lym.srt) <- 'Cell_state'

PlotPDF('1.2.dot.lym_cell_state3_marker_select', 7, 4)
DotPlot2(lym.srt, group.by = 'Cell_state',
         features = c('CD2', 'CD5', 'BCL11B',
                      'CD4', 'CD40LG', 'IL7R',
                      'CD8B', 'TOX', 'THEMIS',
                      'KLRC2', 'FCRL3', 'KLRC4',
                      'GNLY', 'PRF1', 'KLRF1',
                      'GZMB', 'IFNG', 'REL'),
         assay = 'CBN', cols = c('grey85', mycol_BuGr[1]),
         idents = c('CD4 TC', 'CD8 TC', 'NK1', 'NK2', 'NK3'),
         col.min = 0)
dev.off()

df <- as.data.frame(CountCellBarPlot(lym.srt, group.var = 'Cell_state', stack.var = 'name2', percentage = T)$data)
rownames(df) <- paste(df$StackVar, df$GroupVar)
df2 <- data.frame(StackVar = rep(c('Control', '5d PT', '15m PT', '12y PT'), each = 5),
                  GroupVar = rep(c('CD4 TC', 'CD8 TC', 'NK1', 'NK2', 'NK3'), 4),
                  Count = rep(0, 20)
)
rownames(df2) <- paste(df2$StackVar, df2$GroupVar)
for(i in rownames(df)){df2[i, 'Count'] <- df[i, 'Count']}
df2$Count <- df2$Count + 1
df2$Fraction <- df2$Count/rep(Table(full.srt$name2), each = LU(df2$GroupVar))
x <- matrix(NA, 3, 5)
rownames(x) <- c('5d PT', '15m PT', '12y PT')
colnames(x) <- levels(lym.srt$Cell_state)
for(i in 1:3){
        group <- c('5d PT', '15m PT', '12y PT')[i]
        for(j in 1:5){
                celltype <- levels(lym.srt$Cell_state)[j]
                x[i, j] <- df2$Fraction[df2$StackVar == group & df2$GroupVar == celltype]/
                        df2$Fraction[df2$StackVar == 'Control' & df2$GroupVar == celltype]
        }

}

mtx_log <- log2(x)
df3 <- melt(mtx_log)
colnames(df3) <- c('Patient', 'Cell.Type', 'Log2FC')
p3 <- ggplot(df3, aes(x = Cell.Type, y = Log2FC, group = Cell.Type, color = Patient)) +
        #geom_boxplot(outlier.shape = NA, width = 0.5) +
        geom_point(aes(color = Patient), position = position_jitterdodge(jitter.width = 0)) +
        geom_hline(yintercept = c(0), color = 'red4', linetype = 'dashed') +
        #scale_color_manual(values = Color_donor) +
        theme_classic() +
        #scale_y_continuous(minor_breaks = seq(-3, 9, 0.5), breaks = seq(-3, 9, 1), limits = c(-3, 9)) +
        theme(aspect.ratio = 1, panel.grid.major.y = element_line()) +
        RotatedAxis() +
        labs(x = '', y = 'Log2 fold change (Transplant vs. 4 donors)',
             title = 'Immune cell composition change', color = 'Transplant vs.')
p3
PlotPDF('02.6.box.immune_cellstate3_foldchange_har_vs_donors', 4, 4)
print(p3)
dev.off()
####--------------------------------------------------------------------------------------------------------------------


####--------------------------------------------------------------------------------------------------------------------
####  NK toxicity genes  ####
####--------------------------------------------------------------------------------------------------------------------
library('ggtern')

grp <- read.table('external/KEGG_NATURAL_KILLER_CELL_MEDIATED_CYTOTOXICITY.v7.5.1.grp', comment.char = '#', header = T)

df <- AverageExpression(lym.srt[, lym.srt$donor_pub == 'HAR'],
                        features = grp$KEGG_NATURAL_KILLER_CELL_MEDIATED_CYTOTOXICITY)
df <- df$alra[, c('NK1', 'NK2', 'NK3')]
df <- as.data.frame(df)
df$gene <- rownames(df)
p <- ggtern(data = df, aes(NK1, NK3, NK2)) +
        geom_mask() +
        geom_point(fill="red",shape=21,size=1) +
        geom_text(data = df, aes(label = gene),
                                 color = 'black', size = 3) +
        theme_bw()
p

PlotPDF('02.7.ternary.nk_har', 10, 10)
print(p)
dev.off()
####--------------------------------------------------------------------------------------------------------------------


####--------------------------------------------------------------------------------------------------------------------
saveRDS(lym.srt, 'integrated/STEP21.lymphoid_cell_alra.srt.rds')
####--------------------------------------------------------------------------------------------------------------------
