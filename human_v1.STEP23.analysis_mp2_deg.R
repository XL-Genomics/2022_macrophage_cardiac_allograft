####--------------------------------------------------------------------------------------------------------------------
####  Pediatric Acute Rejection
####  2022-08-03 by Xiao LI (Texas Heart Institute, US)
####--------------------------------------------------------------------------------------------------------------------


####--------------------------------------------------------------------------------------------------------------------
####  Initiate  ####
####--------------------------------------------------------------------------------------------------------------------
Ver <- '1'
Step <- 'STEP23_Analysis_MP2_DEG'
Project <- '2022_acute_rejection_dturaga'
Code_dir <- paste0('~/Documents/Bioinformatics/project/', Project, '/code/human_v', Ver, '/')

source(Sys.readlink(paste0(Code_dir, 'src/bioinformatics.R')))
source(Sys.readlink(paste0(Code_dir, 'src/scRNAseq.R')))
source(Sys.readlink(paste0(Code_dir, 'src/geneset_human.R')))
source(paste0(Code_dir, 'human_v', Ver, '.helper_functions.R'))

InitiateProject('Rivendell', Ver, Step, 'human', Project, 'ithil')

library('ggupset')

mycol_genotype <- c('grey85', 'cyan3', 'deeppink2' ,'lightskyblue4')
mycol_group <- c('deeppink2', 'cyan3', 'lightskyblue4')
####--------------------------------------------------------------------------------------------------------------------


####--------------------------------------------------------------------------------------------------------------------
full.srt <- readRDS('integrated/STEP19.annotated_alra.srt.rds')
full.srt$Origin <- factor(full.srt$Fmux_genotype_clean, levels = c('Donor', 'Recipient', 'Ambiguous', 'Control'))
full.srt$Origin[is.na(full.srt$Origin)] <- 'Control'

mye.srt <- full.srt[, full.srt$Cell_type == 'Mye' &
                            full.srt$Origin != 'Ambiguous' &
                            full.srt$Cell_state != 'Mast' &
                            full.srt$Cell_type_non_ambig &
                            full.srt$name2 != 'Control']
mye.srt$Cell_state <- droplevels(mye.srt$Cell_state)
mye.srt$Origin <- droplevels(mye.srt$Origin)
mye.srt$Cell_state_orig <- as.vector(mye.srt$Cell_state)
mye.srt$Cell_state_orig[mye.srt$Origin == 'Donor' & mye.srt$Cell_state == 'MP2'] <- 'ResMP'
mye.srt$Cell_state_orig[mye.srt$Origin == 'Recipient' & mye.srt$Cell_state == 'MP2'] <- 'MoMP2'
mye.srt$Cell_state_orig <- factor(mye.srt$Cell_state_orig, levels = c(
        'ResMP',
        'MoMP2',
        'Mono',
        'MP1'))
Table(mye.srt$Cell_state_orig, mye.srt$name2)

mp2.srt <- mye.srt[, mye.srt$Cell_state == 'MP2' &
                           ! paste(mye.srt$Cell_state_orig, mye.srt$name2) %in% c('ResMP 15m PT',  'ResMP 12y PT')]
Table(mp2.srt$Cell_state_orig, mp2.srt$name2)
mp2.srt$Group <- factor(paste(mp2.srt$Cell_state_orig, mp2.srt$name2),
                        levels = c('ResMP 5d PT', 'MoMP2 5d PT', 'MoMP2 15m PT', 'MoMP2 12y PT'))
Table(mp2.srt$Group)
Idents(mp2.srt) <- 'Group'
####--------------------------------------------------------------------------------------------------------------------


####--------------------------------------------------------------------------------------------------------------------
####  4 way comparison ####
####--------------------------------------------------------------------------------------------------------------------
mp2.srt <- FindVariableFeatures(mp2.srt, nfeatures = 3000)
mp2.srt <- ScaleData(mp2.srt)
mp2.srt <- RunPCA(mp2.srt) |> RunHarmony(group.by.vars = 'name', assay.use = 'CBN')
p1 <- DimPlot2(mp2.srt, reduction = 'sub_umap', cols = mycol_10) /
        DimPlot2(mp2.srt, reduction = 'harmony', cols = mycol_10, dims = c(1, 3))
PlotPDF('1.1.umap.mp2_pc', 5, 5)
p1
dev.off()

df <- data.frame(mp2.srt@reductions$harmony@cell.embeddings[, c(1, 3)])
df$Group <- mp2.srt$Group
df <- df |> group_by(Group) |> summarize(mean_x = mean(harmony_1), mean_y = mean(harmony_3))
PlotPDF('1.2.tree.mp2_pc', 5, 5)
plot(hclust(dist(as.matrix(df[, 2:3]))))
dev.off()
####--------------------------------------------------------------------------------------------------------------------
saveRDS(mp2.srt@reductions$harmony, 'analysis/STEP23.mp2_harmonized_pc.srt_redu.rds')
####--------------------------------------------------------------------------------------------------------------------

deg <- FindAllMarkers(mp2.srt, assay = 'CBN', only.pos = T, test.use = 'MAST')
deg <- deg[deg$p_val_adj < 0.05, ]
deg_hi <- deg[deg$avg_log2FC >= 1, ]
Table(deg_hi$cluster)


deg_1 <- FindMarkers(mp2.srt, assay = 'CBN', only.pos = F, test.use = 'MAST',
                     ident.1 = 'MoMP2 5d PT', ident.2 = 'ResMP 5d PT')
deg_1$cluster <- 'MoMP2 5d PT'
deg_1$gene <- rownames(deg_1)

deg_2 <- FindMarkers(mp2.srt, assay = 'CBN', only.pos = F, test.use = 'MAST',
                     ident.1 = 'MoMP2 15m PT', ident.2 = 'ResMP 5d PT')
deg_2$cluster <- 'MoMP2 15m PT'
deg_2$gene <- rownames(deg_2)

deg_3 <- FindMarkers(mp2.srt, assay = 'CBN', only.pos = F, test.use = 'MAST',
                     ident.1 = 'MoMP2 12y PT', ident.2 = 'ResMP 5d PT')
deg_3$cluster <- 'MoMP2 12y PT'
deg_3$gene <- rownames(deg_3)

deg_all <- rbind(deg_1, deg_2, deg_3)
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
ggplot(deg_df, aes(x=Cluster)) +
        geom_bar() +
        scale_x_upset(order_by = 'degree')
####--------------------------------------------------------------------------------------------------------------------
saveRDS(deg, 'analysis/STEP23.mp2_4_way_compare.srt_mk.rds')
saveRDS(deg_all, 'analysis/STEP23.mp2_3_way_compare_with_resmp.srt_mk.rds')
####--------------------------------------------------------------------------------------------------------------------


####--------------------------------------------------------------------------------------------------------------------
####  2 way comparison ####
####--------------------------------------------------------------------------------------------------------------------
mp2_sub.srt <- mp2.srt[, mp2.srt$Group != 'MoMP2 5d PT']
mp2_sub.srt$Group <- revalue(mp2_sub.srt$Group, replace = c(
        'ResMP 5d PT' = 'ResMP',
        'MoMP2 15m PT' = 'Late MoMP2',
        'MoMP2 12y PT' = 'Late MoMP2'
))
mp2_sub.srt$Group <- factor(mp2_sub.srt$Group, levels = c('ResMP', 'Late MoMP2'))
Idents(mp2_sub.srt) <- 'Group'

deg <- FindMarkers(mp2_sub.srt, assay = 'CBN', test.use = 'MAST', ident.1 = 'Late MoMP2', logfc.threshold = 0)

PlotPDF('2.1.vol.late_momp2_vs_resmp', 8, 8)
MarkerVolcano(deg, label = T, min.log2FC = 0.25)
dev.off()

deg_hi <- deg[deg$p_val_adj < 0.05, ]
deg_hi$gene <- rownames(deg_hi)
deg_hi$cluster <- 'Late MoMP2'
deg_hi$cluster[deg_hi$avg_log2FC < 0] <- 'ResMP'

genelist <- split(deg_hi$gene, deg_hi$cluster)
str(genelist)

enrich <- ModuleEnrichment(module_list = genelist, human_or_mouse = 'human')
momp2_go <- enrich$GO$`GO_Late MoMP2`
resmp_go <- enrich$GO$GO_ResMP
momp2_pw <- enrich$Reactome$`Reactome_Late MoMP2`
resmp_pw <- enrich$Reactome$Reactome_ResMP

module_geneID <- lapply(genelist, function(gr) as.numeric(bitr(gr,
                                                               fromType = "SYMBOL",
                                                               toType = "ENTREZID",
                                                               OrgDb = "org.Hs.eg.db")$ENTREZID))
modules.kegg <- compareCluster(module_geneID, pvalueCutoff = 0.05, fun = 'enrichKEGG')
modules.kegg@compareClusterResult <- group_by(modules.kegg@compareClusterResult, Cluster)
modules.kegg@compareClusterResult <- filter(modules.kegg@compareClusterResult, p.adjust < 0.05)
PlotPDF('2.2.dot.late_momp2_vs_resmp_kegg', 8, 16)
dotplot(modules.kegg, showCategory = 50) &
        scale_color_distiller(palette = 'YlGnBu', direction = -1) &
        RotatedAxis() &
        theme(aspect.ratio = 2)
dev.off()

modules.pathway <- compareCluster(module_geneID, pvalueCutoff = 0.05, fun = 'enrichPathway')
modules.pathway@compareClusterResult <- group_by(modules.pathway@compareClusterResult, Cluster)
modules.pathway@compareClusterResult <- filter(modules.pathway@compareClusterResult, p.adjust < 0.05)
PlotPDF('2.3.dot.late_momp2_vs_resmp_reactome', 8, 16)
dotplot(modules.pathway, showCategory = 30) &
        scale_color_distiller(palette = 'YlGnBu', direction = -1) &
        RotatedAxis() &
        theme(aspect.ratio = 2)
dev.off()

####--------------------------------------------------------------------------------------------------------------------
saveRDS(deg, 'analysis/STEP23.mp2_resmp_late_momp2_deg.mk.rds')
saveRDS(enrich, 'analysis/STEP23.mp2_resmp_late_momp2_deg.enrich.rds')
####--------------------------------------------------------------------------------------------------------------------


deg <- FindMarkers(mp2_sub.srt, assay = 'alra', ident.1 = 'Late MoMP2', logfc.threshold = 0.1)

PlotPDF('2.2.vol.late_momp2_vs_resmp_alra', 8, 8)
MarkerVolcano(deg, label = T, min.log2FC = 0.25)
dev.off()

deg_hi <- deg[deg$p_val_adj < 1e-20 & abs(deg$avg_log2FC) > 1.5, ]
deg_hi$gene <- rownames(deg_hi)
deg_hi$cluster <- 'Late MoMP2'
deg_hi$cluster[deg_hi$avg_log2FC < 0] <- 'ResMP'

genelist <- split(deg_hi$gene, deg_hi$cluster)
str(genelist)

enrich <- ModuleEnrichment(module_list = genelist, human_or_mouse = 'human')
momp2_go <- enrich$GO$`GO_Late MoMP2`
resmp_go <- enrich$GO$GO_ResMP
momp2_pw <- enrich$Reactome$`Reactome_Late MoMP2`
resmp_pw <- enrich$Reactome$Reactome_ResMP
####--------------------------------------------------------------------------------------------------------------------
saveRDS(deg, 'analysis/STEP23.mp2_resmp_late_momp2_deg_alra.mk.rds')
saveRDS(enrich, 'analysis/STEP23.mp2_resmp_late_momp2_deg_alra.enrich.rds')
####--------------------------------------------------------------------------------------------------------------------



endoc <- read.table('external/KEGG_ENDOCYTOSIS.v7.5.1.grp',comment.char = '#', header = T)[,1]
phago <- read.table('external/GOBP_PHAGOCYTOSIS.v7.5.1.grp',comment.char = '#', header = T)[,1]
ttol <- read.table('external/GOBP_T_CELL_TOLERANCE_INDUCTION.v7.5.1.grp',comment.char = '#', header = T)[,1]
mp2_sub.srt <- AddModuleScore2(mp2_sub.srt,
                               features = c(cc.genes.updated.2019, mhc_genes, list(phago, endoc, ttol)), assay = 'CBN',
                        names = c('Score_S', 'Score_G2M', 'Score_MHC1', 'Score_MHC2',
                                  'Score_phago', 'Score_endoc', 'Score_ttol'),
                        return_z = T)

VlnPlot2(mp2_sub.srt,
         features = c('Score_S', 'Score_G2M', 'Score_MHC1', 'Score_MHC2', 'Score_phago', 'Score_endoc', 'Score_ttol'))


