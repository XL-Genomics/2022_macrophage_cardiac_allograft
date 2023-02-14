####--------------------------------------------------------------------------------------------------------------------
####  Pediatric Acute Rejection
####  2022-08-03 by Xiao LI (Texas Heart Institute, US)
####--------------------------------------------------------------------------------------------------------------------


####--------------------------------------------------------------------------------------------------------------------
####  Initiate  ####
####--------------------------------------------------------------------------------------------------------------------
Ver <- '1'
Step <- 'STEP24_Analysis_ResMP_MoMP2_Cell_Communication'
Project <- '2022_acute_rejection_dturaga'
Code_dir <- paste0('~/Documents/Bioinformatics/project/', Project, '/code/human_v', Ver, '/')

source(Sys.readlink(paste0(Code_dir, 'src/bioinformatics.R')))
source(Sys.readlink(paste0(Code_dir, 'src/scRNAseq.R')))
source(Sys.readlink(paste0(Code_dir, 'src/geneset_human.R')))
source(paste0(Code_dir, 'human_v', Ver, '.helper_functions.R'))

InitiateProject('Rivendell', Ver, Step, 'human', Project, 'ithil')

suppressMessages(library('CellChat'))
# library('nichenetr')
# select <- dplyr::select

Color_cell_type <- c(mycol_40[c(3, 6, 11, 14, 35)], 'grey')
####--------------------------------------------------------------------------------------------------------------------


####--------------------------------------------------------------------------------------------------------------------
####  Prepare seurat  ####
####--------------------------------------------------------------------------------------------------------------------
full.srt <- readRDS('integrated/STEP19.annotated_alra.srt.rds')
full.srt$Origin <- factor(full.srt$Fmux_genotype_clean, levels = c('Donor', 'Recipient', 'Ambiguous', 'Control cells'))
full.srt$Origin[is.na(full.srt$Origin)] <- 'Control cells'

srt <- full.srt[, full.srt$name2 != 'Control' & full.srt$Cell_type_non_ambig & full.srt$Origin != 'Ambiguous']
srt <- srt[, srt$Cell_type %in% c('CM', 'FB', 'EC', 'PC', 'SMC', 'Mye')]
srt <- srt[, srt$Cell_state != 'Mast']
srt$Cell_type <- droplevels(srt$Cell_type)
srt$Cell_state <- droplevels(srt$Cell_state)

srt$Group1 <- as.vector(srt$Cell_type)
srt$Group1[srt$Cell_type == 'Mye'] <- paste(as.vector(srt$Cell_state),
                                            srt$Origin,
                                            srt$name)[srt$Cell_type == 'Mye']
# mye.srt <- srt[, ! srt$Group1 %in% c('CM', 'EC', 'FB', 'PC', 'SMC')]
# mye.srt <- mye.srt[, mye.srt$Group1 %in% names(table(mye.srt$Group1))[table(mye.srt$Group1)>100]]
# Idents(mye.srt) <- 'Group1'
# mk <- FindAllMarkers(mye.srt, only.pos = T)
# mk <- mk[order(mk$avg_log2FC, decreasing = T),]
# ResMP_marker <- mk$gene[mk$cluster == 'MP2 Donor 5d PT' & mk$p_val_adj <= 0.01 & mk$avg_log2FC >= 0.5]
# saveRDS(ResMP_marker, file = '~/Documents/Bioinformatics/r/db/ResMP_transplant_marker.rds')

srt <- srt[, srt$Group1 %in% c('CM', 'FB', 'EC', 'PC', 'SMC',
                               'MP2 Donor 5d PT', 'MP2 Recipient 12y PT', 'MP2 Recipient 15m PT')]
srt$Group1 <- revalue(srt$Group1, replace = c('MP2 Donor 5d PT' = 'MP2',
                                              'MP2 Recipient 12y PT' = 'MP2',
                                              'MP2 Recipient 15m PT' = 'MP2'))
srt$Group1 <- factor(srt$Group1, levels = c('CM', 'FB', 'EC', 'PC', 'SMC', 'MP2'))
Table(srt$Group1, srt$Cell_type)

p136.srt <- srt[, srt$name == '5d PT']
p136.srt$Group1 <- droplevels(p136.srt$Group1)
Table(p136.srt$Group1, p136.srt$Cell_type)
p182.srt <- srt[, srt$name == '15m PT']
p182.srt$Group1 <- droplevels(p182.srt$Group1)
Table(p182.srt$Group1, p182.srt$Cell_type)
p170.srt <- srt[, srt$name == '12y PT']
p170.srt$Group1 <- droplevels(p170.srt$Group1)
Table(p170.srt$Group1, p170.srt$Cell_type)
####--------------------------------------------------------------------------------------------------------------------


####--------------------------------------------------------------------------------------------------------------------
####  Prepare CellChat database  ####
####--------------------------------------------------------------------------------------------------------------------
## Adding FGF13 interactions per NicheNet Database
CellChatDB.my <- CellChatDB.human
last <- nrow(CellChatDB.my$interaction)
CellChatDB.my$interaction[(last+1):(last+4), ] <- NA
CellChatDB.my$interaction$interaction_name[(last+1):(last+4)] <- c('FGF13_FGFR1',
                                                                   'FGF13_FGFR2',
                                                                   'FGF13_FGFR3',
                                                                   'FGF13_FGFR4')
rownames(CellChatDB.my$interaction)[(last+1):(last+4)] <- c('FGF13_FGFR1',
                                                            'FGF13_FGFR2',
                                                            'FGF13_FGFR3',
                                                            'FGF13_FGFR4')
CellChatDB.my$interaction$pathway_name[(last+1):(last+4)] <- 'FGF'
CellChatDB.my$interaction$ligand[(last+1):(last+4)] <- 'FGF13'
CellChatDB.my$interaction$receptor[(last+1):(last+4)] <- c('FGFR1','FGFR2','FGFR3','FGFR4')
CellChatDB.my$interaction$agonist[(last+1):(last+4)] <- ''
CellChatDB.my$interaction$antagonist[(last+1):(last+4)] <- ''
CellChatDB.my$interaction$co_A_receptor[(last+1):(last+4)] <- ''
CellChatDB.my$interaction$co_I_receptor[(last+1):(last+4)] <- ''
CellChatDB.my$interaction$evidence[(last+1):(last+4)] <- ''
CellChatDB.my$interaction$annotation[(last+1):(last+4)] <- 'Secreted Signaling'
CellChatDB.my$interaction$interaction_name_2[(last+1):(last+4)] <- c('FGF13 - FGFR1',
                                                                     'FGF13 - FGFR2',
                                                                     'FGF13 - FGFR3',
                                                                     'FGF13 - FGFR4')
DoCellChat <- function(srt, group.by = "Group", assay){
        CellChatDB <- CellChatDB.my
        CellChatDB <- subsetDB(CellChatDB, search = "Secreted Signaling", key = "annotation")
        ####  Build CellChat Object  ####
        cch <- createCellChat(object = srt, group.by = group.by, assay = assay)
        cch <- addMeta(cch, meta = srt@meta.data)
        groupSize <- as.numeric(table(cch@idents))
        cch@DB <- CellChatDB
        ####  Preprocessing the expression data for cell-cell communication analysis  ####
        cch <- subsetData(cch) # subset the expression data of signaling genes for saving computation cost
        cch <- identifyOverExpressedGenes(cch)
        cch <- identifyOverExpressedInteractions(cch)
        cch <- projectData(cch, PPI.human)
        ####  Inference of cell-cell communication network  ####
        cch <- computeCommunProb(cch)
        # Filter out the cell-cell communication if there are only few number of cells in certain cell groups
        cch <- filterCommunication(cch)
        # CellChat computes the communication probability on signaling pathway level
        cch <- computeCommunProbPathway(cch)
        cch <- aggregateNet(cch)
        return(cch)
}
####--------------------------------------------------------------------------------------------------------------------


# ####--------------------------------------------------------------------------------------------------------------------
# ####  CellChat workflow for CBN  ####
# ####--------------------------------------------------------------------------------------------------------------------
# p136.cch <- DoCellChat(p136.srt, group.by = 'Group1', assay = 'CBN')
# p182.cch <- DoCellChat(p182.srt, group.by = 'Group1', assay = 'CBN')
# p170.cch <- DoCellChat(p170.srt, group.by = 'Group1', assay = 'CBN')
#
# cch13.list <- list(p136 = p136.cch, p170 = p170.cch)
# merge13.cch <- mergeCellChat(cch13.list, add.names = names(cch13.list))
# PlotPDF('01.1.dot.cbn_cellchat_mp2_to_cardiac_p136_vs_p170', 8, 10)
# for(i in 1:2){
#         ct <-  c('PC', 'SMC', 'MP2')[i]
#         p1 <- netVisual_bubble(merge13.cch,
#                          comparison = c(1, 2),
#                          sources.use = 'MP2',
#                          targets.use = ct,
#                          max.dataset = 2,
#                          title.name = "Increased signaling in CR",
#                          angle.x = 45,
#                          remove.isolate = T)
#         p2 <- netVisual_bubble(merge13.cch,
#                                comparison = c(1, 2),
#                                sources.use = 'MP2',
#                                targets.use = ct,
#                                max.dataset = 1,
#                                title.name = "Increased signaling in PGF",
#                                angle.x = 45,
#                                remove.isolate = T)
#         print(p1 | p2)
# }
# dev.off()
#
# cch12.list <- list(p136 = p136.cch, p182 = p182.cch)
# merge12.cch <- mergeCellChat(cch12.list, add.names = names(cch12.list))
# PlotPDF('01.2.dot.cbn_cellchat_mp2_to_cardiac_p136_vs_p182', 8, 10)
# for(i in 1:5){
#         ct <-  c('FB', 'MP2')[i]
#         p1 <- netVisual_bubble(merge12.cch,
#                                comparison = c(1, 2),
#                                sources.use = 'MP2',
#                                targets.use = ct,
#                                max.dataset = 2,
#                                title.name = "Increased signaling in CR",
#                                angle.x = 45,
#                                remove.isolate = T)
#         p2 <- netVisual_bubble(merge12.cch,
#                                comparison = c(1, 2),
#                                sources.use = 'MP2',
#                                targets.use = ct,
#                                max.dataset = 1,
#                                title.name = "Increased signaling in PGF",
#                                angle.x = 45,
#                                remove.isolate = T)
#         print(p1 | p2)
# }
# dev.off()
# ####--------------------------------------------------------------------------------------------------------------------
# saveRDS(p136.cch, 'analysis/STEP24.cbn_cellchat.p136_cch.rds')
# saveRDS(p182.cch, 'analysis/STEP24.cbn_cellchat.p182_cch.rds')
# saveRDS(p170.cch, 'analysis/STEP24.cbn_cellchat.p170_cch.rds')
# ####--------------------------------------------------------------------------------------------------------------------


####--------------------------------------------------------------------------------------------------------------------
####  CellChat workflow for ALRA imputed  ####
####--------------------------------------------------------------------------------------------------------------------
p136_alra.cch <- DoCellChat(p136.srt, group.by = 'Group1', assay = 'alra')
p170_alra.cch <- DoCellChat(p170.srt, group.by = 'Group1', assay = 'alra')
p182_alra.cch <- DoCellChat(p182.srt, group.by = 'Group1', assay = 'alra')

cch12.list <- list(p136 = p136_alra.cch, p182 = p182_alra.cch)
cch13.list <- list(p136 = p136_alra.cch, p170 = p170_alra.cch)
merge12.cch <- mergeCellChat(cch12.list, add.names = names(cch12.list))
merge13.cch <- mergeCellChat(cch13.list, add.names = names(cch13.list))
####--------------------------------------------------------------------------------------------------------------------
saveRDS(p136_alra.cch, 'analysis/STEP24.p136_alra_cellchat.cch.rds')
saveRDS(p182_alra.cch, 'analysis/STEP24.p182_alra_cellchat.cch.rds')
saveRDS(p170_alra.cch, 'analysis/STEP24.p170_alra_cellchat.cch.rds')
####--------------------------------------------------------------------------------------------------------------------
## Filter interactions by DE ligands
pos.dataset = "p170"
features.name = pos.dataset
merge13.cch <- identifyOverExpressedGenes(merge13.cch,
                                        group.dataset = "datasets",
                                        pos.dataset = pos.dataset,
                                        features.name = features.name,
                                        only.pos = FALSE,
                                        thresh.pc = 0.1,
                                        thresh.fc = 0.1,
                                        thresh.p = 0.05)
## map the results of differential expression analysis onto the inferred cell-cell communications to
## easily manage/subset the ligand-receptor pairs of interest
net13 <- netMappingDEG(merge13.cch, features.name = features.name)
net13.up <- subsetCommunication(merge13.cch, net = net13, datasets = "p170", ligand.logFC = 0.1,  receptor.logFC = NULL)
net13.dn <- subsetCommunication(merge13.cch, net = net13, datasets = "p136", ligand.logFC = -0.1, receptor.logFC = NULL)

pos.dataset = "p182"
features.name = pos.dataset
merge12.cch <- identifyOverExpressedGenes(merge12.cch,
                                          group.dataset = "datasets",
                                          pos.dataset = pos.dataset,
                                          features.name = features.name,
                                          only.pos = FALSE,
                                          thresh.pc = 0.1,
                                          thresh.fc = 0.1,
                                          thresh.p = 0.05)
net12 <- netMappingDEG(merge12.cch, features.name = features.name)
net12.up <- subsetCommunication(merge12.cch, net = net12, datasets = "p182", ligand.logFC = 0.1,  receptor.logFC = NULL)
net12.dn <- subsetCommunication(merge12.cch, net = net12, datasets = "p136", ligand.logFC = -0.1, receptor.logFC = NULL)

## 2nd round of MAST DEG filtering
tmp <- srt[, srt$Cell_type == 'Mye']
Idents(tmp) <- 'name'
deg <- FindMarkers(tmp, ident.1 = '5d PT', test.use = 'MAST', assay = 'alra')
deg_up <- rownames(deg)[deg$avg_log2FC < 0.1 & deg$p_val_adj < 0.05]
L(deg_up) ## 4045
deg_dn <- rownames(deg)[deg$avg_log2FC > -0.1 & deg$p_val_adj < 0.05]
L(deg_dn) ## 4781

gene12.up <- extractGeneSubsetFromPair(net12.up, merge12.cch)
L(intersect(gene12.up, deg_up))  ## 44
gene12.dn <- extractGeneSubsetFromPair(net12.dn, merge12.cch)
L(intersect(gene12.dn, deg_dn))  ## 39
net12.up <- net12.up[net12.up$ligand %in% intersect(gene12.up, deg_up), ]
net12.dn <- net12.dn[net12.dn$ligand %in% intersect(gene12.dn, deg_dn), ]
pairLR.use.up_12 = net12.up[, "interaction_name", drop = F]
nrow(pairLR.use.up_12) ## 207
pairLR.use.dn_12 = net12.dn[, "interaction_name", drop = F]
nrow(pairLR.use.dn_12) ## 59

gene13.up <- extractGeneSubsetFromPair(net13.up, merge13.cch)
L(intersect(gene13.up, deg_up))  ## 33
gene13.dn <- extractGeneSubsetFromPair(net13.dn, merge13.cch)
L(intersect(gene13.dn, deg_dn))  ## 42
net13.up <- net13.up[net13.up$ligand %in% intersect(gene13.up, deg_up), ]
net13.dn <- net13.dn[net13.dn$ligand %in% intersect(gene13.dn, deg_dn), ]
pairLR.use.up_13 = net13.up[, "interaction_name", drop = F]
nrow(pairLR.use.up_13) ## 123
pairLR.use.dn_13 = net13.dn[, "interaction_name", drop = F]
nrow(pairLR.use.dn_13) ## 75

PlotPDF('02.1.dot.alra_cellchat_mp2_to_cardiac_p182_v_p136', 8, 10)
for(i in 1:5){
        ct <-  c('CM', 'FB', 'EC', 'PC', 'SMC', 'MP2')[i]
        p1 <- netVisual_bubble(merge12.cch,
                               comparison = c(1, 2),
                               sources.use = 'MP2',
                               targets.use = ct,
                               pairLR.use = pairLR.use.up_12,
                               max.dataset = 2,
                               title.name = "Increased signaling in p182",
                               angle.x = 45,
                               remove.isolate = T)
        p2 <- netVisual_bubble(merge12.cch,
                               comparison = c(1, 2),
                               sources.use = 'MP2',
                               targets.use = ct,
                               pairLR.use = pairLR.use.dn_12,
                               max.dataset = 1,
                               title.name = "Increased signaling in p136",
                               angle.x = 45,
                               remove.isolate = T)
        print(p1 | p2)
}
dev.off()

PlotPDF('02.2.dot.alra_cellchat_mp2_to_cardiac_p170_v_p136', 8, 10)
for(i in 1:5){
        ct <-  c('CM', 'FB', 'EC', 'PC', 'SMC', 'MP2')[i]
        p1 <- netVisual_bubble(merge13.cch,
                               comparison = c(1, 2),
                               sources.use = 'MP2',
                               targets.use = ct,
                               pairLR.use = pairLR.use.up_13,
                               max.dataset = 2,
                               title.name = "Increased signaling in p170",
                               angle.x = 45,
                               remove.isolate = T)
        p2 <- netVisual_bubble(merge13.cch,
                               comparison = c(1, 2),
                               sources.use = 'MP2',
                               targets.use = ct,
                               pairLR.use = pairLR.use.dn_13,
                               max.dataset = 1,
                               title.name = "Increased signaling in p136",
                               angle.x = 45,
                               remove.isolate = T)
        print(p1 | p2)
}
dev.off()

pairLR.use.up_12_13 <- data.frame(interaction_name = intersect(pairLR.use.up_12$interaction_name,
                                                               pairLR.use.up_13$interaction_name))
nrow(pairLR.use.up_12_13) ## 29
pairLR.use.dn_12_13 <- data.frame(interaction_name = intersect(pairLR.use.dn_12$interaction_name,
                                                               pairLR.use.dn_13$interaction_name))
nrow(pairLR.use.dn_12_13) ## 24

PlotPDF('02.3.chord.alra_cellchat_mp2_to_all_cardiac_p182_v_p136_common', 8, 6)
netVisual_chord_gene(cch12.list[[2]], sources.use = 'MP2', targets.use = c('CM', 'EC', 'FB', 'PC', 'SMC'),
                     slot.name = 'netP', pairLR.use = pairLR.use.up_12_13, color.use = Color_cell_type,
                     scale = T, lab.cex = 0.8, small.gap = 3.5, reduce = 0.005,
                     title.name = paste0("Up-regulated signaling in ", names(cch12.list)[2]))
netVisual_chord_gene(cch12.list[[1]], sources.use = 'MP2', targets.use = c('CM', 'EC', 'FB', 'PC', 'SMC'),
                     slot.name = 'netP', pairLR.use = pairLR.use.dn_12_13, color.use = Color_cell_type,
                     scale = T, lab.cex = 0.8, small.gap = 3.5, reduce = 0.005,
                     title.name = paste0("Up-regulated signaling in ", names(cch12.list)[1]))
dev.off()

PlotPDF('02.4.chord.alra_cellchat_mp2_to_all_cardiac_p170_v_p136_common', 8, 6)
netVisual_chord_gene(cch13.list[[2]], sources.use = 'MP2', targets.use = c('CM', 'EC', 'FB', 'PC', 'SMC'),
                     slot.name = 'netP', pairLR.use = pairLR.use.up_12_13, color.use = Color_cell_type,
                     scale = T, lab.cex = 0.8, small.gap = 3.5, reduce = 0.005,
                     title.name = paste0("Up-regulated signaling in ", names(cch13.list)[2]))
netVisual_chord_gene(cch13.list[[1]], sources.use = 'MP2', targets.use = c('CM', 'EC', 'FB', 'PC', 'SMC'),
                     slot.name = 'netP', pairLR.use = pairLR.use.dn_12_13, color.use = Color_cell_type,
                     scale = T, lab.cex = 0.8, small.gap = 3.5, reduce = 0.005,
                     title.name = paste0("Up-regulated signaling in ", names(cch13.list)[1]))
dev.off()


hi_light_PGF <- list(
        'FGF13_FGFR1' = 'FGF',
        'HBEGF_EGFR_ERBB2' = 'EGF',
        'VEGFB_VEGFR1' = 'VEGF'
)

PlotPDF('02.5.chord.cellchat_pgf_resmp_to_all_cardiac_highligh', 4, 4)
for(i in 1:3){
        print(netVisual_individual(p136_alra.cch,
                                   sources.use = 'MP2',
                                   targets.use = c('CM', 'EC', 'FB', 'PC', 'SMC'),
                                   pairLR.use = names(hi_light_PGF)[i],
                                   signaling = hi_light_PGF[[i]],
                                   color.use = Color_cell_type,
                                   layout = "circle"))
}
dev.off()

hi_light_CR <- list(
        'PDGFB_PDGFRA' = 'PDGF',
        'PDGFB_PDGFRB' = 'PDGF',
        'WNT5B_FZD6' = 'ncWNT',
        'TGFB2_TGFBR1_TGFBR2' = 'TGFb',
        'IGF1_IGF1R' = 'IGF',
        'TGFA_EGFR' = 'EGF'
        #'C3_C3AR1' = 'COMPLEMENT'
)
PlotPDF('02.6.chord.cellchat_p182_resmp_to_all_cardiac_highligh', 4, 4)
for(i in 1:6){
        print(netVisual_individual(p182_alra.cch,
                                   sources.use = 'MP2',
                                   targets.use = c('CM', 'EC', 'FB', 'PC', 'SMC'),
                                   pairLR.use = names(hi_light_CR)[i],
                                   signaling = hi_light_CR[[i]],
                                   color.use = Color_cell_type,
                                   layout = "circle"))
}
dev.off()

PlotPDF('02.7.chord.cellchat_p170_resmp_to_all_cardiac_highligh', 4, 4)
for(i in 1:6){
        print(netVisual_individual(p170_alra.cch,
                                   sources.use = 'MP2',
                                   targets.use = c('CM', 'EC', 'FB', 'PC', 'SMC'),
                                   pairLR.use = names(hi_light_CR)[i],
                                   signaling = hi_light_CR[[i]],
                                   color.use = Color_cell_type,
                                   layout = "circle"))
}
dev.off()

pair_use <- list(res_pairLR = pairLR.use.dn_12_13,
                 momp2_pairLR = pairLR.use.up_12_13)
####--------------------------------------------------------------------------------------------------------------------
saveRDS(merge12.cch, 'analysis/STEP24.p136_p182_alra_merge_cellchat.cch.rds')
saveRDS(merge13.cch, 'analysis/STEP24.p136_p170_alra_merge_cellchat.cch.rds')
saveRDS(pair_use, 'analysis/STEP24.final_cellchat_lr_pairs.list.rds')
####--------------------------------------------------------------------------------------------------------------------


####--------------------------------------------------------------------------------------------------------------------
####  Check lig rec expression in scRNA-seq  ####
####--------------------------------------------------------------------------------------------------------------------
p182_lig <- c()
p182_rec <- c()
for(i in c('CM', 'EC', 'FB', 'PC', 'SMC')){
        rst <- netVisual_bubble(merge12.cch,
                                comparison = c(1, 2),
                                sources.use = 'MP2',
                                targets.use = i,
                                pairLR.use = pairLR.use.up_12_13,
                                return.data = T,
                                max.dataset = 2,
                                remove.isolate = T)
        p182_lig <- c(p182_lig, rst$communication$ligand)
        p182_rec <- c(p182_rec, rst$communication$receptor)
}
p182_lig <- U(p182_lig)
p182_rec <- U(as.vector(str_split(U(p182_rec), pattern = '_', simplify = T)))
p182_rec <- revalue(p182_rec, replace = c('TGFbR1' = 'TGFBR1', 'TGFbR2' = 'TGFBR2', 'TGFbR' = 'TGFBR'))
p182_rec <- intersect(rownames(p182.srt), p182_rec)

p170_lig <- c()
p170_rec <- c()
for(i in c('CM', 'EC', 'FB', 'PC', 'SMC')){
        rst <- netVisual_bubble(merge13.cch,
                                comparison = c(1, 2),
                                sources.use = 'MP2',
                                targets.use = i,
                                pairLR.use = pairLR.use.up_12_13,
                                return.data = T,
                                max.dataset = 2,
                                remove.isolate = T)
        p170_lig <- c(p170_lig, rst$communication$ligand)
        p170_rec <- c(p170_rec, rst$communication$receptor)
}
p170_lig <- U(p170_lig)
p170_rec <- U(as.vector(str_split(U(p170_rec), pattern = '_', simplify = T)))
p170_rec <- revalue(p170_rec, replace = c('TGFbR1' = 'TGFBR1', 'TGFbR2' = 'TGFBR2', 'TGFbR' = 'TGFBR'))
p170_rec <- intersect(rownames(p170.srt), p170_rec)

common_lig <- intersect(p170_lig, p182_lig)
common_rec <- intersect(p170_rec, p182_rec)

p136_lig <- c()
p136_rec <- c()
for(i in c('CM', 'EC', 'FB', 'PC', 'SMC')){
        rst <- netVisual_bubble(merge12.cch,
                                comparison = c(1, 2),
                                sources.use = 'MP2',
                                targets.use = i,
                                pairLR.use = pairLR.use.dn_12_13,
                                return.data = T,
                                max.dataset = 1,
                                remove.isolate = T)
        p136_lig <- c(p136_lig, rst$communication$ligand)
        p136_rec <- c(p136_rec, rst$communication$receptor)
}
p136_lig <- U(p136_lig)
p136_rec <- U(as.vector(str_split(U(p136_rec), pattern = '_', simplify = T)))
p136_rec <- intersect(rownames(p136.srt), p136_rec)

## Check ligands expression
srt$Group2 <- paste(srt$Group1, srt$name2)
PlotPDF('03.1.dot.cellchat_mp2_to_cardiac_ligand', 24, 8)
DotPlot2(srt, features = c(p136_lig), group.by = 'Group2', assay = 'alra') |
        DotPlot2(srt, features = c(common_lig), group.by = 'Group2', assay = 'alra') |
        DotPlot2(srt, features = c(p136_lig), group.by = 'Group2', assay = 'CBN') |
        DotPlot2(srt, features = c(common_lig), group.by = 'Group2', assay = 'CBN')
dev.off()

PlotPDF('03.2.dot.cellchat_mp2_to_cardiac_receptor', 40, 8)
DotPlot2(srt, features = c(p136_rec), group.by = 'Group2', assay = 'alra') |
        DotPlot2(srt, features = c(common_rec), group.by = 'Group2', assay = 'alra') |
        DotPlot2(srt, features = c(p136_rec), group.by = 'Group2', assay = 'CBN') |
        DotPlot2(srt, features = c(common_rec), group.by = 'Group2', assay = 'CBN')
dev.off()
####--------------------------------------------------------------------------------------------------------------------
