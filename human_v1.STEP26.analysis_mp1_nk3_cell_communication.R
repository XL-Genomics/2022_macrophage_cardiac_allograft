####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
####  Pediatric Acute Rejection
####  2022-08-03 by Xiao LI (Texas Heart Institute, US)
####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
####  Initiate  ####
####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~----
Ver <- '1'
Step <- 'STEP26_Analysis_MP1_NK3_Cell_Communication'
Project <- '2022_acute_rejection_dturaga'
Code_dir <- paste0('~/Documents/Bioinformatics/project/', Project, '/code/human_v', Ver, '/')

source(Sys.readlink(paste0(Code_dir, 'src/bioinformatics.R')))
source(Sys.readlink(paste0(Code_dir, 'src/scRNAseq.R')))
source(Sys.readlink(paste0(Code_dir, 'src/geneset_human.R')))
source(paste0(Code_dir, 'human_v', Ver, '.helper_functions.R'))

InitiateProject('Rivendell', Ver, Step, 'human', Project, 'ithil')

suppressMessages(library('CellChat'))
####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
####  Load data ####
####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~----
full.srt <- readRDS('integrated/STEP19.annotated_alra.srt.rds')
Idents(full.srt) <- 'Cell_type'
p136.srt <- full.srt[, full.srt$name == '5d PT' &
                             full.srt$Cell_type_non_ambig &
                             full.srt$Fmux_genotype != 'Ambiguous']
cell_state_discard <- names(Table(p136.srt$Cell_state))[Table(p136.srt$Cell_state) < 20]
p136.srt <- p136.srt[, ! p136.srt$Cell_state %in% cell_state_discard]
p136.srt$Cell_type <- droplevels(p136.srt$Cell_type)
p136.srt$Cell_state <- droplevels(p136.srt$Cell_state)
Table(p136.srt$Cell_state, p136.srt$Cell_type)
####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
####  Prepare CellChat database  ####
####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~----
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
####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~----


####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
####  All cell states with T NK subtype ####
####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~----
p136.srt$Cell_state_simp <- factor(p136.srt$Cell_type,
                                   levels = c(levels(p136.srt$Cell_type),
                                              'BC', 'CD4 TC', 'CD8 TC', 'NK1', 'NK2', 'NK3', 'MP1', 'MP2', 'Mono'))
p136.srt$Cell_state_simp[p136.srt$Cell_type %in% c('Mye', 'Lym')] <-
        as.vector(p136.srt$Cell_state[p136.srt$Cell_type %in% c('Mye', 'Lym')])
p136.srt$Cell_state_simp <- droplevels(p136.srt$Cell_state_simp)
Table(p136.srt$Cell_state_simp)

p136.cch <- DoCellChat(p136.srt, group.by = 'Cell_state_simp', assay = 'alra')

PlotPDF('1.1.chrod.cellchat_imputated_by_cell_state', 20, 20)
for(i in 1:9){
        sender <- c('BC', 'CD4 TC', 'CD8 TC', 'NK1', 'NK2', 'NK3', 'MP1', 'MP2', 'Mono')
        message(sender[i])
        print(netVisual_chord_gene(p136.cch,
                             sources.use = sender[i],
                             targets.use = c('CM', 'FB', 'EC', 'PC', 'SMC', 'Neu',
                                             'BC', 'CD4 TC', 'CD8 TC', 'NK1', 'NK2', 'NK3', 'MP1', 'MP2', 'Mono'),
                             slot.name = "net",
                             scale = T,
                             link.border= T))
}
dev.off()
####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
saveRDS(p136.cch, 'analysis/STEP26.p136_alra_group_by_all_cell_states.cellchat.rds')
####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
####  All cell states with no T NK subtypes  ####
####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~----
p136.srt$Cell_state_simp2 <- factor(p136.srt$Cell_type,
                                    levels = c(levels(p136.srt$Cell_type), 'BC', 'TC', 'NK', 'MP1', 'MP2', 'Mono'))
p136.srt$Cell_state_simp2[p136.srt$Cell_type %in% c('Mye')] <-
        p136.srt$Cell_state[p136.srt$Cell_type %in% c('Mye')]
p136.srt$Cell_state_simp2[p136.srt$Cell_state %in% c('CD4 TC', 'CD8 TC')] <- 'TC'
p136.srt$Cell_state_simp2[p136.srt$Cell_state %in% c('NK1', 'NK2', 'NK3')] <- 'NK'
p136.srt$Cell_state_simp2[p136.srt$Cell_state %in% c('BC')] <- 'BC'
p136.srt$Cell_state_simp2 <- droplevels(p136.srt$Cell_state_simp2)
Table(p136.srt$Cell_state_simp2, p136.srt$Cell_state_simp)

p136_simp.cch <- DoCellChat(p136.srt, group.by = 'Cell_state_simp2', assay = 'alra')

PlotPDF('03.1.chrod.cellchat_imputated_by_cell_state', 20, 20)
for(i in 1:3){
        sender <- c('MP1', 'MP2', 'Mono')
        message(sender[i])
        netVisual_chord_gene(p136_simp.cch,
                             sources.use = sender[i],
                             targets.use = c('CM', 'FB', 'EC', 'PC', 'SMC', 'Neu',
                                             'BC', 'TC', 'NK', 'MP1', 'MP2', 'Mono'),
                             slot.name = "net",
                             scale = T,
                             link.border= T)
}
dev.off()
####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
saveRDS(p136_simp.cch, 'analysis/STEP26.p136_alra_group_by_all_cell_states_simplified.cellchat.rds')
####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
