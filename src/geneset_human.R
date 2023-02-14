####--------------------------------------------------------------------------------------------------------------------
####  General genesets  ####
####--------------------------------------------------------------------------------------------------------------------
mhc_genes <- list(MHC1 = as.vector(read.csv('~/Documents/Bioinformatics/r/db/human_mhc_gene.csv', header = T)[1:20, 1]),
                  MHC2 = as.vector(read.csv('~/Documents/Bioinformatics/r/db/human_mhc_gene.csv', header = T)[, 2]))
chrY_genes <- readRDS('~/Documents/Bioinformatics/r/db/human_chrY_genes.rds')
chrX_genes <- readRDS('~/Documents/Bioinformatics/r/db/human_chrX_genes.rds')
tf_genes <- read.csv('~/Documents/Bioinformatics/r/db/tf/human_tf.csv', header = T)[, 1]

####--------------------------------------------------------------------------------------------------------------------
####  Markersets from publication  ####
####--------------------------------------------------------------------------------------------------------------------
# Yap_target <- list(
#         Yap_target_cm_dev = as.vector(
#                 read.table('/Volumes/shire/data/other/yap_target_cm_atac_2019_DevCell.csv', header = T)[, 1]),
#         Yap_target_cm = as.vector(
#                 read.table('/Volumes/shire/data/other/yap_target_cm_scenic_Yuka_Unpub.txt', header = F)[, 1]),
#         Yap_target_cf = as.vector(
#                 read.table('/Volumes/shire/data/other/yap_target_cf_scenic_2019_GenesDev.txt', header = T)[, 1]),
#         Yap_target_yap5sa = as.vector(
#                 readRDS('~/Documents/Bioinformatics/r/db/mouse_yap_target_yap5sa.rds')),
#         Yap_target_yap5sa_union = as.vector(
#                 readRDS('~/Documents/Bioinformatics/r/db/mouse_yap_target_yap5sa_union.rds'))
# )
# Yap_target$Yap_target_cm[Yap_target$Yap_target_cm == 'Ptrf'] <- 'Cavin1' ## update symbol
# Yap_target <- list(
#         Yap_target_cm_dev = ConvertGeneSpecies(Yap_target$Yap_target_cm_dev, from = 'mouse', to = 'human'),
#         Yap_target_cm = ConvertGeneSpecies(Yap_target$Yap_target_cm, from = 'mouse', to = 'human'),
#         Yap_target_cf = ConvertGeneSpecies(Yap_target$Yap_target_cf, from = 'mouse', to = 'human'),
#         Yap_target_yap5sa = ConvertGeneSpecies(Yap_target$Yap_target_yap5sa, from = 'mouse', to = 'human'),
#         Yap_target_yap5sa_union = ConvertGeneSpecies(Yap_target$Yap_target_yap5sa_union, from = 'mouse', to = 'human')
# )
# saveRDS(Yap_target, file = '~/Documents/Bioinformatics/r/db/yap_target_human.list.rds')
Yap_target <- readRDS('~/Documents/Bioinformatics/r/db/yap_target_human.list.rds')
Yap_target$Yap_target_cf_cutrun <- readRDS(paste0('~/Documents/Bioinformatics/r/db/',
                                                  'cf_yap_target.latscko_yap1cutrun_based_on_genedev_paper.rds'))[[2]]


#### Neonatal MI induced regenerative CM markers, 2020 DevCell EOlson ####
# Neo_MI_CM_marker_mouse <- list(
#         CM1 = c('Gm30382', 'Fhl2', 'Fgf13', 'Mhrt', 'Ank2'),
#         CM2 = c('Casc5', 'Lockd', 'Top2a', 'Gm14091'),
#         CM3 = c('Mid1', 'Ddc', '1700042O10Rik'),
#         CM4 = c('Mb', 'Sod2', 'Myl3', 'Mdh2', 'Atp5b'),
#         CM5 = c('Enah', 'Ankrd1', 'Cd44', 'Gm13601', 'Xirp2')
# )
# Neo_MI_CM_marker <- list(
#         CM4 = ConvertGeneList(Neo_MI_CM_marker_mouse$CM4, from = 'mouse', to = 'human'),
#         CM5 = ConvertGeneList(Neo_MI_CM_marker_mouse$CM5, from = 'mouse', to = 'human')
# )
# saveRDS(Neo_MI_CM_marker, file = '~/Documents/Bioinformatics/r/db/Neo_MI_CM_marker_human.list.rds')
Neo_MI_CM_marker <- readRDS('~/Documents/Bioinformatics/r/db/Neo_MI_CM_marker_human.list.rds')


#### Ventricular YAP5SA CM C3ar1+Mac, C3+CF markers ####
# Yap5sa_CM_Mac_CF_marker_mouse <- list(
#         CM1 = c("Atp5j","Cox6b1","Atp5g1","Chchd10","Atp5c1","Cox7b","Ndufb10","Mpc2","Sdhb","Uqcrb","Ndufv2",
#                 "Etfa","Acadm","Ckmt2","Mgst3","Cox7a1","Mrpl42","Eci1","Phyh","Fabp3-ps1"),
#         CM2 = c("Flnc","Sorbs2","Itga7","Ahnak","Tgm2","Prnp","Tns1","Ddb1","Eef2","Rtn4","Serpinh1","Ctgf",
#                 "Jam2","Rras2","Ptrf","Lmcd1","Rock2","Acta2","Parm1"),
#         C3ar1_Mac = c('C3ar1', 'Cx3cr1', 'C5ar1', 'Apoe', 'Rbpj', 'F13a1'),
#         C3_FB = c('C3', 'Fstl1', 'Serpina3n', 'Cxcl14', 'Serping1')
# )
# Yap5sa_CM_Mac_CF_marker <- list(
#         CM1 = ConvertGeneList(Yap5sa_CM_Mac_CF_marker_mouse$CM1, from = 'mouse', to = 'human'),
#         CM2 = ConvertGeneList(Yap5sa_CM_Mac_CF_marker_mouse$CM2, from = 'mouse', to = 'human'),
#         C3ar1_Mac = ConvertGeneList(Yap5sa_CM_Mac_CF_marker_mouse$C3ar1_Mac, from = 'mouse', to = 'human'),
#         C3_FB = ConvertGeneList(Yap5sa_CM_Mac_CF_marker_mouse$C3_FB, from = 'mouse', to = 'human')
# )
# saveRDS(Yap5sa_CM_Mac_CF_marker, file = '~/Documents/Bioinformatics/r/db/Yap5sa_CM_Mac_CF_marker_human.list.rds')
Yap5sa_CM_Mac_CF_marker <- readRDS('~/Documents/Bioinformatics/r/db/Yap5sa_CM_Mac_CF_marker_human.list.rds')

#### CM Stress markers ####
CM_stress_genes <- c('NPPA', 'NPPB', 'MYH7', 'MYH7B', 'XIRP2', 'CMYA5', 'ANKRD1', 'TNNI3', 'ACTA1', 'PFKP')


#### ResMP MoMP markers -- K. Lavine's Nature Medicine Paper ####
Mac_marker <- list(
        TRMP_gene = c('CD163L1', 'MRC1', 'MAF', 'SIGLEC1', 'CCL8', 'CCL14', 'LILRB5', 'LYVE1', 'IL2RA',
                     'PDGFC', 'WLS', 'DAB2', 'NRP1', 'SCN9A', 'FGF13', 'GDF15', 'IGF1', 'FMOD', 'SLIT3',
                     'EGFL7', 'ECM1', 'SDC3'),
        MoMP_gene = c('CCL17', 'CCR5', 'CCR2', 'CXCL9', 'CXCL3', 'CXCL2', 'CXCL10', 'CSF2RA', 'CCL5', 'CCL20',
                       'TREM1', 'TET2', 'NFKBIE', 'IL1R1', 'IL1R2', 'NFKB1', 'REL', 'MAP2K3', 'IL1A', 'NLRP3',
                       'TRAF1', 'MAPK6', 'SRC', 'IL1B', 'NOD2', 'JAK3', 'IRAK2', 'MAP3K8', 'RELB', 'SOCS3', 'MYD88',
                       'MMP9', 'IL20', 'AREG', 'PTX3', 'IL23A', 'IL10', 'OSM', 'TIMP1', 'EREG', 'IL27')
)

#### ResMP MoMP markers -- S. Epelman's Science Immunology Paper ####
genes <- read_excel('~/Documents/Bioinformatics/r/db/sciimmunol.abf7777_table_s7.xlsx',
                    sheet = 3, col_names = T, skip = 1) ## These genes are from the E. Slava's Science Immun Paper
genes <- genes[genes$avg_logFC > 0.75 & genes$p_val_adj < 0.001, ]
Mac_marker_2 <- split(genes$gene, genes$cluster)

#### ResMP markers -- 5d postHTx sample ####
ResMac_marker <- list('ResMac_HTx' = readRDS(file = '~/Documents/Bioinformatics/r/db/ResMP_transplant_marker.rds'))


#### Human MI Fibrosis Markers -- Kuppe_et_al_Nature_MI_Fibro_marker_SuppData.xlsx  ####
Fibrosis_marker <- list(c(
        'POSTN', 'FN1', 'TNC', 'COL1A1', 'COL1A2', 'COL3A1', 'DEC1', 'RUNX1', 'ADAM12', 'KIF26B',
        'THBS4', 'FAP', 'COL5A1', 'SERPINE1', 'SLC20A1', 'KALRN', 'PRICKLE1', 'CDH11', 'FGF14', 'THBS2'
))
