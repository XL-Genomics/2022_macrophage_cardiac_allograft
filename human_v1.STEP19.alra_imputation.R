####--------------------------------------------------------------------------------------------------------------------
####  Pediatric Acute Rejection
####  2022-08-03 by Xiao LI (Texas Heart Institute, US)
####--------------------------------------------------------------------------------------------------------------------


####--------------------------------------------------------------------------------------------------------------------
####  Initiate  ####
####--------------------------------------------------------------------------------------------------------------------
Ver <- '1'
Step <- 'STEP19_ALRA_Imputation'
Project <- '2022_acute_rejection_dturaga'
Code_dir <- paste0('~/Documents/Bioinformatics/project/', Project, '/code/human_v', Ver, '/')

source(Sys.readlink(paste0(Code_dir, 'src/bioinformatics.R')))
source(Sys.readlink(paste0(Code_dir, 'src/scRNAseq.R')))
source(Sys.readlink(paste0(Code_dir, 'src/geneset_human.R')))
source(paste0(Code_dir, 'human_v', Ver, '.helper_functions.R'))

InitiateProject('Rivendell', Ver, Step, 'human', Project, 'ithil')
####--------------------------------------------------------------------------------------------------------------------


####--------------------------------------------------------------------------------------------------------------------
####  Impute  ####
####--------------------------------------------------------------------------------------------------------------------
srt <- readRDS('integrated/STEP18.annotated.srt.rds')
srt.list <- SplitObject(srt, split.by = 'sample')

## Impute for each sample
for(i in 1:L(srt.list)) {
        print(i)
        srt.list[[i]] <- RunALRA(srt.list[[i]], assay = 'CBN', genes.use = rownames(srt.list[[i]]))
}
for(i in 1:L(srt.list)) {
        srt.list[[i]] <- DietSeurat(srt.list[[i]])
}

tmp <- merge(srt.list[[1]], srt.list[2:L(srt.list)])
srt[['alra']] <- tmp[['alra']]

DefaultAssay(srt) <- 'CBN'
####--------------------------------------------------------------------------------------------------------------------
saveRDS(srt, 'integrated/STEP19.annotated_alra.srt.rds')
####--------------------------------------------------------------------------------------------------------------------
