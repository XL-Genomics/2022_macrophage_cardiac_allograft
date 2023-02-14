####--------------------------------------------------------------------------------------------------------------------
####  Pediatric Acute Rejection
####  2022-08-03 by Xiao LI (Texas Heart Institute, US)
####--------------------------------------------------------------------------------------------------------------------


####--------------------------------------------------------------------------------------------------------------------
####  Initiate  ####
####--------------------------------------------------------------------------------------------------------------------
Ver <- '1'
Step <- 'STEP01_Data_Collection'
Project <- '2022_acute_rejection_dturaga'
Code_dir <- paste0('~/Documents/Bioinformatics/project/', Project, '/code/human_v', Ver, '/')

source(Sys.readlink(paste0(Code_dir, 'src/bioinformatics.R')))
source(Sys.readlink(paste0(Code_dir, 'src/scRNAseq.R')))
source(Sys.readlink(paste0(Code_dir, 'src/geneset_human.R')))
source(paste0(Code_dir, 'human_v', Ver, '.helper_functions.R'))

InitiateProject('Rivendell', Ver, Step, 'human', Project, 'ithil')
####--------------------------------------------------------------------------------------------------------------------


####--------------------------------------------------------------------------------------------------------------------
####  Global Functions  ####
####--------------------------------------------------------------------------------------------------------------------
MakeSrt <- function(mode, matrix_dir, sample, study, method, platform, protocol, data_process, tissue, enrichment,
                    preparation, procedure, diagnosis, recipient, sex_recipient, age_recipient, sex_donor, age_donor,
                    duration, replicate) {
        if (mode == '10x') {
                matrix <- paste0(matrix_dir, '/outs/filtered_feature_bc_matrix/')
                srt <- CreateSeuratObject(counts = Read10X(data.dir = matrix),
                                          min.cells = 1, min.features = 1, project = study)
        }
        else if (mode == 'cellbender') {
                matrix <- paste0(matrix_dir, '/cellbender_filtered.h5')
                srt <- CreateSeuratObject(counts = ReadCB_h5(matrix), min.cells = 1, min.features = 1, project = study)
        }
        else if (mode == 'matrix') {
                matrix <- read.table(gzfile(paste0(matrix_dir[1], '.matrix.csv.gz')), header = T, sep = ',')
                srt <- CreateSeuratObject(counts = matrix, min.cells = 1, min.features = 1, project = study)
        }
        srt$sample <- sample
        srt$orig.name <- Cells(srt)
        srt$study <- study
        srt$method <- method
        srt$platform <- platform
        srt$protocol <- protocol
        srt$data_process <- data_process
        srt$tissue <- tissue
        srt$enrichment <- enrichment
        srt$preparation <- preparation
        srt$diagnosis <- diagnosis
        srt$recipient <- recipient
        srt$sex_recipient <- sex_recipient
        srt$age_recipient <- age_recipient
        srt$sex_donor <- sex_donor
        srt$age_donor <- age_donor
        srt$duration <- duration
        srt$replicate <- replicate
        srt <- RenameCells(srt, new.names = paste(srt$study,
                                                  srt$sample,
                                                  srt$orig.name,
                                                  sep = ':'), for.merge = F)
        srt <- PercentageFeatureSet(srt, pattern = '^MT-', col.name = 'pct_mito', assay = 'RNA')
        srt$pct_mito[is.nan(srt$pct_mito)] <- 0
        return(srt)
}
MakeDataset <- function(study, study_id, sample_name, mode, matrix_dir, starting_sample = 1){
        srt.list <- list()
        for(i in 1:L(sample_name)) {
                sample_id = paste0(study_id, '_S', str_pad(starting_sample - 1 + i, 3, pad = '0'))
                message('Processing sample:', sample_id)
                sample_meta_sub.df <- sample_meta.df[sample_meta.df$Study == study &
                                                             sample_meta.df$Sample_id == sample_id &
                                                             sample_meta.df$Name_on_disk == sample_name[i], ]
                message('Sample metadata found')
                srt.list[[i]] <- MakeSrt(mode = mode,
                                         matrix_dir = matrix_dir[i],
                                         study = study,
                                         sample = sample_id,
                                         method = sample_meta_sub.df$Method,
                                         platform = sample_meta_sub.df$Platform,
                                         protocol = sample_meta_sub.df$Protocol,
                                         data_process = sample_meta_sub.df$Data_process,
                                         tissue = sample_meta_sub.df$Tissue,
                                         enrichment = sample_meta_sub.df$Enrichment,
                                         preparation = sample_meta_sub.df$Preparation,
                                         procedure = sample_meta_sub.df$Procedure,
                                         diagnosis = sample_meta_sub.df$Diagnosis,
                                         recipient = sample_meta_sub.df$Recipient,
                                         sex_recipient = sample_meta_sub.df$Sex_recipient,
                                         age_recipient = sample_meta_sub.df$Age_recipient,
                                         sex_donor = sample_meta_sub.df$Sex_donor,
                                         age_donor = sample_meta_sub.df$Age_donor,
                                         duration = sample_meta_sub.df$Duration,
                                         replicate = sample_meta_sub.df$Replicate
                )
                print('Seurat generated...')
                # print(srt.list[[i]])
                # cat('\n_____________________________________________________\n')
        }
        if(L(srt.list) > 1) {
                merge.srt <- merge(srt.list[[1]], srt.list[2:L(srt.list)])
        } else {
                merge.srt <- srt.list[[1]]
        }
        return(merge.srt)
}
MakeRawDataset <- function(study, sample_name, raw_matrix_type, starting_sample = 1){
        no. <- names(studies[studies==study])
        message('Collecting Raw Data...')
        merge.srt <- MakeDataset(study = study,
                                 study_id = no.,
                                 sample_name = sample_name,
                                 mode = raw_matrix_type,
                                 matrix_dir = paste0('/Volumes/shire/data/scrnaseq/',
                                                     study, '/matrix/', sample_name),
                                 starting_sample = starting_sample)
        message('Processing Raw Seurat Object...')
        # merge.srt <- Process(merge.srt, assay = 'RNA')
        message('Saving Raw Seurat Object...')
        saveRDS(merge.srt, paste0('individual/', no., '.', study, '.raw.srt.rds'))
        # SaveH5ad(merge.srt, path = 'individual/', name = paste0(no., '.', study, '.raw'),
        #          assay = 'RNA', raw_count_only = F, verbose = F)
        rm(merge.srt)
        gc()
}
MakeCbnDataset <- function(study, sample_name, raw_matrix_type, starting_sample = 1, cb_folder = 'cellbender_v1'){
        no. <- names(studies[studies==study])
        message('Collecting CellBender Data...')
        merge.srt <- MakeDataset(study = study,
                                 study_id = no.,
                                 sample_name = sample_name,
                                 mode = 'cellbender',
                                 matrix_dir = paste0('/Volumes/shire/data/scrnaseq/',
                                                     study, '/matrix/', cb_folder, '/', sample_name),
                                 starting_sample = starting_sample)
        message('Processing CellBender Seurat Object...')
        # merge.srt <- Process(merge.srt, assay = 'RNA')
        message('Saving CellBender Seurat Object...')
        saveRDS(merge.srt, paste0('individual/', no., '.', study, '.cbn.srt.rds'))
        # SaveH5ad(merge.srt, path = 'individual/', name = paste0(no., '.', study, '.cbn'),
        #          assay = 'RNA', raw_count_only = F, verbose = F)
        rm(merge.srt)
        gc()
}
####--------------------------------------------------------------------------------------------------------------------


####--------------------------------------------------------------------------------------------------------------------
####  Load sample metadata  ####
####--------------------------------------------------------------------------------------------------------------------
sample_meta.df <- read.csv(paste0(Docu_dir, 'pediatric_sample_meta.csv'))
studies <- U(sample_meta.df$Study)
names(studies) <- U(sample_meta.df$Study_id)
studies_cellbender <- studies[studies %in% sample_meta.df$Study[sample_meta.df$Platform %in% c('10X', 'Drop-seq')]]
####--------------------------------------------------------------------------------------------------------------------


####--------------------------------------------------------------------------------------------------------------------
####  Dataset #1: 2019_Unpub7_JMartin  ####
####--------------------------------------------------------------------------------------------------------------------
study <- '2019_Unpub7_JMartin'
sample_name <- c(
        '2019_Unpub7_JMartin_ctrl_s4',
        '2019_Unpub7_JMartin_ctrl_s5',
        '2019_Unpub7_JMartin_ctrl_s11',
        '2019_Unpub7_JMartin_ctrl_s12',
        '2019_Unpub7_JMartin_ctrl_s6',
        '2019_Unpub7_JMartin_ctrl_s7'
        )
matrix_type <- '10x'
MakeRawDataset(study = study, sample_name = sample_name, raw_matrix_type = matrix_type)
MakeCbnDataset(study = study, sample_name = sample_name, raw_matrix_type = matrix_type, cb_folder = 'cellbender')
####--------------------------------------------------------------------------------------------------------------------


####--------------------------------------------------------------------------------------------------------------------
####  Dataset #2: 2021_Circulation_EPorrello  ####
####--------------------------------------------------------------------------------------------------------------------
study <- '2021_Circulation_EPorrello'
sample_name <- c(
        '2021_Circulation_EPorrello_Y2',
        '2021_Circulation_EPorrello_Y3'
)
matrix_type <- '10x'
MakeRawDataset(study = study, sample_name = sample_name, raw_matrix_type = matrix_type)
MakeCbnDataset(study = study, sample_name = sample_name, raw_matrix_type = matrix_type, cb_folder = 'cellbender')
####--------------------------------------------------------------------------------------------------------------------


####--------------------------------------------------------------------------------------------------------------------
####  Dataset #3: 2022_Unpub2_DTuraga  ####
####--------------------------------------------------------------------------------------------------------------------
study <- '2022_Unpub2_DTuraga'
sample_name <- c(
        '2022_Unpub2_DTuraga_P136_s1',
        '2022_Unpub2_DTuraga_P170_s1'
)
matrix_type <- '10x'
MakeRawDataset(study = study, sample_name = sample_name, raw_matrix_type = matrix_type)
MakeCbnDataset(study = study, sample_name = sample_name, raw_matrix_type = matrix_type)
####--------------------------------------------------------------------------------------------------------------------


####--------------------------------------------------------------------------------------------------------------------
####  Dataset #4: 2022_Unpub3_DTuraga  ####
####--------------------------------------------------------------------------------------------------------------------
study <- '2022_Unpub3_DTuraga'
sample_name <- c(
        '2022_Unpub3_DTuraga_P182'
)
matrix_type <- '10x'
MakeRawDataset(study = study, sample_name = sample_name, raw_matrix_type = matrix_type)
MakeCbnDataset(study = study, sample_name = sample_name, raw_matrix_type = matrix_type)
####--------------------------------------------------------------------------------------------------------------------
