####--------------------------------------------------------------------------------------------------------------------
####  Pediatric Acute Rejection
####  2022-08-03 by Xiao LI (Texas Heart Institute, US)
####--------------------------------------------------------------------------------------------------------------------


####--------------------------------------------------------------------------------------------------------------------
####  Initiate  ####
####--------------------------------------------------------------------------------------------------------------------
Ver <- '1'
Step <- 'STEP04_Doublet'
Project <- '2022_acute_rejection_dturaga'
Code_dir <- paste0('~/Documents/Bioinformatics/project/', Project, '/code/human_v', Ver, '/')

source(Sys.readlink(paste0(Code_dir, 'src/bioinformatics.R')))
source(Sys.readlink(paste0(Code_dir, 'src/scRNAseq.R')))
source(Sys.readlink(paste0(Code_dir, 'src/geneset_human.R')))
source(paste0(Code_dir, 'human_v', Ver, '.helper_functions.R'))

InitiateProject('Rivendell', Ver, Step, 'human', Project, 'ithil')

suppressMessages(library('DoubletFinder'))
suppressMessages(library('reticulate'))
scr <- import('scrublet')
plt <- import("matplotlib.pyplot")
####--------------------------------------------------------------------------------------------------------------------


####--------------------------------------------------------------------------------------------------------------------
####  Load sample metadata  ####
####--------------------------------------------------------------------------------------------------------------------
merged.dlt.srt <- readRDS('integrated/STEP03.merged.flt.srt.rds')
####--------------------------------------------------------------------------------------------------------------------


####--------------------------------------------------------------------------------------------------------------------
####  Global Functions  ####
####--------------------------------------------------------------------------------------------------------------------
GetDoublet <- function(srt_obj, doublet_rate, dimN.var.toal){
        ## Scrublet (run via reticulate)
        mtx <- srt_obj@assays$SCT@counts
        mtx <- t(mtx)
        scrub_model <- scr$Scrublet(mtx, expected_doublet_rate = doublet_rate)
        rst <- scrub_model$scrub_doublets(min_gene_variability_pctl = dimN.var.toal*100,
                                          n_prin_comps = 30L,
                                          min_counts = 2, min_cells = 3)
        rst[[2]] <- scrub_model$call_doublets(threshold = 0.25) ## adjusted based on histogram
        sc_doublets <- Cells(srt_obj)[rst[[2]]]
        sc_singlets <- Cells(srt_obj)[!rst[[2]]]
        srt_obj$Scrublet_doublet <- 'Singlet'
        srt_obj$Scrublet_doublet[rst[[2]]] <- 'Doublet'
        Scrublet <- rst[[1]]
        names(Scrublet) <- Cells(srt_obj)

        p2 <- DimPlotSplit(srt_obj, split_by = 'Scrublet_doublet', split_order = c('Singlet', 'Doublet'),
                           cols.highlight = mycol_14[c(2, 1)], ncol = 2)
        p2[[1]] <- p2[[1]] + labs(title = paste0('Srub Singlet: ', L(sc_singlets), ' Cells'))
        p2[[2]] <- p2[[2]] + labs(title = paste0('Srub Doublet: ', L(sc_doublets), ' Cells'))
        p <- wrap_plots(
                p2[[1]],
                p2[[2]],
                ncol = 2)
        return(list(
                sc_doublets,
                p,
                Scrublet,
                scrub_model
        ))
}
####--------------------------------------------------------------------------------------------------------------------


###--------------------------------------------------------------------------------------------------------------------
###  Identify doublets for each dataset  (linear processing) ####
###--------------------------------------------------------------------------------------------------------------------
system(paste0('mkdir ', Plot_dir, '/umap_per_sample/'))

doublet_rate <- 0.1 ## Assuming 10% doublet formation rate
all_Doublet_SC <- c()
all_Scrublet <- c()

all_samples <- levels(merged.dlt.srt$sample)
for(i in 1:L(all_samples)){
        gc()
        print(paste0('Processing ', all_samples[i], ' ...'))
        tmp.srt <- merged.dlt.srt[, merged.dlt.srt$sample == all_samples[i]]
        results <- GetDoublet(srt_obj = tmp.srt, doublet_rate = doublet_rate, dimN.var.toal = 0.85)
        all_Doublet_SC <- c(all_Doublet_SC, results[[1]])
        ## plot umap
        PlotPDF(paste0('umap_per_sample/', all_samples[i], '.doublets_found'), 10, 5)
        print(results[[2]])
        dev.off()
        all_Scrublet <- c(all_Scrublet, results[[3]])
        ## plot scrublet histogram
        PlotPDF(paste0('umap_per_sample/', all_samples[i], '.scrublet_hist'), 8, 4)
        print(plt$show(results[[4]]$plot_histogram()[[1]]))
        dev.off()
        ####------------------------------------------------------------------------------------------------------------
        saveRDS(all_Doublet_SC, 'analysis/STEP04.cell_filtered.scrublet_doublets.rds')
        saveRDS(all_Scrublet,   'analysis/STEP04.cell_filtered.scrublet_score.rds')
        ####------------------------------------------------------------------------------------------------------------
}

####--------------------------------------------------------------------------------------------------------------------
####  Evaluate doublets  ####
####--------------------------------------------------------------------------------------------------------------------
## save doublets to the main Seurat
merged.dlt.srt$Doublet_SC <- F
merged.dlt.srt$Doublet_SC[all_Doublet_SC] <- T
merged.dlt.srt$Doublet_SC_score <- NA
merged.dlt.srt$Doublet_SC_score[names(all_Scrublet)] <- all_Scrublet


p <- DimPlot(merged.dlt.srt, cols = 'grey75', reduction = 'hmn_umap', raster = T, pt.size = 0.001,
             cells.highlight = all_Doublet_SC, cols.highlight = 'red', sizes.highlight = 0.001) +
        labs(title = paste0('Total cells: ', ncol(merged.dlt.srt), '  Scrublet Doublets found: ', L(all_Doublet_SC))) +
        NoLegend() +
        theme(aspect.ratio = 1, axis.text = element_blank(), axis.ticks = element_blank())
PlotPDF('01.merged_filtered.scrublet_doublets', 10, 10)
p
dev.off()

p <- FeaturePlot2(merged.dlt.srt, reduction = 'hmn_umap', raster = T, features = 'Doublet_SC_score') +
        NoLegend() +
        theme(aspect.ratio = 1, axis.text = element_blank(), axis.ticks = element_blank())
PlotPDF('02.merged_filtered.scrublet_score', 10, 10)
p
dev.off()

p <- VlnPlot2(merged.dlt.srt, features = 'Doublet_SC_score', group.by = 'sample') +
        NoLegend()
PlotPDF('03.merged_filtered.scrublet_score_per_sample', 5, 5)
p
dev.off()

p1 <- CountCellBarPlot(merged.dlt.srt, group.var = 'study', stack.var = 'Doublet_SC',
                       stack.color = mycol_10)$plot
p2 <- CountCellBarPlot(merged.dlt.srt, group.var = 'study', stack.var = 'Doublet_SC',
                       stack.color = mycol_10, percentage = T)$plot
PlotPDF('04.merged_filtered.pct_sc_dlt_study', 10, 15, onefile = T)
wrap_plots(p1, p2, ncol = 1)
dev.off()

p1 <- CountCellBarPlot(merged.dlt.srt, group.var = 'sample', stack.var = 'Doublet_SC',
                       stack.color = mycol_10)$plot
p2 <- CountCellBarPlot(merged.dlt.srt, group.var = 'sample', stack.var = 'Doublet_SC',
                       stack.color = mycol_10, percentage = T)$plot
PlotPDF('05.merged_filtered.pct_sc_dlt_sample', 20, 15, onefile = T)
wrap_plots(p1, p2, ncol = 1)
dev.off()
####--------------------------------------------------------------------------------------------------------------------
saveRDS(merged.dlt.srt@meta.data, 'integrated/STEP04.merged.dlt.srt_meta.rds') ## Create new
####--------------------------------------------------------------------------------------------------------------------
