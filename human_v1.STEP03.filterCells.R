####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
####  Pediatric Acute Rejection
####  2022-08-03 by Xiao LI (Texas Heart Institute, US)
####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
####  Initiate  ####
####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
Ver <- '1'
Step <- 'STEP03_Filter_Cells'
Project <- '2022_acute_rejection_dturaga'
Code_dir <- paste0('~/Documents/Bioinformatics/project/', Project, '/code/human_v', Ver, '/')

source(Sys.readlink(paste0(Code_dir, 'src/bioinformatics.R')))
source(Sys.readlink(paste0(Code_dir, 'src/scRNAseq.R')))
source(Sys.readlink(paste0(Code_dir, 'src/geneset_human.R')))
source(paste0(Code_dir, 'human_v', Ver, '.helper_functions.R'))

InitiateProject('Rivendell', Ver, Step, 'human', Project, 'bree')
####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
####  Load data  ####
####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
merged.flt.srt <- readRDS('integrated/STEP02.merged.cbn.srt.rds')
studies <- levels(merged.flt.srt$study)
####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
####  Set hard quality filters  ####
####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Filter cells with less than 200 umi(counts)/cell
cutoff1 <- 200
pct_filtered <- rep(NA, L(studies))
names(pct_filtered) <- studies
p1 <- VlnPlot(merged.flt.srt,
              features = c('nCount_CBN'), log = T,
              group.by = 'sample', pt.size = -1, ncol = 1) +
        geom_hline(yintercept = cutoff1) +
        NoLegend()
for(i in 1:L(studies)){
        n_filtered <- sum(merged.flt.srt$nCount_CBN < cutoff1 & merged.flt.srt$study == studies[i])
        pct_filtered[i] <- n_filtered*100/sum(merged.flt.srt$study == studies[i])
}
## Filter cells with less than 150 genes/cell
cutoff2 <- 150
pct_filtered <- rep(NA, L(studies))
names(pct_filtered) <- studies
p2 <- VlnPlot(merged.flt.srt,
              features = c('nFeature_CBN'), log = T,
              group.by = 'sample', pt.size = -1, ncol = 1) +
        geom_hline(yintercept = cutoff2) +
        NoLegend()
for(i in 1:L(studies)){
        n_filtered <- sum(merged.flt.srt$nFeature_CBN < cutoff2 & merged.flt.srt$study == studies[i])
        pct_filtered[i] <- n_filtered*100/sum(merged.flt.srt$study == studies[i])
}
## Filter cells with higher than 5% mitochondrial counts/cell
cutoff3 <- 5
pct_filtered <- rep(NA, L(studies))
names(pct_filtered) <- studies
p3 <- VlnPlot(merged.flt.srt,
              features = c('pct_mito_CBN'),
              group.by = 'sample', pt.size = -1, ncol = 1) +
        geom_hline(yintercept = cutoff3) +
        NoLegend()
for(i in 1:L(studies)){
        n_filtered <- sum(merged.flt.srt$pct_mito_CBN > cutoff3 & merged.flt.srt$study == studies[i])
        pct_filtered[i] <- n_filtered*100/sum(merged.flt.srt$study == studies[i])
}

cells_toss_1 <- Cells(merged.flt.srt)[
        merged.flt.srt$nCount_CBN < cutoff1 |
                merged.flt.srt$nFeature_CBN < cutoff2 |
                merged.flt.srt$pct_mito_CBN > cutoff3
]
L(cells_toss_1) ## 10189
L(cells_toss_1)*100/ncol(merged.flt.srt) ## 7.612423%

p4 <- DimPlot(merged.flt.srt, cols = 'grey75', reduction = 'hmn_umap', raster = T,
              cells.highlight = cells_toss_1, cols.highlight = 'red', sizes.highlight = 0.001, pt.size = 0.001) +
        labs(title = paste0('Total cells: ', ncol(merged.flt.srt), '  Cells filtered (red): ', L(cells_toss_1))) +
        NoLegend() +
        theme(aspect.ratio = 1, axis.text = element_blank(), axis.ticks = element_blank())

PlotPDF('01.1.merged_consensus.hard_cutoff.vln', 12, 6)
p1 | p2 | p3
dev.off()
PlotPDF('01.2.merged_consensus.hard_cutoff.dim', 10, 10)
p4
dev.off()
####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
####  Set dynamic quality filters  ####
####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
L(cells_toss_1) ## 10189

## Remove bad cells that did not pass hard filters before setting dynamic filters
merged.flt.srt2 <- merged.flt.srt[, !Cells(merged.flt.srt) %in% cells_toss_1]
all_samples <- U(merged.flt.srt2$sample)
L(all_samples) ## 11

cutoffs <- data.frame(name = all_samples)

## Set filter for max nCount
cutoff_upper <- c()
cutoff_lower <- c()
cells_toss_2 <- c()
for(i in 1:L(all_samples)){
        print(i)
        cells_in_group <- Cells(merged.flt.srt2)[merged.flt.srt2$sample == all_samples[i]]
        GetOutlier(merged.flt.srt2$nCount_CBN[cells_in_group])
        cutoff_upper[i] <- GetOutlier(merged.flt.srt2$nCount_CBN[cells_in_group])[2]
        cutoff_lower[i] <- GetOutlier(merged.flt.srt2$nCount_CBN[cells_in_group])[1]
        if(cutoff_lower[i] < cutoff1){cutoff_lower[i] <- cutoff1}
        cells_toss_2 <- c(cells_toss_2,
                          cells_in_group[merged.flt.srt2$nCount_CBN[cells_in_group] > cutoff_upper[i] |
                                                 merged.flt.srt2$nCount_CBN[cells_in_group] < cutoff_lower[i]])
}
L(cells_toss_2) ## 6110
p <- list()
for(i in 1:L(all_samples)){
        p[[i]] <- ggplot(merged.flt.srt2@meta.data[merged.flt.srt2$sample == all_samples[i], ],
                         aes(x = sample, y = nCount_CBN)) +
                geom_violin(fill = mycol_30[3]) +
                geom_hline(yintercept = c(cutoff_upper[i], cutoff_lower[i]), color = mycol_14[1], size = 1) +
                labs(title = paste0('Sample # ', all_samples[i]),
                     subtitle = paste0(round(cutoff_lower[i], 0) , ' < UMI < ', round(cutoff_upper[i], 0))) +
                scale_y_continuous(limits = c(0, 5e4))
}
p2 <- wrap_plots(p[order(cutoff_upper)], ncol = 5) &
        NoLegend() &
        WhiteBackground() &
        theme(aspect.ratio = 2,
              axis.title = element_blank(),
              axis.text = element_blank(),
              axis.ticks = element_blank(),
              axis.line = element_line())
p3 <- DimPlot(merged.flt.srt, cols = 'grey75', reduction = 'hmn_umap', raster = T,
              pt.size = 0.001, sizes.highlight = 0.001,
              cells.highlight = cells_toss_2, cols.highlight = 'red') +
        labs(title = paste0('Total cells: ', ncol(merged.flt.srt), '  Cells filtered (red): ', L(cells_toss_2))) +
        NoLegend() +
        theme(aspect.ratio = 1, axis.text = element_blank(), axis.ticks = element_blank())
PlotPDF('02.1.merged_consensus.dynamic_cutoff_nCount.vln', 16, 8)
p2
dev.off()
PlotPDF('02.2.merged_consensus.dynamic_cutoff_nCount.dim', 10, 10)
p3
dev.off()

cutoffs$'nCount_max' <- cutoff_upper
cutoffs$'nCount_min' <- cutoff_lower

## Set filter for max nFeature
cutoff_upper <- c()
cutoff_lower <- c()
cells_toss_3 <- c()
for(i in 1:L(all_samples)){
        print(i)
        cells_in_group <- Cells(merged.flt.srt2)[merged.flt.srt2$sample == all_samples[i]]
        GetOutlier(merged.flt.srt2$nCount_CBN[cells_in_group])
        cutoff_upper[i] <- GetOutlier(merged.flt.srt2$nFeature_CBN[cells_in_group])[2]
        cutoff_lower[i] <- GetOutlier(merged.flt.srt2$nFeature_CBN[cells_in_group])[1]
        if(cutoff_lower[i] < cutoff2){cutoff_lower[i] <- cutoff2}
        cells_toss_3 <- c(cells_toss_3,
                          cells_in_group[merged.flt.srt2$nFeature_CBN[cells_in_group] > cutoff_upper[i] |
                                                 merged.flt.srt2$nFeature_CBN[cells_in_group] < cutoff_lower[i]])
}
L(cells_toss_3) ## 2184
p <- list()
for(i in 1:L(all_samples)){
        p[[i]] <- ggplot(merged.flt.srt2@meta.data[merged.flt.srt2$sample == all_samples[i], ],
                         aes(x = sample, y = nFeature_CBN)) +
                geom_violin(fill = mycol_30[3]) +
                geom_hline(yintercept = c(cutoff_upper[i], cutoff_lower[i]), color = mycol_14[1], size = 1) +
                labs(title = paste0('Sample # ', all_samples[i]),
                     subtitle = paste0(round(cutoff_lower[i], 0) , ' < Gene < ', round(cutoff_upper[i], 0))) +
                scale_y_continuous(limits = c(0, 1.5e4))
}
p2 <- wrap_plots(p[order(cutoff_upper)], ncol = 5) &
        NoLegend() &
        WhiteBackground() &
        theme(aspect.ratio = 2,
              axis.title = element_blank(),
              axis.text = element_blank(),
              axis.ticks = element_blank(),
              axis.line = element_line())
p3 <- DimPlot(merged.flt.srt, cols = 'grey75', reduction = 'hmn_umap', raster = T,
              pt.size = 0.001, sizes.highlight = 0.001,
              cells.highlight = cells_toss_3, cols.highlight = 'red') +
        labs(title = paste0('Total cells: ', ncol(merged.flt.srt), '  Cells filtered (red): ', L(cells_toss_3))) +
        NoLegend() +
        theme(aspect.ratio = 1, axis.text = element_blank(), axis.ticks = element_blank())
PlotPDF('03.1.merged_consensus.dynamic_cutoff_nFeature.vln', 16, 10)
p2
dev.off()
PlotPDF('03.2.merged_consensus.dynamic_cutoff_nFeature.dim', 10, 10)
p3
dev.off()

cutoffs$'nFeature_max' <- cutoff_upper
cutoffs$'nFeature_min' <- cutoff_lower

# ## Set filter for max mito percentage
# cutoff_upper <- c()
# cells_toss_4 <- c()
# for(i in 1:L(all_samples)){
#         print(i)
#         cells_in_group <- Cells(merged.flt.srt2)[merged.flt.srt2$sample == all_samples[i]]
#         cutoff_upper[i] <- GetOutlier(merged.flt.srt2$pct_mito_CBN[cells_in_group])[2]
#         if(cutoff_upper[i] < 1){cutoff_upper[i] <- 1}
#         if(cutoff_upper[i] > cutoff3){cutoff_upper[i] <- cutoff3}
#         cells_toss_4 <- c(cells_toss_4,
#                           cells_in_group[merged.flt.srt2$pct_mito_CBN[cells_in_group] > cutoff_upper[i]])
# }
# L(cells_toss_4) ## 10554
# p <- list()
# for(i in 1:L(all_samples)){
#         p[[i]] <- ggplot(merged.flt.srt2@meta.data[merged.flt.srt2$sample == all_samples[i], ],
#                          aes(x = sample, y = pct_mito_CBN)) +
#                 geom_violin(fill = mycol_30[3]) +
#                 geom_hline(yintercept = c(cutoff_upper[i], cutoff_lower[i]), color = mycol_14[1], size = 1) +
#                 labs(title = paste0('Sample # ', all_samples[i]),
#                      subtitle = paste0('mito% < ', round(cutoff_upper[i], 2))) +
#                 scale_y_continuous(limits = c(0, 10))
# }
# p2 <- wrap_plots(p[order(cutoff_upper)], ncol = 5) &
#         NoLegend() &
#         WhiteBackground() &
#         theme(aspect.ratio = 2,
#               axis.title = element_blank(),
#               axis.text = element_blank(),
#               axis.ticks = element_blank(),
#               axis.line = element_line())
# p3 <- DimPlot(merged.flt.srt, cols = 'grey75', reduction = 'hmn_umap', raster = T,
#               pt.size = 0.001, sizes.highlight = 0.001,
#               cells.highlight = cells_toss_4, cols.highlight = 'red') +
#         labs(title = paste0('Total cells: ', ncol(merged.flt.srt), '  Cells filtered (red): ', L(cells_toss_4))) +
#         NoLegend() +
#         theme(aspect.ratio = 1, axis.text = element_blank(), axis.ticks = element_blank())
# PlotPDF('04.1.merged_consensus.dynamic_cutoff_mito.vln', 16, 12)
# p2
# dev.off()
# PlotPDF('04.2.merged_consensus.dynamic_cutoff_mito.dim', 10, 10)
# p3
# dev.off()
####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
####  Filter cells  ####
####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
cells_toss <- U(c(cells_toss_1, cells_toss_2, cells_toss_3))
L(cells_toss) ## 16347
p1 <- DimPlot(merged.flt.srt, cols = 'grey75', reduction = 'hmn_umap', raster = T, pt.size = 0.001,
              cells.highlight = cells_toss, cols.highlight = 'red', sizes.highlight = 0.001) +
        labs(title = paste0('Total cells: ', ncol(merged.flt.srt),
                            '  Cells filtered (red): ', L(cells_toss),
                            ' (', round(L(cells_toss)*100/ncol(merged.flt.srt)),'%)')) +
        NoLegend() +
        theme(aspect.ratio = 1, axis.text = element_blank(), axis.ticks = element_blank())
PlotPDF('05.1.merged_consensus.all_cutoff', 10, 10)
p1
dev.off()

merged.flt.srt$LowQual <- F
merged.flt.srt$LowQual[cells_toss_3] <- 'cells_toss_3'
merged.flt.srt$LowQual[cells_toss_2] <- 'cells_toss_2'
merged.flt.srt$LowQual[cells_toss_1] <- 'cells_toss_1'

p1 <- CountCellBarPlot(merged.flt.srt, group.var = 'study', stack.var = 'LowQual', stack.color = mycol_10)$plot
p2 <- CountCellBarPlot(merged.flt.srt, group.var = 'study', stack.var = 'LowQual', stack.color = mycol_10, percentage = T)$plot
PlotPDF('05.2.merged_consensus.pct_toss_study', 10, 15, onefile = T)
wrap_plots(p1, p2, ncol = 1)
dev.off()

p1 <- CountCellBarPlot(merged.flt.srt, group.var = 'sample', stack.var = 'LowQual', stack.color = mycol_10)$plot
p2 <- CountCellBarPlot(merged.flt.srt, group.var = 'sample', stack.var = 'LowQual', stack.color = mycol_10, percentage = T)$plot
PlotPDF('05.3.merged_consensus.pct_toss_group', 20, 15, onefile = T)
wrap_plots(p1, p2, ncol = 1)
dev.off()

cell_toss.list <- list(cells_toss_1, cells_toss_2, cells_toss_3)

rm(merged.flt.srt2, p, p1, p2, p3, p4)
gc()
merged.flt.srt <- merged.flt.srt[, !Cells(merged.flt.srt) %in% cells_toss]
####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
saveRDS(cell_toss.list, 'analysis/STEP03.cell_filtered.low_quality_cells.rds')
saveRDS(cutoffs, 'analysis/PART03.cell_filtered.low_quality_cutoffs.df.rds')
saveRDS(merged.flt.srt@meta.data, 'integrated/STEP03.merged.flt.srt_meta.rds')
saveRDS(merged.flt.srt, 'integrated/STEP03.merged.flt.srt.rds')
####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
