####--------------------------------------------------------------------------------------------------------------------
####  Pediatric Acute Rejection
####  2022-02-24 by Xiao LI (Texas Heart Institute, US)
####
####  Helper functions for data processing
####--------------------------------------------------------------------------------------------------------------------


####--------------------------------------------------------------------------------------------------------------------
####  Load Libraries  ####
####--------------------------------------------------------------------------------------------------------------------
suppressMessages(library('hdf5r'))
suppressMessages(library('data.table'))
####--------------------------------------------------------------------------------------------------------------------


####--------------------------------------------------------------------------------------------------------------------
####  Initiate global variables  ####
####--------------------------------------------------------------------------------------------------------------------
s.genes <- cc.genes.updated.2019$s.genes
g2m.genes <- cc.genes.updated.2019$g2m.genes

# markers_lvl1 <- c('Tnnt2',  'Nppa',           # CM
#                   'Col1a1', 'Tcf21',          # CF
#                   'Pecam1',  'Cdh5',          # EC
#                   'Prox1',   'Flt4',          # LEC
#                   'Acta2',   'Rgs5',          # SMC
#                   'Pdgfrb', 'Vtn',            # Pericyte
#                   'Upk3b',  'Wt1',            # Epicardium
#                   'Npr3',   'H19',            # Endocardium
#                   'Plp1',   'Nrn1',           # Schwann
#                   'Nrsn1',  'Npy',            # Neuronal
#                   'Tenm4',  'Car3',           # Adipocyte
#                   'Ptprc',  'H2-D1',          # Immune
#                   'Hba-a1', 'Hbb-bh1',        # RBC
#                   'Mki67',  'Top2a')          # Mitotic
# markers_lvl1_min <- markers_lvl1[seq(1, L(markers_lvl1), 2)]
# markers_lvl2_early <- c('Epcam',  'Krt8',     # Mesoderm progen
#                         'Nkx2-5', 'Isl1',     # Mesoderm
#                         'Foxa2', 'Pax9',      # Endoderm
#                         'Sox2',  'Wnt6')      # Ectoderm
# markers_lvl2_immune <- c('Adgre1', 'Fcgr1',   # Mf, Mono
#                          'Xcr1',   'Cd209a',  # Dc
#                          'Cd79a',  'Ms4a1',   # B
#                          'Cd3e',   'Nkg7',    # T+NK
#                          'S100a9', 'S100a8')  # Granulocyte
# markers_lvl2_bm <- c('Ccl5',   # ILC
#                      'Vpreb1', # B
#                      'Pf4',  # megakaryocytes
#                      'Car1',   # erythrocytes
#                      'Prss34',    # basophils
#                      'Prg2')  # Geosinophils
# markers_lvl2_tc <- c('Cd4',
#                      'Trdc', #
#                      'Cd8a', #
#                      'Foxp3',  # Treg
#                      'Ifngr1', # Treg, Th1
#                      'Slamf6',   # Tfh
#                      'Pdcd1',    # Cd8_exhausted
#                      'Ccr7', # Cd8_mem (ccr7+Pdcd1+)
#                      'Ctla4',
#                      'Xcl1',
#                      'Gzmb',
#                      'Tox')  # Tfh
# markers_lvl2_cm <- c('Nr2f1', 'Cav1', # Atrium
#                      'Isl1', 'Tcn', # Outflow tract
#                      'Myl2', 'Mpped2', # Ventrical
#                      'Shox2', 'Pitx2'
# )
markers_lvl1 <- c('MYL7',   'NPPA',     # Atrial CM
                  'MYL2',   'MYH7',     # Ventricular CM
                  'DCN',    'GSN',      # CF
                  'PECAM1', 'VWF',      # EC
                  'MYH11',  'TAGLN',    # SMC
                  'KCNJ8',  'PDGFRB',   # Pericyte
                  'WT1',    'MSLN',     # Methothelial
                  'NRXN1',  'NRXN3',    # Neuronal
                  'GPAM',   'FASN',     # Adipocyte
                  'PTPRC',  'HLA-DRA',  # Immune
                  'HBB',    'HBA1',     # RBC
                  'MKI67',  'TOP2A')    # Mitotic
markers_lvl1_min <- markers_lvl1[seq(1, L(markers_lvl1), 2)]
markers_lvl2_immune <- c('ADGRE1', 'MRC1',  # Mf, Mono
                         'XCR1',   'CD1A', # Dc
                         'CD79A',  'MS4A1',  # B
                         'CD3E',   'NKG7',   # T+NK
                         'S100A9', 'S100A8') # Granulocyte

## Colors
# study_color <- clr_desaturate(shift = 0,
#                               c(
#                                       colorRampPalette(c("#2077b4", 'white'))(7)[1:5],
#                                       colorRampPalette(c("#d62628", 'white'))(7)[1:5],
#                                       colorRampPalette(c("#FF8C00", 'white'))(7)[1:5],
#                                       colorRampPalette(c("#2ba32b", 'white'))(7)[1:5],
#                                       colorRampPalette(c("#9466bd", 'white'))(7)[1:5],
#                                       colorRampPalette(c("#e377c1", 'white'))(7)[1:5],
#                                       colorRampPalette(c("#18bdcf", 'white'))(7)[1:5],
#                                       colorRampPalette(c("#bcbd21", 'white'))(7)[1:5],
#                                       colorRampPalette(c("#8b564c", 'white'))(7)[c(1,5)]
#                               )
# )
# age_group_color <- ScaleColor(mycol_RYB, 11)[c(11:8, 4:1)]
# condition_color <- c(paste0('grey', round(seq(40, 80, length.out = 10))), ScaleColor(mycol_BuGr, 11), mycol_40[5:8])
# genotype_color <- mycol_20[c(3,4,6,14)]

options(scipen = 999)
options(future.globals.maxSize= 2000*1024^2)
set.seed(505)
####--------------------------------------------------------------------------------------------------------------------


####--------------------------------------------------------------------------------------------------------------------
####  Initiate global functions  ####
####--------------------------------------------------------------------------------------------------------------------
InitiateProject <- function(Machine = 'Rivendell', Ver, Part, Catagory = 'mouse', Project_dir, Data_drive){
        if (Machine == 'Rivendell'){
                Data_dir <<- paste0('/Volumes/', Data_drive, '/project/', Project_dir, '/')
                File_dir <<- paste0('~/Documents/Bioinformatics/project/', Project_dir, '/')
        } else if (Machine == 'Moria') {
                Data_dir <<- paste0('/moria/', Project_dir, '/')
                File_dir <<- Data_dir
        } else if (Machine == 'bobbyd') {
                Data_dir <<- paste0('~/work/', Project_dir, '/')
                File_dir <<- Data_dir
        } else {stop('Machine not defined')}
        ####------------------------------------------------------------------------------------------------------------
        Rdata_dir <<- paste0(Data_dir, 'rdata/', Catagory, '_v', Ver)
        Docu_dir <<- paste0(File_dir, 'doc/', Catagory, '_v', Ver, '/')
        Plot_dir <<- paste0(File_dir, 'plot/', Catagory, '_v', Ver, '/', Part, '/')
        Meta_dir <<- paste0(File_dir, 'meta/', Catagory, '_v', Ver, '/', Part, '/')
        dir.create(file.path(Rdata_dir), showWarnings = F, recursive = T)
        setwd(Rdata_dir)
        dir.create(file.path('./individual'), showWarnings = F, recursive = T)
        dir.create(file.path('./integrated'), showWarnings = F, recursive = T)
        dir.create(file.path('./analysis'), showWarnings = F, recursive = T)
        dir.create(file.path('./external'), showWarnings = F, recursive = T)
        dir.create(file.path('./tmp'), showWarnings = F, recursive = T)
        dir.create(file.path(Plot_dir), showWarnings = F, recursive = T)
        dir.create(file.path(Meta_dir), showWarnings = F, recursive = T)
        dir.create(file.path(Docu_dir), showWarnings = F, recursive = T)
}
ProcessSrt_std <- function(srt_obj, var.toal = 0.75, assay = 'RNA', do.umap = T, npcs = 50, ...) {
        srt.out <- srt_obj %>%
                NormalizeData() %>%
                CellCycleScoring(s.features = s.genes,
                                 g2m.features = g2m.genes,
                                 set.ident = F) %>%
                PercentageFeatureSet(pattern = '^MT-', col.name = paste0('pct_mito_', assay), assay = assay)
        srt.out@meta.data[, paste0('pct_mito_', assay)][is.nan(srt.out@meta.data[, paste0('pct_mito_', assay)])] <- 0
        srt.out <- srt.out %>%
                FindVariableFeatures() %>%
                ScaleData(vars.to.regress = c(paste0('nFeature_', assay),
                                              paste0('pct_mito_', assay),
                                              'S.Score',
                                              'G2M.Score'), ...) %>%
                RunPCA(verbose = F, seed.use = 505, npcs = npcs)
        if(do.umap){
                dimN <- FindDimNumber(srt_obj = srt.out, var.toal = var.toal, reduction = 'pca')
                srt.out <- RunUMAP(srt.out, dims = 1:dimN, seed.use = 505)
        }
        return(srt.out)
}

ProcessSrt_sct <- function(srt_obj, var.toal = 0.75, do.umap = T, assay = 'RNA', regress_cc = T) {
        srt.out <- srt_obj %>%
                NormalizeData() %>%
                CellCycleScoring(s.features = s.genes,
                                 g2m.features = g2m.genes,
                                 set.ident = F) %>%
                PercentageFeatureSet(pattern = '^MT-', col.name = paste0('pct_mito_', assay), assay = assay)
                if(regress_cc){
                        srt.out <- SCTransform(srt.out,
                                               assay = assay,
                                               method = "glmGamPoi",
                                               seed.use = 505,
                                               return.only.var.genes = F,
                                               vars.to.regress = c(paste0('nFeature_', assay),
                                                                   paste0('pct_mito_', assay),
                                                                   'S.Score', 'G2M.Score'))
                } else {
                        srt.out <- SCTransform(srt.out,
                                               assay = assay,
                                               method = "glmGamPoi",
                                               seed.use = 505,
                                               return.only.var.genes = F,
                                               vars.to.regress = c(paste0('nFeature_', assay),
                                                                   paste0('pct_mito_', assay)))
                }
        srt.out <- RunPCA(srt.out, verbose = F, seed.use = 505)
        dimN <- FindDimNumber(srt_obj = srt.out, var.toal = var.toal, reduction = 'pca')
        if(do.umap){srt.out <- RunUMAP(srt.out, dims = 1:dimN, seed.use = 505)}
        return(srt.out)
}

ProcessSrt_hmn <- function(srt_obj, haromize.by = 'sample', harmony_max_iter = 30, var.toal = 0.75, assay = 'RNA',
                           do.umap = T, harmony_early_stop = T, dimN = NULL, ...) {
        if(harmony_early_stop){
                srt.out <- srt_obj %>%
                        RunHarmony(group.by.vars = haromize.by,
                                   max.iter.harmony = harmony_max_iter,
                                   plot_convergence = T,
                                   assay.use = assay)
        } else {
                srt.out <- srt_obj %>%
                        RunHarmony(group.by.vars = haromize.by,
                                   epsilon.cluster = -Inf, epsilon.harmony = -Inf,
                                   plot_convergence = T,
                                   max.iter.harmony = harmony_max_iter,
                                   assay.use = assay)
        }
        if(is.null(dimN)){
                dimN <- FindDimNumber(srt_obj = srt.out, var.toal = var.toal, reduction = 'harmony')
        }
        if(do.umap){srt.out <- RunUMAP(srt.out, dims = 1:dimN, seed.use = 505, reduction = 'harmony',
                                       reduction.name = 'hmn_umap', reduction.key = 'hmnumap_', ...)
        }
        return(srt.out)
}

ProcessSrt_clust <- function(srt_obj, resolution = 0.2, reduction = 'harmony', var.toal = 0.75, dimN = NULL, ...){
        if(is.null(dimN)){
                dimN <- dimN <- FindDimNumber(srt_obj = srt_obj, var.toal = var.toal, reduction = reduction)
        }
        srt.out <- srt_obj %>%
                FindNeighbors(dims = 1:dimN, reduction = reduction, force.recalc = T) %>%
                FindClusters(resolution = resolution, ...)
}

IntersectSrts <- function(ref.srt, que.srt, ref_data.name, que_data.name) {
        ## Identify shared cells
        shared_cells <- intersect(Cells(ref.srt), Cells(que.srt))
        title <- paste0(ref_data.name, ' has: ', ncol(ref.srt), '... ',
                        que_data.name, ' has: ', ncol(que.srt), '... ',
                        'Shared cells: ',        L(shared_cells))
        ref.srt <- ref.srt %>%
                ProcessSrt_std() %>%
                ProcessSrt_hmn()
        p1 <- DimPlot(ref.srt,
                      cells.highlight = shared_cells,
                      cols.highlight = 'grey60', cols = 'red', pt.size = 0.1, sizes.highlight = 0.1, raster = F) +
                labs(title = paste0(ref_data.name, ' has: ', ncol(ref.srt), '... Shared cells: ', L(shared_cells))) +
                theme(aspect.ratio = 1, axis.text = element_blank(), axis.ticks = element_blank(),
                      legend.position = 'bottom')
        que.srt <- que.srt %>%
                ProcessSrt_std() %>%
                ProcessSrt_hmn()
        p2 <- DimPlot(que.srt,
                      cells.highlight = shared_cells,
                      cols.highlight = 'grey60', cols = 'red', pt.size = 0.1, sizes.highlight = 0.1, raster = F) +
                labs(title = paste0(que_data.name, ' has: ', ncol(que.srt), '... Shared cells: ', L(shared_cells))) +
                theme(aspect.ratio = 1, axis.text = element_blank(), axis.ticks = element_blank(),
                      legend.position = 'bottom')
        print(title)
        ## Create new intersect srt from query data
        que.srt <- que.srt[, shared_cells]
        print('New seurat object generated:')
        print(que.srt)

        ## Match meta data
        ref.meta <- ref.srt@meta.data
        que.meta <- ref.meta[Cells(que.srt),
                             c('pub', 'sample', 'orig.name', 'method', 'platform', 'protocol', 'processed',
                               'tissue', 'enrichment', 'preparation', 'sex', 'age', 'genotype_s', 'genotype_l',
                               'condition', 'strain', 'replicate', 'group', 'batch')]
        print(paste0(nrow(que.meta), ' cells found in reference metadata'))
        que.srt@meta.data <- cbind(que.srt@meta.data, que.meta)
        que.srt$processed <- 'CellBender'
        print(head(que.srt@meta.data, n = 3))
        return(list(p1, p2, que.srt))
}

## A Cellbender bug workaround:
ReadCB_h5 <- function(filename, use.names = TRUE, unique.features = TRUE) {
        if (!requireNamespace('hdf5r', quietly = TRUE)) {
                stop("Please install hdf5r to read HDF5 files")
        }
        if (!file.exists(filename)) {
                stop("File not found")
        }
        infile <- hdf5r::H5File$new(filename = filename, mode = 'r')
        genomes <- names(x = infile)
        output <- list()
        if (hdf5r::existsGroup(infile, 'matrix')) {
                # cellranger version 3
                message('CellRanger version 3+ format H5')
                if (use.names) {
                        feature_slot <- 'features/name'
                } else {
                        feature_slot <- 'features/id'
                }
        } else {
                message('CellRanger version 2 format H5')
                if (use.names) {
                        feature_slot <- 'gene_names'
                } else {
                        feature_slot <- 'genes'
                }
        }
        for (genome in genomes) {
                counts <- infile[[paste0(genome, '/data')]]
                indices <- infile[[paste0(genome, '/indices')]]
                indptr <- infile[[paste0(genome, '/indptr')]]
                shp <- infile[[paste0(genome, '/shape')]]
                features <- infile[[paste0(genome, '/', feature_slot)]][]
                barcodes <- infile[[paste0(genome, '/barcodes')]]
                sparse.mat <- sparseMatrix(
                        i = indices[] + 1,
                        p = indptr[],
                        x = as.numeric(x = counts[]),
                        dims = shp[],
                        giveCsparse = FALSE
                )
                if (unique.features) {
                        features <- make.unique(names = features)
                }
                rownames(x = sparse.mat) <- features
                colnames(x = sparse.mat) <- barcodes[]
                sparse.mat <- as(object = sparse.mat, Class = 'dgCMatrix')
                # Split v3 multimodal
                if (infile$exists(name = paste0(genome, '/features'))) {
                        types <- infile[[paste0(genome, '/features/feature_type')]][]
                        types.unique <- unique(x = types)
                        if (length(x = types.unique) > 1) {
                                message("Genome ",
                                        genome, "
                                        has multiple modalities, returning a list of matrices for this genome")
                                sparse.mat <- sapply(
                                        X = types.unique,
                                        FUN = function(x) {
                                                return(sparse.mat[which(x = types == x), ])
                                        },
                                        simplify = FALSE,
                                        USE.NAMES = TRUE
                                )
                        }
                }
                output[[genome]] <- sparse.mat
        }
        infile$close_all()
        if (length(x = output) == 1) {
                return(output[[genome]])
        } else{
                return(output)
        }
}
QuickCheck <- function(srt_obj, markers, ...) {
        p1 <- VlnPlot(srt_obj,
                      features = c('nFeature_RNA', 'nCount_RNA', 'pct_mito'),
                      group.by = 'sample', ncol = 3, pt.size = -1) &
                theme(aspect.ratio = 0.5)
        p1 <- wrap_plots(list((
                p1[[1]] + theme(axis.text.x = element_blank())),
                (p1[[2]] + theme(axis.text.x = element_blank())),
                p1[[3]]), ncol = 1)
        p2.1 <- DimPlot(srt_obj, group.by = "sample", reduction = 'hmn_umap', raster = F) +
                labs(x = "UMAP 1", y = "UMAP 2", title = paste0(ncol(srt_obj), " Cells x ", nrow(srt_obj), " Genes")) +
                guides(col = guide_legend(ncol = 1)) +
                theme(aspect.ratio = 1,
                      axis.text = element_blank(), axis.ticks = element_blank(), legend.position = "bottom")
        p2.2 <- DimPlot(srt_obj, group.by = "sample", reduction = 'umap', raster = F) +
                labs(x = "UMAP 1", y = "UMAP 2", title = paste0(ncol(srt_obj), " Cells x ", nrow(srt_obj), " Genes")) +
                guides(col = guide_legend(ncol = 1)) +
                theme(aspect.ratio = 1,
                      axis.text = element_blank(), axis.ticks = element_blank(), legend.position = "bottom")
        p3 <- FeaturePlot2(srt_obj,
                           reduction = 'hmn_umap',
                           features = intersect(markers, rownames(srt_obj)),
                           ncol = ceiling(L(intersect(markers, rownames(srt_obj))) / 4))
        return(list(p2.1, p2.2, p1, p3))
}
SaveH5ad <- function(srt, path, name, assay, raw_count_only = F, h5seurat_keep = F, verbose = T){
        for(i in 1:ncol(srt@meta.data)){srt@meta.data[, i] <- as.vector(srt@meta.data[, i])} ## convert factor to vectors
        srt <- DietSeurat(srt, scale.data = F,
                          assays = assay,
                          dimreducs = names(srt@reductions),
                          graphs = names(srt@graphs))
        if(raw_count_only){srt <- SetAssayData(srt, slot = 'data', new.data = GetAssayData(srt, slot = 'counts'))}
        if(verbose){message('Raw matrix:')
                print(GetAssayData(srt, slot = 'counts')[1:20, 1:10])
                message('Data matrix:')
                print(GetAssayData(srt, slot = 'data')[1:20, 1:10])
                message('Scaled Data matrix:')}
        if(sum(dim(GetAssayData(srt, assay = assay, slot = 'scale.data')))==0){message('No scaled data slot')} else{
                print(GetAssayData(srt, assay = assay, slot = 'scale.data')[1:20, 1:10])}
        SaveH5Seurat(object = srt, filename = paste0(path, '/', name, '.h5Seurat'), overwrite = T, verbose = verbose)
        Convert(paste0(path, '/', name, '.h5Seurat'), dest = "h5ad", assay = assay, overwrite = T, verbose = verbose)
        if(!h5seurat_keep){system(paste0('rm ', path, '/', name, '.h5Seurat'))}
}

Subcluster <- function(plot_num, srt, celltype, dimN, resolution, features = NULL, markerset,
                       regress = c('nCount_DCX', 'nFeature_DCX', 'pct_mito_DCX')){
        ## Seurat processing
        message('\n', 'Processing Seurat...')
        srt <- srt |>
                FindVariableFeatures(verbose = F) |>
                ScaleData(verbose = F, vars.to.regress = regress) |>
                RunPCA(seed.use = 505, verbose = F) |>
                RunHarmony(group.by.vars = 'sample', max.iter.harmony = 50, assay.use = 'DCX', verbose = F)
        gc()
        STEP <- str_split(Step, '_')[[1]][1]
        ad_name <- paste0(STEP, '.', celltype, '.ann')
        message('\n', 'Writing H5AD...')
        SaveH5ad(srt, path = 'integrated/', name = ad_name, assay = 'DCX', raw_count_only = F, verbose = F)
        message('\n', 'Reading H5AD...')
        ad <- read_h5ad(paste0('integrated/', ad_name, '.h5ad'))
        ## Scanpy processing
        message('\n', 'Reducing with ', dimN, ' components...')
        ad$obsm[['X_harmony_subset']] <- ad$obsm[['X_harmony']][, 1:dimN]
        sc$pp$neighbors(ad, use_rep = 'X_harmony_subset')
        sc$tl$umap(ad, neighbors_key = 'neighbors', min_dist = 0.3)
        ## Transfer reductions to Seurat
        srt@reductions$hmn_sc_umap2 <- srt@reductions$hmn_umap
        srt@reductions$hmn_sc_umap2@cell.embeddings[, 1] <- ad$obsm[['X_umap']][, 1]
        srt@reductions$hmn_sc_umap2@cell.embeddings[, 2] <- ad$obsm[['X_umap']][, 2]
        srt@reductions$hmn_sc_umap2@key <- 'hmnscumap2_'
        colnames(srt@reductions$hmn_sc_umap2@cell.embeddings) <- c("hmnscumap2_1", "hmnscumap2_2")
        srt <- DietSeurat(srt, dimreducs = c('hmn_sc_umap2', 'harmony'))
        gc()
        p <- wrap_plots(list(DimPlot2(srt, group.by = 'study', raster = F, cols = study_color) + NoLegend(),
                             DimPlot2(srt, group.by = 'method',  raster = F, cols = mycol_10),
                             DimPlot2(srt, group.by = 'platform', raster = F, cols = mycol_10),
                             DimPlot2(srt, group.by = 'age_group', cols = age_group_color, raster = F),
                             DimPlot2(srt, group.by = 'condition', cols = condition_color, raster = F),
                             DimPlot2(srt, group.by = 'genotype', cols = mycol_10, raster = F)
        ), ncol = 3)
        PlotPDF(paste0(plot_num, '.', celltype, '.', dimN, 'dims.meta'), 30, 20)
        print(p)
        dev.off()
        if(is.null(features)){ features <- markers_lvl1}
        p <- FeaturePlot2(srt, features = U(features), ncol = ceiling(LU(features)^0.5), raster = T)
        PlotPDF(paste0(plot_num, '.', celltype, '.', dimN, 'dims.feature'),
                ceiling(LU(features)^0.5)*2,
                floor(LU(features)^0.5)*2)
        print(p)
        dev.off()
        ## Subcluster evaluation
        message('\n', 'Clustering with ', dimN, ' components, at ', resolution, ' resolution...')
        sc$tl$leiden(ad, resolution = resolution, key_added = paste0('leiden_', resolution))
        srt@meta.data[[paste0('Subcluster_leiden_', resolution)]] <- ad$obs[paste0('leiden_', resolution)]
        Idents(srt) <- paste0('Subcluster_leiden_', resolution)
        levels(srt) <- str_sort(levels(srt), numeric = T)
        marker <- FindAllMarkers(srt, assay = 'DCX', only.pos = T, return.thresh = 0.0001)
        srt@misc$marker[[paste0('Subcluster_leiden_', resolution)]] <- marker
        cols <- mycol_20
        if(LU(srt@active.ident)>20){cols <- mycol_40}
        p <- DimPlot2(srt, label = T, repel = T, raster = F, cols = cols) + labs(title = paste0('leiden_', resolution))
        PlotPDF(paste0(plot_num, '.', celltype, '.', dimN, 'dims.', resolution, 'res.umap'), 10, 10)
        print(p)
        dev.off()
        p <- MarkerHeatmap(srt, marker.df = marker, top = 20)
        PlotPDF(paste0(plot_num, '.', celltype, '.', dimN, 'dims.', resolution, 'res.heatmap'),
                20,
                L(levels(srt))*20/8)
        print(p)
        dev.off()
        gc()
        message('\n', 'Updating H5AD...')
        ad$write_h5ad(filename = paste0('integrated/', ad_name, '.h5ad'))

        for(i in 1:L(markerset)){
                srt <- AddModuleScore(srt, features = markerset[[i]],
                                      name = paste0('Markerset__', names(markerset)[i], '_'))
        }
        plist <- list()
        for(i in 1:L(markerset)){
                set_names <- grep(paste0('^Markerset__', names(markerset)[i], '_'), colnames(srt@meta.data), value = T)
                plist[[i]] <- DotPlot2(srt, features = set_names, col.min = 0) + labs(title = names(markerset)[i])
        }
        PlotPDF(paste0(plot_num, '.', celltype, '.', dimN, 'dims.', resolution, '.dot_markerset'),
                15, 15, onefile = T)
        print(plist)
        dev.off()
        srt <- AddModuleScore(srt, features = core_markers, name = paste0('Core_marker__'))
        for(i in 1:L(core_markers)){
                srt@meta.data[paste0('Core_marker__', names(core_markers)[i])] <- srt@meta.data[paste0('Core_marker__', i)]
                srt@meta.data[paste0('Core_marker__', i)] <- NULL
        }
        p <- DotPlot2(srt, features = grep('^Core_marker', colnames(srt@meta.data), value = T), col.min = 0) +
                labs(title = 'Core_markers')
        PlotPDF(paste0(plot_num, '.', celltype, '.', dimN, 'dims.', resolution, '.dot_core_marker'), 8, 8)
        print(p)
        dev.off()
        return(srt)
}

####--------------------------------------------------------------------------------------------------------------------
