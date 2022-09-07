# library("Signac")
# library("MultiAssayExperiment")
# wt211 <- readRDS("signac/wt_frozen2021_signac.rds")
# wt212 <- readRDS("signac/wt_frozen2021_2_signac.rds")
# wt172 <- readRDS("signac/wt_frozen2017_2_signac.rds")
#
# rna1 <- wt211@assays$RNA@counts
# rna2 <- wt212@assays$RNA@counts
# rna0 <- wt172@assays$RNA@counts
# # RNA <- cbind(rna0, rna1, rna2)
#
# atc0 <- wt172@assays$ATAC@counts
# atc1 <- wt211@assays$ATAC@counts
# atc2 <- wt212@assays$ATAC@counts
# # ATAC <- cbind(atc0, atc1, atc2)
# #
#
# sampleMapExps <- DataFrame(
#     assay=c(rep("RNA172", dim(rna0)[2]),
#           rep("ATAC172", dim(atc0)[2]),
#           rep("RNA211", dim(rna1)[2]),
#           rep("ATAC211", dim(atc1)[2]),
#           rep("RNA212", dim(rna2)[2]),
#           rep("ATAC212", dim(atc2)[2])
#     ),
#     primary=c(paste0(colnames(rna0), "_172"),
#               paste0(colnames(atc0), "_172"),
#               paste0(colnames(rna1), "_211"),
#               paste0(colnames(atc1), "_211"),
#               paste0(colnames(rna2), "_212"),
#               paste0(colnames(atc2), "_212")),
#     colname=c(colnames(rna0), colnames(atc0), colnames(rna1),
#               colnames(atc1), colnames(rna2), colnames(atc2))
# )
# mae1 <- MultiAssayExperiment(experiments=ExperimentList(RNA172=rna0,
#                                                         ATAC172=atc0,
#                                                         RNA211=rna1,
#                                                         ATAC211=atc1,
#                                                         RNA212=rna2,
#                                                         ATAC212=atc2),
#                              sampleMap=sampleMapExps)
#
# ## problema con nn.dist, impossibilità di salvare (in maniera intuitiva) le distanze in sampleMap perchè
# ## per wt211 sono 3440 e nella sampleMap abbiamo i colnames associati ad atac e rna per ogni esperimento,
# ## quindi sono raddoppiati.
# ## servirebbe un posto specifico per salvare i risultati integrati in ogni esperimento.
# ## Stesso problema per UMAP ecc
# sampleMap(maetest)$nn.dist <- NA
# sampleMap(maetest)$nn.dist[which(sampleMap(maetest)$colname[grep("211", sampleMap(maetest)$primary)] %in%
#     wt211@neighbors$weighted.nn@cell.names)] <- wt211@neighbors$weighted.nn@nn.dist
#
# length(which(sampleMap(maetest)$colname %in%
#         wt211@neighbors$weighted.nn@cell.names))
# #
# # BiocManager::install("LiNk-NY/MultiAssayExtra")
# # library(MultiAssayExtra)
# # mae1 <- MultiAssayExtra(experiments=ExperimentList(RNA=wt211@assays$RNA@counts,
# #                                                     ATAC=wt211@assays$ATAC@counts),
# #                         reducedDims=c())
#
#
#
#
# library(SingleCellMultiModal)
# library(MultiAssayExperiment)
# library(SingleCellExperiment)
# mome <- SingleCellMultiModal::scMultiome(dry.run=FALSE)
#
# colData(mome)
# colData(experiments(mome)$rna)
#
#
