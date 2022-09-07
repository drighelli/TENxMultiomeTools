# library("Signac")
# library("SingleCellExperiment")
# # wt211 <- readRDS("signac/wt_frozen2021_signac.rds")
# # wt212 <- readRDS("signac/wt_frozen2021_2_signac.rds")
# # wt172 <- readRDS("signac/wt_frozen2017_2_signac.rds")
# #
# # addSignacToSce <- function(seu)
# # {
# #     stopifnot(is(seu, "Seurat"))
# #     rna <- seu@assays$RNA@counts
# #     atc <- seu@assays$ATAC@counts
# #     sce <- SingleCellExperiment(assays=list(RNA=rna),
# #         altExps=list(ATAC=SingleCellExperiment(assays=list(ATAC=atc))),
# #         colData= DataFrame(seu@meta.data))
# #     return(sce)
# # }
# #
# # files <- list.files(pattern="*_signac.rds")
# #
# # lapply(files, function(f)
# # {
# #     print(f)
# #     s <- readRDS(f)
# #     print(s)
# #     rna <- s@assays$RNA@counts
# #     atc <- s@assays$ATAC@counts
# #     sce <- SingleCellExperiment(assays=list(RNA=rna), altExps=list(ATAC=SingleCellExperiment(assays=list(ATAC=atc))))
# #     print(sce)
# #     saveRDS(sce, file=gsub("_signac", "", f))
# # })
# #
# # rna1 <- wt211@assays$RNA@counts
# # rna2 <- wt212@assays$RNA@counts
# # rna0 <- wt172@assays$RNA@counts
# #
# # atc0 <- wt172@assays$ATAC@counts
# # atc1 <- wt211@assays$ATAC@counts
# # atc2 <- wt212@assays$ATAC@counts
# #
# # sce211 <- SingleCellExperiment(assays=list(RNA=rna1), altExps=list(ATAC=SingleCellExperiment(assays=list(ATAC=atc1))))
# # sce212 <- SingleCellExperiment(assays=list(RNA=rna2), altExps=list(ATAC=SingleCellExperiment(assays=list(ATAC=atc2))))
# # sce172 <- SingleCellExperiment(assays=list(RNA=rna0), altExps=list(ATAC=SingleCellExperiment(assays=list(ATAC=atc0))))
# #
# # saveRDS(sce211, file="SCE/wt_frozen2021_1_sce.rds")
# # saveRDS(sce212, file="SCE/wt_frozen2021_2_sce.rds")
# # saveRDS(sce172, file="SCE/wt_frozen2017_2_sce.rds")
# #
# # sce172 <- SingleCellExperiment(assays=list(RNA=rna0))
# # sce211 <- SingleCellExperiment(assays=list(RNA=rna1))
# # sce212 <- SingleCellExperiment(assays=list(RNA=rna2))
# # sce_all <- cbind(sce172, sce211, sce212)
# #
# # colData(sce_all)$sample_id <- c(rep("wtfr172", dim(sce172)[2]), rep("wtfr211", dim(sce211)[2]),
# #                                 rep("wtfr212", dim(sce212)[2]))
# #
# # saveRDS(sce_all, file="SCE/wt_frozen_all_sce.rds")
