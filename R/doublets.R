
#' computeDoublets
#' @description computes doublets on a SingleCellExperiment object by aid of
#' scDblFinder package
#' @param sce a SingleCellExperiment object
#' @param method one of "auto", "optim", "dbr", "griffiths", see
#' \link[scDblFinder]{doubletThresholding}
#'
#' @return a SingleCellExperiment object with two additional columns dbl.calls
#' and dbl.scores
#' @export
#' @importFrom scDblFinder computeDoubletDensity doubletThresholding
#' @importFrom S4Vectors cbind.DataFrame
#' @importFrom scran getTopHVGs
#' @examples
#' # TBD
#'
computeDoublets <- function(sce, method=c("auto", "optim", "dbr", "griffiths"))
{

    ## needs to include methods for ATACseq data
    topgs <- getTopHVGs(sce, prop=0.1)
    scores <- computeDoubletDensity(sce, subset.row=topgs)
    dbl.calls <- doubletThresholding(data.frame(score=scores),
        method=method, returnType="call")
    colData(sce)$dbl.calls <- dbl.calls
    colData(sce)$dbl.scores <- scores
    return(sce)
}
#
# scelistdblrm <- lapply(seq_along(scelist), function(i)
# {
#     sce <- scelist[[i]]
#     topgs <- getTopHVGs(sce, prop=0.1)
#     scores <- computeDoubletDensity(sce, subset.row=topgs)
#     dbl.calls <- doubletThresholding(data.frame(score=scores),
#         method="griffiths", returnType="call")
#     print(paste0("------------- ", names(scelist)[i], " -------------"))
#     print(dim(sce))
#     print(summary(dbl.calls))
#     # sce <- sce[,dbl.calls=="singlet"]
#     colData(sce)$dbl.calls <- dbl.calls
#     colData(sce)$dbl.scores <- scores
#     sce
# })
# #
# saveRDS(scelistdblrm, "/mnt/callisto/shared/multiome_2022/scelist_filt_anno_dblts.rds")
#
#
# scef <- lapply(scelistdblrm, function(sce)
# {
#     sce <- sce[,colData(sce)$in_signac]
#     sce
# })
#
#
#
# scelist <- lapply(scelistdblrm, function(sce)
# {
#     set.seed(422)
#     sce <- runPCA(sce, ntop=dim(sce)[1])
#
#     set.seed(422)
#     sce <- runTSNE(sce, dimred="PCA")
#     sce
# })
#
# scelist <- lapply(seq_along(scelist), function(i)
# {
#     sce <- scelist[[i]]
#     topgs <- getTopHVGs(sce, prop=0.1)
#     scores <- computeDoubletDensity(sce, subset.row=topgs)
#     print(paste0("------------- ", names(scelist)[i], " -------------"))
#     # sce <- sce[,dbl.calls=="singlet"]
#     colData(sce) <- cbind.DataFrame(colData(sce), scores)
#     sce
# })
#
# saveRDS(scelist, "/mnt/callisto/shared/multiome_2022/scelist_filt_anno_dblts_tsne.rds")
#
# setwd("/mnt/callisto/shared/multiome_2022/plots/doublets")
#
#
# subclusters <- unique(unlist(lapply(scelist, function(sce) {unique(colData(sce)$subcluster)})))
#
#
#
# rnalist <- lapply(scelist, function(sce){
#     altExp(sce) <- NULL
#     sce
# })
# mice <- do.call(cbind, rnalist)
#
# lapply(subclusters, function(cl) {
#     if( cl %in% unique(colData(mice)$subcluster) )
#     {
#         mice <- mice[, colData(mice)$subcluster==cl]
#         tsne <- plotTSNE(mice, colour_by="dbl.calls") + theme_bw()
#         ggsave(tsne, file=paste("mice/",cl, "_mice_binary.png",sep = ""), width = 14, height = 10, units = "cm")
#         tsne <- plotTSNE(mice, colour_by="scores") + theme_bw()
#         ggsave(tsne, file=paste("mice/",cl, "_mice_scores.png",sep = ""), width = 14, height = 10, units = "cm")
#     }
# })


#' plotDoublets
#' @description plot a tsne and saves a file with the specified filename
#' @param sce a SingleCellExperiment object
#' @param doubl_by the colData column for the doublets to use for colouring the
#' t-SNE
#' @param filename the plot filename with extension
#'
#' @return a ggplot object
#' @importFrom scater plotTSNE
#' @keywords internal
#'
#' @examples
#' #TBD
#'
.plotDoublets <- function(sce, doubl_by, filename)
{
    tsne <- plotTSNE(sce, colour_by=doubl_by) + theme_bw()
    ggsave(tsne, file=filename, width = 14, height = 10, units = "cm")
    return(tsne)
}


#' plotDoublets
#'
#' @param sce a SingleCellExperiment object
#' @param ctCol the colData column with cell types
#' @param doubl_by the colData column for the doublets to use for colouring the
#' t-SNE
#' @param cellTypes a vector of cell types
#' @param path the path to save the plot
#' @param postfix an optional postfix for the filename
#'
#' @return none
#' @export
#'
#' @examples
#' #TBD
#'
plotDoubletsTSNEbyCT <- function(sce, ctCol="subcluster", doubl_by, cellTypes,
    path=".", postfix=NULL)
{
    stopifnot(is(sce, "SingleCellExperiment"))
    if ( !dir.exists(path) ) dir.create(path)

        postfix <- ifelse(is.null(postfix), doubl_by, paste0(postfix, "_", doubl_by))

    lapply(cellTypes, function(cl)
    {
        if( cl %in% unique(colData(sce)[[ctCol]]) )
        {
            filename <- file.path(path, paste0(cl, "_", postfix, ".png"))
            sce <- sce[, colData(sce)[[ctCol]]==cl]
            .plotDoublets(sce, doubl_by=doubl_by, filename)
        }
    })
}

plotDoubletsTSNEbyCTList <- function(scelist, ctCol="subcluster", doubl_by,
    cellTypes, path=".", postfix=NULL)
{
    stopifnot(is(scelist, "list"))
    lapply(seq_along(scelist), function(i)
    {
        sce <- scelist[[i]]
        postfix <- ifelse(is.null(postfix), names(scelist)[i],
            paste0(postfix, "_", names(scelist)[i]))
        plotDoubletsTSNEbyCT(sce, ctCol, doubl_by, cellTypes, path, postfix)
    })
}

#############################
# plotDoubletsTSNEbyCTList(scelist, doubl_by="scores", cellTypes=subclusters)
# plotDoubletsTSNEbyCTList(scelist, doubl_by="dbl.calls", cellTypes=subclusters)
#
# plotDoubletsTSNEbyCT(mice, doubl_by="scores", cellTypes=subclusters, path="mice")
# plotDoubletsTSNEbyCT(mice, doubl_by="dbl.calls", cellTypes=subclusters, path="mice")
#
# tsne <- plotTSNE(mice, colour_by="scores") + theme_bw()
# ggsave(tsne, file="mice/all_CT_tSNE.png", width = 14, height = 10, units = "cm")
#
# library(batchelor)
# mnn <- fastMNN(mice, batch=colData(mice)$Sample ,BSPARAM=BiocSingular::RandomParam(deferred=TRUE))
#
# library(BiocNeighbors)
# mnn <- runTSNE(mnn, dimred="corrected", external_neighbors = TRUE, BNPARAM=AnnoyParam())
# saveRDS(mnn, "/mnt/callisto/shared/multiome_2022/mnn_tsne.rds")
# saveRDS(mice, "/mnt/callisto/shared/multiome_2022/mice.rds")
#
#
# lapply(seq_along(scelist), function(i)
# {
#     filename <- paste0("/mnt/callisto/shared/multiome_2022/plots/doublets/",
#         names(scelist)[i], "_ct.png")
#     tsne <- plotTSNE(scelist[[i]], colour_by="subcluster") + theme_bw()
#     ggsave(tsne, file=, width = 14, height = 10, units = "cm")
#
#     filename <- paste0("/mnt/callisto/shared/multiome_2022/plots/doublets/",
#                        names(scelist)[i], "_.png")
#     tsne <- plotTSNE(scelist[[i]], colour_by="subcluster") + theme_bw()
#     ggsave(tsne, file=, width = 14, height = 10, units = "cm")
# })
# ##################

#
#
#
# dbllist <- lapply(seq_along(scelist), function(i)
# {
#     sce <- scelist[[i]]
#     topgs <- getTopHVGs(sce, prop=0.1)
#     scores <- computeDoubletDensity(sce, subset.row=topgs)
#     dbl.calls <- doubletThresholding(data.frame(score=scores),
#         method="griffiths", returnType="call")
#     return(summary(dbl.calls))
#
# })
# names(dbllist) <- names(scelist)
#
# df <- matrix(unlist(dbllist), byrow=TRUE, ncol=2)
# colnames(df) <- c("singlets", "doublets")
# rownames(df) <- names(dbllist)
# df <- as.data.frame(df)
# df$total <- df$singlets+df$doublets
# df <- df[,c(3,1,2)]
#
#
# sce <- scelist$REP1_SD
#
# dbl.out <- findDoubletClusters(sce, clusters=sce$subcluster)
# rownames(dbl.out)[isOutlier(dbl.out$num.de, type="lower", log=TRUE)]
#
# topgs <- getTopHVGs(sce, prop=0.1)
# scores <- computeDoubletDensity(sce)
# dbl.calls <- doubletThresholding(data.frame(score=scores),
#     method="griffiths", returnType="call", stringency=1)
# table(dbl.calls)
# sce$DoubletScore <- scores
#
# plotTSNE(sce, colour_by="DoubletScore")
#
# sce$DblCalls <- dbl.calls
# table(sce$subcluster,sce$DblCalls)
#
#
# sce1atac <- altExp(scelist$REP1_SD)
#
#
# sce1atacdbl <- scDblFinder(sce1atac, artificialDoublets=20, aggregateFeatures=TRUE, nfeatures=100, processing="default")
#
#
#
#










