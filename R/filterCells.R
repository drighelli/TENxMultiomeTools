#' computeFilterCellsMetrics
#'
#' @description computes metrics on the measure previously computed with
#' computeSignacMetrics creating TRUE/FALSE columns in the sce colData
#' indicating if the cell passes or not the checks
#' @param sce a SingleCellExperiment object
#' @param metric one of quantile/isOutlier/both
#' @param lowQuantThr lower threshold for quantiles to label the cells to
#' filter out (default is 5%)
#' @param highQuantThr higher threshold for quantiles to label the cells to
#' filter out (default is 90%)
#' @param nucleosomeThr threshold for the nucleosome signal (default is 2)
#' @param TSSThr threshold for TSS enrichment (default is 1)
#'
#' @return a SingleCellExperiment with one ligical column in colData per each
#' metric
#' @importFrom scuttle isOutlier
#'
#' @export
#'
#' @examples
#' #TBD
computeFilterCellsMetrics <- function(sce,
    metric=c("quantile", "isOutlier", "both"),
    lowQuantThr="5%", highQuantThr="90%",
    nucleosomeThr=2, TSSThr=1)
{
    stopifnot( all( c("nCount_ATAC", "nCount_RNA", "nFeature_RNA",
                      "nucleosome_signal", "TSS.enrichment") %in%
                        colnames(colData(sce))
        )
    )
    if(is.null(rownames(colData(sce))))
    {
        rownames(colData(sce)) <- colData(sce)$Barcode
    }
    metric=match.arg(metric)

    nsth <- nucleosomeThr
    tssth <- TSSThr

    if( (metric=="quantile") || (metric=="both") )
    {
        qnca <- quantile(sce$nCount_ATAC, c(.01, .05, .10, .90, .95), na.rm=TRUE)
        qncr <- quantile(sce$nCount_RNA, c(.01, .05, .10, .90, .95), na.rm=TRUE)
        qnfr <- quantile(sce$nFeature_RNA, c(.01, .05, .10, .90, .95), na.rm=TRUE)
        qns <- quantile(sce$nucleosome_signal, c(.01, .05, .10, .90, .95, .99), na.rm=TRUE)
        qtsse <- quantile(sce$TSS.enrichment, c(.01, .05, .10, .90, .95), na.rm=TRUE)

        sce$nCount_ATAC_passed <- FALSE
        sce$nCount_RNA_passed <- FALSE
        sce$nFeature_RNA_passed <- FALSE
        sce$nucleosome_passed <- FALSE
        sce$TSS_passed <- FALSE

        filtered <- subset(
            x = colData(sce),
            subset = nCount_ATAC < qnca[highQuantThr] &
                nCount_ATAC >  qnca[lowQuantThr]
        )

        idx <- which(rownames(colData(sce)) %in% rownames(filtered))
        sce$nCount_ATAC_passed[idx] <- TRUE

        filtered <- subset(
            x = colData(sce),
            subset =    nCount_RNA < qncr[highQuantThr] &
                nCount_RNA >  qncr[lowQuantThr]
        )

        idx <- which(rownames(colData(sce)) %in% rownames(filtered))
        sce$nCount_RNA_passed[idx] <- TRUE

        filtered <- subset(
            x = colData(sce),
            subset = nFeature_RNA > qnfr[lowQuantThr]
        )

        idx <- which(rownames(colData(sce)) %in% rownames(filtered))
        sce$nFeature_RNA_passed[idx] <- TRUE

        filtered <- subset(
            x = colData(sce),
            subset = nucleosome_signal < nsth
        )
        idx <- which(rownames(colData(sce)) %in% rownames(filtered))
        sce$nucleosome_passed[idx] <- TRUE

        filtered <- subset(
            x = colData(sce),
            subset = TSS.enrichment > tssth
        )
        idx <- which(rownames(colData(sce)) %in% rownames(filtered))
        sce$TSS_passed[idx] <- TRUE
        sce$in_signac <- (sce$nCount_ATAC_passed &
                              sce$nCount_RNA_passed &
                              sce$nFeature_RNA_passed &
                              sce$nucleosome_passed &
                              sce$TSS_passed)
    }
    if( (metric=="isOutlier") || (metric=="both") )
    {
        require("scuttle")
        sce$nCount_ATAC_passed_io <- !isOutlier(sce$nCount_ATAC)
        sce$nCount_RNA_passed_io <- !isOutlier(sce$nCount_RNA)
        sce$nFeature_RNA_passed_io <- !isOutlier(sce$nFeature_RNA)
        sce$nucleosome_passed_io <- !isOutlier(sce$nucleosome_signal)
        sce$TSS_passed_io <- !isOutlier(sce$TSS.enrichment)
        sce$passed_io <- (sce$nCount_ATAC_passed_io &
                              sce$nCount_RNA_passed_io &
                              sce$nFeature_RNA_passed_io &
                              sce$nucleosome_passed_io &
                              sce$TSS_passed_io)
    }
    if(metric=="both")
    {
        sce$passed_both <- (sce$in_signac & sce$passed_io)
    }
    return(sce)
}

#' filterDoubtfulCells
#' @description
#' It filters out all the cells imputed to be doublets. It takes for granted that
#' some labels have been assigned to the cells and a kind of `delta.next` score
#' (see `dnext` parameter for info) is available in `colData(sce)`.
#' Specifically, it filters out all the cells with a score in the `dnext.col`
#' less than `dnext`.
#' In the same manner it filters out the imputed doublets with a value in
#' the `dbl.col` equal to `dbl.val`. In case the `dbl.val` is a numeric value,
#' it filters out all the cells with a value higher than it.
#'
#' @param sce SingleCellExperiment with assigned labels (it needs the
#' `dnext` column in the `colData`) and, optionally, computed doublets
#' @param dnext.fl logical indicating if to filter for `dnext`
#' @param dnext the threshold (default is 0.03) for the `dnext.col` score to
#' keep the doubtful cells (see \lin[SingleR]{SingleR} `delta.next` column in
#' its returned predictions)
#' @param doubls logical for filtering doublets cells
#' @param dbl.val the value to use for the cells, if character (default is `singlet`)
#' it keeps the cells with this value from the `dnext.col`.
#' If numeric, it keeps the cells with a value lesser that this value.
#' @param dnext.col the `colData(sce)` column to check the `dnext`
#' @param dbl.col the `colData(sce)` column to check the `dbl.val`
#'
#' @return a SingleCellExperiment with filtered cells
#' @export
#'
#' @examples
#' TBD
filterDoubtfulCells <- function(sce, dnext.fl=TRUE, dnext=0.03,
                            doubls=TRUE, dbl.val="singlet",
                            dnext.col="SR.deltanext",
                            dbl.col="dbl.calls")
{
    stopifnot(all(
        is(sce, "SingleCellExperiment"),
        (dnext.col %in% colnames(colData(sce))),
        ifelse(doubls, (dbl.col %in% colnames(colData(sce))), TRUE)
    ))

    if(dnext.fl)
    {
        sce <- sce[, sce[[dnext.col]] > dnext]
    }

    if (doubls)
    {
        if (!is.numeric(dbl.val)) {
            sce <- sce[[dbl.col==dbl.val]]
        } else {
            sce <- sce[[dbl.col<dbl.val]]
        }
    }
    return(sce)
}

# computeFilterCellsMetricsOnList <- function(scelist,
#     metric=c("quantile", "is.outlier", "both"),
#     lowQuantThr="5%", highQuantThr="90%",
#     nucleosomeThr=2, TSSThr=1)
# {
#     stopifnot( all( c("nCount_ATAC", "nCount_RNA", "nFeature_RNA",
#                     "nucleosome_signal", "TSS.enrichment") %in%
#                     colnames(colData(scelist[[1]]))
#                 )
#             )
#
#     metric=match.arg(metric)
#
#     nsth <- nucleosomeThr
#     tssth <- TSSThr
#     nms <- names(scelist)
#     scelist <- lapply(seq_along(scelist), function(i)
#     {
#         computeFilterCellsMetrics(scelist[[i]])
#     })
#     names(scelist) <- nms
#     return(scelist)
# }

#' plotFilteredCells
#' @description creates violin plots per each metric computed with
#' computeFilterCellsMetrics
#' @param sce A SingleCellExperiment Object
#' @param filterCellsBy one of all/custom, in case of custom specify a custom
#' column name with customColName
#' @param inout one of in/out to consider the kept/removed cells
#' @param customColName the colData column to use to subset the SCE object to
#' consider
#' @param name the name to print on top of the plot
#'
#' @return a ggplot2 object
#' @importFrom scater plotColData
#' @importFrom ggpubr ggarrange annotate_figure text_grob
#' @importFrom ggplot2 ggtitle theme element_text
#' @export
#'
#' @examples
#' # TBD
#'
plotFilteredCells <- function(sce,
    filterCellsBy=c("all", "custom"),
    inout=c("out", "in"),
    customColName=c("in_signac", "passed_both"), name)
{
    inout=match.arg(inout)
    filterCellsBy=match.arg(filterCellsBy)
    sce <- .logicalToFactor(sce)
    switch(filterCellsBy,
        "all"={
            mask=rep(TRUE, dim(sce)[2])
        },
        "custom"={
            stopifnot( all(!is.null(customColName),
                (customColName %in% colnames(colData(sce))) ) )
            col <- colData(sce)[[customColName]]
        }
    )
    sceapp <- .subsetcolData(sce, "nCount_ATAC_passed")
    gg1 <- plotColData(sceapp[,mask], y=c("nCount_ATAC"),
        colour_by="filtered_in") +
        ggtitle("nCount_ATAC") +
        theme(plot.title = element_text(size = 8, face = "bold"))

    sceapp <- .subsetcolData(sce, "nCount_RNA_passed")
    gg2 <- plotColData(sceapp[,mask], y=c("nCount_RNA"),
        colour_by="filtered_in") +
        ggtitle("nCount_RNA")+
        theme(plot.title = element_text(size = 8, face = "bold"))

    sceapp <- .subsetcolData(sce, "nCount_ATAC_passed")
    gg3 <- plotColData(sceapp[,mask], y=c("nFeature_RNA"),
        colour_by="filtered_in") +
        ggtitle("nFeature_RNA")+
        theme(plot.title = element_text(size = 8, face = "bold"))

    sceapp <- .subsetcolData(sce, "nucleosome_passed")
    gg4 <- plotColData(sceapp[,mask], y=c("nucleosome_signal"),
        colour_by="filtered_in") +
        ggtitle("nucleosome_signal")+
        theme(plot.title = element_text(size = 8, face = "bold"))

    sceapp <- .subsetcolData(sce, "TSS_passed")
    gg5 <- plotColData(sceapp[,mask], y=c("TSS.enrichment"),
        colour_by="filtered_in") +
        ggtitle("TSS.enrichment")+
        theme(plot.title = element_text(size = 8, face = "bold"))

    gg6 <- ggarrange(gg1, gg2, gg3, gg4, gg5, nrow=1, common.legend=TRUE,
        legend="bottom")

    gg6 <- annotate_figure(gg6, top = text_grob(paste0(name," ", inout," cells"),
            color = "red", face = "bold", size = 14))

    # ggsave(filename=paste0("../qc_filtered_in/", names(scelist)[i], "combined_coloured.png"), plot=gg6)
    return(gg6)
}

.subsetcolData <- function(sce, colname)
{
    colData(sce)$filtered_in <- colData(sce)[,colname]
    return(sce)
}

.logicalToFactor <- function(sce)
{
    for(i in 1:dim(colData(sce))[2])
    {
        if( is.logical(colData(sce)[,i]) )
        {
            colData(sce)[,i] <- as.factor(colData(sce)[,i])
            colData(sce)[,i] <- relevel(colData(sce)[,i], ref="TRUE")
        }
    }
    return(sce)
}

#' @importFrom ggplot2 ggplot_gtable
.get_legend <- function(a.gplot)
{
    tmp <- ggplot_gtable(ggplot_build(a.gplot))
    leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
    legend <- tmp$grobs[[leg]]
    return(legend)
}
