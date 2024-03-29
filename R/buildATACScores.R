#' buildATACScores
#' @description Computes a `score` for the peaks stored in the ATAC exp and
#' for the cell type indicated as input.
#' In case of multiple columns it simply computes a sum of the values for the
#' counts/logNormalizesCounts and stores the results in a new score column in
#' the rowRanges of the ATAC experiment.
#'
#' @param sce object with cellTypesCol column in the
#' colData
#' @param cellTypesCol character indicating the colData column name to use for
#' selecting the barcodes
#' @param cellType character indicating the cell to retrieve from the
#' cellTypesCol
#' @param assayName name of the ATAC assay (default is ATAC)
#' @param lncounts logical for computing the score on the counts or the
#' lognormalized counts (logNormCouts is used), default is FALSE
#'
#' @return a SCE with a `score` column in the rowRanges for the atac exp
#' @importFrom SingleCellExperiment mainExpName altExpNames swapAltExp logcounts
#' counts
#' @importFrom SummarizedExperiment rowRanges
#' @export
#'
#' @examples
#' TBD
buildATACScores <- function(sce, cellTypesCol="SingleR", cellType,
                            assayName="ATAC", lncounts=FALSE)
{
    stopifnot(all(is(sce, "SingleCellExperiment"),
            any(c(mainExpName(sce), altExpNames(sce)) %in% assayName),
            cellTypesCol %in% colnames(colData(sce))))
    ct0<-FALSE
    if(!cellType %in% unique(colData(sce)[[cellTypesCol]]))
    {
        ct0 <- TRUE
        warning("The ", cellType, " is not available for the experiment\n",
                "A 0 score is provided in this case!")
    }

    idx <- which(c(mainExpName(sce), altExpNames(sce)) %in% assayName)
    assayNow <- mainExpName(sce)
    if(idx==2) {sce <- swapAltExp(sce, name=assayName)}

    if (lncounts) sce <- logNormCounts(sce)

    if(!ct0){
        ## This needs to be changed with another wrapper external function
        ## skipping this case
        scect <- sce[,colData(sce)[[cellTypesCol]]==cellType]

        if(dim(scect)[2]>1)
        {
            if(lncounts)
            {
                rowRanges(scect)$score <- rowSums(as.matrix(logcounts(scect)))
            } else {
                rowRanges(scect)$score <- rowSums(as.matrix(counts(scect)))
            }

        } else
        {
            rowRanges(scect)$score <- counts(scect)[,1]
        }
    }else{
        rowRanges(scect)$score <- 0
    }

    rowRanges(sce) <- rowRanges(scect)
    if(idx==2) {sce <- swapAltExp(sce, name=assayNow)}
    return(sce)
}
