#' pseudoBulkRNA
#' @description this function computes pseudobulk counts based on the cell types
#' indicated in the colData ctCol
#' @param sce a SingleCellExperiment object
#' @param ctCol a string indicating the name of the column with the cell types
#' to use for the bulk transformation
#'
#' @return a SingleCellExperiment object with pseudobulk counts
#' @importFrom scuttle aggregateAcrossCells
#' @export
#'
#' @examples
#' # TBD
pseudoBulkRNA <- function(sce, ctCol)
{
    stopifnot(all(
        is(sce, "SingleCellExperiment"),
        (ctCol %in% colData(sce))
    ))
    sce <- aggregateAcrossCells(sce, ids=sce[[ctCol]])

    return(sce)
}

