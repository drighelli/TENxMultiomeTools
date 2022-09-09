#' computePseudoBulk
#' @description this function computes pseudobulk counts based on the cell types
#' indicated in the colData ctCol and using the assayName assay
#'
#' @param sce a SingleCellExperiment object
#' @param ctCol a string indicating the name of the column with the cell types
#' to use for the bulk transformation
#' @param assayName default is counts
#' @param ... additional arguments to pass to \link[scuttle]{aggregateAcrossCells}
#'
#' @return a SingleCellExperiment object with pseudobulk counts
#' @importFrom scuttle aggregateAcrossCells
#' @export
#'
#' @examples
#' # TBD
computePseudoBulk <- function(sce, assayName="counts", ctCol, ...)
{
    stopifnot(all(
        is(sce, "SingleCellExperiment"),
        (ctCol %in% colData(sce))
    ))

    ## ... to be processed and used as input to the following function

    sce_pseudo <- applySCE(sce, FUN="aggregateAcrossCells", ids=sce[[ctCol]],
        use.assay.type=assayName)
    return(sce_pseudo)
}

