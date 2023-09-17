#' splitSCE
#'
#' @param sce a SingleCellExperiment
#' @param by the colname of the rowData or colData to refer to when splitting
#' the `sce`
#' @param col boolean (default `TRUE`) indicating if the `by` column is on the
#' `colData` or `rowData`
#'
#' @return a named list of SingleCellExpreriment objects
#' @export
#'
#' @examples
#' TBD
splitSCE <- function(sce, by, col=TRUE)
{
    stopifnot(any(by %in% colnames(colData(sce)),
                  by %in% colnames(rowData(sce))))

    if (col)
    {
        ap <- colData(sce)[[by]]
    } else {
        ap <- rowData(sce)[[by]]
    }
    iter <- unique(ap)
    scelist <- lapply(iter, function(it){
        ifelse(col, return(sce[,ap==it]), return(sce[ap==it,]))
    })
    names(scelist) <- iter
    return(scelist)
}
