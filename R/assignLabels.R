#' assignLabels
#' @description assigns cell type labels by using a label assignment
#' approach (default is SingleR)
#' @param sce a SingleCellExperiment object
#' @param reference a SingleCellExperiment object for the cell types
#' reference
#' @param refColLab the column in the reference colData to use for the
#' cell types assignment
#' @param ... additional arguments to pass to the method
#'
#' @return a SingleCellExperiment object with a new column in the colData
#' by the name of the method applied (default is SingleR)
#' @importFrom SingleR SingleR
#' @export
#'
#' @examples
#' # TBD
#'
assignLabels <- function(sce, reference, refColLab="SingleR", ...)
{
    stopifnot(all(is(sce, "SingleCellExperiment"),
        is(reference, "SingleCellExperiment"),
        (refColLab %in% colnames(colData(reference)))))
    preds <- SingleR(sce, ref=reference,
        labels=reference[[refColLab]], ...)
    sce$SingleR <- preds$labels
    return(sce)
}
