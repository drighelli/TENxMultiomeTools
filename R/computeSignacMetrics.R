#' computeSignacMetrics
#'
#' @param sce a SingleCellExperiment object created with read10xMultiome
#' @param annotation a GRange object with gene annotations
#'
#' @return a SingleCellExperiment object
#' @importFrom Signac NucleosomeSignal TSSEnrichment CreateChromatinAssay
#' @importFrom Seurat CreateSeuratObject
#' @export
#'
#' @examples
#' #TBD
computeSignacMetrics <- function(sce, annotation)
{
    if(is.null(.getFragpath(altExp(sce)))) stop("Please add a fragpath to the sce object")
    stopifnot(all(is(sce, "SingleCellExperiment"), is(annotation, "GRanges")))
    # create a Seurat object containing the RNA data
    rna <- counts(sce)
    if (is(rna, "DelayedArray")) rna <- Matrix(rna)
    if(is.null(rownames(rna))) rownames(rna) <- rowData(sce)$ID
    if(is.null(colnames(rna))) colnames(rna) <- colData(sce)$Barcode
    seu <- CreateSeuratObject(
        counts = rna,
        assay = "RNA"
    )
    rm(rna)
    # create ATAC assay and add it to the object
    atac <- counts(altExp(sce))
    if(is.null(rownames(atac))) rownames(atac) <- rowData(altExp(sce))$Symbol
    if(is.null(colnames(atac))) colnames(atac) <- colData(sce)$Barcode
    seu[["ATAC"]] <- CreateChromatinAssay(
        counts = atac,
        sep = c(":", "-"),
        # fragments = .getFragpath(sce),
        fragments =.getFragpath(altExp(sce)),

        annotation = annotation
    )
    rm(atac)
    DefaultAssay(seu) <- "ATAC"
    seu <- NucleosomeSignal(seu)
    seu <- TSSEnrichment(seu)

    colData(sce) <- cbind.DataFrame(colData(sce), DataFrame(seu@meta.data))
    return(sce)
}
