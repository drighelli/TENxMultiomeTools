#' createBamCt
#' @description starting from a cell type annotated on the SingleCellExperiment
#' it constructs bam files with only the barcodes of the cells annotated
#' for the specified cell type.
#' Names of the original bam files are: "gex_possorted_bam.bam$",
#' "atac_possorted_bam.bam$". At the moment hard-coded.
#' @param sce a SingleCellExperiment object with cellTypesCol column in the
#' colData
#' @param cellTypesCol character indicating the colData column name to use for
#' selecting the barcodes
#' @param cellType character indicating the cell to retrieve from the
#' cellTypesCol
#' @param bamdir path to the directory where to find the original bam files
#' @param bamType one of "GEX", "ATAC", "both" indicating if the user wants to
#' produce bams for both omics or just GEX or ATAC (default is both)
#' @param outdir the path to output directory
#' @param ncores number of cores to use
#'
#' @return none
#' @export
#'
#' @examples
#' TBD
createBamCt <- function(sce, cellTypesCol="SingleR", cellType, bamdir,
                        bamType=c("both", "GEX", "ATAC"), outdir, ncores=1)
{
    stopifnot( all( is(sce, "SingleCellExperiment"),
                 (cellTypesCol %in% colnames(colData(sce))),
                 ("Barcode" %in% colnames(colData(sce))),
                 (cellType %in% unique(colData(sce)[[cellTypesCol]]))
        )
    )
    bamType <- match.arg(bamType)

    bc <- colData(sce)$Barcode[colData(sce)[[cellTypesCol]] == cellType]

    message("Writing ", cellType, " barcodes on file for sinto usage")
    id <- basename(unique(sce$Sample))
    bcfn <- paste0(outdir, "/", id, "_", cellType,
                   "_barcodes.tsv")
    write.table(x=data.frame(bc, paste0(id,"_", cellType)), file=bcfn,
                quote=FALSE, sep="\t", row.names=FALSE, col.names=FALSE )
    bam <- list.files(bamdir, pattern="*.bam$", recursive=TRUE, full.names=TRUE)

    bamfiles <- switch(bamType,
        "both"=c("gex_possorted_bam.bam$", "atac_possorted_bam.bam$"),
        "GEX"="gex_possorted_bam.bam$",
        "ATAC"="atac_possorted_bam.bam$")

    for( bamfile in bamfiles )
    {
        bami <- bam[grep(bamfile, bam)]
        prefix <- ifelse( length(grep("gex", bami)) != 0, "GEX", "ATAC" )
        bamiout <- paste0(outdir, prefix, "_", id, "_", cellType)
        dir.create(bamiout)
        cmd <- paste0("sinto filterbarcodes -b ", bami, " -c ", bcfn, " --outdir ",
                      bamiout, " -p ", ncores)
        message("executing sinto to create ", id, " ", cellType, " bam file")
        message(cmd)
        system(cmd)
        bamiout <- list.files(path=bamiout, pattern=paste0(id,"_", cellType),
                              full.names=TRUE)
        bamiout <- bamiout[grep("*.bam$",bamiout)]
        cmd <- paste0("samtools index -b ", bamiout, " -@ ", ncores)
        message("executing samtools index to sort ", cellType, " bam file")
        message(cmd)
        system(cmd)
    }
}
