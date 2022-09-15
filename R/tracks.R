#' createTracks
#'
#' @param sce
#' @param cellType
#' @param bamdir
#' @param outdir
#' @param cellTypesCol
#' @param ncores
#'
#' @return
#' @importFrom SingleCellExperiment colData
#' @importFrom
#' @export
#'
#' @examples
createTracks <- function(sce, cellTypesCol="singleR", cellType, bamdir, outdir,
    ncores=1)
{
    stopifnot( all( is(sce, "SingleCellExperiment"),
        (cellTypesCol %in% colnames(colData(sce))),
        ("Barcode" %in% colnames(colData(sce))),
        (cellType %in% unique(colData(sce)[[cellTypesCol]]))
        )
    )


    bc <- colData(sce)$Barcode[colData(sce)[[cellTypesCol]] == cellType]

    message("Writing ", cellType, " barcodes on file for sinto usage")
    bcfn <- paste0(outdir, "/", cellType, "_barcodes.tsv")
    write.table(bcfn, x=data.frame(bc, cellType),
        file=, quote=FALSE, sep="\t",
        row.names=FALSE, col.names=FALSE )
    bam <- list.files(bamdir, pattern="*.bam$", recursive=TRUE, full.names=TRUE)
    for( bamfile in c("gex_possorted_bam.bam$", "atac_possorted_bam.bam$"))
    {
        bami <- bam[grep(bamfile, bam)]
        prefix <- ifelse( length(grep("gex", bami)) != 0, "GEX_", "ATAC_" )
        bamiout <- paste0(outdir, prefix, cellType)
        cmd <- paste0("sinto filterbarcodes -b ", bami, " -c ", bcfn, " --outdir ",
                      bamiout)
        message("executing sinto to create ", cellType, " bam file")
        message(cmd)
        system(cmd)
        bamiout <- list.files(path=bamiout, pattern=cellType, full.names=TRUE)
        bamiout <- bamiout[grep("*.bam$",bamiout)]
        cmd <- paste0("samtools index -b ", bamiout, " -@ ", ncores)
        message("executing samtools index to sort ", cellType, " bam file")
        message(cmd)
        system(cmd)
        cmd <- paste0("bamCoverage -b ", bamiout, " -o ",
            gsub(".bam", paste0("_",prefix,".bw"), bamiout), " -p ",
            ncores)
        message("executing bamCoverage to create ", cellType, " bigwig track")
        message(cmd)
        system(cmd)
    }

}
