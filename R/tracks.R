#' createTracks
#'
#' @param sce
#' @param cellType
#' @param bamdir
#' @param outdir
#' @param cellTypesCol
#' @param bcNorm normalization to use with bamCoverage
#' @param ncores
#'
#' @return
#' @importFrom SingleCellExperiment colData
#' @importFrom
#' @export
#'
#' @examples
createTracks <- function(sce, cellTypesCol="SingleR", cellType, bamdir, outdir,
    bcNorm="CPM", ncores=1)
{
    stopifnot( all( is(sce, "SingleCellExperiment"),
        (cellTypesCol %in% colnames(colData(sce))),
        ("Barcode" %in% colnames(colData(sce))),
        (cellType %in% unique(colData(sce)[[cellTypesCol]]))
        )
    )


    bc <- colData(sce)$Barcode[colData(sce)[[cellTypesCol]] == cellType]

    message("Writing ", cellType, " barcodes on file for sinto usage")
    id <- basename(unique(sce$Sample))
    if (!dir.exists(paste0(outdir,"/bc/"))) dir.create(outdir, recursive=TRUE)
    bcfn <- paste0(outdir, "/bc/", id, "_", cellType,
        "_barcodes.tsv")
    write.table(x=data.frame(bc, paste0(id,"_", cellType)), file=bcfn,
        quote=FALSE, sep="\t", row.names=FALSE, col.names=FALSE )
    bam <- list.files(bamdir, pattern="*.bam$", recursive=TRUE, full.names=TRUE)
    for( bamfile in c("gex_possorted_bam.bam$", "atac_possorted_bam.bam$"))
    {
        bami <- bam[grep(bamfile, bam)]
        prefix <- ifelse( length(grep("gex", bami)) != 0, "GEX", "ATAC" )
        bamiout <- paste0(outdir, "/", prefix, "_", id, "_", cellType)
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
        cmd <- paste0("bamCoverage --normalizeUsing ", bcNorm, " -b ", bamiout,
            " -o ", gsub(".bam", paste0("_", prefix,".bw"), bamiout), " -p ",
            ncores)
        message("executing bamCoverage to create ", cellType, " bigwig track")
        message(cmd)
        system(cmd)
    }
}
