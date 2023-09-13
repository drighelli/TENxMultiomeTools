#' createBamCt
#' @description starting from a cell type annotated on the SingleCellExperiment
#' it constructs bam files with only the barcodes of the cells annotated
#' for the specified cell type.
#' Names of the original bam files are: "gex_possorted_bam.bam$",
#' "atac_possorted_bam.bam$". At the moment hard-coded.
#'
#' @param sce a SingleCellExperiment object with cellTypesCol column in the
#' colData
#' @param cellTypesCol character indicating the colData column name to use for
#' selecting the barcodes
#' @param cellType character indicating the cell to retrieve from the
#' cellTypesCol
#' @param sampleName the name of the sample to use for the output files, if NULL
#' (default) if gets the sample present in the `sampleCol` `colData` column.
#' @param bamdir path to the directory where to find the original bam files
#' @param bamType one of "GEX", "ATAC", "both" indicating if the user wants to
#' produce bams for both omics or just GEX or ATAC (default is both)
#' @param outdir the path to output directory
#' @param sort logical for automatically sort the produced bam files (default TRUE)
#' @param sampleCol name of the colData column to get the name of the sample to
#' use for the output files (default is `Sample`)
#' @param bcCol name of the colData column to get the barcodes (default is `Barcode`)
#' @param ncores number of cores to use
#' @param verbose logical to print more informative messages (default FALSE)
#'
#' @return a `SingleCellExperiment` with the `outdir` stored in `metadata$ct_bams`
#' @export
#'
#'
#' @examples
#' TBD
createBamCt <- function(sce, cellTypesCol="SingleR", cellType, sampleName=NULL,
                        bamdir, bamType=c("both", "GEX", "ATAC"), outdir,
                        sort=TRUE, sampleCol="Sample", bcCol="Barcode",
                        ncores=1, verbose=FALSE)
{
    ### NOTES REMOVE HARDCODING BARCODES AND SAMPLE COLDATA NAMES
    stopifnot( all( is(sce, "SingleCellExperiment"),
        all( c(cellTypesCol, bcCol, sampleCol) %in% colnames(colData(sce))),
        (cellType %in% unique(colData(sce)[[cellTypesCol]]))))
    bamType <- match.arg(bamType)

    bc <- colData(sce)[[bcCol]][colData(sce)[[cellTypesCol]] == cellType]
    ctstr <- gsub("/","_",gsub(" ", "_", cellType))
    if(verbose) message("Writing ", cellType, " barcodes on file for sinto usage")
    if(is.null(sampleName)) id <- basename(unique(sce[[sampleCol]]))

    if (!dir.exists(paste0(outdir,"/bc/")))
        dir.create(paste0(outdir,"/bc/"), recursive=TRUE)
    bcfn <- paste0(outdir, "/bc/", id, "_", ctstr, "_barcodes.tsv")
    write.table(x=data.frame(bc, paste0(id,"_", ctstr)), file=bcfn,
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
        # bamiout <- paste0(outdir, "/", prefix, "_", id, "_", cellType)
        # dir.create(bamiout)
        cmd <- paste0("sinto filterbarcodes -b ", bami, " -c ", bcfn, " --outdir ",
                      outdir, " -p ", ncores)
        if (verbose) message("executing sinto to create ", id, " ", cellType,
                        " bam file")
        message(cmd)
        system(cmd)
        bamiout <- list.files(path=outdir, pattern=paste0(id,"_", ctstr),
                              full.names=TRUE)
        bamiout <- bamiout[grep("*.bam$",bamiout)]
        cmd <- paste0("samtools index -b ", bamiout, " -@ ", ncores)
        if(verbose) message("executing samtools index for ", cellType,
                        " bam file")
        message(cmd)
        system(cmd)
        if(sort)
        {
            if(verbose) message("executing samtools sort and index for ",
                            bamiout, " bam file")
            outsb <- paste0(gsub(".bam$", "", bamiout), "_sorted.bam" )
            cmd <- paste0("samtools sort -@ ", ncores, " ", bamiout," -o ",
                    outsb)
            message(cmd)
            system(cmd)
            cmd <- paste0("samtools index -b ", outsb, " -@ 10 ")
            message(cmd)
            system(cmd)
        }
    }
    metadata(sce)$ct_bams <- outdir
    return(sce)
}
