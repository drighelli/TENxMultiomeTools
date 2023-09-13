
#' wrapMacs2
#'
#' @param sce A SingleCellExperiment object to populate with macs2 peaks
#' @param outdir character for the path where to store the peak files
#' @param genome the genome to use (default is `mm`)
#' @param extsize number of bp to extend the reads (see macs2 callpeak for more
#' info)
#' @param shift it is used by macs2 to cut the ends (5') towards 5'->3'
#' direction then apply EXTSIZE to extend them to fragments.
#' @param broadCall uses the --broad option, for broad peaks.
#' @param bampath the path to the bamfiles, in case of NULL it attempts to
#' retrieve the path from the `metadata(sce)$ct_bams` (default is NULL)
#'
#' @return a SingleCellExperiment populated with the new peaks
#' @export
#'
#' @examples
#' TBD
wrapMacs2 <- function(sce, outdir, genome="mm", extsize=200, shift=-(extsize/2),
        broadCall="--broad", bampath=NULL)
{
    stopifnot(is(sce, "SingleCellExperiment"))

    if(is.null(bampath)) bampath <- metadata(sce)$ct_bams
    if(S4Vectors::isEmpty(bampath))
        stop("please provide a path for locating bam files!")

    bamfiles <- list.files(bampath, pattern="*_sorted.bam$", full.names=TRUE)
    bamfiles <- bamfiles[grep("ATAC", bamfiles)]
    if(isEmpty(bamfiles))
        stop("No BAM files detected in ", bampath)
    if(!dir.exists(outdir)) dir.create(outdir, recursive=TRUE)


    for(bfile in bamfiles)
    {
        cmd <- paste0("macs2 callpeak --treatment ", bfile,
                      " --name ", basename(bfile),
                      " --outdir ", outdir,
                      " --format BAM --gsize ", genome,
                      " ", nomodel,
                      " --extsize ", extsize,
                      " --shift ", shift,
                      " ", broadCall)

        print(cmd)
        system(cmd)
    }

    ## read peaks for all the produced files, store the peaks in the rowRanges of sce
    ## add an extra column for the cell type
}
