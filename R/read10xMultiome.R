#' Load data from a 10X Genomics Multiome experiment
#'
#' @description Creates a \linkS4class{SingleCellExperiment} from the CellRanger
#' ARC output directories for 10X Genomics data.
#'
#' @param sample.path string indicating the path to the experiment outs data.
#' If type is sparse, it must indicate the path where the sparse data are stored.
#' If type is HDF5, it must indicate the name of the HDF5 file.
#' Alternatively, the string may contain a prefix of the name for the three-file
#' system described above, where the rest of the name of each file follows
#' the standard 10X output.
#' (see \link[DropletUtils]{read10xCounts} for additional information)
#' @param sample.name the name of the sample, when NULL is the root of the
#' sample.path (default is NULL)
#' @param type String specifying the type of 10X format to read data from.
#' @param data string specifying if to load filtered or raw data
#' @param compressed Logical scalar indicating whether the text files are
#' compressed for type="sparse" or "prefix".
#' @param col.names logical indicating if to add ID to the cols of the assays.
#' @param annotation boolean specifying if to load the 10x peak annotations
#' NB this process is very slow (default FALSE)
#' @param addfrags boolean indicating to store the fragments file path
#' the path will be stored as metadata of the ATAC assay (default FALSE)
#' @param reference string indicating the main assay to use for the
#' SingleCellExperiment object
#'
#' @return a SingleCellExperiment object with the RNA (or ATAC see
#' \code{reference} parameter) as main assay and ATAC (or RNA see reference
#' parameter) as \code{altExp} assay.
#'
#' @export
#' @import SingleCellExperiment
#' @importFrom S4Vectors cbind.DataFrame
#' @importFrom DropletUtils read10xCounts
#' @importFrom rhdf5 h5read
#' @examples
#' # TBD see data dir
#'
# my.counts <- matrix(rpois(1000, lambda=5), ncol=10, nrow=100)
# my.counts <- as(my.counts, "dgCMatrix")
# cell.ids <- paste0("BARCODE-", seq_len(ncol(my.counts)))
# nfeat <- nrow(my.counts)
# gene.ids <- paste0("ENSG0000", seq_len(nfeat/2))
# gene.symb <- paste0("GENE", seq_len(nfeat/2))
# coord <- GRanges(seqnames="chr1", IRanges(seq_len(nfeat/2), seq_len(nfeat/2)+100))
# peaks <- paste0("chr1:", seq_len(nfeat/2), "-", seq_len(nfeat/2)+100)
# ids <- c(gene.ids, peaks)
# symb <- c(gene.symb, peaks)
# type <- c(rep("Gene Expression", nfeat/2), rep("Peaks", nfeat/2))
# features <- cbind(ids, symb, type, "chr1", rep(start(coord),2), rep(end(coord),2))
# write.table(x=features, "/Users/inzirio/Desktop/gDrive/works/coding/read10xMultiome/exampledata/features.tsv", sep="\t", row.names=FALSE, col.names=FALSE, quote=FALSE)
# # Writing this to file:
# tmpdir <- tempfile()
# write10xCounts(tmpdir, my.counts, gene.id=gene.ids,
#                gene.symbol=gene.symb, barcodes=cell.ids)
#
read10xMultiome <- function(
    sample.path,
    sample.name=NULL,
    type=c("auto", "sparse", "HDF5", "prefix"),
    data=c("filtered", "raw"),
    compressed = NULL,
    col.names=FALSE,
    annotation=FALSE,
    addfrags=FALSE,
    reference=c("RNA", "ATAC"))
{
    type <- match.arg(type)
    data <- match.arg(data)
    reference <- match.arg(reference)
    if ( length(grep("outs", sample.path)) == 0 )
    {
        message(paste0("No \"outs\" folder found, automatically setting",
            " 10x default outs path"))
        sample.path <- file.path(sample.path, "outs")
    }
    if(!is.null(sample.name)) data <- paste0(sample.name, data)
    filenm <- paste0(data, "_feature_bc_matrix", switch(type, HDF5 = ".h5", ""))
    filepath <- file.path(sample.path, filenm)
    if(type!="sparse")
    {
        sce <- read10xCounts(filepath, sample.names=sample.name,
            type=type, compressed=compressed, col.names=col.names)
    } else {
        sce <- .readSparse(filepath, compressed)
    }

    if ( length(grep("h5", filenm)) != 0 )
    {
        interval <- as.character(h5read(filepath,
            paste0("matrix/features/interval")))
        cse <- strsplit(interval, ":", fixed=TRUE)
        cse <- lapply(cse, function(x) { c(x[1], strsplit(x[2], "-",
            fixed=TRUE)[[1]])} )
        pdf <- do.call(rbind.data.frame, cse)
        rowData(sce) <- cbind.DataFrame(rowData(sce), pdf)
    }
    colnames(rowData(sce)) <- c("ID", "Symbol", "Type", "Chr", "Start", "End")
    rowData(sce)$Start <- as.numeric(rowData(sce)$Start)
    rowData(sce)$End <- as.numeric(rowData(sce)$End)
    sce <- splitAltExps(sce, rowData(sce)$Type, ref="Gene Expression")
    altExp(sce) <- .setRowRanges(altExp(sce))
    if (addfrags) altExp(sce) <- .setFragpath(altExp(sce), sample.path)
    metadata(sce)$path <- sample.path
    if ( is.null(sample.name) )
    {
        if ( length(grep("outs/filtered_feature_bc_matrix",
            colData(sce)$Sample)) != 0 )
        {
            colData(sce)$Sample <- gsub("/outs/filtered_feature_bc_matrix", "",
                colData(sce)$Sample, fixed=TRUE)
        }
    }
    colData(altExp(sce)) <- colData(sce)
    mainExpName(sce) <- "RNA"
    altExpNames(sce) <- "ATAC"
    if (annotation)
    {
        annofile <- file.path(sample.path, "atac_peak_annotation.tsv")
        anno <- read.table(annofile, sep="\t", header=TRUE)
        anno <- .processAnnotations(anno)
        map <- match(rowRanges(altExp(sce))$ID, paste0(anno$chrom, ":",
            anno$start, "-", anno$end))
        # mcols(rowRanges(altExp(sce)))$check <- paste0(anno$chrom, ":",
        #     anno$start, "-", anno$end)[map]
        mcols(rowRanges(altExp(sce)))$gene <- anno$gene[map]
        mcols(rowRanges(altExp(sce)))$distance <- anno$distance[map]
        mcols(rowRanges(altExp(sce)))$peak_type <- anno$peak_type[map]
    }

    sce <- .setRowNames(sce)
    altExp(sce) <- .setRowNames(altExp(sce))
    if (reference == "ATAC") sce <- swapAltExp(sce, "ATAC", withColData=FALSE)
    return(sce)
}

.setRowNames <- function(sce)
{
    if(!("ID" %in% colnames(rowData(sce))))
    {
        warning("Impossible to auto-set the rownames to the SingleCellExperiment")
    } else {
        rownames(sce) <- rowData(sce)$ID
    }
    return(sce)
}
#' @importFrom Matrix readMM
.readSparse <- function(filepath, compressed)
{
    m <- readMM(.check_compressed(paste0(filepath, "/matrix.mtx"), compressed))
    bc <- read.table(.check_compressed(paste0(filepath, "/barcodes.tsv"),
        compressed), sep="\t", col.names="Barcode")
    bc <- cbind.DataFrame("Sample"=filepath, bc)
    f <- read.table(.check_compressed(paste0(filepath, "/features.tsv"),
        compressed), sep="\t")
    sce <- SingleCellExperiment(list(counts=m),
        colData=DataFrame(bc),
        rowData=DataFrame(f))
}

.check_compressed <- function(filepath, compressed)
{
    return(paste0(filepath, ifelse(compressed, ".gz", "")))
}

.processAnnotations <- function(anno)
{
    stopifnot(any(c("chrom", "start", "end", "gene", "distance", "peak_type")
        %in% colnames(anno)))
    anstr <- paste0(anno$chrom, ":", anno$start, "-", anno$end)
    dup <- unique(anstr[which(duplicated(anstr))])
    anno1 <- anno
    for( i in 1:length(dup) )
    {
        # if((i %% 100) == 0) print(i)
        idx <- which(anstr %in% dup[i])
        genes <- paste(anno1$gene[idx], collapse=";")
        dists <- paste(anno1$distance[idx], collapse=";")
        type <- paste(anno1$peak_type[idx], collapse=";")
        df <- cbind(anno1$chrom[idx[1]], anno1$start[idx[1]], anno1$end[idx[1]],
                    genes, dists, type)
        anno1[idx[1], ] <- df
        idx <- idx[-1]
        anno1 <- anno1[-idx,]
        anstr <- anstr[-idx]
    }
    return(anno1)
}

.setFragpath <- function(sce, sample.path)
{
    metadata(sce)$fragpath <- file.path(sample.path, "atac_fragments.tsv.gz")
    return(sce)
}

.getFragpath <- function(sce)
{
    return(metadata(sce)$fragpath)
}

.setRowRanges <- function(sce)
{
    rowRanges(sce) <- .createRowRanges(rowData(sce))
    colnames(mcols(rowRanges(sce))) <- c("ID", "Symbol", "Type")
    return(sce)
}
.createRowRanges <- function(rowdata)
{
    stopifnot(any(c("ID", "Symbol", "Type", "Chr", "Start", "End") %in%
        colnames(rowdata)))
    return(GRanges(seqnames=rowdata$Chr,
        ranges=IRanges(start=rowdata$Start,
        end=rowdata$End), mcols=rowdata[,1:3]))
}

