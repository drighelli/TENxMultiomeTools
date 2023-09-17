
#' plotPeaksSamplesDEScan2
#'
#' @param grl a GenomicRangesList with `grcol` stored in each GenomicRanges metadata
#' (default is `k-carriers`, as from `DEScan2::finaRegions` function)
#' @param grcol the column where to take the number of peaks in each GenomicRanges
#' @param title the title of the final plot
#'
#' @return a ggplot object
#' @import ggplot2
#' @import GenomicRanges
#' @import ggpubr
#' @export
#'
#' @examples
#' TBD
plotPeaksSamplesDEScan2 <- function(grl, grcol="k-carriers", title="Peaks-samples")
{
    smax <- max(unlist(lapply(grl, function(gr) {return(max(gr$`k-carriers`))})))
    kdflist <- lapply(grl, function(gr)
    {
        kk <- vector()
        for( k in c(0:1) )
        {
            kk[k+1] <- sum(mcols(gr)[[grcol]]>k)
        }

        kdf <- data.frame(k=c(1:smax), peaks=kk)
        kdf
    })
    kunl <- unlist(kdflist)
    p <- kunl[grep("peaks",names(kunl))]
    ymax <- max(p)+1000

    gglist <- list()
    for(i in seq_along(kdflist))
    {
        kdf <- kdflist[[i]]
        gglist[[i]] <-
            ggplot(data=kdf, aes(x=k, y=peaks)) +
            ylab("peaks") +
            xlab("samples") +
            # xlim(1, 12)+
            #ylim(-100,30000)+
            # geom_smooth(se=FALSE) +
            geom_point() +
            geom_line() +
            scale_x_continuous(breaks = seq(0, smax, by = 1)) +
            scale_y_continuous(breaks = seq(0, ymax, by = 2000))+
            ggtitle(paste0(names(kdflist)[i])
            )
    }
    annotate_figure(ggarrange(plotlist=gglist), top=title)
}
