# library(SingleCellExperiment)
# library(scRNAseq)
# library(scran)
# library(scater)
# library(flexmix)
# library(BiocParallel)
# library(miQC)
# library(batchelor)
# library(SingleR)
# library(pheatmap)
# library(EnsDb.Mmusculus.v79)
# library(biomaRt)
# library(BiocNeighbors)
# library(HDF5Array)
#
# ## reference dataset: 10x genomics aggiornato####
# allen.10xgenomics <- loadHDF5SummarizedExperiment(dir="/mnt/callisto/Zuin", prefix="Allen_mm_21")
# rownames(allen.10xgenomics) <- rowData(allen.10xgenomics)$X
# allen.10xgenomics <- as(allen.10xgenomics, "SingleCellExperiment")
# names(assays(allen.10xgenomics)) <- c("counts")
# colnames(allen.10xgenomics) <- allen.10xgenomics$sample_name
#
# # select labels: hippocampus, non neuronal
# allen.20.hippo <- allen.10xgenomics[,!is.na(allen.10xgenomics$subclass_label) & allen.10xgenomics$subclass_label!=""
#                                     & allen.10xgenomics$subclass_label!="V3d"& allen.10xgenomics$subclass_label!="Meis2"
#                                     & allen.10xgenomics$subclass_label!="CR"& allen.10xgenomics$subclass_label!="Car3"]
# allen.20.hippo <- allen.20.hippo[,-which(grepl("L2", allen.20.hippo$subclass_label))]
# allen.20.hippo <- allen.20.hippo[,-which(grepl("L3", allen.20.hippo$subclass_label))]
# allen.20.hippo <- allen.20.hippo[,-which(grepl("L5", allen.20.hippo$subclass_label))]
# allen.20.hippo <- allen.20.hippo[,-which(grepl("L4", allen.20.hippo$subclass_label))]
# allen.20.hippo <- allen.20.hippo[,-which(grepl("L6", allen.20.hippo$subclass_label))]
# allen.20.hippo <- allen.20.hippo[,-which(grepl("NP", allen.20.hippo$subclass_label))]
# allen.20.hippo <- allen.20.hippo[,-which(grepl("CT", allen.20.hippo$subclass_label))]
#
# symbol_allen <- rownames(allen.20.hippo)
#
# map <- mapIds(EnsDb.Mmusculus.v79, keys= symbol_allen, keytype = "SYMBOL", column = "GENEID")
# stopifnot(length(map) == nrow(allen.20.hippo))
# rowData(allen.20.hippo)$symbol <- symbol_allen
# rowData(allen.20.hippo)$ens <- map
# allen.20.hippo<- allen.20.hippo[!is.na(rowData(allen.20.hippo)$ens),]
# rownames(allen.20.hippo) <- rowData(allen.20.hippo)$ens
#
# ## aggregation subclass####
# allen.subclass <- aggregateAcrossCells(allen.20.hippo, use.assay.type = "counts",
#                                        id=DataFrame(label=allen.20.hippo$subclass_label))
# allen.subclass <- logNormCounts(allen.subclass)
#
# ## aggregation CA1-ProS cluster####
# allen.CA1.ProS <- allen.20.hippo[,allen.20.hippo$subclass_label=="CA1-ProS"]
#
# allen.CA1.ProS.cluster <- aggregateAcrossCells(allen.CA1.ProS, use.assay.type = "counts",
#                                                id=DataFrame(label=allen.CA1.ProS$cluster_label))
# allen.CA1.ProS.cluster <- logNormCounts(allen.CA1.ProS.cluster)
#
# save(allen.subclass, allen.CA1.ProS.cluster, file = "allen_10xgenomics_hippo_aggiornato.RData")
#
# ## multiome dataset####
# multiome <- readRDS("scelist_allexps.RDS")
#
# for (n in seq_along(multiome)) {
#     multiome[[n]] <- logNormCounts(multiome[[n]])
# }
#
# for (n in seq_along(multiome)) {
#     print(n)
#     set.seed(422)
#     multiome[[n]] <- runPCA(multiome[[n]], ntop=length(rownames(multiome[[n]])))
#
#     set.seed(422)
#     multiome[[n]] <- runTSNE(multiome[[n]], dimred="PCA")
# }
#
# ## SingleR ####
# pred.subclass <- list()
# for (n in seq_along(multiome)) {
#     print(n)
#     pred.subclass[[n]] <- SingleR(multiome[[n]], ref=allen.subclass,
#                                   labels=allen.subclass$subclass_label)
#     multiome[[n]]$subclass <- pred.subclass[[n]]$labels
# }
#
# ## MNN correction ####
# same.genes <- Reduce(intersect, list(rownames(multiome[[1]]),rownames(multiome[[2]]),
#                                      rownames(multiome[[3]]),rownames(multiome[[4]]),
#                                      rownames(multiome[[5]]),rownames(multiome[[6]]),
#                                      rownames(multiome[[7]]),rownames(multiome[[8]]),
#                                      rownames(multiome[[9]]),rownames(multiome[[10]]),
#                                      rownames(multiome[[11]]),rownames(multiome[[12]])))
#
# mnn <- fastMNN(multiome[[1]], multiome[[2]], multiome[[3]], multiome[[4]],
#                multiome[[5]], multiome[[6]],multiome[[7]], multiome[[8]],multiome[[9]],
#                multiome[[10]],multiome[[11]],multiome[[12]], d=50, k=20, subset.row=same.genes,
#                BSPARAM=BiocSingular::RandomParam(deferred=TRUE))
#
# mnn.multiome <- list()
# for (n in seq_along(multiome)) {
#     mnn.multiome[[n]] <- mnn[,mnn$batch==n]
#     mnn.multiome[[n]]$subclass <- multiome[[n]]$subclass
# }
#
# uncorrected <- cbind(mnn.multiome[[1]],mnn.multiome[[2]],mnn.multiome[[3]],mnn.multiome[[4]],
#                      mnn.multiome[[5]],mnn.multiome[[6]],mnn.multiome[[7]], mnn.multiome[[8]],mnn.multiome[[9]],
#                      mnn.multiome[[10]],mnn.multiome[[11]],mnn.multiome[[12]])
#
# mnn$subclass <- uncorrected$subclass
#
# set.seed(422)
# mnn <- runTSNE(mnn, dimred="corrected", external_neighbors = TRUE, BNPARAM=AnnoyParam())
#
# save(multiome, mnn, file = "multiome_10xgenomics_new.RData")
#
# ## t-SNE plot subclass ####
# subclass.10xgenomics <- data.frame(allen.subclass$subclass_label, allen.subclass$subclass_color)
# colnames(subclass.10xgenomics) <- c("Label", "color")
#
# # change Astro color
# subclass.10xgenomics$color[subclass.10xgenomics$color=="#665C47"] <- "#957b46"
# # change Oligo color
# subclass.10xgenomics$color[subclass.10xgenomics$color=="#53776C"] <- "#744700"
# # change SMC-Peri color
# subclass.10xgenomics$color[subclass.10xgenomics$color=="#807059"] <- "#4d7647"
# # change VLMC color
# subclass.10xgenomics$color[subclass.10xgenomics$color=="#697255"] <- "#54a04d"
# # change Endo color
# subclass.10xgenomics$color[subclass.10xgenomics$color=="#8D6C62"] <- "#c95f3f"
# # change Sncg color
# subclass.10xgenomics$color[subclass.10xgenomics$color=="#D3408D"] <- "#ffff00"
#
# subclass.10xgenomics.color <- subclass.10xgenomics$color
# names(subclass.10xgenomics.color) <- subclass.10xgenomics$Label
#
# for (n in seq_along(multiome)) {
#     tsne <- plotTSNE(multiome[[n]],colour_by="subclass")+
#         scale_color_manual(values=subclass.10xgenomics.color)
#     ggsave(tsne, file=paste("tSNEplot_multiome_mouse",n,"_10xgenomics.png",sep = ""), width = 14, height = 10, units = "cm")
# }
#
# tsne_mnn <- plotTSNE(mnn,colour_by="subclass")+
#     scale_color_manual(values=subclass.10xgenomics.color)
# ggsave(tsne_mnn, file=paste("tSNEplot_multiome_mnn_10xgenomics.png",sep = ""), width = 14, height = 10, units = "cm")
#
# ## select CA1-ProS cells ####
# ca1.pros <- list()
# for (i in seq_along(multiome)) {
#     ca1.pros[[i]] <- multiome[[i]][,multiome[[i]]$subclass=="CA1-ProS"]
# }
#
# for (n in seq_along(ca1.pros)) {
#     print(n)
#     set.seed(422)
#     ca1.pros[[n]] <- runPCA(ca1.pros[[n]], ntop=length(rownames(ca1.pros[[n]])))
#
#     set.seed(422)
#     ca1.pros[[n]] <- runTSNE(ca1.pros[[n]], dimred="PCA")
# }
#
# ## SingleR cluster ####
# pred.cluster <- list()
# for (n in seq_along(ca1.pros)) {
#     print(n)
#     pred.cluster[[n]] <- SingleR(ca1.pros[[n]], ref=allen.CA1.ProS.cluster,
#                                  labels=allen.CA1.ProS.cluster$cluster_label)
#     ca1.pros[[n]]$cluster_label <- pred.cluster[[n]]$labels
# }
#
# for (n in seq_along(ca1.pros)) {
#     ca1.pros[[n]]$cluster_label <- pred.cluster[[n]]$labels
#     ca1.pros[[n]]$cluster <- ca1.pros[[n]]$cluster_label
#     ca1.pros[[n]]$cluster[which(grepl("-ve", ca1.pros[[n]]$cluster))] <- "CA1-ve"
#     ca1.pros[[n]]$cluster[which(grepl("-do", ca1.pros[[n]]$cluster))] <- "CA1-do"
#     ca1.pros[[n]]$cluster[which(grepl("-ProS", ca1.pros[[n]]$cluster))] <- "CA1-ProS"
#     ca1.pros[[n]]$cluster[which(grepl("_CA1", ca1.pros[[n]]$cluster))] <- "CA1"
# }
#
# ## MNN correction CA1-ProS####
# same.genes <- Reduce(intersect, list(rownames(ca1.pros[[1]]),rownames(ca1.pros[[2]]),
#                                      rownames(ca1.pros[[3]]),rownames(ca1.pros[[4]]),
#                                      rownames(ca1.pros[[5]]),rownames(ca1.pros[[6]]),
#                                      rownames(ca1.pros[[7]]),rownames(ca1.pros[[8]]),
#                                      rownames(ca1.pros[[9]]),rownames(ca1.pros[[10]]),
#                                      rownames(ca1.pros[[11]]),rownames(ca1.pros[[12]])))
#
# mnn.CA1.ProS <- fastMNN(ca1.pros[[1]], ca1.pros[[2]], ca1.pros[[3]], ca1.pros[[4]],
#                         ca1.pros[[5]], ca1.pros[[6]],ca1.pros[[7]], ca1.pros[[8]],ca1.pros[[9]],
#                         ca1.pros[[10]],ca1.pros[[11]],ca1.pros[[12]], d=50, k=20, subset.row=same.genes,
#                         BSPARAM=BiocSingular::RandomParam(deferred=TRUE))
#
# mnn.replic.CA1 <- list()
# for (n in seq_along(ca1.pros)) {
#     mnn.replic.CA1[[n]] <- mnn.CA1.ProS[,mnn.CA1.ProS$batch==n]
#     mnn.replic.CA1[[n]]$cluster_label <- ca1.pros[[n]]$cluster_label
#     mnn.replic.CA1[[n]]$cluster <- ca1.pros[[n]]$cluster
# }
#
# uncorrected <- cbind(mnn.replic.CA1[[1]],mnn.replic.CA1[[2]],mnn.replic.CA1[[3]],mnn.replic.CA1[[4]],
#                      mnn.replic.CA1[[5]],mnn.replic.CA1[[6]],mnn.replic.CA1[[7]], mnn.replic.CA1[[8]],mnn.replic.CA1[[9]],
#                      mnn.replic.CA1[[10]],mnn.replic.CA1[[11]],mnn.replic.CA1[[12]])
#
# mnn.CA1.ProS$cluster_label <- uncorrected$cluster_label
# mnn.CA1.ProS$cluster <- uncorrected$cluster
#
# set.seed(422)
# mnn.CA1.ProS <- runTSNE(mnn.CA1.ProS, dimred="corrected", external_neighbors = TRUE, BNPARAM=AnnoyParam())
#
# ## t-SNE plot CA1-ProS cluster ####
# label <- c("CA1","CA1-do", "CA1-ProS","CA1-ve")
# color <- c("#da0c0c","#ffff00","#11d611","#0f0fcd")
#
# cluster.color <- data.frame(label,color)
# colnames(cluster.color) <- c("Label", "color")
#
# clust.color <- cluster.color$color
# names(clust.color) <- cluster.color$Label
#
# for (n in seq_along(ca1.pros)) {
#     tsne <- plotTSNE(ca1.pros[[n]],colour_by="cluster")+
#         scale_color_manual(values=clust.color)
#     ggsave(tsne, file=paste("tSNEplot_multiome_CA1_mouse",n,".png",sep = ""), width = 14, height = 10, units = "cm")
# }
#
# # Wfs1 gene: ENSMUSG00000039474####
# for (n in seq_along(ca1.pros)) {
#     tsne <- plotTSNE(ca1.pros[[n]],colour_by="ENSMUSG00000039474")+
#         scale_colour_gradient(low = "grey", high = "red")
#     ggsave(tsne, file=paste("tSNEplot_multiome_Wfs1_CA1_mouse",n,".png",sep = ""), width = 14, height = 10, units = "cm")
# }
#
# ## Violin plot Wfs1 on cell type#####
# assays(uncorrected)$logcounts <- cbind(logcounts(multiome[[1]]),logcounts(multiome[[2]]),
#                                        logcounts(multiome[[3]]),logcounts(multiome[[4]]),
#                                        logcounts(multiome[[5]]),logcounts(multiome[[6]]),
#                                        logcounts(multiome[[7]]),logcounts(multiome[[8]]),
#                                        logcounts(multiome[[9]]),logcounts(multiome[[10]]),
#                                        logcounts(multiome[[11]]),logcounts(multiome[[12]]))
#
# #rownames(multiome[[1]])[which(grepl("ENSMUSG00000039474", rownames(multiome[[1]])))]
# violin.plot <- plotExpression(uncorrected, rownames(uncorrected)[7894],
#                               x = "subclass")+theme(axis.text=element_text(size=5))
# ggsave(violin.plot, file=paste("violin_plot_multiome_Wfs1.png",sep = ""), width = 14, height = 10, units = "cm")
#
# ## edgeR S3-WT####
# # select S3 and WT samples
# mnn.WT.SD <- mnn[,(mnn$batch=="3"|mnn$batch=="4"|mnn$batch=="7"|mnn$batch=="8"|
#                        mnn$batch=="11"|mnn$batch=="12")]
#
# # remove Non-Neuronal labels and labels with less 500 cells
# mnn.WT.SD.clean <- mnn.WT.SD[,!(mnn.WT.SD$subcluster=="Endo"|mnn.WT.SD$subcluster=="SMC-Peri"|
#                                     mnn.WT.SD$subcluster=="Micro-PVM"|mnn.WT.SD$subcluster=="VLMC"|
#                                     mnn.WT.SD$subcluster=="Sst Chodl")]
#
# summed <- aggregateAcrossCells(mnn.WT.SD.clean,
#             ids=DataFrame(label=mnn.WT.SD.clean$subcluster, sample=mnn.WT.SD.clean$batch),
#             use.assay.type = "counts")
# summed <- scater::logNormCounts(summed)
#
# summed$Sleep <- summed$batch
# summed$Sleep[summed$Sleep=="4"|summed$Sleep=="8"|summed$Sleep=="12"] <- "WT"
# summed$Sleep[summed$Sleep=="3"|summed$Sleep=="7"|summed$Sleep=="11"] <- "S3"
# colnames(summed) <- paste0(summed$subcluster, summed$batch)
#
# df <- data.frame(colData(summed))[,c("label", "Sleep")]
#
# current <- y <- discarded <- keep <- design <- fit <- res <- res_edgeR <- FDR <- contrast <- list()
# for (i in 1:length(levels(factor(summed$subcluster)))) {
#     print(levels(factor(summed$subcluster))[i])
#     current[[i]] <- summed[,summed$subcluster==levels(factor(summed$subcluster))[i]]
#
#     y[[i]] <- DGEList(counts(current[[i]]), samples=colData(current[[i]]))
#
#     discarded[[i]] <- isOutlier(y[[i]]$samples$lib.size, log=TRUE, type="lower")
#
#     y[[i]] <- y[[i]][,!discarded[[i]]]
#
#     #keep[[i]] <- filterByExpr(y[[i]], group=current[[i]]$Sleep)
#
#     y[[i]] <- calcNormFactors(y[[i]])
#
#     #design[[i]] <- model.matrix(~y[[i]]$samples$Sleep, y[[i]]$samples)
#     design[[i]] <- model.matrix(~0+y[[i]]$samples$Sleep, y[[i]]$samples)
#     colnames(design[[i]]) <- c("S3", "WT")
#
#     contrast[[i]] <- makeContrasts(S3-WT,levels=design[[i]])
#
#     y[[i]] <- estimateDisp(y[[i]], design[[i]])
#
#     #fit[[i]] <- glmFit(y[[i]], design[[i]])
#     fit[[i]] <- glmQLFit(y[[i]], design[[i]])
#
#     #res[[i]] <- glmLRT(fit[[i]], contrast =contrast[[i]])
#     res[[i]] <- glmQLFTest(fit[[i]], contrast=contrast[[i]])
#
#     res_edgeR[[i]] <- as.data.frame(res[[i]]$table)
#
#     #res_edgeR[[i]]$Symbol <- rowData(current[[i]])$symbol
#
#     res_edgeR[[i]] <- as.data.frame(res_edgeR[[i]])[order(as.data.frame(res_edgeR[[i]])$PValue),]
#
#     FDR[[i]] <- as.data.frame(topTags(res[[i]], n=length(rownames(mnn.WT.SD.clean)))[,5])
#
#     res_edgeR[[i]] <- cbind(res_edgeR[[i]], FDR[[i]])
#     res_edgeR[[i]] <- res_edgeR[[i]][order(res_edgeR[[i]]$FDR),]
# }
