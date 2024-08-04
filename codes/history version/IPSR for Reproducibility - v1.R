

## install IPs detection package
## devtools::install_github("ShuoStat/findIPs", force = T)

library(findIPs)

################################################################################
#- Figure 1
################################################################################

tiff("../output/Figure 1.tiff", width = 8, height = 4, units = "in",
     res = 300, compression = "lzw")

par(family = "serif", mar = c(4, 4, 2, 0.5), mfrow = c(1, 2))

lty <- c("solid", "dashed", "dotdash")
x = 1:100
m_set <- c(0.03, 0.05, 0.07)
max_m <- max(m_set)

plot.new()
plot.window(xlim = c(0, 100), ylim = c(0,  max_m * exp(max_m / 2 - max_m * 1)))

axis(1)
axis(2)
box()

n = 100
for (i in 1:length(m_set)) {

  m = m_set[i]
  y =  exp(-m * x) * (1 - exp(m)) / (exp(-m * n) - 1)
  lines(x, y, lty = lty[i])
}

legend("topright", legend = m_set, inset = 0.02, lty = lty,
       bty = "n", title = expression(kappa))

mtext(text = "Ranks", side = 1, line = 2.3)
mtext(text = "Weights", side = 2, line = 2.3)
mtext(text = "A", side = 3, line = .3, adj = 0)

#- B, length n ---
m = 0.03
n_set = c(30, 50, 100)
x = 1
ymax <- exp(-m * x) * (1 - exp(m)) / (exp(-m * min(n_set)) - 1)

plot.new()
plot.window(xlim = c(0, max(n_set)), ylim = c(0, ymax))

axis(1)
axis(2)
box()
for (i in n_set) {
  n = i
  x = 1:n
  y =  exp(-m * x) * (1 - exp(m)) / (exp(-m * n) - 1)
  lines(x, y)
  text(max(x), min(y), labels = max(x), pos = 3, cex = 0.8)
}

legend("topright", legend = bquote(kappa ~ "= 0.03"), inset = 0.02, bty = "n")
mtext(text = "Ranks", side = 1, line = 2.3)
mtext(text = "Weights", side = 2, line = 2.3)
mtext(text = "B", side = 3, line = .3, adj = 0)

dev.off()

################################################################################
#- Using findIPs to get the leave-1-out rankings
################################################################################

#- import data, get data name
s <- dir("../data")

# get ranks, weighted ranks, and scores
getIPsFun <- function(filename, topN = 200){

  load(paste0("../data/", filename))
  name = gsub(".RData", "", filename)

  # t.testFun <- function(x, y)
  #   t.test(x ~ y, var.equal = TRUE)$p.value

  drop1ranks <- getdrop1ranks(X, y,
                              fun = "t.test",
                              topN = topN,
                              decreasing = FALSE,
                              nCores = 10)

  weightedRanks <- sumRanks(drop1ranks$origRank,
                            drop1ranks$drop1Rank,
                            method = "adaptive",
                            topN = topN)

  return(weightedRanks)
}

# obj = getIPsFun(filename = s[7], topN = 200)

drop1out.ranks <- mapply(getIPsFun, filename = s, SIMPLIFY = F)
names(drop1out.ranks) <- substr(names(drop1out.ranks), 6, 9)
# save the output
save(list = "drop1out.ranks", file = "../output/drop1out.ranks.RData")

################################################################################
#- LUSC example
#- output: Table 2, 3 and Figure 2
################################################################################

# load data
load("../output/drop1out.ranks.RData")
obj <- drop1out.ranks[["LUSC"]]

#- Table 2 and 3 ---------------------------------------------------------------
#- Table 2, original ranks

orig <- obj$drop1Rank[1:20, 1:10]
orig <- cbind(orig = 1:20, orig)
colnames(orig) <- paste0("r[", 0:10, ",j]")
rownames(orig) <- obj$origRank[1:20]
write.csv(orig, "../output/Table2.csv", row.names = F)

#- Table 3, weighted ranks

weighted <- cbind(obj$origRankWeighted[1:20],
                  obj$drop1RankWeighted[1:20, c(1:10)])
weighted <- format(round(weighted, 4), nsmall = 4)

colnames(weighted) <- paste0("f(r[", 0:10, ",j])")
rownames(weighted) <- obj$origRank[1:20]
write.csv(weighted, "../output/Table3.csv", row.names = F)

#- Figure 2 --------------------------------------------------------------------

plot.kappa.extension <- function(obj, points = T, dist = T, ylim = NULL){

  n <- nrow(obj$drop1Rank)
  x <- 1:n
  kappa = obj$kappa

  plotAdaptiveWeights(obj$kappa, n, type = "line")

  if (points == T) {
    or <- obj$drop1Rank  #or, original ranks
    all_r <- reshape2::melt(or)
    sel_r <- all_r[all_r[,1] != all_r[,3],]
    x_p = sel_r[,3]
    y_p = exp(-kappa * x_p) * (1 - exp(kappa)) / (exp(-kappa * n) - 1)
    points(sel_r[,1], y_p, cex = 0.7, col = "grey60", pch = 19)
  }

  if (dist == T){
    top = order(obj$score, decreasing = T)[1]
    y0 = obj$origRankWeighted
    y1 = obj$drop1RankWeighted[,top]

    for(i in seq_along(y0)){
      lines(c(x[i], x[i]), c(y0[i], y1[i]), lwd = 0.7)
    }
    points(x, y1, cex = 0.7, pch = 19)
  }

  y <- exp(-kappa * x) * (1 - exp(kappa))/(exp(-kappa * n) - 1)
  lines(x, y)
}

tiff("../output/Figure 2.tiff", width = 12, height = 3.5, compression = "lzw",
     units = "in", res = 600)
#
# png("../output/Figure2.png", width = 12, height = 3.5,
#      units = "in", res = 900)

par(mar = c(4, 4, 2, 1), mfrow = c(1, 3))

plotRankScatters(obj)
mtext("A", side = 3, adj = 0, line = 0.2)

plot.kappa.extension(obj)
mtext("B", side = 3, adj = 0, line = 0.2)

plotIPs(obj,  ylim = c(0, 13))
mtext("C", side = 3, adj = 0, line = 0.2)

dev.off()

################################################################################
# Clustering to explore Obs51, Figure 3, Figure A1
################################################################################

library(dendextend)
## load LUSC data

load(paste0("../data/TCGA-LUSC.RData"))
colnames(X) <- seq_len(ncol(X))

# case control without obs51
pvalues <- apply(X[,-51], 1,
                 function(x) t.test(x ~ y[-51], var.equal = TRUE)$p.value)
top200 <- rownames(X)[order(pvalues, decreasing = F)][1:200]
Z <- X[top200, ]
Z <- t(apply(Z[top200, ], 1, scale))

d <- stats::dist(t(Z))
h <- hclust(d)

col <- ifelse(y == "01", "darkcyan", "darkgoldenrod1")
labels <- rep("", length(y))
labels[51] <- "51"

tiff("../output/Figure A1.tiff", width = 6, height = 3, compression = "lzw",
     units = "in", res = 600)

par(mar = c(1, 1, 1, 1))
dplot <- as.dendrogram(h) %>%
    set("labels", labels[h$order]) %>%
    set("labels_cex", 0.7) %>%
    set("branches_k_color", col[h$order]) %>%
    # set("leaves_pch", 20) %>%
    # set("leaves_col", leaves_col[h$order]) %>%
    plot(yaxt  = "n")

dev.off()

# import oncoGenes
# https://www.oncokb.org/cancer-genes
geneList <- read.csv("./cancerGeneList.tsv", sep = "\t")[,1]
topOncoGenes <- intersect(top200, geneList)
Dat <- t(apply(X[topOncoGenes, ], 1, scale))

d <- stats::dist(t(Dat))
h <- hclust(d)

# define labels
labels <- rep("", length(y))
labels[51] <- "51"

library(dendextend)
col_branch <- ifelse(y == "01", "darkcyan", "darkgoldenrod1")

dplot <- as.dendrogram(h) %>%
    set("labels", labels[h$order]) %>%
    # set("labels_cex", 0.3) %>%
    set("branches_k_color", col_branch[h$order])

col_fun  <- circlize::colorRamp2(quantile(Dat, c(0.01, 0.5, 0.99)),
                                 c("green", "white", "red"))

tiff("../output/Figure 3.tiff", width = 8, height = 4, compression = "lzw",
     units = "in", res = 600)

library(ComplexHeatmap)
ComplexHeatmap::Heatmap(Dat,
                        col = col_fun,
                        name = "mat",
                        # border = T,
                        rect_gp = gpar(col = "gray"),
                        cluster_columns = dplot,
                        column_dend_height = unit(45, "mm"),
                        cluster_rows = T,
                        show_row_dend = F,
                        column_labels = labels,
                        row_names_gp = gpar(fontsize = 6),
                        column_names_gp = gpar(fontsize = 10, col = "black"),
                        show_heatmap_legend = FALSE)
dev.off()


#- cluster plots for more data

clusterPlot <- function(dataName, obs) {

    library(dendextend)
    ## load LUSC data

    load(paste0("../data/TCGA-", dataName, ".RData"))
    colnames(X) <- seq_len(ncol(X))

    # case control without obs51
    pvalues <- apply(X[,-obs], 1,
                     function(x) t.test(x ~ y[-obs], var.equal = TRUE)$p.value)

    top200 <- rownames(X)[order(pvalues, decreasing = F)][1:200]

    # geneList <- read.csv("./cancerGeneList.tsv", sep = "\t")[,1]
    # topOncoGenes <- intersect(top200, geneList)

    Z <- X[top200, ]
    Z <- t(apply(Z, 1, scale))

    d <- stats::dist(t(Z))
    h <- hclust(d)

    col <- ifelse(y == "01", "darkcyan", "darkgoldenrod1")
    labels <- rep("", length(y))
    labels[obs] <- paste0(obs)

    dplot <- as.dendrogram(h) %>%
        set("labels", labels[h$order]) %>%
        set("labels_cex", 0.7) %>%
        set("branches_k_color", col[h$order]) %>%
        # set("leaves_pch", 20) %>%
        # set("leaves_col", leaves_col[h$order]) %>%
        plot(yaxt  = "n")

    return(dplot)
}

#-------------------------------------------------------------------------------

tiff("../output/Figure A10.tiff", width = 8, height = 7, units = "in", res = 600,
     compression = "lzw")

par(mfrow = c(4, 2), mar = c(1, 1, 1, 1), oma = c(1, 1, 1, 1))

clusterPlot("LUSC", 51)
mtext("A, LUSC", side = 3, adj = 0)

clusterPlot("STAD", 22)
mtext("B, STAD", side = 3, adj = 0)

clusterPlot("BRCA", 338)
mtext("C, BRCA", side = 3, adj = 0)

clusterPlot("COAD", 73)
mtext("D, COAD", side = 3, adj = 0)

clusterPlot("HNSC", 93)
mtext("E, HNSC", side = 3, adj = 0)

clusterPlot("KIRC", 133)
mtext("F, KIRC", side = 3, adj = 0)

c6 <- clusterPlot("PRAD", 107)
mtext("G, PRAD", side = 3, adj = 0)

c7 <- clusterPlot("THCA", 134)
mtext("H, THCA", side = 3, adj = 0)

dev.off()

#- cancer purity

# PCA plot in R ----------------------------------------------------------------

# load("../output/drop1out.ranks.RData")
# library(dplyr)
# library(ggplot2)
#
# pcaPlot <- function(datName, IPs) {
#
#     load(paste0("../data/TCGA-", datName, ".RData"))
#     colnames(X) <- seq_len(ncol(X))
#
#     top200 <- drop1out.ranks[[datName]]$origRank
#
#     Z <- X[top200, ]
#     Z <- t(apply(Z, 1, scale))
#
#     load(paste0("../data/TCGA-", datName, ".RData"))
#
#     PCs <- prcomp(t(Z))
#     PCs  <- t(Z) %*% PCs$rotation
#
#     PCs <- PCs %>%
#         as.data.frame() %>%
#         dplyr::select(1:2) %>%
#         mutate(obs = seq_len(nrow(PCs)),
#                IPs = ifelse(obs == IPs, IPs, ""),
#                group = as.factor(y))
#
#     ggplot(PCs, aes(x = PC1, y = PC2, colour = group, label = IPs)) +
#         geom_point() +
#         geom_text(vjust = 1.5) +
#         theme_bw() +
#         theme(legend.position = "none")
# }
#
#
# pcaPlot("LUSC", 51)
# pcaPlot("STAD", 22)
# pcaPlot("BRCA", 338)
# pcaPlot("HNSC", 93)
# pcaPlot("KIRC", 133)
# pcaPlot("PRAD", 107)
# pcaPlot("THCA", 134)
#
# ggsave(filename = "../output/PCA.tiff", plot = pcaPlot, width = 5,
#        height = 4, compression = "lzw",
#        units = "in", dpi = 600)

################################################################################
#- Weights and distribution
#- Output Figure 4
################################################################################

load("../output/drop1out.ranks.RData")
plot.kappa <- function(kappa, col, n, name) {

  x <- 1:n
  plot.new()
  plot.window(xlim = c(1, n), ylim = c(0, max(kappa)))
  axis(1)
  axis(2)
  box()

  for(i in seq_along(kappa)){
    y <- exp(-kappa[i] * x) * (1 - exp(kappa[i])) / (exp(-kappa[i] * n) - 1)
    lines(x, y, col = col[i])
  }

  legend = paste0(round(kappa, 3), "(", name, ")")
  legend("topright", legend = legend, col = col,
         lty = c("solid", "solid"))
  mtext("Original rankings", side = 1, line = 2.2)
  mtext("Weighted rankings", side = 2, line = 2.2)
}

s <- c("LIHC", "PRAD")

tiff("../output/Figure 4.tiff", width = 10, height = 3, units = "in", res = 600,
     compression = "lzw")

# png("../output/Figure3.png", width = 10, height = 3, units = "in", res = 600)

par(mfrow = c(1, 3), mar = c(4, 4, 2, 2))

plotRankScatters(drop1out.ranks[[s[1]]], top = F,
                 points.arg = list(col = "grey"))
mtext(paste0("A,", s[1]), side = 3, adj = 0, line = 0.2)

plotRankScatters(drop1out.ranks[[s[2]]], top = F,
                 points.arg = list(col = "black"))
mtext(paste0("B,", s[2]), side = 3, adj = 0, line = 0.2)

kappa <- c(drop1out.ranks[[s[1]]]$kappa, drop1out.ranks[[s[2]]]$kappa)
plot.kappa(kappa, col = c("grey", "black"), n = 200, s)
mtext("C", side = 3, adj = 0, line = 0.2)

dev.off()

################################################################################
#- Exhibit Influential Observations
#- Output: Figure 5 and Figure A2
################################################################################

load("../output/drop1out.ranks.RData")

#- Figure 5 --------------------------------------------------------------------

tiff("../output/Figure 5.tiff", width = 10, height = 3 * 3, units = "in",
     res = 600, compression = "lzw")

# png("../output/Figure4.png", width = 10, height = 3 * 3, units = "in",
#      res = 600)

par(mar = c(4, 4, 1, 1), mfrow = c(3, 3), oma = c(1, 3, 1, 1))

for(i in c("STAD", "LIHC", "BRCA")){

  obj <- drop1out.ranks[[i]]
  kappa <- obj$kappa
  n <- nrow(obj$drop1Rank)

  plotRankScatters(obj)
  mtext(i, side = 2, line = 4)
  plotAdaptiveWeights(kappa, n, type = "line", ylim = c(0, 0.02))
  plotIPs(obj, ylim = c(0, 12))
}

dev.off()

#- Figure A1--------------------------------------------------------------------

tiff("../output/Figure A2.tiff", width = 10, height = 6 * 3, units = "in",
     res = 300, compression = "lzw")

# png("../output/Figure A1.png", width = 10, height = 6 * 3, units = "in",
#      res = 600)

par(mar = c(4, 4, 1, 1), mfrow = c(6, 3), oma = c(1, 3, 1, 1))

for(i in c("COAD", "HNSC", "KIRC", "LUAD", "PRAD", "THCA")){

  obj <- drop1out.ranks[[i]]
  kappa <- obj$kappa
  n <- length(obj$score)

  plotRankScatters(obj)
  mtext(i, side = 2, line = 4)
  plotAdaptiveWeights(kappa, n, type = "line", ylim = c(0, 0.025))
  plotIPs(obj, ylim = c(0, 12))
}

dev.off()

################################################################################
#- Comparison: unweighted correlations, weighted correlations, and adaptive
#- Output: Figure 6, 7 and Figure A3
################################################################################

source("./helper.R")

dist <- function(obj, method = c("spearman", "kendall"), weights = FALSE) {

  orig <- obj$drop1Rank
  n = nrow(orig)
  m <- 1:n

  D <- function(x)
    (x - m)^2

  W <- function(x)
    (x - m)^2 * (2 * n + 2 - x - m)

  K <- function(x, weights = NULL) {
      calcTopTau(m, x, posWeights = weights)
  }


  if (method == "spearman"){
      if (weights) {
          score <- apply(orig, 2, W)
          score <- colSums(score)
          score <- score / sd(score)
      } else {
          score <- apply(orig, 2, D)
          score <- colSums(score)
          score <- score / sd(score)
      }
  }

  if (method == "kendall") {
      if (weights) {
          score <- apply(orig, 2, K, weights = "DCG")
          score <- score / sd(score)
      } else {
          score <- apply(orig, 2, K, weights = NULL)
          score <- score / sd(score)
      }
  }

  return(score)
}

lollipop <- function(x, y, topn = 5, ylim = NULL){

  plot.new()
  if (is.null(ylim)) ylim = max(y)
  plot.window(xlim = c(min(x), max(x)), ylim = c(0, ylim))
  axis(1)
  axis(2)
  box()

  points(x, y, pch = 19, cex = 0.8)
  for(i in seq_along(x)){
    lines(c(x[i], x[i]), c(0, y[i]),  lwd = 0.7)
  }

  topn = order(y, decreasing = T)[1:topn]
  text(x[topn], y[topn], pos = 3, cex = 0.9, labels = topn)
  mtext("Observations", side = 1, line = 2.2)
  mtext("Scores", side = 2, line = 2.2)
}

#- Figure 6 --------------------------------------------------------------------

load("../output/drop1out.ranks.RData")

tiff("../output/Figure 6.tiff", width = 12, height = 3 * 2.5,
     units = "in", res = 600, compression = "lzw")

# png("../output/Figure 5.png", width = 10, height = 3 * 3,
#      units = "in", res = 600)

par(mfrow = c(3, 5), mar = c(3.5, 3.5, 1.5, 1.5), oma = c(1, 2, 1, 1))

for(i in c("STAD", "LIHC", "BRCA")){

  obj <- drop1out.ranks[[i]]
  lollipop(1:ncol(obj$drop1Rank), dist(obj, "spearman", FALSE), ylim = 13)
  mtext("Unweighted Spearman", side = 3, adj = 0, line = 0.2)
  mtext(i, side = 2, line = 4)

  lollipop(1:ncol(obj$drop1Rank), dist(obj, "spearman", TRUE), ylim = 13)
  mtext("Weighted Spearman", side = 3, adj = 0, line = 0.2)

  lollipop(1:ncol(obj$drop1Rank), dist(obj, "kendall", FALSE), ylim = 13)
  mtext("Unweighted Kendall", side = 3, adj = 0, line = 0.2)

  lollipop(1:ncol(obj$drop1Rank), dist(obj, "kendall", TRUE), ylim = 13)
  mtext("Weighted Kendall", side = 3, adj = 0, line = 0.2)

  plotIPs(obj, ylim = c(0, 13))
  mtext("Adaptive Weights", side = 3, adj = 0, line = 0.2)

}

dev.off()

#- Figure A4 -------------------------------------------------------------------

tiff("../output/Figure A5.tiff", width = 12, height = 6 * 2, units = "in",
     res = 300, compression = "lzw")

par(mfrow = c(6, 5), mar = c(3.5, 3.5, 1.5, 1.5), oma = c(1, 2, 1, 1))

for(i in c("COAD", "HNSC", "KIRC", "LUAD", "PRAD", "THCA")){

    obj <- drop1out.ranks[[i]]
    lollipop(1:ncol(obj$drop1Rank), dist(obj, "spearman", FALSE), ylim = 13)
    mtext("Unweighted Spearman", side = 3, adj = 0, line = 0.2)
    mtext(i, side = 2, line = 4)

    lollipop(1:ncol(obj$drop1Rank), dist(obj, "spearman", TRUE), ylim = 13)
    mtext("Weighted Spearman", side = 3, adj = 0, line = 0.2)

    lollipop(1:ncol(obj$drop1Rank), dist(obj, "kendall", FALSE), ylim = 13)
    mtext("Unweighted Kendall", side = 3, adj = 0, line = 0.2)

    lollipop(1:ncol(obj$drop1Rank), dist(obj, "kendall", TRUE), ylim = 13)
    mtext("Weighted Kendall", side = 3, adj = 0, line = 0.2)

    plotIPs(obj, ylim = c(0, 13))
    mtext("Adaptive Weights", side = 3, adj = 0, line = 0.2)
}

dev.off()

#- Further check on BRCA -------------------------------------------------------
#- Figure 7 -------------------------------------------------------------------

top.plot <- function(obj, obs, topn = 10){

  orirank <- obj$drop1Rank
  meltr <- reshape2::melt(orirank)
  plot(meltr[,1], meltr[,3], cex = 0.8, xlab = "", ylab = "",
       col = "grey", pch = 19)

  n <- nrow(obj$drop1Rank)
  points(1:n, obj$drop1Rank[,obs], cex = 0.7, pch = 19)

  top <- order(abs((1:n) - obj$drop1Rank[,obs]), decreasing = T)[1:topn]
  points(top, obj$drop1Rank[top, obs], pch = 15, cex = 1.2)
  abline(v = c(50, 100, 150), lty = "dashed")

  legend("topleft", pch = c(19, 15, 19), col = c("black", "black", "grey"),
         legend = c(paste0("obs", obs), paste0("top10 of obs", obs),"The others"))

  mtext("Original ranking", side = 1, line = 2.2)
  mtext("Leave-one-out ranking", side = 2, line = 2.2)
}

obj <- drop1out.ranks$BRCA

tiff("../output/Figure 7.tiff", width = 10, height = 2.5, units = "in",
     res = 600, compression = "lzw")
# png("../output/Figure 6.png", width = 10, height = 3, units = "in",
#      res = 600)
par(mfrow = c(1, 4), mar = c(4, 4, 1.5, 1.5))

top.plot(obj, 175, topn = 10)
top.plot(obj, 53, topn = 10)
top.plot(obj, 2, topn = 10)
top.plot(obj, 338, topn = 10)

dev.off()


################################################################################
# Figure 8, A4: Effects of IPs on GSEA
################################################################################

#- get IPs
load("../output/drop1out.ranks.RData")
ips <- lapply(drop1out.ranks, function(x) order(x$score, decreasing = T)[1])

rks <- function(X, y, ips, topn = 200){

  f <- function(x, y) t.test(x ~ y, var.equal = TRUE)$p.value
  xp <- apply(X, 1, f, y = y)

  gns <- rownames(X)[order(xp, decreasing = F)[1:topn]] #- gns, genes
  Xt <- X[gns, ]

  genelist <- list(so = sort(xp, decreasing = F)[1:topn])
  for(i in ips) {
    tmp <- sort(apply(Xt[,-i], 1, f, y = y[-i]), decreasing = F)
    genelist[[paste0("s", i)]] <- tmp
  }
  return(genelist)
}

library(clusterProfiler) # version 4.4.1
library(org.Hs.eg.db) # version 3.18.0

getgse <- function(nam){

  d <- dir("../data/")
  s <- grepl(nam, d)
  load(paste0("../data/", d[s]))
  grks <- rks(X, y, ips = ips[[nam]])

  #- to ENTREZID
  grks <- lapply(grks, function(x){
    tmp <- select(org.Hs.eg.db,
                  keys = names(x),
                  columns = c("ENTREZID"),
                  keytype = "SYMBOL")$ENTREZID
    setNames(x, tmp)
  })

  gses <- lapply(grks, function(x)
    gseGO(rev(x), OrgDb = 'org.Hs.eg.db', pvalueCutoff = 1)@result)

  list(s0 = gses[[1]][1:10,1], s1 = gses[[2]][1:10,1])
}

plot.fun <- function(s0, s1, ip) {

  plot.new()
  plot.window(xlim = c(0.3, 3.7), ylim = c(0.7, 10.2))

  for(i in 1:10){
    rect(0.5, i - 0.4, 1.5, i + 0.4)
    rect(2.5, i - 0.4, 3.5, i + 0.4)

    text(1, i, labels = s0[i])
    text(3, i, labels = s1[i])
  }
  #- add lines
  d <- na.omit(cbind(1:10, match(s0, s1)))
  if(nrow(d) > 0){
    for(i in 1:nrow(d)){
      lines(c(1.5, 2.5), d[i,], lty = "dashed")
    }
  }

  mtext("Original", side = 3, line = 0, at = 1, cex = 1.2)
  mtext(paste0("-Obs ", ip), side = 3, line = 0, at = 3, cex = 1.2)
}

dat <- c("LUSC", "STAD", "LIHC", "BRCA")

tiff("../output/Figure 8.tiff", width = 8, height = 6, units = "in", res = 600,
     compression = "lzw")

par(mfrow = c(2, 2), mar = c(1, 0, 3, 2), oma = c(1, 1, 1, 1))
j = 1
set.seed(7651189) #-2022, 1037641
for(i in dat){

  s <- getgse(i)
  s0 <- s[[1]]
  s1 <- s[[2]]

  plot.fun(s0, s1, ip = ips[[i]])
  mtext(paste0(LETTERS[j],", ", i), side = 3, line = 1,
        cex = 1.2, at = 2)
  j = j + 1
}
dev.off()

#- Figure A4 -------------------------------------------------------------------

dat <- c("COAD", "HNSC", "KIRC", "LUAD", "PRAD", "THCA")
tiff("../output/Figure A4.tiff", width = 8, height = 9, units = "in", res = 600,
     compression = "lzw")
# png("../output/Figure A3.png", width = 8, height = 9, units = "in", res = 600)
par(mfrow = c(3, 2), mar = c(1, 0, 3, 2), oma = c(1, 1, 1, 1))
j = 1
set.seed(9852174)
for(i in dat){

  s <- getgse(i)
  s0 <- s[[1]]
  s1 <- s[[2]]

  plot.fun(s0, s1, ip = ips[[i]])
  mtext(paste0(LETTERS[j],", ", i), side = 3, line = 1,
        cex = 1.2, at = 2)
  j = j + 1
}
dev.off()

#-------------------------------------------------------------------------------
# revison V1, sensitivity analysis on kappa selection
#-------------------------------------------------------------------------------

load("../output/drop1out.ranks.RData")
# tiff("../output/Figure V1.tiff", width = 10, height = 3 * 3, units = "in",
#      res = 600, compression = "lzw")

# png("../output/Figure4.png", width = 10, height = 3 * 3, units = "in",
#      res = 600)

# par(mar = c(4, 4, 1, 1), mfrow = c(3, 3), oma = c(1, 3, 1, 1))


# for(i in c("STAD", "LIHC", "BRCA")){
# i = "BRCA"
# obj <- drop1out.ranks[[i]]
# kappa <- obj$kappa
# n <- nrow(obj$drop1Rank)
#
# plotRankScatters(obj)
# mtext(i, side = 2, line = 4)
# plotAdaptiveWeights(kappa, n, type = "line", ylim = c(0, 0.02))
# plotIPs(obj, ylim = c(0, 12))

# }

# dev.off()

wfun <- function(x, kappa, n) {
    exp(-kappa * x) * (1 - exp(kappa))/(exp(-kappa * n) - 1)
}


tiff("../output/FigureNew5.tiff", width = 9, height = 7, units = "in", res = 600,
     compression = "lzw")
par(mfrow = c(3, 4), oma = c(1, 1, 1, 1), mar = c(4, 4, 2, 1))
dat <- names(drop1out.ranks)
use <- c("LIHC", "PRAD", "BRCA")
for(i in use){

    obj <- drop1out.ranks[[i]]
    drop1Rank <- obj$drop1Rank
    kappa = c(obj$kappa, 0.01, 0.03, 0.05)

    # get score
    for(j in kappa) {
        drop1RankWeighted <- apply(drop1Rank, 2, wfun, kappa = j, n = 200)
        origRankWeighted <- wfun(seq_len(200), kappa = j, n = 200)
        score <- apply(drop1RankWeighted, 2,
                       function(x) sqrt(sum((origRankWeighted - x)^2)))
        score <- score/sd(score)
        findIPs:::lollipop(seq_along(score), score, topn = 5, ylim = c(0, 10))
        if (j == kappa[1]) {
            mtext(bquote(kappa ~ "=" ~ .(round(j, 3)) ~ "(optimized)"),
                  side = 3, line = 0.5, adj = 0)
            mtext(i, side = 2, line = 3.5)
        }else{
            mtext(bquote(kappa ~ "=" ~ .(round(j, 3))), side = 3, line = 0.5, adj = 0)
        }
    }
}

dev.off()

# "BRCA" "COAD" "HNSC" "KIRC" "LIHC" "LUAD" "LUSC" "PRAD" "STAD" "THCA"

#-------------------------------------------------------------------------------
# Revision, using limma for feature ranking
#-------------------------------------------------------------------------------

# Limma method
s <- dir("../data")

rankFun <- function(X, y){

    design <- model.matrix( ~ as.factor(y))
    fit <- limma::lmFit(X, design)
    limma::eBayes(fit)$p.value[,2]
}

# get ranks, weighted ranks, and scores
getIPsFun <- function(filename, topN = 200){

    load(paste0("../data/", filename))
    name = gsub(".RData", "", filename)

    library(limma)
    # t.testFun <- function(x, y)
    #   t.test(x ~ y, var.equal = TRUE)$p.value

    drop1ranks <- getdrop1ranks(X, y,
                                fun = rankFun,
                                topN = topN,
                                decreasing = FALSE,
                                nCores = 25)

    weightedRanks <- sumRanks(drop1ranks$origRank,
                              drop1ranks$drop1Rank,
                              method = "adaptive",
                              topN = topN)

    return(weightedRanks)
}

# obj = getIPsFun(filename = s[7], topN = 200)

drop1out.ranks <- mapply(getIPsFun, filename = s, SIMPLIFY = F)
names(drop1out.ranks) <- substr(names(drop1out.ranks), 6, 9)

# save the output
# save(list = "drop1out.ranks", file = "../output/drop1out.ranks_limma.RData")

# output influential plot

tiff("../output/Figure A3.tiff", width = 12, height = 5, units = "in",
     res = 600, compression = "lzw")

# png("../output/Figure4.png", width = 10, height = 3 * 3, units = "in",
#      res = 600)

par(mar = c(4, 3.5, 1, 1), mfrow = c(2, 5), oma = c(1, 3, 1, 1))

for(i in names(drop1out.ranks)){

    obj <- drop1out.ranks[[i]]
    plotIPs(obj)
    mtext(i, side = 3, line = 0.2)
}

dev.off()


## kruskal test

s <- dir("../data")
# get ranks, weighted ranks, and scores
getIPsFun <- function(filename, topN = 200){

    load(paste0("../data/", filename))
    name = gsub(".RData", "", filename)

    library(limma)
    # t.testFun <- function(x, y)
    #   t.test(x ~ y, var.equal = TRUE)$p.value

    drop1ranks <- getdrop1ranks(X, y,
                                fun = "kruskal.test",
                                topN = topN,
                                decreasing = FALSE,
                                nCores = 25)

    weightedRanks <- sumRanks(drop1ranks$origRank,
                              drop1ranks$drop1Rank,
                              method = "adaptive",
                              topN = topN)

    return(weightedRanks)
}

# obj = getIPsFun(filename = s[7], topN = 200)

drop1out.ranks <- mapply(getIPsFun, filename = s, SIMPLIFY = F)
names(drop1out.ranks) <- substr(names(drop1out.ranks), 6, 9)

# save the output
# save(list = "drop1out.ranks", file = "../output/drop1out.ranks_kruskal.RData")

# output influential plot

tiff("../output/Figure A4.tiff", width = 12, height = 5, units = "in",
     res = 600, compression = "lzw")

# png("../output/Figure4.png", width = 10, height = 3 * 3, units = "in",
#      res = 600)

par(mar = c(4, 3.5, 1, 1), mfrow = c(2, 5), oma = c(1, 3, 1, 1))

for(i in names(drop1out.ranks)){

    obj <- drop1out.ranks[[i]]
    plotIPs(obj)
    mtext(i, side = 3, line = 0.2)
}

dev.off()
























