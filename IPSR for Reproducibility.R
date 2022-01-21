
#- install package
devtools::install_github("ShuoStat/IPSR")
library(IPSR)

################################################################################
#- Using IPSR to get the leave-1-out rankings 
################################################################################

#- import data, get data name
s <- dir("./data")
# get ranks, weighted ranks, and scores
getIPSfun <- function(filename, topn = 200){
  
  load(paste0("./data/", filename))
  name = gsub(".RData", "", filename)
  X <- t(X)
  
  drop1ranks <- getdrop1ranks(X, y, fun = "t.test", topn = topn, 
                            decreasing = FALSE, ncores = 6)
  weightedranks <- rank.compare(drop1ranks$orig, drop1ranks$drop1rank,
                                dist = "L2")
  return(weightedranks)
}

drop1out.ranks <- mapply(getIPSfun, filename = s, SIMPLIFY = F)
names(drop1out.ranks) <- substr(names(drop1out.ranks), 6, 9)
save(list = "drop1out.ranks", file = "./output/drop1out.ranks.RData")

################################################################################
#- LUSC example
#- output: Table 2, 3 and Figure 2
################################################################################

load("./output/drop1out.ranks.RData")
obj <- drop1out.ranks[["LUSC"]]

#- Table 2 and 3 ---------------------------------------------------------------

#- Table 2, original ranks

orig <- obj$orig.ranks[1:20, c(1:10)]
colnames(orig) <- paste0("obs ", 1:10)
orig <- cbind(orig = 1:20, orig)
write.csv(orig, "./output/Table2.csv", row.names = F)

#- Table 3, weighted ranks 

weighted <- obj$weighted.ranks[1:20, c(1:11)]
weighted <- format(round(weighted, 4), nsmall = 4)
colnames(weighted) <- c("orig", paste0("obs ", 1:10))
write.csv(weighted, "./output/Table3.csv", row.names = F)

#- Figure 2 --------------------------------------------------------------------

plot.kappa.extension <- function(obj, points = T, dist = T, ylim = NULL){
  
  n <- nrow(obj$orig.ranks)
  x <- 1:n
  kappa = obj$kappa
  
  plot_kappa_weight(obj$kappa, n, type = "line")
  
  if (points == T) {
    or <- obj$orig.ranks  #or, original ranks
    all_r <- reshape2::melt(or)
    sel_r <- all_r[all_r[,1] != all_r[,3],]
    x_p = sel_r[,3]
    y_p = exp(-kappa * x_p) * (1 - exp(kappa)) / (exp(-kappa * n) - 1)
    points(sel_r[,1], y_p, cex = 0.7, col = "grey60", pch = 19)
  } 
  
  if (dist == T){
    top = order(obj$score, decreasing = T)[1]
    y0 = obj$weighted.ranks[,1]
    y1 = obj$weighted.ranks[,top + 1]
    
    for(i in seq_along(y0)){
      lines(c(x[i], x[i]), c(y0[i], y1[i]), lwd = 0.7)
    }
    points(x, y1, cex = 0.7, pch = 19)
  }
  
  y <- exp(-kappa * x) * (1 - exp(kappa))/(exp(-kappa * n) - 1)
  lines(x, y)
}

jpeg("./output/Figure2.jpeg", width = 12, height = 3.5, units = "in", res = 400)
par(mar = c(4, 4, 2, 1), mfrow = c(1, 3))

plot_unweighted_ranks(obj)
mtext("A", side = 3, adj = 0, line = 0.2)

plot.kappa.extension(obj)
mtext("B", side = 3, adj = 0, line = 0.2)

plot(obj,  ylim = c(0, 13))
mtext("C", side = 3, adj = 0, line = 0.2)

dev.off()


################################################################################
#- weights and distribution
#- Output Figure 3
################################################################################

load("./output/allweightedranks.RData")
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

jpeg("./output/Figure3.jpeg", width = 10, height = 3, res = 300, units = "in")
par(mfrow = c(1, 3), mar = c(4, 4, 2, 2))

plot_unweighted_ranks(drop1out.ranks[[s[1]]], top = F, 
                      points.arg = list(col = "grey"))
mtext(paste0("A,", s[1]), side = 3, adj = 0, line = 0.2)

plot_unweighted_ranks(drop1out.ranks[[s[2]]], top = F, 
                      points.arg = list(col = "black"))
mtext(paste0("B,", s[2]), side = 3, adj = 0, line = 0.2)

kappa <- c(drop1out.ranks[[s[1]]]$kappa, drop1out.ranks[[s[2]]]$kappa)
plot.kappa(kappa, col = c("grey", "black"), n = 200, s)
mtext("C", side = 3, adj = 0, line = 0.2)

dev.off()

################################################################################
#- Exhibit Influential Observations
#- Output: Figure 4 and Figure A1
################################################################################

load("./output/drop1out.ranks.RData")

#- Figure 4 --------------------------------------------------------------------

jpeg("./output/Figure4.jpeg", width = 10, height = 3 * 3, units = "in", res = 300)
par(mar = c(4, 4, 1, 1), mfrow = c(3, 3), oma = c(1, 3, 1, 1))

for(i in c("STAD", "LIHC", "BRCA")){
  
  obj <- drop1out.ranks[[i]]
  kappa <- obj$kappa
  n <- length(obj$score)
  
  plot_unweighted_ranks(obj)
  mtext(i, side = 2, line = 4)
  plot_kappa_weight(kappa, n, type = "line", ylim = c(0, 0.02))
  plot(obj, ylim = c(0, 12))
}  

dev.off()

#- Figure A1--------------------------------------------------------------------

jpeg("./output/FigureA1.jpeg", width = 10, height = 6 * 3, units = "in", res = 300)
par(mar = c(4, 4, 1, 1), mfrow = c(6, 3), oma = c(1, 3, 1, 1))

for(i in c("COAD", "HNSC", "KIRC", "LUAD", "PRAD", "THCA")){
  
  obj <- drop1out.ranks[[i]]
  kappa <- obj$kappa
  n <- length(obj$score)
  
  plot_unweighted_ranks(obj)
  mtext(i, side = 2, line = 4)
  plot_kappa_weight(kappa, n, type = "line", ylim = c(0, 0.025))
  plot(obj, ylim = c(0, 12))
}  

dev.off()

################################################################################
#- Comparison: unweighted correlations, weighted correlations, and adaptive 
#- Output: Figure 5,6 and Figure A2
################################################################################

load("./output/drop1out.ranks.RData")

dist <- function(obj, method = c("D", "W")) {
  
  orig <- obj$orig.ranks
  n = nrow(orig)
  m <- 1:n
  
  D <- function(x)
    (x - m)^2 
  
  W <- function(x) 
    (x - m)^2 * (2 * n + 2 - x - m)
  
  if (method == "D"){
    score <- apply(orig, 2, D)
    score <- colSums(score)
    score <- score / sd(score)
  }
  
  if (method == "W"){
    score <- apply(orig, 2, W)
    score <- colSums(score)
    score <- score / sd(score)
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

#- Figure 5 --------------------------------------------------------------------

jpeg(paste0("./output/Figure5.jpeg"), width = 10, height = 3 * 3, units = "in", 
     res = 300)
par(mfrow = c(3, 3), mar = c(3.5, 3.5, 1.5, 1.5), oma = c(1, 2, 1, 1))

for(i in c("STAD", "LIHC", "BRCA")){
  
  obj <- drop1out.ranks[[i]]
  lollipop(1:ncol(obj$orig.ranks), dist(obj, "D"), ylim = 13)
  mtext("Unweighted", side = 3, adj = 0, line = 0.2)
  mtext(nam, side = 2, line = 4)
  
  lollipop(1:ncol(obj$orig.ranks), dist(obj, "W"), ylim = 13)
  mtext("Weighted", side = 3, adj = 0, line = 0.2)
  
  plot(obj, ylim = c(0, 13))
  mtext("Adaptive Weights", side = 3, adj = 0, line = 0.2)
}

dev.off()

#- Figure A2 -------------------------------------------------------------------

jpeg(paste0("./output/FigureA2.jpeg"), width = 10, height = 6 * 3, units = "in", 
     res = 300)
par(mfrow = c(6, 3), mar = c(3.5, 3.5, 1.5, 1.5), oma = c(1, 2, 1, 1))

for(i in c("COAD", "HNSC", "KIRC", "LUAD", "PRAD", "THCA")){
  
  obj <- drop1out.ranks[[i]]
  
  lollipop(1:ncol(obj$orig.ranks), dist(obj, "D"), ylim = 13)
  mtext("Unweighted", side = 3, adj = 0, line = 0.2)
  mtext(nam, side = 2, line = 4)
  
  lollipop(1:ncol(obj$orig.ranks), dist(obj, "W"), ylim = 13)
  mtext("Weighted", side = 3, adj = 0, line = 0.2)
  
  plot(obj, ylim = c(0, 13))
  mtext("Adaptive Weights", side = 3, adj = 0, line = 0.2)
  
}
dev.off() 

#- Further check on BRCA -------------------------------------------------------
#- Figure A6 -------------------------------------------------------------------

top.plot <- function(obj, obs, topn = 10){
  
  orirank <- obj$orig.ranks
  meltr <- reshape2::melt(orirank)
  plot(meltr[,1], meltr[,3], cex = 0.8, xlab = "", ylab = "", 
       col = "grey", pch = 19)
  
  n <- nrow(obj$orig.ranks)
  points(1:n, obj$orig.ranks[,obs], cex = 0.7, pch = 19)
  
  top <- order(abs((1:n) - obj$orig.ranks[,obs]), decreasing = T)[1:topn]
  points(top, obj$orig.ranks[top, obs], pch = 15, cex = 1.2)
  abline(v = c(50, 100, 150), lty = "dashed")
  
  legend("topleft", pch = c(19, 15, 19), col = c("black", "black", "grey"), 
         legend = c(paste0("obs", obs), paste0("top10 of obs", obs),"The others"))
  
  mtext("Original ranking", side = 1, line = 2.2)
  mtext("Leave-one-out ranking", side = 2, line = 2.2)
}

obj <- drop1out.ranks$BRCA

jpeg("./output/Figure6.jpeg", width = 10, height = 3, units = "in",
     res = 300)
par(mfrow = c(1, 3), mar = c(4, 4, 1.5, 1.5))

top.plot(obj, 175, topn = 10)
top.plot(obj, 53, topn = 10)
top.plot(obj, 338, topn = 10)

dev.off()

#-------------------------------------------------------------------------------












