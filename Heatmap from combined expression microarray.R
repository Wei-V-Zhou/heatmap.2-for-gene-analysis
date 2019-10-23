##################
# Load libraries #
##################
# Clear and garbage
rm(list = ls())
gc()
options(stringsAsFactors = FALSE)
# Load packages
{
  library(openxlsx)
  library(gdata)
  library(gplots)
  library(scater)
  library(export)
  library(ggplot2)
  library(pheatmap)
  library(geneplotter)
}

# Load functions
cor.dist.spearman <- function (x) 
{
  as.dist(1 - cor(t(x), method="spearman"))
}
hclust.ward.d <- function(x){hclust(x, method="complete")}

# Load variation-stablized RNA-seq data
if(!file.exists("combined expression matrix.Rdata")){
  expr1 <- read.xlsx("41467_2015_BFncomms7351_MOESM956_ESM.xlsx", sheet = 1)
  row.names(expr1) <- expr1[ , 1]
  expr1 <- expr1[ , -1]
  for (i in 1 : 5){
    for (i in 1 : nrow(expr1)){
      if (is.na(expr1[i, 1])){
        expr1 <- expr1[-i, ]
      }
    }
    print(dim(expr1))
  }
  
  expr2 <- read.xlsx("1-s2.0-S2211124716314784-mmc3.xlsx", sheet = 1)
  expr2 <- expr2[ , c(-6, -7, -8, -9, -10, -11)]
  row.names(expr2) <- expr2[ , 1]
  expr2 <- expr2[ , -1]
  expr2 <- expr2[rownames(expr1), ]
  for (i in 1: 5){
    for (i in 1 : nrow(expr2)){
      if (is.na(expr2[i, 1])){
        expr2 <- expr2[-i, ]
      }
    }
    print(dim(expr2))
  }

  expr1 <- expr1[rownames(expr2), ]
  exprs <- cbind(expr1, expr2)
  expr <- exprs[!duplicated(exprs[ , 1]), ]
  rownames(expr) <- expr[ , 1]
  expr <- expr[ , -1]
  save(expr, file = "combined expression matrix.Rdata")
}
load("combined expression matrix.Rdata")


markergene <- c("NOTCH1", "NOTCH2", "NOTCH3", "NOTCH4", "JAG1", "JAG2", "DLL1", "DLL3", "DLL4")
exprs <- expr[markergene, ]
n <- t(scale(t(exprs)))
n[n > 2] <- 2
n[n < -2] <- -2

heatmap.2(t(n), dendrogram = "row", trace = "none", scale = "column", col=bluered(100), 
          distfun = cor.dist.spearman, hclust=hclust.ward.d, 
          srtCol = 45, Colv=FALSE, margins=c(8,10))


#============================#
#       Musician: Resonance  #
#           Date: 2019/10/23 #
# Revised author: Resonance  #
#           Time: 2019/10/23 #
#============================#