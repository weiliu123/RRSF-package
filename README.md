# RRSF-package
Integration of gene interaction information into a reweighted random survival forest approach for accurate survival prediction and survival biomarker discovery.
# Details

Package:    RRSF

Type:    Package

Title:    A reweighted random survival forest approach by integrating gene interactionn information

Version:    1.0

Date:    2018-06-18

Author;    Wei Liu

Maintainer:    Wei Liu <freelw@qq.com>

Depends:    R (>= 2.10), igraph, Matrix, survival, randomForestSRC

Description:    Integration of gene interaction information into a reweighted random survival forest approach for accurate survival prediction and survival biomarker discovery.

License:     GPL(>=2)
# Index of help topics:

CreatLadder:     Create a decreasing sequence

DRW:    Directed Random Walk

RRSF:    Reweighted random survival forest

RRSF-package:    A reweighted random survival forest approach by integrating gene interactionn information

dGMMirGraph:    The global pathway graph

getW:     Calculating the weights of genes

mRNA_matrix:     The expression data

predict.RRSF:     Make predictions from a "RRSF" object

survData:     Survival data
# Examples

data(dGMMirGraph)

data(mRNA_matrix)

data(survData)

trainSmpl.Idx <- sample(1:dim(mRNA_matrix)[1], floor(1/2*dim(mRNA_matrix)[1]))

testSmpl.Idx <- setdiff(1:dim(mRNA_matrix)[1], trainSmpl.Idx)

trainSmpl <- mRNA_matrix[trainSmpl.Idx ,]

testSmpl <- mRNA_matrix[testSmpl.Idx ,]


res <- RRSF(x=trainSmpl, y=survData[trainSmpl.Idx ,], geneSel = TRUE, DEBUG=TRUE, standardize=TRUE,globalGraph = dGMMirGraph, topo.Wt = NULL, Gamma=0.3, ntree=100)

lp <- predict.RRSF(res$RRSFmodel,testSmpl)
