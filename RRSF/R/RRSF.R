RRSF <-
function(x, y, geneSel=FALSE, DEBUG=FALSE, standardize=TRUE, globalGraph = NULL, topo.Wt = NULL, Gamma=0.3, ...){
if(standardize){
x <- scale(x, center = TRUE, scale = TRUE)
}
commonmRNA <- intersect(colnames(x),V(globalGraph)$name)
x <- x[, commonmRNA]

if(is.null(topo.Wt)){
# calculate the p-value of each gene using Cox PH model
if(DEBUG) cat('Calculating Cox score...')

colnames(y) <- c("status", "time")
Survdata.mRNA <- data.frame(x, y,check.names = TRUE)
geneCoxZP <- matrix(NA,nrow=ncol(x),ncol=2)
rownames(geneCoxZP) <- colnames(x)
colnames(geneCoxZP) <- c("ZScore","p-value")
for(i in 1 : ncol(x)){
res.coxph <- coxph(as.formula(paste("Surv(time, status)~", colnames(Survdata.mRNA)[i])), Survdata.mRNA)
geneCoxZP[i,] <- summary(res.coxph)$coefficients[c(4,5)]
}
geneWeight <- -log(geneCoxZP[ ,2]+2.2e-16)
geneWeight[which(is.na(geneWeight))] <- 0
geneWeight <- (geneWeight - min(geneWeight)) / (max(geneWeight) - min(geneWeight))

W <- getW(geneWeight=geneWeight, globalGraph=globalGraph)

if(DEBUG) cat('Done\n')

if(DEBUG) cat('Performing directed random walk...')
vertexWeight <- DRW(igraphM = globalGraph, p0 = W, gamma = Gamma)
if(DEBUG) cat('Done\n')

topo.Wt <- vertexWeight[commonmRNA]

}
# browser()
topo.Wt.temp <- topo.Wt
names(topo.Wt.temp) <- colnames(Survdata.mRNA)[1:(ncol(Survdata.mRNA)-2)]

res <- rfsrc(Surv(time, status) ~ ., Survdata.mRNA, xvar.wt = topo.Wt.temp, ...)
class(res) <- c("rfsrc","grow","surv")

# gene selection
if(geneSel){
ladder <- CreatLadder(Ntotal = ncol(x), pRatio = 0.9, Nmin = 2)
geneSelected <- matrix(0, nrow = ncol(x), ncol = length(ladder)-1)  # 记录基因被选择次数
rownames(geneSelected) <- colnames(Survdata.mRNA)[1:(ncol(Survdata.mRNA)-2)]
colnames(geneSelected) <- ladder[-1]

genes <- colnames(Survdata.mRNA)[1:(ncol(Survdata.mRNA)-2)]
for(i in 1:length(ladder)){
x.temp <- Survdata.mRNA[,c(genes,"status", "time")]
x.grow <- rfsrc(Surv(time, status) ~ ., x.temp, xvar.wt = topo.Wt.temp[genes], ...)
vimp.i <- sort(abs(vimp(x.grow)$importance),decreasing = TRUE)   # 变量重要性排序

if(i < length(ladder)){
genes <- names(vimp.i)[1:ladder[i+1]]
geneSelected[genes, i] <- geneSelected[genes, i] + 1
}

}
rownames(geneSelected) <- colnames(x)

}else{
geneSelected = NULL
}

return(list(RRSFmodel = res, topo.Wt = topo.Wt, geneSelected = geneSelected))


}
