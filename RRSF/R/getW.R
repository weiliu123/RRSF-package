getW <-
function(geneWeight, globalGraph)
{
Vertexs <- V(globalGraph)
W <- rep(0, length(Vertexs))
names(W) <- Vertexs$name
for(i in 1 : length(W)){
idx <- which(names(geneWeight) == names(W[i]))
if(length(idx) > 0){
W[i] <- geneWeight[idx]
}
}
W <- W/sum(W)
}
