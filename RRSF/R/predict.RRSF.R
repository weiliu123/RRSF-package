predict.RRSF <-
function(object, newx, ...){
newx <- data.frame(newx)
predict.rfsrc(object, newx, ...)
}
