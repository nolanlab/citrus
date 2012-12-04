

pamr.plotfdr.new = function(fdrfit, threshold, call.win.metafile = FALSE,errorRate=NULL,thresholdHighlight=NULL){
  om = fdrfit$results[, "Number of significant genes"] > 0
  na.min = function(x) {
    min(x[!is.na(x)])
  }
  na.max = function(x) {
    max(x[!is.na(x)])
  }
  par(oma=c(0,0,0,0),mar=c(5,5,5,5))
  plot(fdrfit$results[,1], fdrfit$results[, "Median FDR"], 
       pch=20,
       xlab = "Regularization Threshold", 
       ylab = "Cross validation error rate", 
       type = "o", ylim = c(0, 0.5),xlim=c(rev(range(fdrfit$results[,1]))))
  #lines(fdrfit$results[,1], fdrfit$results[, "90th percentile of FDR"],col="grey",lty=2) 
  #points(fdrfit$results[,1], fdrfit$results[, "90th percentile of FDR"],col="grey") 
  #nog = which(fdrfit$results[,1]==threshold)
  #points(x=fdrfit$results[nog,2],y=fdrfit$results[nog,4],col="red",pch=20,cex=3)
  #points(x=fdrfit$results[nog,2],y=fdrfit$results[nog,4],col="black",pch=20,cex=1)
  if (!is.null(errorRate)){
    lines(x=fdrfit$results[,1],y=errorRate,type='o',col="blue",pch=20)
  }
  x = fdrfit$results[, "Threshold"]
  upper = fdrfit$results[, "90th percentile of FDR"]
  lower = fdrfit$results[, "Median FDR"]
  segments(x, upper, x, lower, lty = 2)
  segments(x - .1, upper, x + .1, upper, lty = 2)
  axis(3, at = fdrfit$results[,1], 
       labels = round(fdrfit$results[,2], 2))
  mtext("Number of Model Features", 3, 3, cex = 1)
  axis(4, at = seq(from=0,to=0.5,by=0.1))
  mtext("Estimated FDR: Median and 90th percentile", 4, 3, cex = 1)
  if (!is.null(threshold)){
    lines(x=c(thresholdHighlight,thresholdHighlight),y=c(0,1),lty=2,lwd=2,col="red")
  }
  #mtext("Number of Model Features", 3, 2, cex = 1)
  if (call.win.metafile) {
    dev.off()
  }
  legend(x=(max(fdrfit$results[,1])),y=0.5,col=c("black","blue"),legend=c("Median Estimated Feature FDR Rate","Cross Validation Error Rate"),fill=c("black","blue"),cex=.75,box.col=rgb(0,0,0,0))
  return()
}

pamr.fdr.new = function (trained.obj, data, nperms = 100, xl.mode = c("regular", 
                                                                      "firsttime", "onetime", "lasttime"), xl.time = NULL, xl.prevfit = NULL) 
{
  this.call <- match.call()
  xl.mode = match.arg(xl.mode)
  if (xl.mode == "regular" | xl.mode == "firsttime") {
    y = data$y
    m = nrow(data$x)
    nclass = length(table(y))
    threshold <- trained.obj$threshold
    n.threshold = length(threshold)
    tt <- scale((trained.obj$centroids - trained.obj$centroid.overall)/trained.obj$sd, 
                FALSE, trained.obj$threshold.scale * trained.obj$se.scale)
    ttstar <- array(NA, c(m, nperms, nclass))
    results = NULL
    pi0 = NULL
  }
  if (xl.mode == "onetime" | xl.mode == "lasttime") {
    y = xl.prevfit$y
    m = xl.prevfit$m
    nclass = xl.prevfit$nclass
    threshold = xl.prevfit$threshold
    n.threshold = xl.prevfit$n.threshold
    tt = xl.prevfit$tt
    ttstar = xl.prevfit$ttstar
    nperms = xl.prevfit$nperms
    results = xl.prevfit$results
    pi0 = xl.prevfit$pi0
  }
  if (xl.mode == "regular") {
    first = 1
    last = nperms
  }
  if (xl.mode == "firsttime") {
    first = 1
    last = 1
  }
  if (xl.mode == "onetime") {
    first = xl.time
    last = xl.time
  }
  if (xl.mode == "lasttime") {
    first = nperms
    last = nperms
  }
  for (i in first:last) {
    cat("", fill = T)
    cat(c("perm=", i), fill = T)
    ystar <- sample(y)
    data2 <- data
    data2$y <- ystar
    foo <- pamr.train(data2, threshold = 0, scale.sd = trained.obj$scale.sd, remove.zeros=F,
                      threshold.scale = trained.obj$threshold.scale, se.scale = trained.obj$se.scale, 
                      offset.percent = 50, hetero = trained.obj$hetero, 
                      prior = trained.obj$prior, sign.contrast = trained.obj$sign.contrast)
    sdstar = foo$sd - foo$offset + trained.obj$offset
    ttstar[, i, ] = scale((foo$centroids - foo$centroid.overall)/sdstar,FALSE, foo$threshold.scale * foo$se.scale)
  }
  if (xl.mode == "regular" | xl.mode == "lasttime") {
    fdr = rep(NA, n.threshold)
    fdr90 = rep(NA, n.threshold)
    ngenes = rep(NA, n.threshold)
    for (j in 1:n.threshold) {
      nobs = sum(drop((abs(tt) - threshold[j] > 0) %*% rep(1, ncol(tt))) > 0)
      temp = abs(ttstar) - threshold[j] > 0
      temp2 = rowSums(temp, dim = 2)
      nnull = colSums(temp2 > 0)
      fdr[j] = median(nnull)/nobs
      fdr90[j] = quantile(nnull, 0.9)/nobs
      ngenes[j] = nobs
    }
    q1 <- quantile(ttstar, 0.25)
    q2 <- quantile(ttstar, 0.75)
    pi0 <- min(sum(tt > q1 & tt < q2)/(0.5 * m * nclass),1)
    fdr <- fdr * pi0
    fdr90 = fdr90 * pi0
    fdr = pmin(fdr, 1)
    fdr90 = pmin(fdr90, 1)
    results <- cbind(threshold, ngenes, fdr * ngenes, fdr, fdr90)
    om = is.na(fdr)
    results[om, 3:5]=0
    dimnames(results) <- list(NULL, c("Threshold", "Number of significant genes", 
                                      "Median number of null genes", "Median FDR", "90th percentile of FDR"))
    y = NULL
    x = NULL
    m = NULL
    threshold = NULL
    n.threshold = NULL
    tt = NULL
    nperms = NULL
    ttstar = NULL
  }
  return(list(results = results, pi0 = pi0, y = y, m = m, threshold = threshold, 
              n.threshold = n.threshold, tt = tt, ttstar = ttstar, 
              nperms = nperms))
}


