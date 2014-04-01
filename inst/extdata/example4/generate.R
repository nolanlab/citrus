library("flowCore")
events = 10000
dims = 2
p1c1 = .90
p2c1 = .09
p3c1 = .01
p1c2 = .89
p2c2 = .03
p3c2 = .08
p1c3 = .86
p2c3 = .08
p3c3 = .08



for (i in 1:10){
  err=runif(n=1,min=-0.01,max=0.01)
  
  p1n = floor(events*(p1c1+err/2))
  p2n = floor(events*(p2c1+err/2))
  p3n = floor(events*(p3c1-err))
  
  p1data = matrix(rnorm(p1n*dims,sd=.85),ncol=dims)
  p2data = matrix(rnorm(p2n*dims,mean=2.75,sd=.5),ncol=dims)
  p3data = cbind(rnorm(p3n,mean=3.2,sd=.5),rnorm(p3n,sd=.5))
  
  d = rbind(p1data,p2data,p3data)
  colnames(d) = c("LineageMarker1","LineageMarker2")
  f = flowFrame(exprs=d)  
  write.FCS(f,filename=paste("~/Desktop/work/citrus/inst/extdata/example4/Patient",sprintf("%02.0f",i),"_Class1.fcs",sep=""))
}

for (i in 11:20){
  err=runif(n=1,min=-0.01,max=0.01)
  
  p1n = floor(events*(p1c2+err/2))
  p2n = floor(events*(p2c2+err/2))
  p3n = floor(events*(p3c2-err))
  
  p1data = matrix(rnorm(p1n*dims,sd=.85),ncol=dims)
  p2data = matrix(rnorm(p2n*dims,mean=2.75,sd=.5),ncol=dims)
  p3data = cbind(rnorm(p3n,mean=3.2,sd=.5),rnorm(p3n,sd=.5))
  
  d = rbind(p1data,p2data,p3data)
  colnames(d) = c("LineageMarker1","LineageMarker2")
  f = flowFrame(exprs=d)  
  write.FCS(f,filename=paste("~/Desktop/work/citrus/inst/extdata/example4/Patient",sprintf("%02.0f",i),"_Class2.fcs",sep=""))
}

for (i in 21:30){
  err=runif(n=1,min=-0.01,max=0.01)
  
  p1n = floor(events*(p1c3+err/2))
  p2n = floor(events*(p2c3+err/2))
  p3n = floor(events*(p3c3-err))
  
  p1data = matrix(rnorm(p1n*dims,sd=.85),ncol=dims)
  p2data = matrix(rnorm(p2n*dims,mean=2.75,sd=.5),ncol=dims)
  p3data = cbind(rnorm(p3n,mean=3.2,sd=.5),rnorm(p3n,sd=.5))
  
  d = rbind(p1data,p2data,p3data)
  colnames(d) = c("LineageMarker1","LineageMarker2")
  f = flowFrame(exprs=d)  
  write.FCS(f,filename=paste("~/Desktop/work/citrus/inst/extdata/example4/Patient",sprintf("%02.0f",i),"_Class3.fcs",sep=""))
}
