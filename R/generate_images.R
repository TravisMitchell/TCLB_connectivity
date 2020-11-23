# This script generates a set of png images of slices through a fracture

library(rfracture)

power.spectrum = function(f) 0.0001/(f^3.5)

refine = 30
  ret = fracture_geom(
    width = 1,
    refine = refine,
    corr.profile = function(lambda) 0.8,
    gap=0.03,
    power.iso=power.spectrum,
    seed=123
  )
  ret = slice(ret,value = "above")

N = 512
for (i in 1:N) {
  ed = slice(ret,by="x", value = "edge", level=i/N)
  png(sprintf("cut_%04d.png",i),width=N,height=N)
  par(mar=c(0,0,0,0))
  plot(NA, xlim=c(0,1),ylim=c(-0.5,0.5), asp=1,xaxt="n",xaxs="i",yaxt="n",yaxs="i",bty="n")
  quads = function(a,b,c,d,e,f,...) { apply(cbind(a,b,c,d,e,f), 1, function(x,...) polygon(c(x[1],x[2],x[2],x[1]),c(x[3],x[4],x[6],x[5]),...), ...) }
  quads(ed$points$y[ed$edge[,1]],ed$points$y[ed$edge[,2]],ed$points$f1[ed$edge[,1]],ed$points$f1[ed$edge[,2]],ed$points$f2[ed$edge[,1]],ed$points$f2[ed$edge[,2]],col=2,border=NA)
  dev.off()
}
