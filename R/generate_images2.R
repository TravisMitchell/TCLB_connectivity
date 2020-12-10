# This script generates a set of png images of slices through a fracture

library(rfracture)

power.spectrum = function(f) 0.0001/(f^3.5)

ret = fracture_matrix(dims = c(200,200,200),power.iso=power.spectrum)

mat = ret$f1<0.05

N = dim(mat)[3]
for (i in 1:N) {
  png(sprintf("cut_%04d.png",i),width=N,height=N)
  par(mar=c(0,0,0,0))
  image(mat[,,i],asp=1,xaxt="n",xaxs="i",yaxt="n",yaxs="i",bty="n",col=c("red","white"))
  dev.off()
}
