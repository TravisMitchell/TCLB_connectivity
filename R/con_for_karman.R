## This script generates the connectivity for the karman_arb.xml example
source("R/lib.R")

d2q9 = matrix(c(
  0,  1,  0, -1,  0,  1, -1, -1,  1,
  0,  0,  1,  0, -1,  1,  1, -1, -1,
  0,  0,  0,  0,  0,  0,  0,  0,  0), 9, 3);

U = d2q9

A = matrix(0,1024,100)
A[] = 1
A[1,] = 3
A[1024,] = 4
A[,1] = 2
A[,100] = 2
image(1:nrow(A),1:ncol(A),A,asp=1)
X = row(A)-0.5
Y = col(A)-0.5

A[abs(Y-50)+abs(X-140)<20] = 2

I = 1:length(A) -1
dim(I) = dim(A)
B = A
ri = 1:nrow(B)
rp = c(2:nrow(B),1)
rm = c(nrow(B),1:(nrow(B)-1))
ci = 1:ncol(B)
cp = c(2:ncol(B),1)
cm = c(ncol(B),1:(ncol(B)-1))
range(A)

tab = data.frame(
  x=as.vector(X),
  y=as.vector(Y),
  z=0.5,
  c1=as.vector(I[,]),
  c2=as.vector(I[rp,]),
  c3=as.vector(I[,cp]),
  c4=as.vector(I[rm,]),
  c5=as.vector(I[,cm]),
  c6=as.vector(I[rp,cp]),
  c7=as.vector(I[rm,cp]),
  c8=as.vector(I[rm,cm]),
  c9=as.vector(I[rp,cm]),
  labels=as.factor(c("Main", "Wall", "Inlet", "Outlet"))[as.vector(A)]
)

write.connectivity(tab,U,"karman.cxn")
