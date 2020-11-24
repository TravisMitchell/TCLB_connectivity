# This scipt generates the connectivity from a set of png images

library(png)
source("R/lib.R")

pngs = sprintf("cut_%04d.png",1:512)
if (! all(file.exists(pngs))) stop("Some of the png's don't exist")

U = as.matrix(expand.grid(x=-1:1,y=-1:1,z=-1:1))

N = length(pngs)

get.slice = function(i) {
  i = (i-1) %% N + 1
  readPNG(pngs[i])[,,2] < 0.5
}

B = get.slice(1)
# image(B,asp=1)

X = row(B)-1
Y = col(B)-1

i = 2

seq_mod = function(n,d) (seq_len(n) + d - 1) %% n + 1
get.mod = function(M, dx=0, dy=0, dz=0) {
  if (length(dim(M)) == 2) {
    M[seq_mod(nrow(M),dx), seq_mod(ncol(M),dy)]
  } else if (length(dim(M)) == 3) {
    M[seq_mod(dim(M)[1],dx), seq_mod(dim(M)[2],dx), seq_mod(dim(M)[3],dz)]
  } else stop("unknown array in get.mod")
}

get.slice.indexes = function(i) {
  A = get.slice(i-1)
  B = get.slice(i+0)
  C = get.slice(i+1)
  BB = A | B | C
  BB = get.mod(BB,-1,0) | get.mod(BB,0,0) | get.mod(BB,1,0)
  BB = get.mod(BB,0,-1) | get.mod(BB,0,0) | get.mod(BB,0,1)
  I = matrix(NA,nrow(B),ncol(B))
  I[BB] = seq_len(sum(BB))-1
  I
}

offsetB = 0

ret = lapply(seq_mod(N,+1), function(i) {
  print(i)
  A = get.slice.indexes(i-1)
  B = get.slice.indexes(i+0)
  C = get.slice.indexes(i+1)
  offsetA = offsetB
  offsetB <<- offsetA + max(A,na.rm=TRUE)
  if (i == 1) offsetB = 0
  offsetC = offsetB + max(B,na.rm=TRUE)
  if (i == N) offsetC = 0
  BIG = c(A+offsetA, B+offsetB, C+offsetC)
  dim(BIG) = c(nrow(B),ncol(B),3)
  
  sel = !is.na(B)
  I1 = B[sel] + offsetB
  TP = get.slice(i)[sel]
  ret = lapply(1:nrow(U), function(k) {
    get.mod(BIG[,,2+U[k,3]], U[k,1], U[k,2])[sel]
  })
  names(ret) = paste0("c",1:nrow(U))
  I2 = do.call(cbind,ret)
  tab = cbind(data.frame(x=X[sel],y=Y[sel],z=i,idx=I1,type=TP),I2)
  tab
})

ret = ret[seq_mod(length(ret),-1)] # first slice was the last
tab = do.call(rbind,ret) #glue everything together
rm(ret)
gc()
range(tab$idx)
range(tab[,paste0("c",1:nrow(U))],na.rm=TRUE)

names(tab)

tab$labels = ifelse(tab$type, "Interior","Wall")
for (i in paste0("c",1:nrow(U))) {
  sel = is.na(tab[,i])
  tab[sel,i] = tab$idx[sel]
}

write.connectivity(tab, U, "frac.cxn.gz")
