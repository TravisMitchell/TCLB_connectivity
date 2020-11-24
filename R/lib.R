
write.connectivity = function(x, U, filename) {
  if (! "data.frame" %in% class(x)) stop("x has to be a data.frame")
  if (length(dim(U)) != 2) stop("U has to be an array")
  if (ncol(U) != 3) stop("U has to have 3 columns (even in 2D)")
  Q = nrow(U)
  C = paste0("c",1:Q)
  need = c("x","y","z",C,"labels")
  if (! all(need %in% names(x))) stop("Missing columnts in tab:", setdiff(need,names(x)))
  
  if (is.character(x$labels)) {
    nlab = 1
    lab = x$labels
  } else if (is.factor(x$labels)) {
    nlab = 1
    lab = as.character(x$labels)
  } else if (is.list(x$labels)) {
    nlab = sapply(x$labels, length)
    lab = sapply(x$labels, paste)
  } else {
    stop("")
  }
  if (any(is.na(x[,C]))) stop("NA in connectivity")
  if (! all(x[,C] < nrow(x))) stop("Elements of connectivity out of range")
  if (! all(x[,C] >= 0)) stop("Elements of connectivity out of range")
  tab = cbind(1:nrow(x)-1, x$x, x$y, x$z, x[,C],nlab,lab)
  if (grepl("[.]gz$",filename)) {
    f = gzfile(filename,open="w")
  } else {
    f = file(filename,open="w")
  }
  cat("LATTICESIZE ",nrow(tab),"\n",file=f,sep="")
  cat("BASE_LATTICE_DIM ",paste(apply( x[,c("x","y","z")], 2, function(x)max(x)-min(x)+1),collapse=" "),"\n",file=f,sep="")
  cat("d ", sum(apply(U,2,function(x) any(x!=0))),"\n",file=f,sep="")
  cat("Q ", nrow(U),"\n",file=f,sep="")
  cat("OFFSET_DIRECTIONS\n",file=f,sep="")
  cat(paste0("[",U[,1],",",U[,2],",",U[,3],"]",collapse=","),"\n",file=f,sep="")
  cat("NODES\n",file=f,sep="")
  options(scipen=10)
  write.table(tab, file=f, row.names = FALSE, col.names=FALSE,quote = FALSE)
  options(scipen=0)
  close(f)
}

