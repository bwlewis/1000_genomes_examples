# Example parallel irlba unit test
# Input: none
# Output: none, but exits with error if correctness check fails
# invoke with mpirun -np 2 R --slave -e 'source("unit.R")'
# adjusting np as required (e.g., 2, 4, 8, 32, 64, ...)
#
# this test scales soft; the problem size increases with the number of nodes

suppressMessages(library(irlba))
suppressMessages(library(doMPI))
suppressMessages(library(methods))
cl = startMPIcluster()
registerDoMPI(cl)
N = clusterSize(cl) # number of nodes
X = 100
y = matrix(runif(N * X * 100), nrow=N * X)
u = qr.Q(qr(y))
v = qr.Q(qr(t(y)))
d = exp(-seq(1, min(100, N * X)))
y = u %*% (d * t(v))  # an (N * X by 100) matrix with known SVD

setClass("fake", contains="list", S3methods=TRUE, slots=c(dims="numeric", M="matrix"))
setMethod("%*%", signature(x="fake", y="numeric"), function(x ,y) 
  {
    foreach(i=x, .combine=c) %dopar%
    {
      x@M[i, ] %*% y
    }
  })
setMethod("%*%", signature(x="numeric", y="fake"), function(x ,y) 
  {
    foreach(i=y, .combine=`+`) %dopar%
    {
      drop(x[i] %*% y@M[i,])
    }
  })
A = new("fake", split(1:nrow(y), rep(1:N, each=X)), dims=c(nrow(y), ncol(y)), M=y)
dim.fake = function(x) x@dims
nrow.fake = function(x) x@dims[1]
ncol.fake = function(x) x@dims[2]

s = irlba(A, nv=2, fastpath=FALSE)
message("matrix products: ", s$mprod)
stopifnot(all.equal(d[1:2], s$d))

closeCluster(cl)
mpi.quit()
