# VCF Principal Components Example
# Assumes that parse.c (https://raw.githubusercontent.com/bwlewis/1000_genomes_examples/master/parse.c)
# is compiled and copied to executable parsevcf program in PATH. This implementation
# assumes MPI and a global shared file system.
#
# Input: One or more variant files in *.vcf.gz
# Optional input:
#   CHUNKSIZE environment variable (number of VCF file rows per chunk)
#   NCOMP environment variable (number of components, defaults to 3)
# Output: pca.rdata R data file
# Example invocation
# CHUNKSIZE=10000000 mpirun -np 4 Rscript --slave pca.R

suppressMessages(library(Matrix))
suppressMessages(library(doMPI))
suppressMessages(library(methods))
suppressMessages(library(irlba))
cl = startMPIcluster()
registerDoMPI(cl)

# 1. VCF parsing
# Parse the VCF input data files in *.vcf.gz into partitioned sparse
# R submatrices with at most CHUNKSIZE nonzer elements per matrix, saving
# the partitions to R data files for use by the PCA computation following
# this section.

t0 = proc.time()
# Establish an approximate chunk size based on 1000 Genomes project estimates
chunksize = as.numeric(Sys.getenv("CHUNKSIZE"))
if(is.na(chunksize))
{
   chunksize = floor(100000000 * as.numeric(system("free -b | sed -n 2p | tr -s ' ' | cut -f 2 -d ' '", intern=TRUE)) / 16e9)
}
if(is.na(chunksize)) stop("error setting chunksize")
message("chunksize: ", chunksize)

meta = foreach(f=dir(pattern="*\\.vcf\\.gz"), .combine=rbind, .packages=c("methods", "Matrix")) %dopar%
{
  name = gsub("\\.gz", "", f); message(name)
  chunk = 1
  p = pipe(sprintf("zcat %s  | sed /^#/d | cut  -f '10-' | parsevcf | cut -f '1-2'", f), open="r")
  meta = data.frame()
  while(chunk > 0)
  {
    x = tryCatch(read.table(p, colClasses=c("integer", "integer"), fill=TRUE, row.names=NULL, nrows=chunksize),
                 error=function(e) data.frame())
    if(nrow(x) < 1) chunk = 0
    else
    {
      x = sparseMatrix(i=x[, 1] - x[1, 1] + 1, j=x[, 2], x=1.0)
      attr(x, "rowmeans") = rowMeans(x)
      cfn = sprintf("%s-%d.pmat.rdata", name, chunk)
      cf = file(cfn, open="wb")
      serialize(x, connection=cf, xdr=FALSE)
      close(cf)
      meta = rbind(meta, data.frame(file=cfn, nrow=nrow(x), ncol=ncol(x), stringsAsFactors=FALSE))
      chunk = chunk + 1
    }
    rm(x)
    gc()
  }
  close(p)
  meta
}
message("parsing time: ", (proc.time() - t0)[[3]])
meta$end = cumsum(meta$nrow)
meta$start = c(1, meta$end[-length(meta$end)] + 1)
meta$file = sprintf("%s/%s", getwd(), meta$file)



# 2. Principal components computation

ncomp = as.numeric(Sys.getenv("NCOMP"))
if(is.na(ncomp)) ncomp = 3

setClass("pmat", contains="list", S3methods=TRUE, slots=c(dims="numeric"))
setMethod("%*%", signature(x="pmat", y="numeric"), function(x ,y)
  {
    ans = rep(0.0, nrow(x))
    p = foreach(i=1:length(x$file), .packages=c("methods", "Matrix")) %dopar%
    {
      f = file(x$file[i], open="rb")
      a = unserialize(f)
      close(f)
      r = attr(a, "rowmeans")
      drop(a %*% y - r * drop(crossprod(rep(1, length(y)), y)))
    }
    i = 1
    for(j in 1:length(p))
    {
      k = length(p[[j]])
      ans[i:(i + k - 1)] = p[[j]]
      i = i + k
    }
    gc()
    ans
  })
setMethod("%*%", signature(x="numeric", y="pmat"), function(x ,y)
  {
    ans = foreach(i=1:length(y$file), .packages=c("methods", "Matrix"), .combine="+") %dopar%
    {
      f = file(y$file[i], open="rb")
      a = unserialize(f)
      close(f)
      j = seq(from=y$start[i], to=y$end[i])
      drop(x[j] %*% a - drop(crossprod(x[j], attr(a, "rowmeans"))))
    }
    gc()
    ans
  })
A = new("pmat", as.list(meta), dims=c(tail(meta$end, 1), meta$ncol[1]))
dim.pmat = function(x) x@dims
nrow.pmat = function(x) x@dims[1]
ncol.pmat = function(x) x@dims[2]

t1 = proc.time()
L  = irlba(A, nv=ncomp, tol=1e-5, right_only=TRUE, work=4)
dt = proc.time() - t1
saveRDS(L, file="pca.rdata")
message("PCA time: ", dt[[3]])

# Remove scratch files
system("rm -f *.pmat.rdata")

# Shutdown the cluster and quit
closeCluster(cl)
mpi.quit()
