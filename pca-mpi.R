# VCF Principal Components Example
# Assumes that parse.c
# (https://raw.githubusercontent.com/bwlewis/1000_genomes_examples/master/parse.c)
# is compiled and copied to executable parsevcf program in PATH on every node.
# For instance:
#
# cc -O2 -o parsevcf parse.c
# sudo mv parsevcf /usr/local/bin
# (alternatively, don't move the parsevcf program but invoke R with PATH=$(pwd):$PATH)
#
# This implementation assumes MPI and that the input *.vcf.gz files are split
# up among the workers. Each worker will process only the *.vcf.gz files
# present in its local working directory. It does not assume a shared file
# system.
#
# NOTE: This version uses MPI to coordinate cross-node work in master/slave
# fashion, but simple forked parallelism within nodes. The number of MPI nodes
# should be set to the number of nodes + 1 (not generally the number of CPU
# cores). The extra MPI node (the + 1) is used by Rmpi as the master
# coordinating process.
#
# Input: One or more variant files in *.vcf.gz
# Optional input:
#   CHUNKSIZE environment variable (number of VCF file rows per chunk)
#   NCOMP environment variable (number of components, defaults to 3)
#   NP environment variable number of R worker CPUs to use per MPI host, default parallel::detectCores()
#   OMP_NUM_THREADS environment variable number of BLAS threads per worker to use, you probably want this to be 1
# Output: pca.rdata R data file
#
# Example invocation across four computers with 8 cores per computer
# NP=8 CHUNKSIZE=10000000 mpirun -np 4 -x NP -x CHUNKSIZE Rscript --slave pca-mpi.R

suppressMessages(
{
  library(Matrix)
  library(doMPI)
  library(methods)
  library(irlba)
  library(parallel)
})
cl = startMPIcluster()
N  = clusterSize(cl)
registerDoMPI(cl)

# Set the environment variable SKIP_PARSE=TRUE to skip parsing step.
SKIP_PARSE = (Sys.getenv("SKIP_PARSE") == "TRUE")

# If your system has hyperthreading/SMT enabled, consider setting NP to the
# number of physical CPU cores per computer.
NP = as.integer(Sys.getenv("NP"))
if(is.na(NP)) NP = detectCores()


if(SKIP_PARSE)
{
  load("meta.rdata")
} else
{
  # 1. VCF parsing
  # Parse the VCF input data files in *.vcf.gz into partitioned sparse
  # R submatrices with at most CHUNKSIZE nonzer elements per matrix, saving
  # the partitions to R data files for use by the PCA computation following
  # this section.

  t0 = proc.time()
  chunksize = as.numeric(Sys.getenv("CHUNKSIZE"))
  if(is.na(chunksize)) chunksize = 1e7  # adjust as needed to fit your memory constraints
  if(is.na(chunksize)) stop("error setting chunksize")
  message("chunksize: ", chunksize)

  meta = foreach(i=1:N, .combine=rbind, .packages=c("methods", "Matrix", "parallel")) %dopar%
  {
    Reduce(rbind, mclapply(dir(pattern="*\\.vcf\\.gz"), function(f)
    {
      name = gsub("\\.gz", "", f); message(name)
      chunk = 1
      p = pipe(sprintf("zcat %s  | cut  -f '10-' | parsevcf | cut -f '1-2'", f), open="r")
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
      meta$nodename = Sys.info()["nodename"]
      meta
    }, mc.cores=NP))
  }
  message("parsing time: ", (proc.time() - t0)[[3]])
  meta$end = cumsum(meta$nrow)
  meta$start = c(1, meta$end[-length(meta$end)] + 1)
  meta$file = sprintf("%s/%s", getwd(), meta$file)
  save(meta, file="meta.rdata")
}

# 2. Principal components computation

ncomp = as.numeric(Sys.getenv("NCOMP"))
if(is.na(ncomp)) ncomp = 3

setClass("pmat", contains="list", S3methods=TRUE, slots=c(dims="numeric"))
setMethod("%*%", signature(x="pmat", y="numeric"), function(x ,y)
  {
    ans = rep(0.0, nrow(x))
    p = foreach(k=1:N, .packages=c("methods", "Matrix", "parallel"), .combine=c) %dopar%
    {
      fidx = which(x$nodename == Sys.info()["nodename"])
      q = mclapply(fidx, function(i)
      {
        f = file(x$file[i], open="rb")
        a = unserialize(f)
        close(f)
        r = attr(a, "rowmeans")
        drop(a %*% y - r * drop(crossprod(rep(1, length(y)), y)))
      }, mc.cores=NP)
      names(q) = as.character(fidx)
      q
    }
    for(j in 1:length(p))
    {
      k = as.integer(names(p)[j])
      ans[x$start[k]:x$end[k]] = p[[j]]
    }
    gc()
    ans
  })
setMethod("%*%", signature(x="numeric", y="pmat"), function(x ,y)
  {
    ans = foreach(k=1:N, .packages=c("methods", "Matrix", "parallel"), .combine="+") %dopar%
    {
      fidx = which(y$nodename == Sys.info()["nodename"])
      q = Reduce(`+`, mclapply(fidx, function(i)
      {
        f = file(y$file[i], open="rb")
        a = unserialize(f)
        close(f)
        j = seq(from=y$start[i], to=y$end[i])
        drop(x[j] %*% a - drop(crossprod(x[j], attr(a, "rowmeans"))))
      }, mc.cores=NP))
    }
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


# Shutdown the cluster and quit
closeCluster(cl)
mpi.quit()
