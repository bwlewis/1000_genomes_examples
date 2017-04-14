# VCF Principal Components Example, SMP variation
# Assumes that parse.c
# (https://raw.githubusercontent.com/bwlewis/1000_genomes_examples/master/parse.c)
# is compiled and copied to executable parsevcf program in PATH on every node.
# For instance:
#
# cc -O2 -o parsevcf parse.c
# sudo mv parsevcf /usr/local/bin
#
# Input: One or more variant files in *.vcf.gz
# Optional input:
#   CHUNKSIZE environment variable number nonzero matrix elements per chunk, default 1e8
#   NP environment variable number of R worker CPUs to use, default parallel::detectCores()
#   OMP_NUM_THREADS environment variable number of BLAS threads per worker to use, you probably want this to be 1
#   NCOMP environment variable number of components, default 3
# Output: pca.rdata R data file, a list with principal component vectors and values,
#         meta.rdata temporary file left from the parsing step for optional re-use.
# Example invocation
# NP=8  Rscript --slave pca-smp.R

suppressMessages(
{
  library(Matrix)
  library(methods)
  library(irlba)
  library(parallel)
})

NP = as.integer(Sys.getenv("NP"))
if(is.na(NP)) NP = detectCores()
message("cores ", NP)

# 1. VCF parsing
# Parse the VCF input data files in *.vcf.gz into partitioned sparse
# R submatrices saving the partitions to R data files for use by the PCA computation following
# this section.

# Set the environment variable SKIP_PARSE=TRUE to skip parsing step
# and use the output from a previous run.
SKIP_PARSE = (Sys.getenv("SKIP_PARSE") == "TRUE")

if(SKIP_PARSE)
{
  message("Loading previously parsed data")
  load("meta.rdata")
} else
{
  t0 = proc.time()
  chunksize = as.numeric(Sys.getenv("CHUNKSIZE"))
  if(is.na(chunksize))
  {
     chunksize = 1e8
  }
  message("chunksize: ", chunksize)
  a = mclapply(dir(pattern="*\\.vcf\\.gz"), function(f)
    {
      name = gsub("\\.gz", "", f); message(name)
# work around for http://r.789695.n4.nabble.com/R-3-0-1-parallel-collection-triggers-quot-long-memory-not-supported-yet-quot-td4668388.html  -- we save temporary parse results to files and then collate instead of using sendMaster() in the mc code, which limits partial results to 2GB :/
      vfd = sprintf("%s.temp.rdata", name)
      mfd = sprintf("%s.meta.temp.rdata", name)
      chunk = 1
      p = pipe(sprintf("zcat %s  | cut  -f '10-' | parsevcf | cut -f '1-2'", f), open="r")
      meta = data.frame()
      values = list()
      while(chunk > 0)
      {
        x = tryCatch(read.table(p, colClasses=c("integer", "integer"), fill=TRUE, row.names=NULL, nrows=chunksize),
                     error=function(e) data.frame())
        if(nrow(x) < 1) chunk = 0
        else
        {
          x = sparseMatrix(i=x[, 1] - x[1, 1] + 1, j=x[, 2], x=1.0)
          attr(x, "rowmeans") = rowMeans(x)
          meta = rbind(meta, data.frame(source_file=f, file_chunk=chunk, nrow=nrow(x), ncol=ncol(x), stringsAsFactors=FALSE))
          values = c(values, x)
          chunk = chunk + 1
        }
      }
      close(p)
      saveRDS(meta, file=mfd)
      saveRDS(values, file=vfd)
      c(mfd, vfd)
  }, mc.cores=NP)
  meta = Reduce(rbind, Map(function(x)
    {
      readRDS(x[1])
    }, a))
  meta$end = cumsum(meta$nrow)
  meta$start = c(1, meta$end[-length(meta$end)] + 1)
  meta = as.list(meta)
  meta$values = Reduce(c, Map(function(x)
    {
      readRDS(x[2])
    }, a))
  message("parsing time: ", (proc.time() - t0)[[3]])
  unlink(dir(pattern=".*.temp.rdata"))
  save(meta, file="meta.rdata", compress=FALSE) # for optional re-use (NOTE: this may take a *long* time without compress=FALSE)
}


# 2. Principal components computation

ncomp = as.numeric(Sys.getenv("NCOMP"))
if(is.na(ncomp)) ncomp = 4

setClass("pmat", contains="list", S3methods=TRUE, slots=c(dims="numeric"))
setMethod("%*%", signature(x="pmat", y="numeric"), function(x ,y)
  {
    ans = rep(0.0, nrow(x))
    p = mcMap(function(i)
    {
      drop(x$values[[i]] %*% y - attr(x$values[[i]], "rowmeans") * drop(crossprod(rep(1, length(y)), y)))
    }, seq(1, length(x$start)), mc.cores=NP)
    for(j in 1:length(p))
    {
      ans[x$start[j]:x$end[j]] = p[[j]]
    }
    ans
  })
setMethod("%*%", signature(x="numeric", y="pmat"), function(x ,y)
  {
    Reduce(`+`, mcMap(function(i)
    {
      j = seq(from=y$start[i], to=y$end[i])
      drop(x[j] %*% y$values[[i]] - drop(crossprod(x[j], attr(y$values[[i]], "rowmeans"))))
    }, seq(1, length(y$start)), mc.cores=NP))
  })

A = new("pmat", meta, dims=c(tail(meta$end, 1), meta$ncol[1]))
dim.pmat = function(x) x@dims
nrow.pmat = function(x) x@dims[1]
ncol.pmat = function(x) x@dims[2]

message("Starting PCA computation")
t1 = proc.time()
L  = irlba(A, nv=ncomp, tol=1e-5, right_only=TRUE, work=4)
dt = proc.time() - t1
saveRDS(L, file="pca.rdata")
message("PCA time: ", dt[[3]])
