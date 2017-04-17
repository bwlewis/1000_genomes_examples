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
# This implementation assumes that the input *.vcf.gz files are in a global
# shared directory available to all workers.
#
# Input: One or more variant files in *.vcf.gz
# Optional input:
#   CHUNKSIZE environment variable (number of VCF file rows per chunk)
#   NCOMP environment variable (number of components, defaults to 3)
#   OMP_NUM_THREADS environment variable number of BLAS threads per worker to use, you probably want this to be 1
# Output: pca.rdata R data file
#
# Example invocation
# 1) start soee doRedis R workers
# 2) CHUNKSIZE=10000000 Rscript --slave pca-doRedis.R

suppressMessages(
{
  library(Matrix)
  library(doRedis)
  library(methods)
  library(irlba)
})

SKIP_PARSE = (Sys.getenv("SKIP_PARSE") == "TRUE")
ONLY_PARSE = (Sys.getenv("ONLY_PARSE") == "TRUE")

QUEUE = Sys.getenv("QUEUE")
if(nchar(QUEUE) < 1) QUEUE = "RJOBS"
registerDoRedis(QUEUE)

PARSEVCF = system("which parsevcf", intern=TRUE)

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

  DIR = getwd()

  meta = foreach(f=dir(pattern="*\\.vcf\\.gz"), .combine=rbind, .packages=c("methods", "Matrix")) %dopar%
  {
    name = gsub("\\.gz", "", f); message(name)
    setwd(DIR)
    chunk = 1
    p = pipe(sprintf("zcat %s  | cut  -f '10-' | %s | cut -f '1-2'", f, PARSEVCF), open="r")
    ans = data.frame(stringsAsFactors=FALSE)
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
        ans = rbind(ans, data.frame(file=cfn, nrow=nrow(x), ncol=ncol(x), stringsAsFactors=FALSE))
        chunk = chunk + 1
      }
      rm(x)
      gc()
    }
    close(p)
    ans$nodename = Sys.info()["nodename"]
    ans
  }
  meta$end = cumsum(meta$nrow)
  meta$start = c(1, meta$end[-length(meta$end)] + 1)
  meta$file = sprintf("%s/%s", getwd(), meta$file)
  save(meta, file="meta.rdata")
  message("parsing time: ", (proc.time() - t0)[[3]])
}

if(ONLY_PARSE) q(save="no")

# 2. Principal components computation

# 2a) Initialize the foreach workers
workers = unique(foreach(j=1:nrow(meta), .packages="Matrix", .combine=c) %dopar%
{
  if(!exists("values", envir=globalenv())) assign("values", list(), envir=globalenv())
  values[sprintf("%d", j)] = readRDS(as.character(meta$file[j]))
  assign("values", values, envir=globalenv())
  Sys.getpid()
})


# 2b) Run the principal components computation
ncomp = as.numeric(Sys.getenv("NCOMP"))
if(is.na(ncomp)) ncomp = 3

setClass("pmat", contains="list", S3methods=TRUE, slots=c(dims="numeric"))
setMethod("%*%", signature(x="pmat", y="numeric"), function(x ,y)
  {
    p = foreach(k=workers, .packages=c("methods", "Matrix"), .combine=c) %dopar%
    {
      q = lapply(values, function(a)
      {
        r = attr(a, "rowmeans")
        drop(a %*% y - r * drop(crossprod(rep(1, length(y)), y)))
      })
      names(q) = names(values)
      q
    }
    ans = rep(0.0, nrow(x))
    for(j in 1:length(p))
    {
      k = as.integer(names(p)[j])
      ans[x$start[k]:x$end[k]] = p[[j]]
    }
    gc()
    ans
  })
setMethod("%*%", signature(x="numeric", y="pmat"), function(x, y)
  {
    ans = foreach(k=workers, .packages=c("methods", "Matrix"), .combine="+") %dopar%
    {
      q = Reduce(`+`, lapply(1:length(values), function(k)
      {
        i = as.integer(names(values)[k])
        a = values[[k]]
        j = seq(from=y$start[i], to=y$end[i])
        drop(x[j] %*% a - drop(crossprod(x[j], attr(a, "rowmeans"))))
      }))
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
