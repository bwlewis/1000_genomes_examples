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
# Optional environment variable input:
#   CHUNKSIZE (number of nonzero matrix elements per chunk)
#   NP (number of worker processes, defaults to CPU cores)
#   NCOMP (number of components, defaults to 3)
#   OMP_NUM_THREADS number of BLAS threads per worker to use, you probably want this to be 1
#   SKIP_PARSE if SKIP_PARSE=TRUE, then skip the VCF data file parsing step
#   ONLY_PARSE if ONLY_PARSE=TRUE, then only parse the VCF data files
#
# Output: pca.rdata R data file, parsed VCF data sparse matrix partitions in *.pmat.rdata
#
# Example invocations
# PATH=$(pwd):$PATH CHUNKSIZE=10000000 R --slave < pca-fork.R
# PATH=$(pwd):$PATH SKIP_PARSE=TRUE R --slave < pca-fork.R

suppressMessages(
{
  library(Matrix)
  library(parallel)
  library(methods)
  library(irlba)
})

SKIP_PARSE = (Sys.getenv("SKIP_PARSE") == "TRUE")
ONLY_PARSE = (Sys.getenv("ONLY_PARSE") == "TRUE")
PARSEVCF = system("which parsevcf", intern=TRUE)

NP = as.integer(Sys.getenv("NP"))
if(is.na(NP)) NP = detectCores()

NCOMP = as.numeric(Sys.getenv("NCOMP"))
if(is.na(NCOMP)) NCOMP = 3

cl = makeForkCluster(NP)

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
  clusterExport(cl, "chunksize", environment())

  meta = Reduce(rbind, parLapplyLB(cl, dir(pattern="*\\.vcf\\.gz"), function(f)
  {
    name = gsub("\\.gz", "", f); message(name)
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
  }))
  meta$end = cumsum(meta$nrow)
  meta$start = c(1, meta$end[-length(meta$end)] + 1)
  meta$file = sprintf("%s/%s", getwd(), meta$file)
  save(meta, file="meta.rdata")
  message("parsing time: ", (proc.time() - t0)[[3]])
}

if(ONLY_PARSE) q(save="no")

clusterExport(cl, "meta", environment())

# 2. Principal components computation

# 2a) Initialize the workers by dividing up the matrix partitions
#     between them and loading them into worker memory spaces.
null = parLapplyLB(cl, 1:nrow(meta), function(j)
{
  if(!exists("values", envir=globalenv())) assign("values", list(), envir=globalenv())
  values[sprintf("%d", j)] = readRDS(as.character(meta$file[j]))
  assign("values", values, envir=globalenv())
  NULL
})


# 2b) Run the principal components computation
setClass("pmat", contains="list", S3methods=TRUE, slots=c(dims="numeric"))
setMethod("%*%", signature(x="pmat", y="numeric"), function(x ,y)
  {
    clusterExport(cl, "y", environment())
    p = Reduce(c, clusterEvalQ(cl,
    {
    if(exists("values"))
     {
      q = lapply(values, function(a)
      {
        r = attr(a, "rowmeans")
        drop(a %*% y - r * drop(crossprod(rep(1, length(y)), y)))
      })
      names(q) = names(values)
     } else q = NULL
     q
    }))
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
    clusterExport(cl, c("x", "y"), environment())
    ans = Reduce(`+`, clusterEvalQ(cl,
    {
     if(exists("values"))
     {
      q = Reduce(`+`, lapply(1:length(values), function(k)
      {
        i = as.integer(names(values)[k])
        a = values[[k]]
        j = seq(from=y$start[i], to=y$end[i])
        drop(x[j] %*% a - drop(crossprod(x[j], attr(a, "rowmeans"))))
      }))
     } else q = 0
    }))
  })

A = new("pmat", as.list(meta), dims=c(tail(meta$end, 1), meta$ncol[1]))
dim.pmat = function(x) x@dims
nrow.pmat = function(x) x@dims[1]
ncol.pmat = function(x) x@dims[2]

t1 = proc.time()
L  = irlba(A, nv=NCOMP, tol=1e-5, right_only=TRUE, work=4)
dt = proc.time() - t1
saveRDS(L, file="pca.rdata")
stopCluster(cl)
message("PCA time: ", dt[[3]])
