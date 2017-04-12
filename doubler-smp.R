# Simulate double the number of genomes by replicating data and randomly
# removing about 10% of the entries in the replications.
# Assumes the chunked 'meta' data list is loaded, for instance from
# the parsing output in meta.rdata. For use with pca-smp.R.


r = meta$values
for(k in 1:length(r))
{
  N = length(r[[k]]@x)
  r[[k]]@x[sample(N, ceiling(N / 10))] = 0
}

for(k in 1:4) meta[[k]] = c(meta[[k]], meta[[k]])
meta$end = cumsum(meta$nrow)
meta$start = c(1, head(meta$end, -1) + 1)

# Number of genomes doubled. Save it with, for instance:
# save(meta, file="meta.rdata", compress=FALSE)
