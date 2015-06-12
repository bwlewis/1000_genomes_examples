# 1000_genomes_examples

The 1000 Genomes project (http://www.1000genomes.org/) is an amazing public
resource coordinated by the NIH. It's a "deep catalog of human genetic
variation" for more than 2500 people, all carefully curated, organized and
published on the web.

Because of the size and complexity of the data, and because of the significant
interest in understanding human genetics, these data are often used as examples
of "big data" analysis. A typical example of this kind of thing can be found
in the Google genomics project https://github.com/googlegenomics.

## That's cool and all, but...

my personal opinion is that folks get too worked up over "big data" solutions
and sometimes too focused on distributed system implementation details (Spark,
etc.) when simple tools like R work just fine.

I am very impressed by a recent article by Michael Isard, Frank McSherry and
Derek Murray (summarized 
http://www.frankmcsherry.org/graph/scalability/cost/2015/01/15/COST.html).  In
that article, McSherry and his coauthors show that careful thinking about some
network analysis problems let them run better on a laptop than on a 128 CPU
core cluster using a few common "big data" computing frameworks.

This note presents a few simple but useful examples using the 1000 genomes data
in the spirit of Isard, McSherry and Murray. The examples are often used to
showcase big data analysis frameworks. This note shows how they can sometimes
be easily and efficiently solved on a decent laptop.

## Variant data

The examples use "variant call format" (VCF) files following an NCBI-specific
format available from links shown in the code snippets below. Loosely,
"variants" are places on the genome that vary from a reference genome in a
cataloged way. Variants include single-nucleotide polymorphisms (a.k.a. SNPs,
basically a single base change along the genome) and larger "structural"
alterations. The 1000 genomes project catalogs about 81 million variants.

Variant data are stored in files by chromosome. Each file contains a set of
header lines beginning with the `#` character, followed by one line per
variant. 
Full file format details are described in
www.1000genomes.org/wiki/analysis/variant%20call%20format/vcf-variant-call-format-version-41.

Variant data files for each chromosome are available from
ftp://ftp-trace.ncbi.nih.gov/1000genomes/ftp/release/20130502/.

Discounting the header lines, a variant line in the data files consists of some
information columns about the variant followed by one column for each sample
(person) that indicates if they exhibit the variant on either or both
chromosomes. For example, part of a typical line (showing only the first 5
columns and 10-15 columns) looks like:

```
zcat ALL.chr20.phase3_shapeit2_mvncall_integrated_v5.20130502.genotypes.vcf.gz | sed /^#/d | cut -f "1-5,10-15" | head -n 1

20      60343   .       G       A   ...    0|0     0|0     0|0     0|0     0|0     0|0
```
This variant is on chromosome 20 at position 60343. The reference nucleotide is G and
the variant is A. Note that in some cases, more than one possibility may be listed
in which case the variants are numbered 1,2,3,...
Columns 10-15 show that none of the first 6 people in the database
have this variant on either strand of DNA `0|0`. Someone with the G to A variant
on the 2nd strand of DNA will display `0|1`, for example.

Numerous full-featured VCF file parsers exist for R, see for example 
the http://bioconductor.org project. But the simple
analyses considered in this project don't need to read VCF files in full
generality, and we can also benefit from the knowledge that the 1000 genomes
project follows a somewhat restricted VCF subset. I wrote the really simple
32-line C parser program
https://github.com/bwlewis/1000_genomes_examples/blob/master/parse.c to take
advantage of these observations and load a subset of VCF data from the 1000
genomes project into R pretty quickly.

The simple parser program turns VCF files into comma-separated output with four
or three columns: variant number (just the line offset in file), sample number
(person), alternate number on first strand or haploid alternate, optional
alternate number on 2nd strand (diploid). For example chromosome 20 again:
```
cc parse.c
zcat ALL.chr20.phase3_shapeit2_mvncall_integrated_v5.20130502.genotypes.vcf.gz  | sed /^#/d  | cut  -f "10-" | ./a.out | head

1,898,0,1
2,1454,0,1
3,462,0,1
3,463,1,0
```
The output says that person 898 has variant number 1 (the first listed) for
chromosome 20 present on the 2nd strand of DNA. ANd person 1454 has variant
number 2 present on the 2nd strand of DNA, and so on.

For our purposes in the following examples, this simple C parser quickly
converts the 1000 genomes VCF data into a variant number by person number
table.

The 1000 genomes project also includes ethnicity and some other phenotypic data
about each of the more than 2500 people whose genomes are represented.
Those data are available from:
ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/working/20130606_sample_info/20130606_g1k.ped
and are easy to directly read into R.


## Principal components analysis

This example is very basic but often used to illustrate analysis of 1000
genomes data. It's popular perhaps because it's very effective at clustering
the people by ethnicity, but otherwise I'm not sure how interesting it really
is!

The example uses principal components analysis to project all of the variant
data for one chromosome into a three-dimensional subspace, and then plot
the result.

The example uses:
- my simple C parsing program to read variant data into an R sparse matrix,
- the irlba package (development version at https://github.com/bwlewis/IRL) to efficiently compute principal components,
- the threejs package to visualize the result.

All of these steps, from reading the data in to visualization, only take a few
minutes on a decent laptop, expressed in just a few lines of R code.

### Reading variant data into an R sparse matrix

This step assumes that you've downloaded and compiled my simple VCF parser and
downloaded at least the chromosome 20 and phenotype data files from the 1000
genomes project, for example (from a Mac or Linux shell):
```
# 1000 genomes example variant data file (chromosome 20)
wget ftp://ftp-trace.ncbi.nih.gov/1000genomes/ftp/release/20130502/ALL.chr20.phase3_shapeit2_mvncall_integrated_v5.20130502.genotypes.vcf.gz

# 1000 genomes phenotype data file
wget ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/working/20130606_sample_info/20130606_g1k.ped

# My simple but fast parser program (after compilation you'll have a program called a.out)
wget https://raw.githubusercontent.com/bwlewis/1000_genomes_examples/master/parse.c
cc -O2 parse.c
```
Note that we _could_ use R alone to read and parse the VCF file, it would just
take a little longer.

All the remaining steps in this document run from R. Let's read the variant
data for chromosome 20 into an R sparse matrix. Note that we only care about
the variant number and sample (person) number in this exercise and ignore
everything else.
```{r}
library(Matrix)
p = pipe("zcat ALL.chr20.phase3_shapeit2_mvncall_integrated_v5.20130502.genotypes.vcf.gz  | sed /^#/d  | cut  -f '10-' | ./a.out | cut -f '1-2'")
x = read.table(p, colClasses=c("integer","integer"), fill=TRUE, row.names=NULL)

# Convert to a sparse matrix of people (rows) x variant (columns)
chr20 = sparseMatrix(i=x[,2], j=x[,1], x=1.0)

# Inspect the dimensions of this matrix
print(dim(chr20))
# [1]    2504 1812841
```
That was pretty easy!  We've loaded a sparse matrix with 2504 rows (people) by
1,812,841 columns (variants).  The next step computes the first three principal
component vectors with the irlba package and plots them with a cool 3d
scatterplot. It's pretty easy too!
```{r}
library(irlba)
cm = colMeans(chr20)
p = irlba(chr20, nv=3, nu=3, tol=0.1, dU=rep(1,nrow(chr20)), ds=1, dV=cm)

library(threejs)
scatterplot3js(p$u)
```
Not bad. We see that there are three obvious groups.
