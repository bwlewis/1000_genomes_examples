# 1000_genomes_examples

The 1000 Genomes project (http://www.1000genomes.org/) is an amazing public
resource coordinated by the NIH. It's a "deep catalog of human genetic
variation" for more than 2500 people, all carefully curated, organized and
published on the web.

Because of the size and complexity of the data, and because of the significant
interest in understanding human genetics, these data are often used as examples
of "big data" analysis. A typical example of this kind of thing can be found
in the Google genomics project https://github.com/googlegenomics.

## That's cool and all, but

my personal opinion is that folks get too worked up over "big data" solutions
and sometimes too focused on distributed system implementation details (Spark,
etc.) when simple tools like R work just fine.

I am very impressed by a recent paper by Michael Isard, Frank McSherry and
Derek Murray (summarized 
http://www.frankmcsherry.org/graph/scalability/cost/2015/01/15/COST.html).  In
that article, McSherry and his coauthors show that careful thinking about some
network analysis problems let them run better on a laptop than on a 128 CPU
core cluster using a few common "big data" computing frameworks.

We present a few simple but useful examples with the 1000 genomes data in the
spirit of Isard, McSheery and Murray. Although our examples are often used to
showcase big data analysis frameworks, we show how they can sometimes be easily
and efficiently solved on a decent laptop.

## Variant data

Our examples use "variant call format" (VCF) files following an NCBI-specific
format available from links shown in the code snippets below. Loosely,
"variants" are places on the genome that vary from a reference genome in a
cataloged way. Variants include single-nucleotide polymorphisms (a.k.a. SNPs,
basically a single base change along the genome) and larger "structural"
alterations. The 1000 genomes project catalogs about 81 million variants.

Variant data are stored in files by chromosome. Each file contains a set of
header lines beginning with the `#` character, followed by one line per
variant.  Some file parsing details are described in the examples below.
Full file format details are described in
www.1000genomes.org/wiki/analysis/variant%20call%20format/vcf-variant-call-format-version-41.

The 1000 genomes project also includes ethnicity and some other specific data
about each of the more than 2500 people whose genomes are represented.


Discounting the header lines, a variant line in the data files consists


## Principal components analysis


