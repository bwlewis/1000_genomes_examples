/*
 * Parse a 1000 genomes vcf file into an output like
 * variant_number\tsample_number
 * that is, it only shows whether or not a variant occured not what it is.
 * This parser requires a filter to cut out comments and info fields.
 * Example use:
 * cc -O2 parse.c
 * zcat ALL.chr20.phase3_shapeit2_mvncall_integrated_v5.20130502.genotypes.vcf.gz  | cut  -f "10-" | ./a.out 
 *
 * NOTE!! This program ignores phasing and assumes specific 1000 genomes vcf format.
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
int main(int argc, char **argv)
{
  int k = 0, t;
  char *entry;
  char *line = NULL;
  size_t size, n;
  while((n = getline(&line, &size, stdin)) != -1)
  {
    if(line[0] == '#') continue;
    k++;
    t = 1;
    entry = strtok(line, "\t");
    while(entry != NULL)
    {
      if(strncmp("0|0", entry, 3) != 0)
        printf("%d\t%d\n", k, t);
      entry = strtok(NULL, "\t");
      t++;
    }
  }
  if(line) free(line);
  return 0;
}
