/*
 * Parse a 1000 genomes vcf file into an output like:
 * variant_number, sample_number, ref_number, alt_number
 * This parser requires a filter to cut out comments and info fields.
 * Example use:
 * cc -O2 parse.c
 * cat ALL.chr20.phase3_shapeit2_mvncall_integrated_v5.20130502.genotypes.vcf.gz  | zcat | sed /^#/d  | cut  -f "10-" | ./a.out 
 *
 * NOTE!! This program ignores phasing and assumes specific 1000 genomes vcf format.
 */
#include <stdio.h>
#include <stdlib.h>
int main()
{
  int k=0,t;
  char *line = NULL;
  size_t size, j, n;
  while((n = getline(&line, &size, stdin)) != -1)
  {
    k++;
    n--;   // trailing delimiter
    j=0;
    t=1;
    while(j<n)
    {
      if((j+1)==n && line[j]!=48){printf("%d,%d,%c\n",k,t,line[j]); break;} // haploid at end of line
      j++;
      if(line[j]=='|' || line[j]=='\\') // note ignore phasing XXX
      { // diploid
        j++;
        if(line[j-2]!=48 || line[j]!=48) printf("%d,%d,%c,%c\n",k,t,line[j-2],line[j]);
        j++;
      } else // haploid
      {
        if(line[j-1]!=48) printf("%d,%d,%c\n",k,t,line[j-1]);
      }
      t++;
      j++;
    }
  }
  if(line) free(line);
  return 0;
}
