#  Data preprocessing

After obtaining the gamete Hi-C library, you can process the data for example as follows:

## 1. Alignment and duplicate detection

- Align to a reference genome using `bwa mem`, importantly with the `-5SP -T0` options [link](https://bio-bwa.sourceforge.net) 

- Sort the alignment with `samtools sort` [link](http://www.htslib.org) 

- Mark duplicates using `picard MarkDuplicates` [link](https://broadinstitute.github.io/picard/) 

## 2. Variant calling/filtering

- Call variants directly using the Hi-C reads with `bcftools` [link](http://www.htslib.org). Example:

```console
$ bcftools mpileup --count-orphans -Ou -f reference.fa DEDUPLICATED_READS.bam | bcftools call --threads 2 -mv -Oz -o UNFILTERED_VARIANTS.vcf.gz 
```

- Filter variants.Variant filtering is a complex topic. However, `Hi-reComb` is relatively robust to errors. For the analyses in the manuscript, we used filering on depth, masking repetitive regions, and a set of hard filters.

 ## 3. Haplotype phasing
