#  Obtaining input data

First you need to make a Hi-C sequencing library from gametes of a single individual (donor). We have been using the Dovetail Omni-C kit, because it's endonuclease-based chromating digestion delivers a more uniform coverage across the genome when compared with restriction-enzyme-based kits. For starting material, we aimed for between 100,000 and 1 million fresh or flash frozen gametes, then we followed the standard protocol, and then aimed fot about 100x coverage paired-end Illumina sequencing. After obtaining the gamete Hi-C library, you can process the data for example as follows:

## 1. Alignment and duplicate detection

- Align to a reference genome using [bwa](https://bio-bwa.sourceforge.net) with `bwa mem` and importantly the `-5SP -T0` options 

- Sort the alignment with [samtools](http://www.htslib.org) using `samtools sort` 

- Mark duplicates with [picard](https://broadinstitute.github.io/picard/) using `MarkDuplicates`

## 2. Variant calling/filtering

- Call variants directly using the Hi-C reads with [bcftools](http://www.htslib.org). Example:

```console
$ bcftools mpileup --count-orphans -Ou -f reference.fa DEDUPLICATED_READS.bam | bcftools call --threads 2 -mv -Oz -o UNFILTERED_VARIANTS.vcf.gz 
```

- Filter variants. Variant filtering is a complex topic; luckily, `Hi-reComb` is reasonably robust to errors as we demostrated in the manuscript. For the analyses in the manuscript, we used filering on depth, masking repetitive regions, and a set of hard filters.

 ## 3. Haplotype phasing

- First divide up the VCF file per-chromosome using [bcftools](http://www.htslib.org) `bcftools view`.

- If using the Hi-C reads, also divide the BAM file per-chromosome with [samtools](http://www.htslib.org) `samtools view` and then run [hapcut2](https://github.com/vibansal/HapCUT2) to obtain phased donor genotypes as follows:

```console
extractHAIRS --hic 1 --VCF chr1.vcf --bam chr1.bam --out hapcutFragments_chr1.txt

hapcut2 --VCF chr1.vcf --hic 1 --fragments hapcutFragments_chr1.txt --output HAPCUT2_PHASE.txt
```

- If using trio phasing, you can go back to instructions for `Hi-reComb TrioPhase`
  
  
