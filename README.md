# Hi-reComb

## Estimating recombination maps from sperm Hi-C data

Author: Milan Malinsky 
Email: millanek@gmail.com

## Quickstart:
```
Commands:
           FindInfoPairs       Find read pairs that would be informative for estimating the recombination
           RecombMap           Estimate the recombination map from informative Hi-C read pairs
Utilities:
           Simulate            Simulate informative read-pairs for a given map to evaluate confidence in reconstruction
           TrioPhase           Generate a phased het file for use with Hi-Recomb from a VCF with trio(s) (mother-father-offspring)

Usage:
1) samtools view ALIGNEMENT.bam | Hi-reComb FindInfoPairs [OPTIONS] HAPCUT2_PHASE.txt > INFORMATIVE_READS.sam
2) Hi-reComb RecombMap [OPTIONS] HAPCUT2_PHASE.txt INFORMATVE_PAIRS.sam

Utilities usage:
Hi-reComb Simulate [OPTIONS] ReferenceMap.txt INFORMATVE_PAIRS.sam
Hi-reComb TrioPhase [OPTIONS] OFFSPRING.vcf(.gz) PARENTS.vcf(.gz) PARENT1,PARENT2
```

## Input files:
### Required files:

All files should be divided by chromosome/scaffold.

1. An alignment file with Hi-C read pairs from gametes (sperm, pollen ...) from a single individual in the [SAM / BAM](https://samtools.github.io/hts-specs/SAMv1.pdf) format.   
2. Phased SNPs at sites that are heterozygous in the donor (the individual from whom the gametes were obtained from). The phasing can be obtained directly from the Hi-C data using the [hapcut2 software](https://github.com/vibansal/HapCUT2) (make sure you use the `--hic 1` option). Alternatively, if you have variant calls for the parents of the donor individual, you can use the `Hi-reComb TrioPhase` utility. In any case, the 'HAPCUT2_PHASE.txt' file should be in the [hapcut2 format](HapCUT2_format.md).

For suggestions/examples on producing the input files see [here](Data_preprocessing.md). 
   
## Installation
To compile you must have a reasonably recent GCC (>=4.9.0) or clang compiler (on mac OS this comes with 'Command Line Tools').

```console
$ git clone https://github.com/millanek/Hi-reComb
$ cd Hi-reComb
$ make
```

The Hi-reComb executable will be in the Build folder, so to run it type e.g. `./Build/Hi-reComb`; this will show the available commands. To execute e.g. the RecombMap command, type `./Build/Hi-reComb RecombMap`.
