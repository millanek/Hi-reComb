### Description of phased haplotype file (block format).

The phased haplotypes should be provided to Hi-reComb in a block format, as output by HapCUT2. The full HapCUT2 format is described below, taken from  https://github.com/vibansal/HapCUT2/blob/master/outputformat.md

Each block starts with a block header with the following format:

BLOCK: offset: \<SNV offset\> len: \<SNV span of block\> phased: \<\# SNVs phased\> SPAN: \<base pair span of block\> fragments \<\# of fragments in block\>

Following the header, there is one line per SNV with the following tab-delimited fields:

1. VCF file index (1-based index of the line in the input VCF describing variant)
2. allele on haploid chromosome copy A (0 means reference allele, 1 means variant allele, - for an unphased variant)
3. allele on haploid chromosome copy B (0 means reference allele, 1 means variant allele, - for an unphased variant)
4. chromosome
5. position
6. reference allele (allele corresponding to 0 in column 2 or 3)
7. variant allele (allele corresponding to 1 in column 2 or 3)
8. VCF genotype field (unedited, directly from original VCF)
9. discrete pruning status (1 means pruned, 0 means phased)
10. switch quality: phred-scaled estimated probability that there is a switch error starting at this SNV (0 means switch error is likely, 100 means switch is unlikely)
11. mismatch quality: phred-scaled estimated probability that there is a mismatch [single SNV] error at this SNV (0 means SNV is low quality, 100 means SNV is high quality)
12. number of haplotype-informative fragments covering this variant in the input

Each block ends with a line with '********'. 

Some variants can be unphased after post-processing. Such variants have the allele as '-'  in columns 2 and 3. 


