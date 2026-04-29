Downloads:
- Reference genomes
    - We download both reference genome versions:
        - mm39.fa (USCS naming, i.e. chr1, chr2, ...)
        - mm39.ensmembl.fa (Ensembl naming, i.e. 1, 2, ... for VCF compatibility)
    - We download the Mouse Genome Project (MGP) VCF. This contains the whole genome collections of SNP and short indel variants for 52
        inbred strains. The variants are derived from Illumina short-read sequencing data generated between 2013-2021 from a mixture of 
        read lengths (100-150bp) (See https://www.mousegenomes.org/snps-indels/ for full description). 
    - We then extract strain-specific SNPs and indels 
        - NOTE: available strains: