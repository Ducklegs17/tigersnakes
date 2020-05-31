## Brief Introduction

Snakemake was used to create the variant calling workflow. A visual representation of the rulegraph for the workflow can be seen in the TigerSnakeRulegraph.png. Descriptions of the purpose of each of the rules can be found below. For further details regarding pre-processing, see the Snakefile where the rules are defined and R/tiger_snake_project.Rmd for a brief description of pre-processing steps. Further analysis was undertaken in R and can again be found in the R/tiger_snake_project.Rmd and tigerSnake_pt2.Rmd files. Result files from this analysis are located in the results folder (see below for a description of each type of result file).
 
The following links lead to Rmd reports detailing the project (also stored in ./R/).

- [01 Identification of variants]( https://htmlpreview.github.io/?https://github.com/Ducklegs17/tigersnakes/blob/master/R/tiger_snake_project.html)
- [02 Identification of genes with highest non-synonymous variant density]( https://htmlpreview.github.io/?https://github.com/Ducklegs17/tigersnakes/blob/master/R/tigerSnake_pt2.html)

Separate files are provided for growth, lipid and pigmentation gene sets. 

**geneset_promoter_variants**

Variants identified by all three variant callers that pass coverage filters and are located within the promoter region of the genes of interest. 

**geneset_nonsynonymous_variants**

All non-synonymous variants located within the genes of interest (includes the reference and alternate allele)

**geneset_genes_without_variants**

Genes with zero non-synonymous variants.

**geneset_variant_density**

Genes with calculated non-synonymous variant density per kbp of coding regions.

## Explanation of rule names:

**bwa_index:**
Indexes reads. Necessary for alignment.

**cutadapt:**
Removes poly-G tails.

**bwa_mem:**
Aligns reads to reference genome.

**name_sort:**
Sorts aligned reads in .bam by name. Necessary for fixmate.

**fixmate:**
Fills in mate coordinates and mate related flags (requires a name-sorted bam). Necessary for compatibility with gatk.

**position_sort:**
Sorts by location (necessary for markduplicates).

**markduplicates:**
Marks duplicate alignments (necessary for Gatk compatibility).

**gunzip_ref:**
Unzips the reference genome.

**freebayes_call:**
Calls variants using freebayes.

**freebayes_filter:**
Filters freebayes variants using vcftools.

**addOrReplaceReadGroups:**
Adds read group tags (required for Gatk).

**bcftools_call:**
Calls variants with bcftools.

**bcftools_bcfToVcf:**
Converts bcf file to vcf.

**bcftools_filter:**
Filters bcftools variants using vcftools.

**join_bcf_vcfs:**
Joins bcftools .vcf files into a single vcf.

**promoter_regions:**
Creates .bed files describing promoter regions of the genes of interest from the .bed files describing the coordinates of the genes of interest creates in R. 

**picard_dict:**
Creates a sequence dictionary. Used to create region files necessary for parallelising bcftools_call. 

**gatk_make_region_files:**
Makes files used as input into gatk_call to parallelise gatk_call.

**bcf_make_region_files:**
Edits gatk region files to make files compatible with bcftools for variant calling.

**gatk_call:**
Call variants with gatk.

**gatk_filter:**
Filter variants called with Gatk using vcftools.

**join_gatk_vcfs:**
Join resulting gatk .vcf files together into a single vcf. 

**bedtools:**
Overlaps vcf files (variants) with .bed files specifying target genes.

