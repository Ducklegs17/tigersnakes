SAMPLES = ["VT211_HG2W3DSXX_GATTCTGC-GACGAGAG_L001",]
MAX_THREADS = 32
NUMBERS = ["1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13", "14", "15", "16", "17", "18", "19", "20", "21", "22", "23", "24", "25", "26", "27", "28", "29", "30", "31", "32", "33", "34", "35", "36", "37", "38", "39", "40", "41", "42", "43", "44", "45", "46", "47", "48", "49", "50", "51", "52", "53", "54", "55", "56", "57", "58", "59", "60", "61", "62", "63", "64", "65", "66", "67", "68", "69", "70", "71", "72", "73", "74", "75", "76", "77", "78", "79", "80", "81", "82", "83", "84", "85", "86", "87", "88", "89", "90", "91", "92", "93", "94", "95", "96", "97", "98", "99", "100", "101", "102", "103", "104", "105", "106", "107", "108", "109", "110", "111", "112", "113", "114", "115", "116", "117", "118", "119", "120", "121", "122", "123", "124", "125", "126", "127", "128", "129", "130", "131", "132", "133", "134", "135", "136", "137", "138", "139", "140", "141", "142", "143", "144", "145", "146", "147", "148", "149", "150", "151", "152", "153", "154", "155", "156", "157", "158", "159", "160", "161", "162", "163", "164", "165", "166", "167", "168", "169", "170", "171", "172", "173", "174", "175", "176", "177", "178", "179", "180", "181", "182", "183", "184", "185", "186", "187", "188", "189", "190", "191", "192", "193", "194", "195", "196", "197", "198", "199", "200", "201", "202", "203", "204", "205", "206", "207", "208", "209", "210", "211", "212", "213", "214", "215", "216", "217", "218", "219", "220", "221", "222", "223", "224", "225", "226", "227", "228", "229", "230", "231", "232", "233", "234", "235", "236", "237", "238", "239", "240", "241", "242", "243", "244", "245", "246", "247", "248", "249", "250", "251", "252", "253", "254", "255", "256", "257", "258", "259", "260", "261", "262", "263", "264", "265", "266", "267", "268", "269", "270", "271", "272", "273", "274", "275", "276", "277", "278", "279", "280", "281", "282", "283", "284", "285", "286", "287", "288", "289", "290", "291", "292", "293", "294", "295", "296", "297", "298", "299", "300", "301", "302", "303", "304", "305", "306", "307", "308", "309", "310", "311", "312", "313", "314", "315", "316", "317", "318", "319", "320", "321", "322", "323", "324", "325", "326", "327", "328", "329", "330", "331", "332", "333", "334", "335", "336", "337", "338", "339", "340", "341", "342", "343", "344", "345", "346", "347", "348", "349", "350", "351", "352", "353", "354", "355", "356", "357", "358", "359", "360", "361", "362", "363", "364", "365", "366", "367", "368", "369", "370", "371", "372", "373", "374", "375", "376", "377", "378", "379", "380", "381", "382", "383", "384", "385", "386", "387", "388", "389", "390", "391", "392", "393", "394", "395", "396", "397", "398", "399", "400", "401", "402", "403", "404", "405", "406", "407", "408", "409", "410", "411", "412", "413", "414", "415", "416", "417", "418", "419", "420", "421", "422", "423", "424", "425", "426", "427", "428", "429", "430", "431", "432", "433", "434", "435", "436", "437", "438", "439", "440", "441", "442", "443", "444", "445", "446", "447", "448", "449", "450", "451", "452", "453", "454", "455", "456", "457", "458", "459", "460", "461", "462", "463", "464", "465", "466", "467", "468", "469", "470", "471", "472", "473", "474", "475", "476", "477", "478", "479", "480", "481", "482", "483", "484", "485", "486", "487", "488", "489", "490", "491", "492", "493", "494", "495", "496", "497", "498", "499", "500", "501", "502", "503", "504", "505", "506", "507", "508", "509", "510", "511", "512", "513", "514", "515", "516", "517", "518", "519", "520", "521", "522", "523", "524", "525", "526", "527", "528", "529", "530", "531", "532", "533",]

singularity:
	"docker://continuumio/miniconda3:4.5.12"

############################
# Include other Snakefiles #
############################
include:
	"rules/misc.smk"

#######################################
# Convenient rules to define targets #
#######################################
localrules:
	all

rule all:
	input:
		expand("0_rawData/FastQC/{SAMPLE}_{read}_fastqc.html", SAMPLE=SAMPLES, read=['R1', 'R2']),
		expand("reference/GCF_900518725.1_TS10Xv2-PRI_genomic.fna.gz.{ext}",ext=['amb', 'ann', 'bwt', 'pac', 'sa']),
		expand("1_trimmedData/fastq/{SAMPLE}_{reads}.fastq.gz",SAMPLE=SAMPLES, reads=['R1', 'R2']),
                expand("1_trimmedData/FastQC/{SAMPLE}_{reads}_fastqc.html", SAMPLE=SAMPLES, reads=['R1', 'R2']),
		expand("2_alignedData/bam/{SAMPLE}.bam", SAMPLE=SAMPLES),
		"3_snps/gatk/VT211_HG2W3DSXX_GATTCTGC-GACGAGAG_L001_gatk_filtered.vcf",
		"reference/GCF_900518725.1_TS10Xv2-PRI_genomic.fna",
		expand("2_alignedData/bam/{SAMPLE}_sorted.bam", SAMPLE=SAMPLES),
		"3_snps/gatk/VT211_HG2W3DSXX_GATTCTGC-GACGAGAG_L001_gatk_calls.vcf",
		"3_snps/freebayes/VT211_HG2W3DSXX_GATTCTGC-GACGAGAG_L001_freebayes_calls.vcf",
		"2_alignedData/bam/VT211_HG2W3DSXX_GATTCTGC-GACGAGAG_L001_namesort.bam",
		"2_alignedData/bam/VT211_HG2W3DSXX_GATTCTGC-GACGAGAG_L001_fixmate.bam",
		"2_alignedData/bam/VT211_HG2W3DSXX_GATTCTGC-GACGAGAG_L001_coord.bam",
		"2_alignedData/bam/VT211_HG2W3DSXX_GATTCTGC-GACGAGAG_L001_markdup.bam",
		"3_snps/freebayes/Q20_DP5_freebayes_variants.recode.vcf",
		"5_subset/pigmentation.vcf",
		"4_subset/lipid.vcf",
		"4_subset/growth.vcf",
		"reference/GCF_900518725.1_TS10Xv2-PRI_genomic.dict",
		"2_alignedData/bam/VT211_HG2W3DSXX_GATTCTGC-GACGAGAG_L001_final.bam",
#		"2_alignedData/bam/VT211_HG2W3DSXX_GATTCTGC-GACGAGAG_L001_temp.bam",
		expand("targets/{TARGETS}_promoter.bed", TARGETS=['pigmentation', 'growth', 'lipid']), 
		expand("3_snps/bcftools/VT211_HG2W3DSXX_GATTCTGC-GACGAGAG_L001_bcftools_calls_{NUMBER}.bcf", NUMBER=NUMBERS),
		"3_snps/bcftools/VT211_HG2W3DSXX_GATTCTGC-GACGAGAG_L001_bcftools_filtered.bcf",
		expand("3_snps/gatk/VT211_HG2W3DSXX_GATTCTGC-GACGAGAG_L001_gatk_calls_{NUMBER}.vcf", NUMBER=NUMBERS),
		"2_alignedData/bam/dummy.txt",
		"reference/chromosome.list",
		expand("reference/chromosome.list_{NUMBER}", NUMBER=NUMBERS),	
		expand("3_snps/gatk/filtered/VT211_HG2W3DSXX_GATTCTGC-GACGAGAG_L001_gatk_filtered_{NUMBER}.recode.vcf", NUMBER=NUMBERS),
		"3_snps/gatk/Q20_DP5_gatk_variants.recode.vcf",
		expand("4_subset/{tool}_{prefix}.vcf", tool=['gatk','freebayes'], prefix=['growthg', 'pigmentationg','lipidg', 'growthp', 'pigmentationp', 'lipidp']),

################
# Rules Proper #
################

rule fastqc_raw:
	input:
		"0_rawData/fastq/{prefix}.fastq.gz",
	output: 
		zip  = "0_rawData/FastQC/{prefix}_fastqc.zip",
		html = "0_rawData/FastQC/{prefix}_fastqc.html",
	resources:
		time_mins = 90,
	conda:
		"envs/default.yml"
	threads:
		MAX_THREADS
	shell:
		"""
		module load FastQC/0.11.7
		fastqc --threads {threads} -o 0_rawData/FastQC/ {input}
		"""

# poly-G artefacts from NovaSeq require removal
rule cutadapt:
	input:
		r1 = "0_rawData/fastq/VT211_HG2W3DSXX_GATTCTGC-GACGAGAG_L001_R1.fastq.gz",
		r2 = "0_rawData/fastq/VT211_HG2W3DSXX_GATTCTGC-GACGAGAG_L001_R2.fastq.gz",
	output:
		r1 = "1_trimmedData/fastq/VT211_HG2W3DSXX_GATTCTGC-GACGAGAG_L001_R1.fastq.gz",
		r2 = "1_trimmedData/fastq/VT211_HG2W3DSXX_GATTCTGC-GACGAGAG_L001_R2.fastq.gz",
		qc = "1_trimmedData/FastQC/VT211_HG2W3DSXX_GATTCTGC-GACGAGAG_L001_cutadapt.txt",
	shell:
		"""
		module load cutadapt/1.14-foss-2016b-Python-2.7.13
		cutadapt --nextseq-trim=20 --minimum-length 50 -o {output.r1} -p {output.r2} {input.r1} {input.r2} > {output.qc}
		"""

rule fastqc_trim:
	input:
		"1_trimmedData/fastq/{prefix}.fastq.gz",
	output:
		zip = "1_trimmedData/FastQC/{prefix}_fastqc.zip",
		html = "1_trimmedData/FastQC/{prefix}_fastqc.html",
	shell:
		"""
		module load FastQC/0.11.7
		fastqc --threads {threads} -o 1_trimmedData/FastQC/ {input}
		"""

rule bwa_index:
	input:
		"reference/{ref}",
	output:
		"reference/{ref}.amb",
		"reference/{ref}.ann",
		"reference/{ref}.bwt",
		"reference/{ref}.pac",
		"reference/{ref}.sa",
	conda:
		"envs/default.yml"
	shell:
		"""
		module load BWA/0.7.15-foss-2017a
		bwa index {input}
		"""


rule bwa_mem:
	input:
		index = expand("reference/GCF_900518725.1_TS10Xv2-PRI_genomic.fna.gz.{ext}", ext=["amb","ann","bwt","pac","sa"]),
		r1 = "1_trimmedData/fastq/VT211_HG2W3DSXX_GATTCTGC-GACGAGAG_L001_R1.fastq.gz",
		r2 = "1_trimmedData/fastq/VT211_HG2W3DSXX_GATTCTGC-GACGAGAG_L001_R2.fastq.gz",
	output:
                "2_alignedData/bam/VT211_HG2W3DSXX_GATTCTGC-GACGAGAG_L001.bam",
	log:
		"logs/bwa_mem/VT.log",
	threads:
		MAX_THREADS
	conda:
                "envs/default.yaml",
	shell:
                """
              	bwa mem -t {threads} reference/GCF_900518725.1_TS10Xv2-PRI_genomic.fna.gz {input.r1} {input.r2} | samtools view -b -@ {threads} -o {output} -  
                """

rule name_sort:
	input:
		"2_alignedData/bam/VT211_HG2W3DSXX_GATTCTGC-GACGAGAG_L001.bam",
	output:
		"2_alignedData/bam/VT211_HG2W3DSXX_GATTCTGC-GACGAGAG_L001_namesort.bam",
	threads:
		MAX_THREADS
	conda:
		"envs/default.yaml",
	shell:
		"""
		samtools sort -n -@ {threads} -o {output} {input}
		"""

rule fixmate:
	input:
		"2_alignedData/bam/VT211_HG2W3DSXX_GATTCTGC-GACGAGAG_L001_namesort.bam",
	output:
		"2_alignedData/bam/VT211_HG2W3DSXX_GATTCTGC-GACGAGAG_L001_fixmate.bam",
	threads:
		MAX_THREADS
	conda:
		"envs/default.yaml",
	shell:
		"""
		samtools fixmate -m {input} {output}
		"""

rule position_sort:
	input: 
		"2_alignedData/bam/VT211_HG2W3DSXX_GATTCTGC-GACGAGAG_L001_fixmate.bam",
	output:
		"2_alignedData/bam/VT211_HG2W3DSXX_GATTCTGC-GACGAGAG_L001_coord.bam",
	threads:
		MAX_THREADS
	conda:
		"envs/default.yaml",
	shell:
		"""
		samtools sort -@ {threads} -o {output} {input}
		"""

rule markduplicates:
	input:
		"2_alignedData/bam/VT211_HG2W3DSXX_GATTCTGC-GACGAGAG_L001_coord.bam",
	output:
		"2_alignedData/bam/VT211_HG2W3DSXX_GATTCTGC-GACGAGAG_L001_markdup.bam",
	threads:
		MAX_THREADS
	conda:
		"envs/default.yaml",
	shell:
		"""
		samtools markdup {input} {output}
		"""

rule split_bam:
	input:
		"2_alignedData/bam/VT211_HG2W3DSXX_GATTCTGC-GACGAGAG_L001_sorted.bam",
	output:
		"2_alignedData/bam/dummy.txt",
	log:
		"logs/split_bam/VT.log",
	conda:
		"envs/bamtools.yaml",
	shell:
		"""
		(bamtools split -in {input} -reference) 2> {log}
		touch {output}
		"""

#rule bcftools_mpileup:
#	input:
#		ref = "reference/GCF_900518725.1_TS10Xv2-PRI_genomic.fna",
#		bam = "2_alignedData/bam/VT211_HG2W3DSXX_GATTCTGC-GACGAGAG_L001_final.bam",
#
#	output:
#		"reference/GCF_900518725.1_TS10Xv2-PRI_genomic.mpileup",
#	log:
#		"logs/bcftools_mpileup/mpileup.log",
#	threads:
#		1
#	conda:
#		"envs/bcftools.yaml",
#	shell:
#		"""
#		(bcftools mpileup -Ou -f {input.ref} {input.bam}) 2> {log}
#		"""

rule bcftools_call:
	input:
		bam = "2_alignedData/bam/VT211_HG2W3DSXX_GATTCTGC-GACGAGAG_L001_final.bam",
		ref = "reference/GCF_900518725.1_TS10Xv2-PRI_genomic.fna",
		chromosomes = "reference/bcf_chromosome_lists/chromosome_{number}.list",
	output:
		"3_snps/bcftools/VT211_HG2W3DSXX_GATTCTGC-GACGAGAG_L001_bcftools_calls_{number}.bcf",
	log:
		"logs/bcftools_snp_call/{number}.log",
	threads:
		MAX_THREADS
	conda:
		"envs/bcftools.yaml",
	shell:
		"""
		(bcftools mpileup -Ou -R {input.chromosomes} -f {input.ref} {input.bam} | bcftools call -mv --variants-only -Ob -o {output}) 2> {log}
		"""

rule bcftools_filter:
	input:
 		"3_snps/bcftools/VT211_HG2W3DSXX_GATTCTGC-GACGAGAG_L001_bcftools_calls.vcf",
	output:
		"3_snps/bcftools/VT211_HG2W3DSXX_GATTCTGC-GACGAGAG_L001_bcftools_filtered.recode.vcf",
	log:
		"logs/bcftools_filter/VT.log",
	threads:
		MAX_THREADS
	conda:
		"envs/default.yaml",
	shell:
		"""
		(vcftools --vcf {input} --minQ 20 --recode --recode-INFO-all --minDP 5 --out 3_snps/bcftools/VT211_HG2W3DSXX_GATTCTGC-GACGAGAG_L001_bcftools_filtered) 2> {log}
		"""

rule gunzip_ref:
	input:
		"reference/GCF_900518725.1_TS10Xv2-PRI_genomic.fna.gz",
	output:
		"reference/GCF_900518725.1_TS10Xv2-PRI_genomic.fna",
	threads:
		1
	conda:
		"envs/default.yaml",
	shell:
		"""
		gunzip -c {input} > {output}
		"""

rule freebayes_call:
	input:
		bam = "2_alignedData/bam/VT211_HG2W3DSXX_GATTCTGC-GACGAGAG_L001_markdup.bam",
		ref = "reference/GCF_900518725.1_TS10Xv2-PRI_genomic.fna",
	output:
		"3_snps/freebayes/VT211_HG2W3DSXX_GATTCTGC-GACGAGAG_L001_freebayes_calls.vcf",
	log:
		"logs/freebayes_call/VT.log",
	threads:
		1
	conda:
		"envs/default.yaml",
	shell:
		"""
		(freebayes -f {input.ref} {input.bam} > {output}) 2> {log}
		"""

rule freebayes_filter:
	input:
		"3_snps/freebayes/VT211_HG2W3DSXX_GATTCTGC-GACGAGAG_L001_freebayes_calls.vcf",
	output:
		"3_snps/freebayes/Q20_DP5_freebayes_variants.recode.vcf",
	log:
		"logs/freebayes_filter/VT.log",
	threads:
		1
	conda:
		"envs/default.yaml",
	shell:
		"""
		(vcftools --vcf {input} --minQ 20 --recode --recode-INFO-all --minDP 5 --out 3_snps/freebayes/Q20_DP5_freebayes_variants) 2> {log}
		"""

rule picard_dict:
	input:
		"reference/GCF_900518725.1_TS10Xv2-PRI_genomic.fna",
	output:
		"reference/GCF_900518725.1_TS10Xv2-PRI_genomic.dict",
	log:
		"logs/picard_dict/GCF.log",
	threads:
		1
	shell:
		"""
		module load picard/2.2.4-Java-1.8.0_121
		(java -jar $EBROOTPICARD/picard.jar CreateSequenceDictionary R= {input} O= {output}) 2> {log}
		"""

rule fix_bam_reads:
	input:
		"2_alignedData/bam/VT211_HG2W3DSXX_GATTCTGC-GACGAGAG_L001_markdup.bam",
	output:
		"2_alignedData/bam/VT211_HG2W3DSXX_GATTCTGC-GACGAGAG_L001_final.bam",
	threads:
		1
	shell:
		"""
		module load picard/2.2.4-Java-1.8.0_121
		java -jar $EBROOTPICARD/picard.jar AddOrReplaceReadGroups I= {input} O= {output} \
			RGID= 1 \
			RGLB= lib1 \
			RGPL=illumina \
			RGPU=unit1 \
			RGSM=sample1 \
			CREATE_INDEX=True
		"""
#extracts all contig names from GCF file and creates files with 500 contigs per file. This allows gatk variant calling to be parallelised.
#Doesn't work.
 
rule makeChromosomeList:
	input:
		"reference/GCF_900518725.1_TS10Xv2-PRI_genomic.fna",
	output:
		"reference/chromosome.list"
	threads:
		1
	shell:
		"""
		set +u; bash chromosome_split.sh; set -u;
#		mkdir -p reference/chromosome_lists
#		grep '^>' {input} | cut -d " " -f1 | cut -c2- > {output}
#		bash chromosome_split.sh

#		split -d -l 500 {output} reference/chromosome_lists/chromosome_
#		rename -v -- 's/$/.list/' reference/chromosome_lists/chromosome_*
	
#		count=1
#
#		for run in {1..66}
#		do
#			sed -n ${run}p {output} > /reference/chromosome_lists/chromosome_${count}.list
#			((count=count+1))
#		done
#
#		for run in {67..52066..500}
#		do
#		end=$((run+499))
#			sed -n "${run},${end}p" {output} > reference/chromosome_lists/chromosome_${count}.list
#			((count=count+1))
#		done
		"""

#call gatk with subsets of contigs to parallelise process. 
rule gatk_call:
	input:
		bam = "2_alignedData/bam/VT211_HG2W3DSXX_GATTCTGC-GACGAGAG_L001_final.bam",
		ref = "reference/GCF_900518725.1_TS10Xv2-PRI_genomic.fna",
		chromosome = "reference/chromosome_lists/chromosome_{number}.list",
		vcf = "3_snps/gatk/VT211_HG2W3DSXX_GATTCTGC-GACGAGAG_L001_gatk_calls_{number}.vcf",
	log:
		"logs/gatk_call/VT211_HG2W3DSXX_GATTCTGC-GACGAGAG_L001_{number}.log",
	conda:
		"envs/gatk.yaml",
	threads:
		1
	shell:
		"""
		(gatk HaplotypeCaller -L {input.chromosome} -R {input.ref} -I {input.bam} -O {output.vcf}) 3> {log}

		"""

#rule picard_gathervcfs:
#	input:
#		expand("3_snps/gatk/VT211_HG2W3DSXX_GATTCTGC-GACGAGAG_L001_gatk_calls_{number}.g.vcf", number=['531', '532', '533']),
#	output:
#		"3_snps/gatk/VT211_HG2W3DSXX_GATTCTGC-GACGAGAG_L001_gatk_calls_total.g.vcf",
#	log:
#		"logs/picard_gathervcfs/VT_final.log",
#	params:
#		files = lambda wildcards, input: " I= ".join(input)
#	threads:
#		MAX_THREADS
#	shell:
#		"""
#		module load picard/2.2.4-Java-1.8.0_121
#		(java -jar $EBROOTPICARD/picard.jar GatherVcfs I= {params.files} O= {output}) 2> {log}
#		"""

rule picard_mergevcfs:
	input:
		expand("3_snps/gatk/filtered/VT211_HG2W3DSXX_GATTCTGC-GACGAGAG_L001_gatk_filtered_{number}.recode.vcf", number=['1', '2', '3']),

	output:
		"3_snps/gatk/filtered/VT211_HG2W3DSXX_GATTCTGC-GACGAGAG_L001_gatk_calls_total.vcf.gz",
	log:
		"logs/picard_gathervcfs/VT_total.log",
	params:
		files = lambda wildcards, input: " I= ".join(input)
	threads:
		MAX_THREADS
	shell:
		"""
		module load picard/2.2.4-Java-1.8.0_121
		(java -jar $EBROOTPICARD/picard.jar MergeVcfs I= {params.files} O= {output}) 2> {log}
		"""

rule gatk_filter:
	input:
		"3_snps/gatk/VT211_HG2W3DSXX_GATTCTGC-GACGAGAG_L001_gatk_calls_{number}.vcf",
	output:
		"3_snps/gatk/filtered/VT211_HG2W3DSXX_GATTCTGC-GACGAGAG_L001_gatk_filtered_{number}.recode.vcf",
	log:
		"logs/gatk/haplotypecaller/VT211_HG2W3DSXX_GATTCTGC-GACGAGAG_L001_{number}.log",
	threads:
		1
	conda:
		"envs/default.yaml",
	shell:
		"""
		(vcftools --vcf {input} --minQ 20 --recode --recode-INFO-all --minDP 5 --out 3_snps/gatk/filtered/VT211_HG2W3DSXX_GATTCTGC-GACGAGAG_L001_gatk_filtered_{wildcards.number}) 2> {log} 
		"""

rule join_gatk_vcfs:
	input:
		expand("3_snps/gatk/filtered/VT211_HG2W3DSXX_GATTCTGC-GACGAGAG_L001_gatk_filtered_{number}.recode.vcf", number=NUMBERS),
	output:
		"3_snps/gatk/Q20_DP5_gatk_variants.recode.vcf",
	log:
		"logs/join_gatk_vcfs/final.log",
	threads:
		1
	shell:
		"""
		sed -n '1,52440p' 3_snps/gatk/filtered/VT211_HG2W3DSXX_GATTCTGC-GACGAGAG_L001_gatk_filtered_1.recode.vcf >> {output}
        
		for i in {{1..533}}; do sed -n '52441,$p' 3_snps/gatk/filtered/VT211_HG2W3DSXX_GATTCTGC-GACGAGAG_L001_gatk_filtered_${{i}}.recode.vcf >> {output}; done
        
		"""


#======================================================================
#Make .bed files for promoter regions for each of the three gene sets
#======================================================================
rule promoter_regions:
	input:
		"targets/{prefix}_genes.bed",
	output:
		"targets/{prefix}_promoter.bed",
	threads:
		1
	shell:
		"""
		awk '{if($2 > 15000) {$3=$2; $2=$2-15000; print;} else {$3=$2; $2=1; print;}}' {input} > {output}
		"""

#===========================================
#intersect .bed files with vcf files
#===========================================
rule bedtools:
	input:
		bed = expand("targets/{{prefix}}.bed", prefix=['pigmentationg','lipidg','growthg','pigmentationp','lipidp','growthp']),
		vcf = expand("3_snps/{{tool}}/Q20_DP5_{{tool}}_variants.recode.vcf", tool=['freebayes','gatk']),
	output:
		"4_subset/{tool}_{prefix}.vcf",
	threads:
		1
	conda:
		"envs/default.yaml",
	log:
		"logs/bedtools/{tool}_{prefix}.log",
	shell:
		"""
		(bedtools intersect -header -a {input.vcf} -b {input.bed} -wa > {output}) 2> {log}
		"""
