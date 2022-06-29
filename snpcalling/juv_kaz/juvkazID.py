#2022-2-2
#HU Xuejing
#snp calling pipeline of juv_kaz trimmed sample use L.sin as outgroup species
##################


#configfile:"config.yaml"


SAMPLES = ["24","23","22","21","20","18","17","15","14","13"]


rule all:
    input:
        "/proj/uppstore2017185/b2014034/nobackup/Xuejing/snpcalling/juv_kaz/ancestor/anc/juv_kaz_ancestor_snp_all.vcf"
rule bwa_map:
    input:
        "/proj/uppstore2017185/b2014034/nobackup/Xuejing/snpcalling/ref/L_sin.fasta",
        "/proj/uppstore2017185/b2014034/nobackup/Xuejing/trimmed_data/juv_kaz/juv_kaz_sample_{sample}_R1_trimmed.fastq.gz",
        "/proj/uppstore2017185/b2014034/nobackup/Xuejing/trimmed_data/juv_kaz/juv_kaz_sample_{sample}_R2_trimmed.fastq.gz"
    params:
        rg=r"@RG\tID:juv_kaz{sample}\tSM:{sample}\tPL:ILLUMINA\tLB:1"
    output:
        "/proj/uppstore2017185/b2014034/nobackup/Xuejing/snpcalling/juv_kaz/ancestor/mapped_reads/juv_kaz_sample_{sample}_trimmed.bam"
    threads:10
    conda:
        "/proj/uppstore2017185/b2014034/nobackup/Xuejing/snpcalling/environment.yaml"
    shell:
        "bwa mem -R '{params.rg}' -t {threads} -M {input} | samtools view -Sb - > {output}"


rule samtools_sort:
    input:
        "/proj/uppstore2017185/b2014034/nobackup/Xuejing/snpcalling/juv_kaz/ancestor/mapped_reads/juv_kaz_sample_{sample}_trimmed.bam"
    output:
        "/proj/uppstore2017185/b2014034/nobackup/Xuejing/snpcalling/juv_kaz/ancestor/sorted_reads/juv_kaz_sample_{sample}_trimmed_sorted.bam"
    threads:8
    conda:
        "/proj/uppstore2017185/b2014034/nobackup/Xuejing/snpcalling/environment.yaml"
    shell:
        "samtools sort -T sorted_reads/{wildcards.sample} "
        "-O bam {input} > {output}"


rule mark_duplicates:
    input:
        "/proj/uppstore2017185/b2014034/nobackup/Xuejing/snpcalling/juv_kaz/ancestor/sorted_reads/juv_kaz_sample_{sample}_trimmed_sorted.bam"
    output:
        bam="/proj/uppstore2017185/b2014034/nobackup/Xuejing/snpcalling/juv_kaz/ancestor/sorted_reads/juv_kaz_sample_{sample}_trimmed_sorted_dedup.bam",
        metrics="/proj/uppstore2017185/b2014034/nobackup/Xuejing/snpcalling/juv_kaz/ancestor/sorted_reads/juv_kaz_sample_{sample}_trimmed_sorted_dedup.metrics.txt"
    threads:4
    resources:
        mem_gb=64
    conda:
        "/proj/uppstore2017185/b2014034/nobackup/Xuejing/snpcalling/environment.yaml"
    shell:
        "module load java/sun_jdk1.8.0_151 bioinfo-tools picard/2.23.4 | "
        "java -Xmx32g -jar /sw/bioinfo/picard/2.23.4/rackham/picard.jar MarkDuplicates "
        "INPUT={input} OUTPUT={output.bam} METRICS_FILE={output.metrics} "
        "REMOVE_DUPLICATES=true CREATE_INDEX=true  ASSUME_SORTED=true"

#samtools merge

rule samtools_index:
    input:
        "/proj/uppstore2017185/b2014034/nobackup/Xuejing/snpcalling/juv_kaz/ancestor/sorted_reads/juv_kaz_sample_{sample}_trimmed_sorted_dedup.bam"
    output:
        "/proj/uppstore2017185/b2014034/nobackup/Xuejing/snpcalling/juv_kaz/ancestor/sorted_reads/juv_kaz_sample_{sample}_trimmed_sorted_dedup.bam.bai"
    threads:4
    conda:
        "/proj/uppstore2017185/b2014034/nobackup/Xuejing/snpcalling/environment.yaml"
    shell:
        "samtools index {input}"

#indel realignment
#1.RealignerTargetCreator
rule realignertargetcreator:
    input:
        bam="/proj/uppstore2017185/b2014034/nobackup/Xuejing/snpcalling/juv_kaz/ancestor/sorted_reads/juv_kaz_sample_{sample}_trimmed_sorted_dedup.bam",
        ref="/proj/uppstore2017185/b2014034/nobackup/Xuejing/snpcalling/ref/L_sin.fasta"
    output:
        intervals="/proj/uppstore2017185/b2014034/nobackup/Xuejing/snpcalling/juv_kaz/ancestor/sorted_reads/juv_kaz_sample_{sample}_trimmed_sorted_dedup.intervals"
    conda:
        "/proj/uppstore2017185/b2014034/nobackup/Xuejing/snpcalling/environment.yaml"
    resources:
        mem_gb=64
    threads: 10
    shell:
        "module load bioinfo-tools GATK/3.8-0 | "
        "java -Xmx32g -jar /sw/apps/bioinfo/GATK/3.8-0/GenomeAnalysisTK.jar \
            -I {input.bam} \
            -R {input.ref}  \
            -T RealignerTargetCreator  \
            -o {output.intervals}"

#2.IndelRealigner
rule indelrealigner:
    input:
        bam="/proj/uppstore2017185/b2014034/nobackup/Xuejing/snpcalling/juv_kaz/ancestor/sorted_reads/juv_kaz_sample_{sample}_trimmed_sorted_dedup.bam",
        ref="/proj/uppstore2017185/b2014034/nobackup/Xuejing/snpcalling/ref/L_sin.fasta",
        target_intervals="/proj/uppstore2017185/b2014034/nobackup/Xuejing/snpcalling/juv_kaz/ancestor/sorted_reads/juv_kaz_sample_{sample}_trimmed_sorted_dedup.intervals"
    output:
        "/proj/uppstore2017185/b2014034/nobackup/Xuejing/snpcalling/juv_kaz/ancestor/realigned/juv_kaz_sample_{sample}_trimmed_sorted_dedup_realign.bam"
    conda:
        "/proj/uppstore2017185/b2014034/nobackup/Xuejing/snpcalling/environment.yaml"
    log:
        "/proj/uppstore2017185/b2014034/nobackup/Xuejing/snpcalling/juv_kaz/ancestor/logs/gatk3/indelrealigner/indelrealigner_{sample}.log"
    params:
        extra=""  # optional
    threads: 10
    resources:
        mem_gb = 64
    shell:
        "module load bioinfo-tools GATK/3.8-0 | "
        "java -jar /sw/apps/bioinfo/GATK/3.8-0/GenomeAnalysisTK.jar \
          -I {input.bam} \
          -R {input.ref} \
          -T IndelRealigner --filter_bases_not_stored \
          -o {output} \
          -targetIntervals {input.target_intervals}"



rule samtools_re_index:
    input:
        "/proj/uppstore2017185/b2014034/nobackup/Xuejing/snpcalling/juv_kaz/ancestor/realigned/juv_kaz_sample_{sample}_trimmed_sorted_dedup_realign.bam"
    output:
        "/proj/uppstore2017185/b2014034/nobackup/Xuejing/snpcalling/juv_kaz/ancestor/realigned/juv_kaz_sample_{sample}_trimmed_sorted_dedup_realign.bam.bai"
    threads:4
    conda:
        "/proj/uppstore2017185/b2014034/nobackup/Xuejing/snpcalling/environment.yaml"
    shell:
        "samtools index {input}"


rule bcftools_call:
    input:
        fa="/proj/uppstore2017185/b2014034/nobackup/Xuejing/snpcalling/ref/L_sin.fasta",
        bam=expand("/proj/uppstore2017185/b2014034/nobackup/Xuejing/snpcalling/juv_kaz/ancestor/realigned/juv_kaz_sample_{sample}_trimmed_sorted_dedup_realign.bam", sample=SAMPLES),
        bai=expand("/proj/uppstore2017185/b2014034/nobackup/Xuejing/snpcalling/juv_kaz/ancestor/realigned/juv_kaz_sample_{sample}_trimmed_sorted_dedup_realign.bam.bai", sample=SAMPLES)
    output:
        "/proj/uppstore2017185/b2014034/nobackup/Xuejing/snpcalling/juv_kaz/ancestor/anc/temp.bcf"
    threads:20
    resources:
        mem_gb = 64
    conda:
        "/proj/uppstore2017185/b2014034/nobackup/Xuejing/snpcalling/environment.yaml"
    shell:
        "bcftools mpileup --skip-indels -Ou -f {input.fa} {input.bam} | "
        "bcftools call -mv - > {output}"


rule trans_call:
    input:
        "/proj/uppstore2017185/b2014034/nobackup/Xuejing/snpcalling/juv_kaz/ancestor/anc/temp.bcf"
    output:
    	"/proj/uppstore2017185/b2014034/nobackup/Xuejing/snpcalling/juv_kaz/ancestor/anc/temp_calls.vcf"
    threads:20
    conda:
        "/proj/uppstore2017185/b2014034/nobackup/Xuejing/snpcalling/environment.yaml"
    shell:
        "bcftools view -i '%QUAL>=40' {input} > {output}"



rule vcftools_call:
    input:
        "/proj/uppstore2017185/b2014034/nobackup/Xuejing/snpcalling/juv_kaz/ancestor/anc/temp_calls.vcf"
    output:
        dbsnp="/proj/uppstore2017185/b2014034/nobackup/Xuejing/snpcalling/juv_kaz/ancestor/anc/dbSNP_calls.recode.vcf"
    threads:20
    resources:
        mem_gb = 64
    conda:
        "/proj/uppstore2017185/b2014034/nobackup/Xuejing/snpcalling/environment.yaml"
    shell:
        "vcftools --vcf {input} \
           --remove-indels \
           --exclude-bed /proj/uppstore2017185/b2014034_nobackup/Xuejing/samples/P14502_103.FINAL-deduped-nuc.filtered.fasta.out.gff.bed \
           --recode \
	       --recode-INFO-all \
           --out /proj/uppstore2017185/b2014034/nobackup/Xuejing/snpcalling/juv_kaz/ancestor/anc/dbSNP_calls"


#3.BaseRecalibrator(skip and redo if don't have dbSNP)
rule baserecalibrator:
    input:
        bam="/proj/uppstore2017185/b2014034/nobackup/Xuejing/snpcalling/juv_kaz/ancestor/realigned/juv_kaz_sample_{sample}_trimmed_sorted_dedup_realign.bam",
        ref="/proj/uppstore2017185/b2014034/nobackup/Xuejing/snpcalling/ref/L_sin.fasta",
        known="/proj/uppstore2017185/b2014034/nobackup/Xuejing/snpcalling/juv_kaz/ancestor/anc/dbSNP_calls.recode.vcf"
    output:
        "/proj/uppstore2017185/b2014034/nobackup/Xuejing/snpcalling/juv_kaz/ancestor/realigned/juv_kaz_sample_{sample}_trimmed_sorted_dedup_realign_calibration.csv"
    conda:
        "/proj/uppstore2017185/b2014034/nobackup/Xuejing/snpcalling/environment.yaml"
    resources:
        mem_gb = 32
    threads: 16
    shell:
        "module load bioinfo-tools GATK/3.8-0 | "
        "java -Xmx30g -jar /sw/apps/bioinfo/GATK/3.8-0/GenomeAnalysisTK.jar \
          -I {input.bam} \
          -R {input.ref} \
          -T BaseRecalibrator  \
          -o {output} \
          -knownSites {input.known}"




##4.PrintReads(optional if don't have dbSNP)
rule printreads:
    input:
        bam="/proj/uppstore2017185/b2014034/nobackup/Xuejing/snpcalling/juv_kaz/ancestor/realigned/juv_kaz_sample_{sample}_trimmed_sorted_dedup_realign.bam",
        ref="/proj/uppstore2017185/b2014034/nobackup/Xuejing/snpcalling/ref/L_sin.fasta",
        recal_data="/proj/uppstore2017185/b2014034/nobackup/Xuejing/snpcalling/juv_kaz/ancestor/realigned/juv_kaz_sample_{sample}_trimmed_sorted_dedup_realign_calibration.csv"
    output:
        "/proj/uppstore2017185/b2014034/nobackup/Xuejing/snpcalling/juv_kaz/ancestor/alignment/juv_kaz_sample_{sample}_trimmed_sorted_dedup_realign_calibration_bqsr.bam"
    conda:
        "/proj/uppstore2017185/b2014034/nobackup/Xuejing/snpcalling/environment.yaml"
    params:
        extra=""  # optional
    threads: 16
    shell:
        "module load bioinfo-tools GATK/3.8-0 | "
        "java -Xmx30g -jar /sw/apps/bioinfo/GATK/3.8-0/GenomeAnalysisTK.jar \
          -I {input.bam} \
          -R {input.ref} \
          -T PrintReads  \
          -o {output} \
          -BQSR {input.recal_data}"


rule samtools_flagstat:
    input:
        "/proj/uppstore2017185/b2014034/nobackup/Xuejing/snpcalling/juv_kaz/ancestor/alignment/juv_kaz_sample_{sample}_trimmed_sorted_dedup_realign_calibration_bqsr.bam"
    output:
        "/proj/uppstore2017185/b2014034/nobackup/Xuejing/snpcalling/juv_kaz/ancestor/alignment/juv_kaz_sample_{sample}_trimmed_sorted_dedup_realign_calibration_bqsr.flagstat"
    conda:
        "/proj/uppstore2017185/b2014034/nobackup/Xuejing/snpcalling/environment.yaml"
    shell:
        "samtools flagstat {input} > {output}"


rule haplotype_caller:
    input:
        # juvgle or list of bam files
        bam="/proj/uppstore2017185/b2014034/nobackup/Xuejing/snpcalling/juv_kaz/ancestor/alignment/juv_kaz_sample_{sample}_trimmed_sorted_dedup_realign_calibration_bqsr.bam",
        ref="/proj/uppstore2017185/b2014034/nobackup/Xuejing/snpcalling/ref/L_sin.fasta",
        known="/proj/uppstore2017185/b2014034/nobackup/Xuejing/snpcalling/juv_kaz/ancestor/anc/dbSNP_calls.recode.vcf"
    output:
       "/proj/uppstore2017185/b2014034/nobackup/Xuejing/snpcalling/juv_kaz/ancestor/alignment/juv_kaz_sample_{sample}_trimmed_sorted_dedup_realign_calibration_bqsr.bam.g.vcf"
    conda:
        "/proj/uppstore2017185/b2014034/nobackup/Xuejing/snpcalling/environment.yaml"
    threads:
        10
    resources:
        mem_gb = 64
    shell:
        "module load bioinfo-tools GATK/3.8-0 | "
    	"java -Xmx30g -jar /sw/apps/bioinfo/GATK/3.8-0/GenomeAnalysisTK.jar \
            -T HaplotypeCaller \
            -R {input.ref} \
            -I {input.bam} \
            --emitRefConfidence GVCF \
            -o {output}"




rule CombineGVCFs:
    input:
        gvcfs=expand("/proj/uppstore2017185/b2014034/nobackup/Xuejing/snpcalling/juv_kaz/ancestor/alignment/juv_kaz_sample_{sample}_trimmed_sorted_dedup_realign_calibration_bqsr.bam.g.vcf", sample=SAMPLES),
        ref="/proj/uppstore2017185/b2014034/nobackup/Xuejing/snpcalling/ref/L_sin.fasta"
    output:
        "/proj/uppstore2017185/b2014034/nobackup/Xuejing/snpcalling/juv_kaz/ancestor/alignment/juv_kaz_trimmed_sorted_dedup_realign_calibration_bqsr.bam.cohort.g.vcf"
    conda:
        "/proj/uppstore2017185/b2014034/nobackup/Xuejing/snpcalling/environment.yaml"
    script:
        "/proj/uppstore2017185/b2014034/nobackup/Xuejing/snpcalling/juv_kaz/scripts/CombineGVCFs.py"



rule GenotypeGVCFs:
    input:
        gvcf="/proj/uppstore2017185/b2014034/nobackup/Xuejing/snpcalling/juv_kaz/ancestor/alignment/juv_kaz_trimmed_sorted_dedup_realign_calibration_bqsr.bam.cohort.g.vcf",
        ref="/proj/uppstore2017185/b2014034/nobackup/Xuejing/snpcalling/ref/L_sin.fasta"
    output:
        "/proj/uppstore2017185/b2014034/nobackup/Xuejing/snpcalling/juv_kaz/ancestor/alignment/juv_kaz_trimmed_sorted_dedup_realign_calibration_bqsr.bam.g.vcf.HC.vcf"
    conda:
        "/proj/uppstore2017185/b2014034/nobackup/Xuejing/snpcalling/environment.yaml"
    shell:
        "module load bioinfo-tools GATK/3.8-0 | "
    	"java -Xmx30g -jar /sw/apps/bioinfo/GATK/3.8-0/GenomeAnalysisTK.jar \
            -T GenotypeGVCFs \
            -R {input.ref} \
            -V {input.gvcf} \
            -o {output}"


rule extract_SNPs:
    input:
        vcf="/proj/uppstore2017185/b2014034/nobackup/Xuejing/snpcalling/juv_kaz/ancestor/alignment/juv_kaz_trimmed_sorted_dedup_realign_calibration_bqsr.bam.g.vcf.HC.vcf",
        ref="/proj/uppstore2017185/b2014034/nobackup/Xuejing/snpcalling/ref/L_sin.fasta"
    output:
        "/proj/uppstore2017185/b2014034/nobackup/Xuejing/snpcalling/juv_kaz/ancestor/anc/juv_kaz_ancestor_snp_all.vcf"
    conda:
        "/proj/uppstore2017185/b2014034/nobackup/Xuejing/snpcalling/environment.yaml"
    threads:
        10
    shell:
        "module load bioinfo-tools GATK/3.8-0 | "
    	"java -Xmx128g -jar /sw/apps/bioinfo/GATK/3.8-0/GenomeAnalysisTK.jar \
            -T SelectVariants \
            -R {input.ref} \
            -V {input.vcf} \
            -selectType SNP \
            -o {output}"
