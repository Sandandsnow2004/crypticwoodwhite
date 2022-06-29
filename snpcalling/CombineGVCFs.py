#CombineGVCFs.py
import os
from snakemake.shell import shell
gvcfs = list(map("-V {}".format, snakemake.input.gvcfs))

shell(
    "module load bioinfo-tools GATK/3.8-0 | "
    "java -Xmx128g -jar /sw/apps/bioinfo/GATK/3.8-0/GenomeAnalysisTK.jar "
    "-T CombineGVCFs "
    "-R {snakemake.input.ref} "
    "{gvcfs} "
    "-o {snakemake.output} "
)

