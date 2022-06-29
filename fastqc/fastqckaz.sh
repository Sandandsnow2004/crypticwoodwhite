# a pipeline for fastqc L.juv.kaz samples 
sample_names=" Sample_14_R2_trimmed.fastq.gz  Sample_17_R2_trimmed.fastq.gz  Sample_20_R2_trimmed.fastq.gz  Sample_22_R2_trimmed.fastq.gz  Sample_24_R2_trimmed.fastq.gz
Sample_13_R1_trimmed.fastq.gz  Sample_15_R1_trimmed.fastq.gz  Sample_18_R1_trimmed.fastq.gz  Sample_21_R1_trimmed.fastq.gz  Sample_23_R1_trimmed.fastq.gz 
Sample_13_R2_trimmed.fastq.gz  Sample_15_R2_trimmed.fastq.gz  Sample_18_R2_trimmed.fastq.gz  Sample_21_R2_trimmed.fastq.gz  Sample_23_R2_trimmed.fastq.gz
Sample_14_R1_trimmed.fastq.gz  Sample_17_R1_trimmed.fastq.gz  Sample_20_R1_trimmed.fastq.gz  Sample_22_R1_trimmed.fastq.gz  Sample_24_R1_trimmed.fastq.gz "
# this loop writes the bash script
for i in $sample_names; do
echo "the sample name is: "$i

echo "#!/bin/bash -l" > $i"_kaz.sh"
echo "#SBATCH -A snic2021-5-20" >> $i"_kaz.sh"
echo "#SBATCH -p core -n 10" >> $i"_kaz.sh"
echo "#SBATCH -J fasqc_$i" >> $i"_kaz.sh"
echo "#SBATCH -t 00:10:00" >> $i"_kaz.sh"
echo "module load bioinfo-tools" >> $i"_kaz.sh"
echo "module load FastQC/0.11.9">> $i"_kaz.sh"
echo "fastqc -o ./fastqc $i ">> $i"_kaz.sh"


done


# this runs/submits the bash script
for i in $sample_names; do
sbatch $i"_kaz.sh"
done

