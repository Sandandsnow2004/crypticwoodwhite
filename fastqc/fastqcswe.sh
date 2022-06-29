# a pipeline for fastqc L.juv.swe samples 
sample_names=" S20_L001_R1_trimmed.fastq.gz  S24_L001_R2_trimmed.fastq.gz  S26_L003_R1_trimmed.fastq.gz	S28_L003_R2_trimmed.fastq.gz  S86_L001_R1_trimmed.fastq.gz
S17_L001_R1_trimmed.fastq.gz  S20_L001_R2_trimmed.fastq.gz  S24_L003_R1_trimmed.fastq.gz  S26_L003_R2_trimmed.fastq.gz	S32_L001_R1_trimmed.fastq.gz  S86_L001_R2_trimmed.fastq.gz
S17_L001_R2_trimmed.fastq.gz  S20_L003_R1_trimmed.fastq.gz  S24_L003_R2_trimmed.fastq.gz  S27_L001_R1_trimmed.fastq.gz	S32_L001_R2_trimmed.fastq.gz  S86_L003_R1_trimmed.fastq.gz
S17_L003_R1_trimmed.fastq.gz  S20_L003_R2_trimmed.fastq.gz  S25_L001_R1_trimmed.fastq.gz  S27_L001_R2_trimmed.fastq.gz	S32_L003_R1_trimmed.fastq.gz  S86_L003_R2_trimmed.fastq.gz
S17_L003_R2_trimmed.fastq.gz  S21_L001_R1_trimmed.fastq.gz  S25_L001_R2_trimmed.fastq.gz  S27_L003_R1_trimmed.fastq.gz	S32_L003_R2_trimmed.fastq.gz  
S19_L001_R1_trimmed.fastq.gz  S21_L001_R2_trimmed.fastq.gz  S25_L003_R1_trimmed.fastq.gz  S27_L003_R2_trimmed.fastq.gz	S33_L001_R1_trimmed.fastq.gz
S19_L001_R2_trimmed.fastq.gz  S21_L003_R1_trimmed.fastq.gz  S25_L003_R2_trimmed.fastq.gz  S28_L001_R1_trimmed.fastq.gz	S33_L001_R2_trimmed.fastq.gz
S19_L003_R1_trimmed.fastq.gz  S21_L003_R2_trimmed.fastq.gz  S26_L001_R1_trimmed.fastq.gz  S28_L001_R2_trimmed.fastq.gz	S33_L003_R1_trimmed.fastq.gz
S19_L003_R2_trimmed.fastq.gz  S24_L001_R1_trimmed.fastq.gz  S26_L001_R2_trimmed.fastq.gz  S28_L003_R1_trimmed.fastq.gz	S33_L003_R2_trimmed.fastq.gz "
# this loop writes the bash script
for i in $sample_names; do
echo "the sample name is: "$i

echo "#!/bin/bash -l" > $i"_swe.sh"
echo "#SBATCH -A snic2021-5-20" >> $i"_swe.sh"
echo "#SBATCH -p core -n 10" >> $i"_swe.sh"
echo "#SBATCH -J fasqc_$i" >> $i"_swe.sh"
echo "#SBATCH -t 00:10:00" >> $i"_swe.sh"
echo "module load bioinfo-tools" >> $i"_swe.sh"
echo "module load FastQC/0.11.9">> $i"_swe.sh"
echo "fastqc -o ./fastqc $i ">> $i"_swe.sh"


done


# this runs/submits the bash script
for i in $sample_names; do
sbatch $i"_swe.sh"
done
