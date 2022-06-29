# a pipeline for fastqc L.juv.ire samples 
sample_names="Ire-juv-21C_R2_trimmed.fastq.gz  Ire-juv-41C_R2_trimmed.fastq.gz  Ire-juv-62C_R2_trimmed.fastq.gz  Ire-juv-22C_R1_trimmed.fastq.gz  Ire-juv-42C_R1_trimmed.fastq.gz  Ire-juv-81C_R1_trimmed.fastq.gz Ire-juv-22C_R2_trimmed.fastq.gz  Ire-juv-42C_R2_trimmed.fastq.gz  Ire-juv-81C_R2_trimmed.fastq.gz Ire-juv-1C_R1_trimmed.fastq.gz	   Ire-juv-2C_R1_trimmed.fastq.gz   Ire-juv-61C_R1_trimmed.fastq.gz  Ire-juv-82C_R1_trimmed.fastq.gz Ire-juv-1C_R2_trimmed.fastq.gz	   Ire-juv-2C_R2_trimmed.fastq.gz   Ire-juv-61C_R2_trimmed.fastq.gz  Ire-juv-82C_R2_trimmed.fastq.gz Ire-juv-21C_R1_trimmed.fastq.gz    Ire-juv-41C_R1_trimmed.fastq.gz  Ire-juv-62C_R1_trimmed.fastq.gz "
# this loop writes the bash script
for i in $sample_names; do
echo "the sample name is: "$i

echo "#!/bin/bash -l" > $i"_ire.sh"
echo "#SBATCH -A snic2021-5-20" >> $i"_ire.sh"
echo "#SBATCH -p core -n 10" >> $i"_ire.sh"
echo "#SBATCH -J fasqc_$i" >> $i"_ire.sh"
echo "#SBATCH -t 00:10:00" >> $i"_ire.sh"
echo "module load bioinfo-tools" >> $i"_ire.sh"
echo "module load FastQC/0.11.9">> $i"_ire.sh"
echo "fastqc -o ./fastqc $i ">> $i"_ire.sh"


done


# this runs/submits the bash script
for i in $sample_names; do
sbatch $i"_ire.sh"
done




