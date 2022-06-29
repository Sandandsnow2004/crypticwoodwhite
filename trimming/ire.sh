# a pipeline for L.juv.ire samples 
sample_names="Ire-juv-1C  Ire-juv-2C  Ire-juv-21C  Ire-juv-22C  Ire-juv-41C  Ire-juv-42C  Ire-juv-61C  Ire-juv-62C  Ire-juv-81C  Ire-juv-82C"


# this loop writes the bash script
for i in $sample_names; do
echo "the sample name is: "$i
read_1=$(less paths_juv_ire.txt | grep "$i" | grep 'R1')
read_2=$(less paths_juv_ire.txt | grep "$i" | grep 'R2')
echo "the sample name is: "$i  "   the read 1 is " $read_1   "   the read 2 is "$read_2

echo "#!/bin/bash -l" > $i"_ire.sh"
echo "#SBATCH -A snic2021-5-20" >> $i"_ire.sh"
echo "#SBATCH -p core -n 10" >> $i"_ire.sh"
echo "#SBATCH -J Trimming_$i" >> $i"_ire.sh"
echo "#SBATCH -t 00:20:00" >> $i"_ire.sh"
echo "module load bioinfo-tools" >> $i"_ire.sh"
echo "module load cutadapt">> $i"_ire.sh"
echo "cutadapt -a AGATCGGAAGAGC -A AGATCGGAAGAGC -j 10 -q 30,30 -o /proj/uppstore2017185/b2014034_nobackup/Xuejing/trimmed_data/juv_ire/"$i"_R1_trimmed.fastq.gz -p /proj/uppstore2017185/b2014034_nobackup/Xuejing/trimmed_data/juv_ire/"$i"_R2_trimmed.fastq.gz  "$read_1" "$read_2 "" >> $i"_ire.sh"
echo " "  >> $i"_ire.sh"
echo "created a script for ${i}_ire.sh"
echo -e "\n"

done


# this runs/submits the bash script
for i in $sample_names; do
sbatch $i"_ire.sh"
done




