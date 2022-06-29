# a pipeline for L.juv.swe samples 
sample_names="S17_L003  S17_L001  S19_L003  S19_L001  S20_L003  S20_L001  S21_L003  S21_L001  S24_L003  S24_L001 S25_L003 S25_L001 S26_L003  S26_L001 S27_L003 S27_L001 S28_L003 S28_L001 S32_L003 S32_L001 S33_L003 S33_L001 S86_L003 S86_L001"


# this loop writes the bash script
for i in $sample_names; do
echo "the sample name is: "$i
read_1=$(less paths_juv_swe.txt | grep "$i" | grep 'R1')
read_2=$(less paths_juv_swe.txt | grep "$i" | grep 'R2')
echo "the sample name is: "$i  "   the read 1 is " $read_1   "   the read 2 is "$read_2

echo "#!/bin/bash -l" > $i"_swe.sh"
echo "#SBATCH -A snic2021-5-20" >> $i"_swe.sh"
echo "#SBATCH -p core -n 10" >> $i"_swe.sh"
echo "#SBATCH -J Trimming_$i" >> $i"_swe.sh"
echo "#SBATCH -t 00:20:00" >> $i"_swe.sh"
echo "module load bioinfo-tools" >> $i"_swe.sh"
echo "module load cutadapt">> $i"_swe.sh"
echo "cutadapt -a AGATCGGAAGAGC -A AGATCGGAAGAGC -j 10 -q 30,30 -o /proj/uppstore2017185/b2014034_nobackup/Xuejing/trimmed_data/juv_swe/"$i"_R1_trimmed.fastq.gz -p /proj/uppstore2017185/b2014034_nobackup/Xuejing/trimmed_data/juv_swe/"$i"_R2_trimmed.fastq.gz  "$read_1" "$read_2 "" >> $i"_swe.sh"
echo " "  >> $i"_swe.sh"
echo "created a script for ${i}_swe.sh"
echo -e "\n"

done


# this runs/submits the bash script
for i in $sample_names; do
sbatch $i"_swe.sh"
done




