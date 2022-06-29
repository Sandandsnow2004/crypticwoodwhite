# a pipeline for L.juv.Kaz samples 
sample_names="Sample_13  Sample_14  Sample_15  Sample_17  Sample_18  Sample_20  Sample_21  Sample_22  Sample_23  Sample_24"


# this loop writes the bash script
for i in $sample_names; do
echo "the sample name is: "$i
read_1=$(less paths_juv_kaz.txt | grep "$i" | grep 'R1')
read_2=$(less paths_juv_kaz.txt | grep "$i" | grep 'R2')
echo "the sample name is: "$i  "   the read 1 is " $read_1   "   the read 2 is "$read_2

echo "#!/bin/bash -l" > $i"_kaz.sh"
echo "#SBATCH -A snic2021-5-20" >> $i"_kaz.sh"
echo "#SBATCH -p core -n 10" >> $i"_kaz.sh"
echo "#SBATCH -J Trimming_$i" >> $i"_kaz.sh"
echo "#SBATCH -t 00:20:00" >> $i"_kaz.sh"
echo "module load bioinfo-tools" >> $i"_kaz.sh"
echo "module load cutadapt">> $i"_kaz.sh"
echo "cutadapt -a AGATCGGAAGAGC -A AGATCGGAAGAGC -j 10 -q 30,30 -o /proj/uppstore2017185/b2014034_nobackup/Xuejing/trimmed_data/juv_kaz/"$i"_R1_trimmed.fastq.gz -p /proj/uppstore2017185/b2014034_nobackup/Xuejing/trimmed_data/juv_kaz/"$i"_R2_trimmed.fastq.gz  "$read_1" "$read_2 "" >> $i"_kaz.sh"
echo " "  >> $i"_kaz.sh"
echo "created a script for ${i}_kaz.sh"
echo -e "\n"

done


# this runs/submits the bash script
for i in $sample_names; do
sbatch $i"_kaz.sh"
done




