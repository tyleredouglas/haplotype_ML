#! /bin/bash
# Job name:
#SBATCH --job-name=simul_100
#
# Account:
#SBATCH --account=fc_poison
#
# Partition:
#SBATCH --partition=savio3
#
# Wall clock limit:
#SBATCH --time=48:00:00
#
## Command(s) to run:
source activate lgbm3.1  
python -c "import pysam"
python -c "import pandas as pd" 
python -c "import numpy as np"
python -c "from Bio import SeqIO"
python simul_100.py

module load bwa
module load samtools
module load bcftools

for file in *_1.fastq; do
    name=$(echo $file | sed 's/_1.fastq//')
    file2="${name}_2.fastq"                
    bwa mem -t 20 -M dmel_r6.fna "$file" "$file2" | samtools sort -@ 20 -o "${name}.sorted.bam"  # Align and sort
    samtools index -@ 20 "${name}.sorted.bam"  
done

bcftools mpileup -I -q 0 -Q 0 --threads 20 -t NT_037436.4 -a "FORMAT/AD,FORMAT/DP" -f dmel_r6.fna -b bam_names.txt -o sim100_mpile.txt
bcftools call --threads 20 sim100_mpile.txt -mv -Ov > sim100_mpile.txt 
