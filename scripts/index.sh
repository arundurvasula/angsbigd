#!/bin/bash -l
#SBATCH -D /home/adurvasu/mergedBams
#SBATCH -J index
#SBATCH -o outs/out-%j.txt
#SBATCH -p bigmem
#SBATCH -e errors/error-%j.txt
#SBATCH -c 8

module load gcc
module load bwa

#samtools faidx data/TRIP.fa.gz
#samtools faidx /home/adurvasu/genomes/Zea_mays.AGPv3.21.dna.genome.fa.gz
for f in *.bam; do
samtools index /home/adurvasu/mergedBams/$f;
done