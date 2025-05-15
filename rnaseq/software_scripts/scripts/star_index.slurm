#!/bin/bash

#SBATCH --job-name=star_index # Job name
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=8
#SBATCH --account=publicgrp
#SBATCH --partition=high
#SBATCH --time=120
#SBATCH --mem=40000 # Memory pool for all cores (see also --mem-per-cpu)
#SBATCH --output=slurmout/star-index_%A.out # File to which STDOUT will be written
#SBATCH --error=slurmout/star-index_%A.err # File to which STDERR will be written

start=`date +%s`
echo $HOSTNAME

outpath="References"
cd ${outpath}

refdir="/quobyte/ikorfgrp/bis180l/mRNAseq/References"
FASTA="$refdir/GRCm39.primary_assembly.genome.fa"

#wget https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/release_M32/gencode.vM32.primary_assembly.annotation.gtf.gz
#gunzip gencode.vM32.primary_assembly.annotation.gtf.gz
#GTF="../gencode.vM32.primary_assembly.annotation.gtf"
GTF="$refdir/gencode.vM37.basic.annotation.gtf"

mkdir star.overlap100.gencode.M35
cd star.overlap100.gencode.M35

module load star/2.7.11a

call="STAR
    --runThreadN ${SLURM_NTASKS} \
    --runMode genomeGenerate \
    --genomeDir . \
    --genomeFastaFiles ${FASTA} \
    --sjdbGTFfile ${GTF} \
    --sjdbOverhang 100"

echo $call
eval $call

end=`date +%s`
runtime=$((end-start))
echo $runtime
