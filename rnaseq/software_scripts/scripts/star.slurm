#!/bin/bash

#SBATCH --job-name=star # Job name
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=8
#SBATCH --account=publicgrp
#SBATCH --partition=low
#SBATCH --time=60
#SBATCH --mem=32000 # Memory pool for all cores (see also --mem-per-cpu)
#SBATCH --array=1-22
#SBATCH --output=slurmout/star_%A_%a.out # File to which STDOUT will be written
#SBATCH --error=slurmout/star_%A_%a.err # File to which STDERR will be written

start=`date +%s`
echo $HOSTNAME
echo "My SLURM_ARRAY_TASK_ID: " $SLURM_ARRAY_TASK_ID

sample=`sed "${SLURM_ARRAY_TASK_ID}q;d" samples.txt`
REF="/quobyte/ikorfgrp/bis180l/mRNAseq/References/star.overlap100.gencode.M35"

outpath='02-STAR_alignment'
[[ -d ${outpath} ]] || mkdir ${outpath}
[[ -d ${outpath}/${sample} ]] || mkdir ${outpath}/${sample}

echo "SAMPLE: ${sample}"

module load star/2.7.11a

call="STAR
     --runThreadN 8 \
     --genomeDir $REF \
     --outSAMtype BAM SortedByCoordinate \
     --readFilesCommand zcat \
     --readFilesIn 01-HTS_Preproc/${sample}/${sample}_R1.fastq.gz 01-HTS_Preproc/${sample}/${sample}_R2.fastq.gz \
     --quantMode GeneCounts \
     --outFileNamePrefix ${outpath}/${sample}/${sample}_ \
     > ${outpath}/${sample}/${sample}-STAR.stdout 2> ${outpath}/${sample}/${sample}-STAR.stderr"

echo $call
eval $call

end=`date +%s`
runtime=$((end-start))
echo $runtime
