#!/usr/bin/env bash 
#SBATCH --time 0:15:00
#SBATCH --partition quick

#module load samtools/1.21

samplefile='samples.txt'
N=${SLURM_ARRAY_TASK_ID}

TRT=$(sed -n ${N}p ${samplefile} | cut -f 1)
SAMPLE=$(sed -n ${N}p ${samplefile} | cut -f 2)

BAM="RNA_alignment/${TRT}/${SAMPLE}.dedup.bam"


GENEGTF='gencode.v45.genes.gtf'
module load rnaseqc
rnaseqc ${GENEGTF} ${BAM} ${BAM%.bam}.rnaseqc.txt

