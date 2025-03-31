#!/usr/bin/env bash 
#SBATCH --time 0:15:00
#SBATCH --partition quick

#module load samtools/1.21
#REFFLAT='/vf/users/CARD_ARDIS/users/wellerca/VR-boyden-kd/refFlat.txt.gz'

samplefile='samples.txt'
N=${SLURM_ARRAY_TASK_ID}

TRT=$(sed -n ${N}p ${samplefile} | cut -f 1)
SAMPLE=$(sed -n ${N}p ${samplefile} | cut -f 2)

BAM="RNA_alignment/${TRT}/${SAMPLE}.dedup.bam"

# samtools stats --threads 4 ${BAM} > ${BAM}.stats

# module load multiqc/1.28

# module load GATK

# gatk CollectRnaSeqMetrics \
#     --INPUT ${BAM} \
#     --OUTPUT ${BAM%.bam}.metrics.txt \
#     --REF_FLAT ${REFFLAT} \
#     --STRAND_SPECIFICITY NONE \
#     --RIBOSOMAL_INTERVALS rRNA_intervals.hg38.txt

GENEGTF='gencode.v45.genes.gtf'
module load rnaseqc
rnaseqc ${GENEGTF} ${BAM} ${BAM%.bam}.rnaseqc.txt

