#!/usr/bin/env bash
#SBATCH --mem 100G
#SBATCH --time 3:59:00
#SBATCH --nodes 1
#SBATCH --ntasks 1
#SBATCH --cpus-per-task 16
#SBATCH --gres=lscratch:400
#SBATCH --partition quick,norm

cd /data/CARD_ARDIS/users/wellerca/VR-boyden-kd
REF='/fdb/STAR_current/GENCODE/Gencode_human/release_45/genes-100'
REF=$(realpath $REF)
ANNO=$(realpath gencode.v45.primary_assembly.annotation.gtf)
N=${SLURM_ARRAY_TASK_ID}
TREATMENT=$(sed -n ${N}p samples.txt | cut -f 1)
SAMPLEID=$(sed -n ${N}p samples.txt | cut -f 2)
OUTDIR=$(realpath .)
OUTDIR="${OUTDIR}/processing/$TREATMENT/"
mkdir -p $OUTDIR

READ1s=$(find -L data/$TREATMENT -type f -regex ".*/${SAMPLEID}_S[0-9]+.*_R1.*.fastq.gz" -exec realpath {} \;)
READ2s=$(find -L data/$TREATMENT -type f -regex ".*/${SAMPLEID}_S[0-9]+.*_R2.*.fastq.gz" -exec realpath {} \;)

TMPDIR="/lscratch/$SLURM_JOB_ID"
mkdir -p $TMPDIR && cd $TMPDIR

# Concatenate reads
cat ${READ1s} > R1.fastq.gz
cat ${READ2s} > R2.fastq.gz

## Deduplication with fastp
# https://github.com/OpenGene/fastp

module load fastp/0.24

if [[ ${TREATMENT} == 'WT' ]]; then
    # Deduplicate based on full sequence
    fastp \
        --disable_adapter_trimming \
        -i  R1.fastq.gz \
        -o R1.dedup.untrimmed.fastq.gz \
        -I R2.fastq.gz \
        -O R2.dedup.untrimmed.fastq.gz \
        --dedup
else
    # Collapse reads based on UMI
    fastp \
        --disable_adapter_trimming \
        -i R1.fastq.gz \
        -o R1.dedup.untrimmed.fastq.gz \
        -I R2.fastq.gz \
        -O R2.dedup.untrimmed.fastq.gz \
        -U \
        --umi_loc=read2 \
        --umi_len=8
fi

# Clean up raw reads
rm R1.fastq.gz
rm R2.fastq.gz

# Trim first 6 nt from front of R2
fastp \
    --disable_adapter_trimming \
    -i R1.dedup.untrimmed.fastq.gz \
    -o R1.final.fastq.gz \
    -I R2.dedup.untrimmed.fastq.gz \
    -O R2.final.fastq.gz \
    --trim_front2 6

# Clean up untrimmed reads

rm R1.dedup.untrimmed.fastq.gz
rm R2.dedup.untrimmed.fastq.gz

module load STAR/2.7.11b

STAR \
    --runThreadN 16 \
    --genomeDir $REF \
    --sjdbOverhang 100 \
    --quantMode GeneCounts \
    --readFilesIn R1.final.fastq.gz R2.final.fastq.gz \
    --readFilesCommand zcat \
    --outSAMtype BAM SortedByCoordinate \
    --outFileNamePrefix "${SAMPLEID}_"

echo "After STAR:"
ls

# STAR for grp 2
mv *.bam star.bam

ls
module load GATK

# Add read groups
gatk AddOrReplaceReadGroups \
       --INPUT star.bam \
       --OUTPUT rg.bam \
       --RGID S1 \
       --RGLB lib1 \
       --RGPL ILLUMINA \
       --RGPU unit1 \
       --RGSM S1

echo "After ReadGroups:"
ls

# Remove duplicate reads
gatk MarkDuplicates \
    --INPUT rg.bam \
    --OUTPUT ${SAMPLEID}.dedup.bam \
    --METRICS_FILE ${SAMPLEID}.dupstats.txt \
    --REMOVE_DUPLICATES


echo "After Duplicates:"
ls

# Run featureCounts
module load subread/2.0.6
featureCounts \
    -p \
    --countReadPairs \
    -M \
    --ignoreDup \
    --fraction \
    -T 16 \
    -t exon \
    -g gene_id \
    -a ${ANNO} \
    -o "${SAMPLEID}.featurecounts.txt" \
    ${SAMPLEID}.dedup.bam

echo "After featureCounts:"
ls

cp ${SAMPLEID}.dedup.bam $OUTDIR
cp ${SAMPLEID}.dupstats.txt $OUTDIR
cp ${SAMPLEID}.featurecounts.txt $OUTDIR
cp ${SAMPLEID}.featurecounts.txt.summary $OUTDIR
