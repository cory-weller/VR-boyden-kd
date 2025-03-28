#!/usr/bin/env bash
#SBATCH

DB='/fdb/kraken/20220803_standard_kraken2/'
TMPDIR="/lscratch/${SLURM_JOB_ID}"

mkdir -p $TMPDIR && cd $TMPDIR
cp -r ${DB} .

DB=$(basename $DB)

module load kraken/2.1.2




SAMPLEID='neurite-4'


R1='/vf/users/CARD_ARDIS/2012_G3boydenRNAseq/raw_data/fastq210201/neurite_4_S1_R1_001.fastq.gz'
R2='/vf/users/CARD_ARDIS/2012_G3boydenRNAseq/raw_data/fastq210201/neurite_4_S1_R2_001.fastq.gz'

cp ${R1} .
cp ${R2} .

R1=$(basename $R1)
R2=$(basename $R2)

SLURM_CPUS_PER_TASK=12


kraken2 \
    --db ${DB} \
    --paired \
    --threads ${SLURM_CPUS_PER_TASK} \
    --output ${SAMPLEID}.kraken \
    --classified-out ${SAMPLEID}#.fq \
    --report ${SAMPLEID}.report \
    ${R1} ${R2}

kraken-report --db $DBNAME kraken.output