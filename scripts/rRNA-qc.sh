#!/usr/bin/env bash


DICT='/fdb/STAR_current/GENCODE/Gencode_human/release_45/ref.dict'

# Start with header, tab-separated
awk -v OFS='\t' 'NR > 1 {print $1','$2','$3}' $DICT > rRNA_intervals.hg38.txt

# Add rRNA transcripts to tmp file
grep -F 'gene_type "rRNA"' ${GTF} | awk -v FS='\t' -v OFS='\t' '$3=="transcript" {print $1,$4,$5,$7,$9}' > tmp1.txt

grep -F 'gene_type "rRNA_pseudogene"' ${GTF} | awk -v FS='\t' -v OFS='\t' '$3=="transcript" {print $1,$4,$5,$7,$9}' > tmp2.txt

cat tmp1.txt tmp2.txt > tmp3.txt


# Reformat and append to final output
paste <(cut -f 1-4 tmp.txt) \
<(cut -f 5 tmp.txt | sed -r 's/.*transcript_id "//g' | sed -r 's/"; .*$//g') > tmp4.txt

module load bedtools
bedtools sort -i tmp4.txt >> rRNA_intervals.hg38.txt


# Remove tmp files
rm tmp{1,2,}.txt

python3 scripts/collapse_annotation.py gencode.v45.primary_assembly.annotation.gtf gencode.v45.genes.gtf