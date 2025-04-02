# VR-kolf-boyden-rna


## rRNA interval list

```bash

wget https://gist.githubusercontent.com/slowkow/b11c28796508f03cdf4b/raw/38d337698ff1e6915578dfa08826c73631c3e0b5/make_rRNA.sh



GTF='/fdb/STAR_current/GENCODE/Gencode_human/release_45/genes.gtf'

gencode.v45.primary_assembly.annotation.gtf



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


rRNA_intervals.hg38.txt

# Remove tmp files
rm tmp{1,2,}.txt

sbatch --array=1-7676 scripts/multiqc.sh

rnaseqc ${GTF} ${BAM} ${BAM}.rnaseqc

python3 scripts/collapse_annotation.py gencode.v45.primary_assembly.annotation.gtf gencode.v45.genes.gtf
```
## Preprocessing sequencing data
Run UMI deduplication, mark optical duplicates, align and count features

```bash
sbatch --array=1-76 scripts/align-and-counts.sh 
```


## Differential RNA Expression Analysis
```bash
ml R/4.3
 # Concatenate raw counts table, normalize with DESeq2, PCA analysis, volcano and MA plots
Rscript scripts/DESeq2.R
```

Note: one sample `TDP-43.20` was excluded due to separation on PC1 (99% var) from all other samples

## Alternative Poly-Adenylation

For REPAC, significance: |compositional Fold Change (cFC)| >= 0.25 AND a FDR corrected p-value < 0.05
For APAlyzer, FDR-corrected p-value < 0.05
```R
Rscript scripts/REPAC.R
Rscript scripts/APAlyzer.R
Rscript scripts/APA-plot-distributions.R
```