# VR-kolf-boyden-rna


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