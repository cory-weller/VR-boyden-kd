# VR-kolf-boyden-rna


# Preprocessing
Run UMI deduplication, mark optical duplicates, align and count features

```bash
sbatch --array=1-76 scripts/align-and-counts.sh 
```



```bash
ml R/4.3
 # Concatenate raw counts table, normalize with DESeq2, PCA analysis, volcano and MA plots
Rscript scripts/DESeq2.R
```

One sample `TDP-43.20` was excluded due to separation on PC1 (99% var) from all other samples