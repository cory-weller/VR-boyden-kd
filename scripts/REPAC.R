# 1 load and install libraries

#BiocManager::install("edgeR")
#BiocManager::install("limma")

library(REPAC)
library(edgeR)
library(limma)
library(SummarizedExperiment)
library(data.table)

# gtf <- 'gencode.v45.primary_assembly.annotation.gtf'

bam_dir <- 'RNA_alignment'
data(hg38_pa)

files <- data.table(path=unlist(lapply(list.files(bam_dir, pattern='soma-.*.bam|neurite-.*.bam', recursive=T, full.names=T), normalizePath)))
files[, SampleID := tstrsplit(basename(path), split='\\.')[1]]
files[, Condition := tstrsplit(SampleID, split='-')[1]]





# 5 External dataset function from REPAC is rse_from_bam
# 6 Analysis
rse <- rse_from_bam(
            bams=files$path, 
            annotation=hg38_pa, 
            colData=files[, .SD, .SDcols=c('SampleID','Condition')], 
            strandSpecific = 0,
            allowMultiOverlap = F, 
            fraction = F,
            countMultiMappingReads=F, 
            isPairedEnd = T,
            nthreads = 12, 
            countReadPairs = T,
            requireBothEndsMapped = F)


# 7 Check analysis and all columns are there
colnames(colData(rse)) # list available metadata columns

# 8 Set group name (should be condition from above)
groups <- factor(rse$Condition, levels=c('soma','neurite'))
table(groups) # Ensure both samples are there 

# 9 Filter any low expression if desired 
keep <- edgeR::filterByExpr(rse, group = groups)
rse <- rse[keep, ]

# 10 Create Design Matrix
dMat <- model.matrix(~ 0 + groups)
colnames(dMat) <- levels(groups)
print(dMat)

cMat <- makeContrasts(
              levels = colnames(dMat),
              NeuritevsSoma = neurite-soma)

# 11 Fit model
fit <- fit_repac(rse,dMat,cMat)
head(fit)

results <- resolve_comparisons(fit)$NeuritevsSoma
setDT(results)

dir.create('APA')
fwrite(results, file="APA/REPAC-results.csv", sep=',', row.names=F, col.names=T, quote=F)



library(data.table)
library(ggplot2)
library(ggthemes)
dat <- fread('APA/REPAC-results.csv')
setnames(dat, 'adjusted.P.Value','p.adj')

p_threshold <- 0.05
min_fc <- 0.5

source('scripts/colors.R')

colors <- c('Whole cell'=soma_color, 'Neurite'=neurite_color, 'NS'=not_signif_color)

dat[p.adj < p_threshold & cFC < 0, Significance := 'Whole cell']
dat[p.adj < p_threshold & cFC > 0, Significance := 'Neurite']
dat[p.adj >= p_threshold, Significance := 'NS']
dat[cFC %between% c(-min_fc, min_fc), Significance := 'NS']


ggplot(dat, aes(x=cFC, y=-1*log10(p.adj), color=Significance)) +
    geom_point() +
    labs(x='Compositional fold-change', y='-log10(P)', title='REPAC Neurite vs Whole Cell') +
    scale_color_manual(values=colors) +
    theme_few() +
    geom_hline(yintercept=-1*log10(p_threshold), linetype='dashed', alpha=0.5) +
    geom_vline(xintercept=c(-min_fc, min_fc), linetype='dashed', alpha=0.5)
    
dat.test <- as.data.table(rse@assays@data@listData, keep.rownames=T)


ZNF644   ZNF644_02   ZNF644_04 -1.4066513 -0.45556850 -10.346881 4.766929e-05 9.533858e-05   Whole cell

quit(status=0)