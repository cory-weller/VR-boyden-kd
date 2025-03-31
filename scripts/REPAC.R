# 1 load and install libraries

#BiocManager::install("edgeR")
#BiocManager::install("limma")

library(REPAC)
library(edgeR)
library(limma)
library(SummarizedExperiment)
library(data.table)
library(ggplot2)
library(ggthemes)

# gtf <- 'gencode.v45.primary_assembly.annotation.gtf'

bam_dir <- 'RNA_alignment'
data(hg38_pa)
rse_file <- 'APA/REPAC-rse.RDS'
results_file <- 'APA/REPAC-results.csv'

files <- data.table(path=unlist(lapply(list.files(bam_dir, pattern='soma-.*.bam|neurite-.*.bam', recursive=T, full.names=T), normalizePath)))
files[, SampleID := tstrsplit(basename(path), split='\\.')[1]]
files[, Condition := tstrsplit(SampleID, split='-')[1]]

tryCatch(dir.create('APA'), warning=function(w) {cat('APA dir already exists\n') } )


p_threshold <- 0.05
min_fc <- 0.25

source('scripts/colors.R')



tryCatch(rse <- readRDS(rse_file), 
        error = function(e) { 
            cat('oh no error\n')
        }, 
        warning = function(w) { 
            cat('Generating rse object\n')
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
        saveRDS(rse, rse_file)
        }, 
        finally = {cat('rse loaded\n')}
)


tryCatch(dat <- fread(results_file),
        error = function(e) { 
            cat('oh no error\n')
        }, 
        warning = function(w) { 
            # 8 Set group name (should be condition from above)
            groups <- factor(rse$Condition, levels=c('soma','neurite'))
            
            # Filter low-expressed genes
            keep <- edgeR::filterByExpr(rse, group = groups)
            rse <- rse[keep, ]
            
            # 10 Create Design Matrix
            dMat <- model.matrix(~ 0 + groups)
            colnames(dMat) <- levels(groups)
            print(dMat)
            
            cMat <- makeContrasts(
                        levels = colnames(dMat),
                        NeuritevsSoma = neurite-soma)
            
            # Generate model fit
            fit <- fit_repac(rse,dMat,cMat)
            head(fit)
            
            # Get contrast table
            dat <- resolve_comparisons(fit)$NeuritevsSoma
            setDT(dat)
            
            setnames(dat, 'adjusted.P.Value','p.adj')
            
            dat[p.adj < p_threshold & cFC < 0, Significance := 'Shortened in Neurites']
            dat[p.adj < p_threshold & cFC > 0, Significance := 'Lengthened in Neurites']
            dat[p.adj >= p_threshold, Significance := 'NS']
            dat[cFC %between% c(-min_fc, min_fc), Significance := 'NS']
            
            fwrite(dat, file=results_file, sep=',', row.names=F, col.names=T, quote=F)
        },
        finally = {cat('results file loaded\n')}
)


# anno <- as.data.table(as.data.frame(rse@rowRanges, keep.rownames=T))

# anno[strand == '-', pos := end]
# anno[strand == '+', pos := start]
# anno_sub <- anno[, .SD, .SDcols=c('geneID','pos','strand')]
# setnames(anno_sub, 'geneID','site')

# dat <- merge(dat, anno_sub, by.x='Site', by.y='site')
# setnames(dat, 'pos', 'site_tested_pos')
# dat <- merge(dat, anno_sub, by.x='Ref', by.y='site')
# setnames(dat, 'pos', 'site_ref_pos')
# dat[, ref_IDN := tstrsplit(Ref, split='_')[2]]
# dat[, ref_IDN := as.numeric(ref_IDN)]
# dat[, test_IDN := tstrsplit(Site, split='_')[2]]
# dat[, test_IDN := as.numeric(test_IDN)]


# dat <- merge(anno, dat, by.x='geneID', by.y='Site')

colors <- c('Shortened in Neurites'=neurite_down_color, 'Lengthened in Neurites'=neurite_color, 'NS'=not_signif_color, 'Gene-level mean'='blue')


g1 <- ggplot(dat, aes(x=cFC, y=-1*log10(p.adj), color=Significance)) +
    geom_point() +
    labs(x='Compositional fold-change', y='-log10(P)', title='REPAC: all transcript comparisons') +
    scale_color_manual(values=colors) +
    theme_few() +
    geom_hline(yintercept=-1*log10(p_threshold), linetype='dashed', alpha=0.5) +
    geom_vline(xintercept=c(-min_fc, min_fc), linetype='dashed', alpha=0.5) +
    theme(legend.position="bottom",legend.direction="vertical")

ggsave(g1, file='APA/REPAC-volcano.png', width=12, height=15, units='cm')
ggsave(g1, file='APA/REPAC-volcano.pdf', width=12, height=15, units='cm')


sizes <- c('comparison'=10, 'mean'=16)
alphas <- c('comparison'=0.8, 'mean'=1)

dat[,  'meancFC' := mean(cFC), by=gene_name]

means <- dat[, list('cFC' = mean(cFC)), by=gene_name][order(cFC)]
gene_order <- dat[, list('cFC' = mean(cFC)), by=gene_name][order(cFC)]$gene_name
means[, gene_name := factor(gene_name, levels=gene_order)]
means[, Significance := 'Gene-level mean']
dat[, gene_name := factor(gene_name, levels=gene_order)]
setkey(dat, gene_name)



g2 <- ggplot(dat, aes(x=gene_name, y=cFC, color=Significance)) +
    geom_point() +
    geom_point(data=dat[Significance != 'NS']) +
    geom_point(data=means) +
    labs(x='Genes ranked by mean cFC across all transcript comparisons (blue)', y='Compositional fold change', title='REPAC: gene-level means') +
    scale_color_manual(values=colors) +
    theme_few() +
    geom_hline(yintercept=0, linetype='dashed', alpha=0.5) +
    theme(legend.position="bottom",legend.direction="vertical") +
    geom_vline(xintercept=c('SNX27', 'PITPNC1', 'CCNI'), linetype='dashed', alpha=0.5) +
    theme(axis.text.x = element_blank(), axis.ticks.x=element_blank()) +
    annotate("text", x = 280, y = -2, label = "mean cFC < -0.25", hjust=0, angle=90, vjust=0) +
    annotate("text", x = 1664, y = -2, label = "mean cFC > +0.25", hjust=0, angle=90, vjust=0)

ggsave(g2, file='APA/REPAC-gene-level-mean.png', width=18, height=15, units='cm')
ggsave(g2, file='APA/REPAC-gene-level-mean.pdf', width=18, height=15, units='cm')

writeLines(as.character(unique(dat$gene_name)), con='APA/REPAC-set.txt')
writeLines(as.character(unique(means[cFC < -0.25]$gene_name)), con='APA/REPAC-neurite-shortened.txt')
writeLines(as.character(unique(means[cFC >  0.25]$gene_name)), con='APA/REPAC-neurite-lengthened.txt')





# Save means file for plotting
means[, 'Significance' := 'NS']
means[cFC < -0.25, 'Significance' := 'Shortened in Neurites']
means[cFC > 0.25, 'Significance' := 'Lengthened in Neurites']
fwrite(means, file='APA/REPAC-gene-level-significance.csv', sep=',', quote=F, row.names=F, col.names=T)


quit(status=0)


# Longest vs other sites
dat <- dat[annotation == '3UTR']