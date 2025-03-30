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

colors <- c('Whole cell'=soma_color, 'Neurite'=neurite_color, 'NS'=not_signif_color)


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
            
            dat[p.adj < p_threshold & cFC < 0, Significance := 'Whole cell']
            dat[p.adj < p_threshold & cFC > 0, Significance := 'Neurite']
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

g <- ggplot(dat, aes(x=cFC, y=-1*log10(p.adj), color=Significance)) +
    geom_point() +
    labs(x='Compositional fold-change', y='-log10(P)', title='REPAC Neurite vs Whole Cell') +
    scale_color_manual(values=colors) +
    theme_few() +
    geom_hline(yintercept=-1*log10(p_threshold), linetype='dashed', alpha=0.5) +
    geom_vline(xintercept=c(-min_fc, min_fc), linetype='dashed', alpha=0.5)

# g <- ggplot(dat[Significance != 'NS'], aes(x=cFC, y=-1*log10(p.adj), color=annotation)) +
#     geom_point() +
#     labs(x='Compositional fold-change', y='-log10(P)', title='REPAC Neurite vs Whole Cell') +
#     theme_few() +
#     geom_hline(yintercept=-1*log10(p_threshold), linetype='dashed', alpha=0.5) +
#     geom_vline(xintercept=c(-min_fc, min_fc), linetype='dashed', alpha=0.5)


ggsave(g, file='APA/REPAC-volcano.png', width=15, height=15, units='cm')
ggsave(g, file='APA/REPAC-volcano.pdf', width=15, height=15, units='cm')


writeLines(as.character(unique(dat$gene_name)), con='APA/REPAC-set.txt')
writeLines(as.character(unique(dat[Significance=='Neurite']$gene_name)), con='APA/REPAC-Neurite.txt')
writeLines(as.character(unique(dat[Significance=='Whole cell']$gene_name)), con='APA/REPAC-Soma.txt')


dat.test <- as.data.table(rse@assays@data@listData, keep.rownames=T)






quit(status=0)


# Longest vs other sites
dat <- dat[annotation == '3UTR']