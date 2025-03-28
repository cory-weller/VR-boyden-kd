#!/usr/bin/env Rscript

library(data.table)
library(foreach)
library(edgeR)
library(limma)
library(ggrepel)
library(ggthemes)

# Import Count Data
if (! file.exists('DE/raw_counts.tsv')) {
    # Get list of all featurecounts files
    files <- list.files('processing',
                    recursive=TRUE, 
                    pattern='*.featurecounts.txt$', 
                    full.names=TRUE
                    )
    # Read featurecounts data into single list
    o <- foreach(x=files) %do% fread(x, select=c(1,7))
    
    # Merge into single table
    countData <- Reduce(function(...) merge(..., all = TRUE, by='Geneid'), o)
    
    # Remove '.dedup.bam' from column names
    setnames(countData, gsub('.dedup.bam', '', colnames(countData)))
    
    # Modify names to match experiment
    setnames(countData, gsub('VR-8153-', 'FUS-', colnames(countData)))
    setnames(countData, gsub('JR-7961-', 'TDP43-', colnames(countData)))
    setnames(countData, gsub('VR-8147-', 'hnRNPA1-', colnames(countData)))
    setnames(countData, gsub('-', '.', colnames(countData)))
    
    # Ensure output directory exists
    dir.create('DE')
    
    # Save merged counts as tabular output
    fwrite(countData, file='DE/raw_counts.tsv', sep='\t', quote=F, row.names=F, col.names=T)
} else {
    # Read pre-saved counts data
    countData <- fread("DE/raw_counts.tsv")
}

# Import metadata
metaData <- fread('metadata/rnaseq-metadata.csv')
metaData[, filestem := gsub('VR-8153-', 'FUS-', filestem)]
metaData[, filestem := gsub('JR-7961-', 'TDP43-', filestem)]
metaData[, filestem := gsub('VR-8147-', 'hnRNPA1-', filestem)]
metaData[, filestem := gsub('-', '.', filestem)]

# get annotations
ANNO <- fread('https://gist.github.com/cory-weller/566a8984b713fc4626ec423913fdab43/raw/90e5eeaace6ee4d6900a404cbd50c952015b40dc/gene_table_GRCh38-2020-A.tsv')
ANNO <- ANNO[, .SD, .SDcols=c('ENSEMBL','SYMBOL')]
ANNO <- ANNO[!duplicated(ENSEMBL)]

# Define function
getContrast <- function(counts, metadata, experiment, fractions, contrast_baseline, annotations, manuallabel, min_pval=0.01, min_fc=1, labelup=10, labeldown=10) {
    
    # Define colors
    neurite_color <- '#d63cc0'
    soma_color <- '#00a347'
    neurite_down_color <- '#6c1760'
    soma_down_color <- '#003f1e'
    not_signif_color <- '#c4c4c4'

    metadata.sub <- metadata[Experiment==experiment & Fraction %in% fractions]
    samples <- metadata.sub$filestem
    counts.sub <- counts[, .SD, .SDcols=c('Geneid',samples)]    # samples should now be in order of metadata
    groups <- as.factor(as.numeric(!(metadata.sub$Contrast == contrast_baseline)) + 1)
    counts.sub <- data.frame(counts.sub[,-1], row.names=counts.sub$Geneid)
    dim(counts.sub)
    dge <- DGEList(counts.sub, group=groups)
    keep <- filterByExpr(dge, group=groups, min.prop=0.75)
    dge <- dge[keep,,keep.lib.sizes=FALSE]
    dge <- DGEList(dge, group=groups)
    
    #dge  <- calcNormFactors(dge)
    v <- voom(dge, plot=FALSE) #check voom plot for curve to see if we need to do more filtering
    vfit <- lmFit(v)
    efit <- eBayes(vfit)
    output <- as.data.table(topTable(efit, sort.by = "P", n = Inf, adjust='BH'), keep.rownames=TRUE)    #output full table
    setnames(output, 'rn', 'Geneid')
    output[, Significance := 'Not Significant']
    output[adj.P.Val < min_pval & logFC > min_fc, Significance := 'Up']
    output[adj.P.Val < min_pval & logFC < -min_fc, Significance := 'Down']
    output[adj.P.Val >= min_pval, Significance := 'Not Significant']
    output[, Geneid := gsub('\\.[0-9]+$','', Geneid)]
    output <- merge(output, annotations, by.x='Geneid', by.y='ENSEMBL', all.x=T)
    output[is.na(SYMBOL), SYMBOL := Geneid]
    output[SYMBOL %in% manuallabel, 'lbl' := SYMBOL]
    output[logFC > 0, up_rank := frank(adj.P.Val, ties.method='first')]
    output[logFC < 0, down_rank := frank(adj.P.Val, ties.method='first')]
    output[up_rank <= labelup, 'lbl' := SYMBOL]
    output[down_rank <= labeldown, 'lbl' := SYMBOL]
    output[Significance == 'Not Significance', 'lbl' := NA]
    
    
    if(experiment == 'WT') {
        # if WT
        output[logFC > 0, 'clr' := neurite_color]
        output[logFC < 0, 'clr' := soma_color]
    } else if (fractions == 'Soma') {
        output[logFC > 0, 'clr' := soma_color]
        output[logFC < 0, 'clr' := soma_down_color]
    } else if (fractions == 'Neurite') {
        output[logFC > 0, 'clr' := neurite_color]
        output[logFC < 0, 'clr' := neurite_down_color]
    }
    
    output[Significance == 'Not Significant', clr := 'gray' ]

    
    setnames(output, 'logFC', 'log2FC')
    if(experiment == 'WT') {
    g.volcano <- ggplot(data=output, aes(x=log2FC, y=-1*log10(adj.P.Val), color=clr)) +
        geom_point() +
        scale_color_identity() +
        geom_label_repel(aes(label=lbl), max.overlaps=20, min.segment.length=0, segment.color='black') +
        theme_few() +
        geom_hline(yintercept=-1*log10(min_pval), linetype='dashed', alpha=0.4) +
        geom_vline(xintercept=c(-min_fc, min_fc), linetype='dashed', alpha=0.4) +
        xlab('Enriched in Whole Cell <--               log2(Fold Change)             --> Enriched in Neurites') +
        ylab('-log10(P)')
    } else {
    g.volcano <- ggplot(data=output, aes(x=log2FC, y=-1*log10(adj.P.Val), color=clr)) +
        geom_point() +
        scale_color_identity() +
        geom_label_repel(aes(label=lbl), max.overlaps=20, min.segment.length=0, segment.color='black') +
        theme_few() +
        geom_hline(yintercept=-1*log10(min_pval), linetype='dashed', alpha=0.4) +
        geom_vline(xintercept=c(-min_fc, min_fc), linetype='dashed', alpha=0.4) +
        xlab('Decreased in knockdown <--               log2(Fold Change)             --> Increased in knockdown') +
        ylab('-log10(P)')
    }

    g.MA <- ggplot(data=output, aes(x=AveExpr, y=log2FC, color=clr)) +
        geom_point() +
        scale_color_identity() +
        geom_label_repel(aes(label=lbl), max.overlaps=20, min.segment.length=0, segment.color='black') +
        theme_few() +
        geom_hline(yintercept=c(-min_fc, min_fc), linetype='dashed', alpha=0.4) +
        xlab('log2(Mean Expression)') +
        ylab('log2(Fold Change)')
    
    if(experiment == 'WT') {
        ggsave(g.volcano, file=paste0('DE/', experiment,'-volcano.png'), width=22, height=22, units='cm')
        ggsave(g.volcano, file=paste0('DE/', experiment,'-volcano.pdf'), width=22, height=22, units='cm')
        ggsave(g.MA, file=paste0('DE/', experiment,'-MA.png'), width=22, height=22, units='cm')
        ggsave(g.MA, file=paste0('DE/', experiment,'-MA.pdf'), width=22, height=22, units='cm')
    } else {
        ggsave(g.volcano, file=paste0('DE/', experiment,'-', fractions, '-volcano.png'), width=22, height=22, units='cm')
        ggsave(g.volcano, file=paste0('DE/', experiment,'-', fractions, '-volcano.pdf'), width=22, height=22, units='cm')
        ggsave(g.MA, file=paste0('DE/', experiment,'-', fractions, '-MA.png'), width=22, height=22, units='cm')
        ggsave(g.MA, file=paste0('DE/', experiment,'-', fractions, '-MA.pdf'), width=22, height=22, units='cm')
    }

    return(list(output[], g.volcano, g.MA))
}



# List of genes to pass as manually labeled
tolabel <- c('HNRNPA1', 'TARDBP', 'FUS','YBX1', 'RPS2','UBQLN1','UBE2V1')

# Wild type neurite vs soma
o <- getContrast(countData, metaData, experiment='WT', fractions=c('Neurite','Soma'), contrast_baseline='Soma', annotations=ANNO, manuallabel=tolabel)


# TDP-43 knockdown vs control
o <- getContrast(countData, metaData, experiment='TDP-43', fractions=c('Neurite'), contrast_baseline='nt_control', annotations=ANNO, manuallabel=tolabel)
o <- getContrast(countData, metaData, experiment='TDP-43', fractions=c('Soma'), contrast_baseline='nt_control', annotations=ANNO, manuallabel=tolabel)

# hnRNPA1 knockdown vs control
o <- getContrast(countData, metaData, experiment='hnRNPA1', fractions=c('Neurite'), contrast_baseline='nt_control', annotations=ANNO, manuallabel=tolabel)
o <- getContrast(countData, metaData, experiment='hnRNPA1', fractions=c('Soma'), contrast_baseline='nt_control', annotations=ANNO, manuallabel=tolabel)

# FUS knockdown vs control
o <- getContrast(countData, metaData, experiment='FUS', fractions=c('Neurite'), contrast_baseline='nt_control', annotations=ANNO, manuallabel=tolabel)
o <- getContrast(countData, metaData, experiment='FUS', fractions=c('Soma'), contrast_baseline='nt_control', annotations=ANNO, manuallabel=tolabel)






























































