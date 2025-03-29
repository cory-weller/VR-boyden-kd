#!/usr/bin/env Rscript

# #load libraries ----
# library(tidyverse) #42; basis for data manipulation
# library(ggpubr) #publication ready plots
# library(pheatmap) #heatmap
# library(GOfuncR) #pathway analysis
# library(EnhancedVolcano)
# library(GGally)
# library(devtools) #load R scripts from github




library(data.table)
library(foreach)
library(ggrepel)
library(ggthemes)
library(DESeq2) #for normalizing counts

excluded_samples <- 'TDP43.20'   # Exclude due to PCA outlier

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
    
    # Exclude PCA outlier TDP43.20
    
    
    # Remove '.dedup.bam' from column names
    setnames(countData, gsub('.dedup.bam', '', colnames(countData)))
    
    # Modify names to match experiment
    setnames(countData, gsub('VR-8153-', 'FUS-', colnames(countData)))
    setnames(countData, gsub('JR-7961-', 'TDP43-', colnames(countData)))
    setnames(countData, gsub('VR-8147-', 'hnRNPA1-', colnames(countData)))
    setnames(countData, gsub('-', '.', colnames(countData)))
    
    # Ensure directories exists
    dir.create("DE")
    
    # Save merged counts as tabular dat
    fwrite(countData, file='DE/raw_counts.tsv', sep='\t', quote=F, row.names=F, col.names=T)
} else {
    # Read pre-saved counts data
    countData <- fread("DE/raw_counts.tsv")
}

# Import metadata
metaData <- fread('metadata/rnaseq-metadata.csv')

# Drop excluded samples
countData[, (excluded_samples) := NULL]
metaData <- metaData[! Sample %in% excluded_samples]

# get annotations
ANNO <- fread('https://gist.github.com/cory-weller/566a8984b713fc4626ec423913fdab43/raw/90e5eeaace6ee4d6900a404cbd50c952015b40dc/gene_table_GRCh38-2020-A.tsv')
ANNO <- ANNO[, .SD, .SDcols=c('ENSEMBL','SYMBOL')]
ANNO <- ANNO[!duplicated(ENSEMBL)]



get_DE_table <- function(countData, metaData, experimentName, fraction, contrast) {
   counts <- copy(countData)
   metadata <- copy(metaData)
   metadata <- metadata[Experiment == experimentName & Fraction %in% fraction & Contrast %in% contrast][order(Contrast)]
   metadata[, Contrast := factor(Contrast, levels=contrast)]
   setkey(metadata, Contrast, Sample)
   counts <- counts[, .SD, .SDcols=c('Geneid', metadata$Sample)]
   if(! identical(metadata$Sample, colnames(counts)[-1])) {
      cat('Metadata and sample columns are not identical. Exiting...\n')
      stop()
   }
   # Build colData df
   colData <- metadata[, .SD, .SDcols=c('Experiment','Fraction','Sample','Contrast')]
   
   smallestGroupSize <- min(rle(as.numeric(colData$Contrast))$lengths)
   samples <- colData$Sample
   colData <- as.data.frame(colData)
   rownames(colData) <- samples
   
   # Build counts df
   genes <- counts$Geneid
   counts <- as.data.frame(counts[, lapply(.SD, round, digits = 0),  .SDcols = samples])
   rownames(counts) <- genes
   
   # Create DESeq object
    dds <- DESeqDataSetFromMatrix(countData = counts, colData = colData, design = ~Contrast)
    
    # Filter out low count genes
    keep <- rowSums(counts(dds) >= 10) >= smallestGroupSize
    dds <- dds[keep,]
    
    # Relevel with respect to first name in the 'contrast' argument
    dds$Contrast <- relevel(dds$Contrast, ref = contrast[1])
    
    # Run DESeq
    dds <- DESeq(dds)
    vsd <- vst(dds)
    g.pca <- DESeq2::plotPCA(vsd, intgroup = c('Contrast'))
    g.pca <- g.pca + theme_few() + theme(aspect.ratio=1)
    
    if(experimentName == 'WT') {
        base_name <- paste0('PCA/', experimentName, '-', paste0(fraction, collapse='-'))
        g.pca <- g.pca + labs(title=experimentName)
    } else {
        base_name <- paste0('PCA/', experimentName, '-', fraction)
        g.pca <- g.pca + labs(title=paste0(experimentName, ' ', fraction))
    }
    dir.create("PCA")

    ggsave(g.pca, file=paste0(base_name, '-PCA.png'), width=15, height=15, units='cm')
    ggsave(g.pca, file=paste0(base_name, '-PCA.pdf'), width=15, height=15, units='cm')


    # Run log(foldChange) shrinkage
    coef_name <- names(coef(dds)[1,])[2]
    res <- lfcShrink(dds, coef=coef_name, type="apeglm")
    
    # Convert to data.table
    dat <- as.data.table(as.data.frame(res), keep.rownames=TRUE)
    setnames(dat, 'rn', 'Geneid')

    dat <- add_Anno(dat)[]
    dat[is.na(padj), padj := 1]
   return(dat)
}

add_Anno <- function(DE.DT) {
    DE.DT[, Geneid := gsub('\\.[0-9]+$','',Geneid)]
    #return(DE.DT)
    out <- merge(DE.DT, ANNO, by.x='Geneid', by.y='ENSEMBL', all.x=T)
    out[is.na(SYMBOL), SYMBOL := Geneid][]
    return(out)
}


plot_DE <- function(DT, down_color, up_color, NS_color='gray', min_fc=1, p_threshold=0.01, label_ENSG=NA) {
    # Mark Non-significant by p-value
    DT[padj >= p_threshold, Significance := 'NS']
    DT[is.na(padj), Significance := 'NS']
    
    # Mark Non-significant by LFC
    DT[log2FoldChange %between% c(-min_fc, min_fc), Significance := 'NS']
    
    DT[padj < p_threshold & log2FoldChange > min_fc, Significance := 'Up']
    DT[padj < p_threshold & log2FoldChange < -min_fc, Significance := 'Down']
    
    # Modify colors
    
    # Add manually labeled genes
    
    
    # If WT, add text labels to manually labeled gens
    if(length(label_ENSG)==4) {
        DT[Geneid %in% label_ENSG, 'text_label' := SYMBOL]
    }
    
    # Set point colors for significant up/down, on middle layer
    DT[Significance == 'Up', 'plot_color' := up_color]
    DT[Significance == 'Up', 'plot_layer' := 'middle']
    DT[Significance == 'Down', 'plot_color' := down_color]
    DT[Significance == 'Down', 'plot_layer' := 'middle']

    # Override significance labeling for manually labeled genes, on top layer
    DT[Geneid %in% label_ENSG, 'plot_color' := manual_color]    
    DT[Geneid %in% label_ENSG, 'plot_layer' := 'top']    # Overrides significance labeling
    
    #  Remaining (insignificant) points are NS_color, on bottom layer
    DT[is.na(plot_color), plot_color := NS_color]
    DT[plot_color == NS_color, 'plot_layer' := 'bottom']
    if(any(is.na(DT$plot_layer))) { cat('NAs in plot_layer !!\n'); stop() }
    if(any(is.na(DT$plot_color))) { cat('NAs in plot_color !!\n'); stop() }
    
    
    # Add top 10 up and 10 down
    DT[Significance=='Down', down_rank := frank(log2FoldChange)]
    DT[Significance=='Up', up_rank := frank(-log2FoldChange)]
    DT[down_rank <= 10 | up_rank <= 10, 'text_label' := SYMBOL]
    
    # Manually label WT genes
    
    g.volcano <- ggplot(DT, aes(x=log2FoldChange, y=-1*log10(padj), color=plot_color, fill=plot_color)) +
        geom_point(data=DT[plot_layer=='bottom'], shape=21, alpha=0.8) +
        geom_point(data=DT[plot_layer=='middle'], shape=21, alpha=0.8) +
        geom_point(data=DT[plot_layer=='top'], shape=21, alpha=0.8) +
        scale_color_identity() +
        scale_fill_identity() +
        geom_label_repel(aes(label=text_label), max.overlaps=9999, min.segment.length=0, segment.color='black', color='black', fill='white') +
        theme_few() +
        geom_hline(yintercept=-1*log10(p_threshold), linetype='dashed', alpha=0.4) +
        geom_vline(xintercept=c(-min_fc, min_fc), linetype='dashed', alpha=0.4) +
        xlab('log2(FoldChange)') +
        ylab('log10(P)')
    
    g.MA <- ggplot(DT, aes(x=log2(baseMean), y=log2FoldChange, color=plot_color, fill=plot_color)) +
        geom_point(data=DT[plot_layer=='bottom'], shape=21, alpha=0.8) +
        geom_point(data=DT[plot_layer=='middle'], shape=21, alpha=0.8) +
        geom_point(data=DT[plot_layer=='top'], shape=21, alpha=0.8) +
        scale_color_identity() +
        scale_fill_identity() +
        geom_label_repel(aes(label=text_label), max.overlaps=9999, min.segment.length=0, segment.color='black', color='black', fill='white') +
        theme_few() +
        geom_hline(yintercept=c(-min_fc, min_fc), linetype='dashed', alpha=0.4) +
        xlab('log2(baseMean)') +
        ylab('log2(FoldChange)')
    
        return(list(dt=DT, volcano=g.volcano, MA=g.MA))
}

neurite_color <- '#d63cc0'      # Green
soma_color <- '#00a347'         # Purple
neurite_down_color <- '#6c1760' # Dark green
soma_down_color <- '#003f1e'    # Dark purple
not_signif_color <- '#c4c4c4'   # Gray
manual_color <- '#2e47ff'      # Blue

save_outputs <- function(O, EXPERIMENT, FRACTION, CONTRAST) {
    # Get file stem for all output
    if(EXPERIMENT == 'WT') {
        base_name <- paste0('DE/', EXPERIMENT, '-', paste0(FRACTION, collapse='-'))
    } else {
        base_name <- paste0('DE/', EXPERIMENT, '-', FRACTION)
    }
    # Save DE table
    fwrite(O$dt, file=paste0(base_name, '.csv'), quote=F, row.names=F, col.names=T, sep=',')
    
    # Save volcano plot as png and pdf
    ggsave(O$volcano, file=paste0(base_name, '-volcano.png'))
    ggsave(O$volcano, file=paste0(base_name, '-volcano.pdf'))
    
    # Save MA plot as png and pdf
    ggsave(O$MA, file=paste0(base_name, '-MA.png'))
    ggsave(O$MA, file=paste0(base_name, '-MA.pdf'))
}

####################################################################################################
# WT, Neurite vs Soma
####################################################################################################
EXPERIMENT <- 'WT'
FRACTION <- c('Soma','Neurite')
CONTRAST <- c('Soma','Neurite')
# YBX1,RPS2,UBQLN1,UBE2V1
MANUALLABEL <- c('ENSG00000065978','ENSG00000140988','ENSG00000135018','ENSG00000244687')

DEstats <- get_DE_table(countData, metaData, experimentName=EXPERIMENT, fraction=FRACTION, contrast=CONTRAST)
o <- plot_DE(DEstats, up_color=neurite_color, down_color=soma_color, NS_color=not_signif_color, label_ENSG=MANUALLABEL)
save_outputs(o, EXPERIMENT, FRACTION, CONTRAST)


####################################################################################################
# TDP-43 knockdown vs Control
####################################################################################################
majiq.tdp <- fread('/data/CARD_ARDIS/2023_07_21_veronica_analysis/rnaseq/data/majiq_tdp_kd_whole_cell.csv', select='gene_id')$gene_id
majiq.tdp <- gsub('\\..+$','', majiq.tdp)
majiq.tdp <- unique(majiq.tdp)
MANUALLABEL <- majiq.tdp

EXPERIMENT <- 'TDP-43'
FRACTION <- c('Soma')
CONTRAST <- c('nt_control','knockdown')
DEstats <- get_DE_table(countData, metaData, experimentName=EXPERIMENT, fraction=FRACTION, contrast=CONTRAST)
o <- plot_DE(DEstats, up_color=soma_color, down_color=soma_down_color, NS_color=not_signif_color, label_ENSG=MANUALLABEL)
save_outputs(o, EXPERIMENT, FRACTION, CONTRAST)

EXPERIMENT <- 'TDP-43'
FRACTION <- c('Neurite')
CONTRAST <- c('nt_control','knockdown')
DEstats <- get_DE_table(countData, metaData, experimentName=EXPERIMENT, fraction=FRACTION, contrast=CONTRAST)
o <- plot_DE(DEstats, up_color=neurite_color, down_color=neurite_down_color, NS_color=not_signif_color, label_ENSG=MANUALLABEL)
save_outputs(o, EXPERIMENT, FRACTION, CONTRAST)


####################################################################################################
# hnRNPA1 knockdown vs Control
####################################################################################################
majiq.hnrnpa1 <- fread('/data/CARD_ARDIS/2023_07_21_veronica_analysis/rnaseq/data/majiq_hnrnpa1_kd_whole_cell.csv', select='gene_id')$gene_id
majiq.hnrnpa1 <- gsub('\\..+$','', majiq.hnrnpa1)
majiq.hnrnpa1 <- unique(majiq.hnrnpa1)
MANUALLABEL <- majiq.hnrnpa1

EXPERIMENT <- 'hnRNPA1'
FRACTION <- c('Soma')
CONTRAST <- c('nt_control','knockdown')
DEstats <- get_DE_table(countData, metaData, experimentName=EXPERIMENT, fraction=FRACTION, contrast=CONTRAST)
o <- plot_DE(DEstats, up_color=soma_color, down_color=soma_down_color, NS_color=not_signif_color, label_ENSG=MANUALLABEL)
save_outputs(o, EXPERIMENT, FRACTION, CONTRAST)

EXPERIMENT <- 'hnRNPA1'
FRACTION <- c('Neurite')
CONTRAST <- c('nt_control','knockdown')
DEstats <- get_DE_table(countData, metaData, experimentName=EXPERIMENT, fraction=FRACTION, contrast=CONTRAST)
o <- plot_DE(DEstats, up_color=neurite_color, down_color=neurite_down_color, NS_color=not_signif_color, label_ENSG=MANUALLABEL)
save_outputs(o, EXPERIMENT, FRACTION, CONTRAST)


quit(status=0)

