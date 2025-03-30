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
library(ggrastr)

excluded_samples <- 'TDP43.20'   # Exclude due to PCA outlier

# Import Count Data
if (! file.exists('RNA_DE/raw_counts.tsv')) {
    # Get list of all featurecounts files
    files <- list.files('RNA_alignment',
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
    fwrite(countData, file='RNA_DE/raw_counts.tsv', sep='\t', quote=F, row.names=F, col.names=T)
} else {
    # Read pre-saved counts data
    countData <- fread("RNA_DE/raw_counts.tsv")
}

# Import metadata
metaData <- fread('metadata/rna-metadata.csv')

# Drop excluded samples
countData[, (excluded_samples) := NULL]
metaData <- metaData[! Sample %in% excluded_samples]

# get annotations
ANNO <- fread('https://gist.github.com/cory-weller/566a8984b713fc4626ec423913fdab43/raw/90e5eeaace6ee4d6900a404cbd50c952015b40dc/gene_table_GRCh38-2020-A.tsv')
ANNO <- ANNO[, .SD, .SDcols=c('ENSEMBL','SYMBOL')]
ANNO <- ANNO[!duplicated(ENSEMBL)]
ANNO[SYMBOL=='', SYMBOL := NA]

plot_PC <- function(countData, metaData, experimentName, fraction, contrast) {
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
    
    # calculate the variance for each gene
    rv <- rowVars(assay(vsd))

    # Top n genes by variance to keep.
    ntop <- 500

    # select the ntop genes by variance
    select <- order(rv, decreasing=TRUE)[seq_len(min(ntop, length(rv)))]

    # perform a PCA on the data in assay(x) for the selected genes
    pca <- prcomp(t(assay(vsd)[select,]), scale. = TRUE)
    
    pc1.var <- round(summary(pca)$importance[2,1] * 100, digits=1)
    pc2.var <- round(summary(pca)$importance[2,2] * 100, digits=1)

    pca <- as.data.table(pca$x, keep.rownames=TRUE)
    setnames(pca, 'rn', 'Sample')
    pca <- pca[, .SD, .SDcols=c('Sample','PC1','PC2')] 
    pca <- merge(pca, metadata, by='Sample')
    pca[Fraction=='Soma', Fraction := 'Whole cell']
    pca[ Contrast == 'knockdown', shape_grp := paste0(Fraction, ' KD')]
    pca[ Contrast == 'nt_control', shape_grp := paste0(Fraction, ' NT')]
    pca[ Contrast == 'Neurite', shape_grp := 'Neurite']
    pca[ Contrast == 'Soma', shape_grp := 'Whole cell']
    shape_mapping <- c('Neurite NT' = 21, 'Neurite KD' = 16, 'Whole cell NT'=22, 'Whole cell KD'=15, 'Whole cell'=15, 'Neurite'=16)
    color_mapping <- c('Neurite'=neurite_color, 'Whole cell'=soma_color)

    g <- ggplot(data=pca, aes(x=PC1, y=PC2, shape=shape_grp, color=Fraction)) +
        geom_point(size=4) +
        scale_shape_manual(values=shape_mapping) +
        scale_color_manual(values=color_mapping) +
        labs(title=experimentName, x=paste0('PC1: ', pc1.var,'%'), y=paste0('PC2: ', pc2.var,'%')) +
        theme_few() + theme(aspect.ratio = 1)
    
    ggsave(g, file=paste0('PCA/', experimentName, '-combined-PCA.png'), width=15, height=15, units='cm')
    ggsave(g, file=paste0('PCA/', experimentName, '-combined-PCA.pdf'), width=15, height=15, units='cm')
}


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
    
    # Disable PCA here, calculate as whole group in plot_PC instead
    if(FALSE) {
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
    }

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
    out[is.na(SYMBOL), SYMBOL := Geneid]
    out[SYMBOL=='', SYMBOL := Geneid]
    return(out)
}


plot_DE <- function(DT, down_color, up_color, NS_color='gray', min_fc=1, p_threshold=0.01, color_ENSG=NA, label_ENSG=NA) {
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
    
    
    DT[Geneid %in% label_ENSG, 'text_label' := SYMBOL]
    
    
    # Set point colors for significant up/down, on middle layer
    DT[Significance == 'Up', 'plot_color' := up_color]
    DT[Significance == 'Up', 'plot_layer' := 'middle']
    DT[Significance == 'Down', 'plot_color' := down_color]
    DT[Significance == 'Down', 'plot_layer' := 'middle']

    # Override significance labeling for manually labeled genes, on top layer
    DT[Geneid %in% color_ENSG, 'plot_color' := manual_color]    
    DT[Geneid %in% color_ENSG, 'plot_layer' := 'top']    # Overrides significance labeling
    
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
        rasterize(geom_point(data=DT[plot_layer=='bottom'], shape=21, alpha=0.8), dpi=300) +
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
        rasterize(geom_point(data=DT[plot_layer=='bottom'], shape=21, alpha=0.8), dpi=300) +
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

source('scripts/colors.R')



save_outputs <- function(O, EXPERIMENT, FRACTION, CONTRAST) {
    # Get file stem for all output
    if(EXPERIMENT == 'WT') {
        base_name <- paste0('RNA_DE/', EXPERIMENT, '-', paste0(FRACTION, collapse='-'))
    } else {
        base_name <- paste0('RNA_DE/', EXPERIMENT, '-', FRACTION)
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

# WT ALL PCA

EXPERIMENT <- 'WT'
FRACTION <- c('Soma','Neurite')
CONTRAST <- c('Soma','Neurite')
plot_PC(countData, metaData, experimentName=EXPERIMENT, fraction=FRACTION, contrast=CONTRAST)


EXPERIMENT <- 'WT'
FRACTION <- c('Soma','Neurite')
CONTRAST <- c('Soma','Neurite')
# YBX1,RPS2,UBQLN1,UBE2V1
LABEL <- c('ENSG00000065978','ENSG00000140988','ENSG00000135018','ENSG00000244687')

DEstats <- get_DE_table(countData, metaData, experimentName=EXPERIMENT, fraction=FRACTION, contrast=CONTRAST)
o <- plot_DE(DEstats, up_color=neurite_color, down_color=soma_color, NS_color=not_signif_color, label_ENSG=LABEL)
save_outputs(o, EXPERIMENT, FRACTION, CONTRAST)

# TDP-43 ALL PCA
EXPERIMENT <- 'TDP-43'
FRACTION <- c('Soma','Neurite')
CONTRAST <- c('nt_control','knockdown')
plot_PC(countData, metaData, experimentName=EXPERIMENT, fraction=FRACTION, contrast=CONTRAST)

####################################################################################################
# TDP-43 knockdown vs Control
####################################################################################################
majiq.tdp <- fread('/data/CARD_ARDIS/2023_07_21_veronica_analysis/rnaseq/data/majiq_tdp_kd_whole_cell.csv', select='gene_id')$gene_id
majiq.tdp <- gsub('\\..+$','', majiq.tdp)
majiq.tdp <- unique(majiq.tdp)
COLOR <- majiq.tdp
LABEL <- c('ENSG00000120948', 'ENSG00000104435')

EXPERIMENT <- 'TDP-43'
FRACTION <- c('Soma')
CONTRAST <- c('nt_control','knockdown')
DEstats <- get_DE_table(countData, metaData, experimentName=EXPERIMENT, fraction=FRACTION, contrast=CONTRAST)
o <- plot_DE(DEstats, up_color=soma_color, down_color=soma_down_color, NS_color=not_signif_color, color_ENSG=COLOR, label_ENSG=LABEL)
save_outputs(o, EXPERIMENT, FRACTION, CONTRAST)

EXPERIMENT <- 'TDP-43'
FRACTION <- c('Neurite')
CONTRAST <- c('nt_control','knockdown')
DEstats <- get_DE_table(countData, metaData, experimentName=EXPERIMENT, fraction=FRACTION, contrast=CONTRAST)
o <- plot_DE(DEstats, up_color=neurite_color, down_color=neurite_down_color, NS_color=not_signif_color, color_ENSG=COLOR, label_ENSG=LABEL)
save_outputs(o, EXPERIMENT, FRACTION, CONTRAST)


# hnRNPA1 ALL PCA
EXPERIMENT <- 'hnRNPA1'
FRACTION <- c('Soma','Neurite')
CONTRAST <- c('nt_control','knockdown')
plot_PC(countData, metaData, experimentName=EXPERIMENT, fraction=FRACTION, contrast=CONTRAST)


####################################################################################################
# hnRNPA1 knockdown vs Control
####################################################################################################
majiq.hnrnpa1 <- fread('/data/CARD_ARDIS/2023_07_21_veronica_analysis/rnaseq/data/majiq_hnrnpa1_kd_whole_cell.csv', select='gene_id')$gene_id
majiq.hnrnpa1 <- gsub('\\..+$','', majiq.hnrnpa1)
majiq.hnrnpa1 <- unique(majiq.hnrnpa1)
COLOR <- majiq.hnrnpa1
LABEL <- c('ENSG00000135486')

EXPERIMENT <- 'hnRNPA1'
FRACTION <- c('Soma')
CONTRAST <- c('nt_control','knockdown')
DEstats <- get_DE_table(countData, metaData, experimentName=EXPERIMENT, fraction=FRACTION, contrast=CONTRAST)
o <- plot_DE(DEstats, up_color=soma_color, down_color=soma_down_color, NS_color=not_signif_color, color_ENSG=COLOR, label_ENSG=LABEL)
save_outputs(o, EXPERIMENT, FRACTION, CONTRAST)

EXPERIMENT <- 'hnRNPA1'
FRACTION <- c('Neurite')
CONTRAST <- c('nt_control','knockdown')
DEstats <- get_DE_table(countData, metaData, experimentName=EXPERIMENT, fraction=FRACTION, contrast=CONTRAST)
o <- plot_DE(DEstats, up_color=neurite_color, down_color=neurite_down_color, NS_color=not_signif_color, color_ENSG=COLOR, label_ENSG=LABEL)
save_outputs(o, EXPERIMENT, FRACTION, CONTRAST)


quit(status=0)

