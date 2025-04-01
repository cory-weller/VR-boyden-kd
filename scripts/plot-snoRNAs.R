#!/usr/bin/env Rscript

library(data.table)
library(ggplot2)
library(ggthemes)
library(viridis)
library(pheatmap)
library(foreach)

counts <- fread('/vf/users/CARD_ARDIS/users/wellerca/VR-boyden-kd/RNA_DE/raw_counts.tsv')
counts[, 'Geneid' := gsub('\\.[0-9]+$', '', Geneid)]
metadata <- fread('/vf/users/CARD_ARDIS/users/wellerca/VR-boyden-kd/metadata/rna-metadata.csv')
samplecols <- colnames(counts)[-1]


convert_log_cpm <- function(v) {
    v <- 1e6 * v/sum(v)
    v <- log10(v + 1)
    return(v)
}

library("AnnotationDbi")
library("org.Hs.eg.db")

ENSG_to_SYMBOL <- as.data.table(as.data.frame(mapIds(org.Hs.eg.db,
                    keys=counts$Geneid, 
                    column="SYMBOL",
                    keytype="ENSEMBL",
                    multiVals="first") ), keep.rownames=T)


setnames(ENSG_to_SYMBOL, c('Geneid','SYMBOL'))
ENSG_to_SYMBOL <- ENSG_to_SYMBOL[!is.na(SYMBOL)]

# Convert count data to log(counts per million + 1)
counts <- cbind(counts[,1], counts[, lapply(.SD, convert_log_cpm), .SDcols=samplecols])
counts <- merge(counts, ENSG_to_SYMBOL, by.x='Geneid', by.y='Geneid')
counts <- counts[!is.na(SYMBOL)]
counts <- counts[, lapply(.SD, sum), by=SYMBOL, .SDcols=samplecols]


counts.long <- melt(counts, measure.vars=samplecols, variable.name = 'Sample', value.name='log_cpm')
counts.long <- merge(counts.long, metadata, by.x='Sample', by.y='Sample')

# GOOD

snoRNAs <- readLines('snoRNA-SYMBOLS.txt')
o <- counts.long[SYMBOL %in% snoRNAs]



o[Fraction %like% 'Neurite', zone := 'Neurite']
o[Fraction %like% 'Soma', zone := 'Whole cell']
o[, Contrast := gsub('nt_control', 'NT', Contrast)]
o[, Contrast := gsub('knockdown', 'KD', Contrast)]
o[Experiment == 'WT', sample := paste0(zone, ' ', Replicate)]
o[Experiment != 'WT', sample := paste0(zone, ' ', Contrast, ' ', Replicate)]


# Ensure all rows have Neurite or Soma labeled zone
stopifnot(nrow(o[is.na(zone)]) == 0)


gene_order <- rev(sort(unique(o$SYMBOL)))
o[, SYMBOL := factor(SYMBOL, levels=gene_order)]


# Exclude sample
o <- o[Sample != 'TDP43.20']

# Add in combinations that don't exist, give NA values, so levels aren't dropped
all_levels <- CJ('SYMBOL'=sort(unique(o$SYMBOL)),
    'Contrast'=sort(unique(o$Contrast)),
    'zone'=sort(unique(o$zone)),
    'Experiment'=sort(unique(o$Experiment)))



# setkey(all_levels, SYMBOL, Contrast, zone, Experiment)
# setkey(o, SYMBOL, Contrast, zone, Experiment)

# hnRNPA1
dt <- o[Experiment=='hnRNPA1']
g1 <- ggplot(dt, aes(x=sample, y=SYMBOL, fill=log_cpm)) +
        geom_tile() +
        scale_fill_viridis(limits=c(0,max(o$log_cpm)), option='magma') +
        theme_few() +
        labs(y='snoRNA') +
        theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

ggsave(g1, file='hnRNPA1-snoRNAs.png', width=16, height=46, units='cm')
ggsave(g1, file='hnRNPA1-snoRNAs.pdf', width=16, height=46, units='cm')



# TDP-43
dt <- o[Experiment=='TDP-43']
g2 <- ggplot(dt, aes(x=sample, y=SYMBOL, fill=log_cpm)) +
        geom_tile() +
        scale_fill_viridis(limits=c(0,max(o$log_cpm)), option='magma') +
        theme_few() +
        labs(y='snoRNA') +
        theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

ggsave(g2, file='TDP-43-snoRNAs.png', width=16, height=46, units='cm')
ggsave(g2, file='TDP-43-snoRNAs.pdf', width=16, height=46, units='cm')


# WT

# TDP-43
dt <- o[Experiment=='WT']
g3 <- ggplot(dt, aes(x=sample, y=SYMBOL, fill=log_cpm)) +
        geom_tile() +
        scale_fill_viridis(limits=c(0,max(o$log_cpm)), option='magma') +
        theme_few() +
        labs(y='snoRNA') +
        theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

ggsave(g3, file='WT-snoRNAs.png', width=10, height=46, units='cm')
ggsave(g3, file='WT-snoRNAs.pdf', width=10, height=46, units='cm')