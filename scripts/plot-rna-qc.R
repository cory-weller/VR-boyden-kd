#!/usr/bin/env Rscript

library(data.table)
dat <- fread('RNA_alignment/multiqc_data/multiqc_rna_seqc.txt')
dat[, 2 := NULL]


mdata <- fread('metadata/rna-metadata.csv')


dat <- merge(mdata, dat, by.x='filestem', by.y='Sample')



dat <- dat[, .SD, .SDcols=c('Sample','Experiment','Fraction','Contrast','Replicate','Genes Detected','rRNA rate','Expression Profiling Efficiency')]

library(ggplot2)
library(ggthemes)

dat[, 'rn' := paste0(Fraction, ' ', Contrast, ' ', Replicate)]
dat[Experiment == 'WT', 'rn' := paste0(Experiment, ' ', Fraction, ' ', Replicate)]
dat[, 'rn' := gsub('nt_control', 'NT', rn)]
dat[, 'rn' := gsub('knockdown', 'KD', rn)]
dat <- dat[Experiment != 'FUS']
setkey(dat, Experiment, Contrast, Fraction)
dat[, rn := gsub('Soma', 'Whole cell', rn)]

dat.long <- melt(dat, id.vars=c('Experiment','rn'), measure.vars=c('rRNA rate','Expression Profiling Efficiency','Genes Detected'))
dat.long[variable  == 'rRNA rate' & rn %like% 'Neurite', 'plot_color' := neurite_down_color]
dat.long[variable  == 'rRNA rate' & rn %like% 'Whole cell', 'plot_color' := soma_down_color]
dat.long[variable  == 'Expression Profiling Efficiency' & rn %like% 'Neurite', 'plot_color' := neurite_color]
dat.long[variable  == 'Expression Profiling Efficiency' & rn %like% 'Whole cell', 'plot_color' := soma_color]
dat.long[variable  == 'Genes Detected' & rn %like% 'Neurite', 'plot_color' := neurite_color]
dat.long[variable  == 'Genes Detected' & rn %like% 'Whole cell', 'plot_color' := soma_color]
# Expression efficiency

source('scripts/colors.R')

plot_qc1 <- function(DT, experiment) {
    dt <- DT[Experiment == experiment]
    sample_order <- dt$rn
    dt <- dt[variable %in% c('rRNA rate','Expression Profiling Efficiency')]
    g <- ggplot(data=dt, aes(x=value*100, y=rn, group=variable, fill=plot_color)) +
        geom_bar(stat='identity', position='dodge') +
        scale_fill_identity() +
        xlim(0,100) +
        labs(title=experiment, x='Percent', y='Sample') +
        theme_few()
    return(g)
}

plot_qc2 <- function(DT, experiment) {
    dt <- DT[Experiment == experiment]
    sample_order <- dt$rn
    dt <- dt[variable=='Genes Detected']
    g <- ggplot(data=dt, aes(x=value, y=rn, fill=plot_color)) +
        geom_bar(stat='identity', position='dodge') +
        scale_fill_identity() +
        xlim(0,30000) +
        labs(title=experiment, x='Genes Detected', y='Sample') +
        theme_few()
    return(g)
}


g1 <- plot_qc1(dat.long, 'TDP-43')
g2 <- plot_qc1(dat.long, 'hnRNPA1')

g3 <- plot_qc2(dat.long, 'TDP-43')
g4 <- plot_qc2(dat.long, 'hnRNPA1')

ggsave(g1, file='RNA_QC/TDP-43-efficiency-rrna.png', width=12, height=12, units='cm')
ggsave(g2, file='RNA_QC/hnRNPA1-efficiency-rrna.png', width=12, height=12, units='cm')
ggsave(g3, file='RNA_QC/TDP-43-genes-detected.png', width=12, height=12, units='cm')
ggsave(g4, file='RNA_QC/hnRNPA1-genes-detected.png', width=12, height=12, units='cm')

fwrite(dat, file='RNA_QC/RNA-metrics-wide.csv', sep=',', quote=F, row.names=F, col.names=T)
fwrite(dat.long, file='RNA_QC/RNA-metrics-long.csv', sep=',', quote=F, row.names=F, col.names=T)
