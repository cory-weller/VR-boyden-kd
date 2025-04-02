#!/usr/bin/env Rscript

library(ggplot2)
library(ggthemes)
library(data.table)
library(ggrepel)

dat <- fread('APA/APA-Scan/output/APA_Scan_soma-links_Vs_neurite-links.csv')

dat <- dat[, .SD, .SDcols=c('Gene','Ratio Difference','p-value')]
setnames(dat, c('Gene','Ratio','p'))
source('scripts/colors.R')

colors <- c('Shortening in Neurites'=neurite_down_color, 
            'Lengthening in Neurites'=neurite_color,
            'NS'=not_signif_color)

dat[p < 0.05 & Ratio > 0, Significance := 'Shortening in Neurites']
dat[p < 0.05 & Ratio < 0, Significance := 'Lengthening in Neurites']
dat[p >= 0.05, Significance := 'NS']
dat[Ratio ==0, Significance := 'NS']



dat[Significance=='Lengthening in Neurites', lbl := Gene]
dat[Significance=='Shortening in Neurites' & Ratio > 0.4, lbl := Gene]
dat[Significance=='Shortening in Neurites' & p < 0.00005, lbl := Gene]
dat[Gene=='HNRNPA1', lbl := Gene]
dat[Gene=='NEFL', lbl := Gene]

dat[, Ratio := -1 * Ratio]
g <- ggplot(dat, aes(x=Ratio, y=-log10(p), color=Significance, label=lbl)) +
    geom_hline(yintercept=-log10(0.05), linetype='dashed', alpha=0.5) +
    geom_point() +
    geom_label_repel(max.overlaps=20, min.segment.length=0, segment.color='black', show.legend=F) +
    scale_color_manual(values=colors) +
    theme_few() +
    theme(legend.position="bottom",legend.direction="vertical")


ggsave(g, file='APA/APA-scan-volcano.png', width=14, height=16, units='cm')
ggsave(g, file='APA/APA-scan-volcano.pdf', width=14, height=16, units='cm')