#!/usr/bin/env Rscript

if (! require('ggplot2'))        { install.packages('ggplot2'); library(ggplot2) }
if (! require('ggthemes'))        { install.packages('ggthemes'); library(ggthemes) }
if (! require('data.table'))        { install.packages('data.table'); library(data.table) }
if (! require('cowplot'))        { install.packages('cowplot'); library(cowplot) }

facet_order <- c(
'REPAC',
"APAlyzer",
"APA-Scan"
)

ipa_plot <- "APAlyzer Intron Polyadenylation"


# APAlyzer 3UTR
dat1 <- fread('APA/APAlyzer-3UTR.csv')
dat1 <- dat1[, .SD, .SDcols=c('gene_symbol','RED','p_adj','Significance')]
dat1[, 'Measure' := facet_order[2]]
dat1[, 'p_adj' := NULL]
setnames(dat1, c('site','Value','Significance', 'Measure'))
dat1[Significance == 'Shortened in Neurites', Significance := 'Shortening in Neurites']
dat1[Significance == 'Lengthened in Neurites', Significance := 'Lengthening in Neurites']

# APAlyzer IPA
dat2 <- fread('APA/APAlyzer-IPA.csv')
dat2 <- dat2[, .SD, .SDcols=c('PASid','RED','p_adj','Significance')]
dat2[, 'Measure' := ipa_plot]
dat2[, p_adj := NULL]
setnames(dat2, c('site','Value','Significance', 'Measure'))


# REPAC
dat3 <- fread('APA/REPAC-gene-level-significance.csv')
setnames(dat3, 'gene_name', 'site')
dat3 <- dat3[, .SD, .SDcols=c('site','cFC','Significance')]
dat3[, 'Measure' := facet_order[1]]
setnames(dat3, c('site','Value','Significance', 'Measure'))

dat3[Significance == 'Shortened in Neurites', Significance := 'Shortening in Neurites']
dat3[Significance == 'Lengthened in Neurites', Significance := 'Lengthening in Neurites']

# APA-Scan
dat4 <- fread('APA/APA-Scan/output/APA_Scan_soma-links_Vs_neurite-links.csv')
dat4 <- dat4[, .SD, .SDcols=c('Gene','Ratio Difference','p-value')]
setnames(dat4, c('Gene','Ratio','p'))

dat4[p < 0.05 & Ratio > 0, Significance := 'Shortening in Neurites']
dat4[p < 0.05 & Ratio < 0, Significance := 'Lengthening in Neurites']
dat4 <- dat4[! is.na(Significance)]
dat4[, p := NULL]
setnames(dat4, 'Gene','site')
setnames(dat4, 'Ratio','Value')
dat4[, Measure := facet_order[3]]

dat <- rbindlist(list(dat1, dat3, dat4))
dat[, Measure := factor(Measure, levels=rev(facet_order))]

dat[, Significance := factor(Significance, levels=rev(c('NS','Shortening in Neurites','Lengthening in Neurites')))]
dat <- dat[Significance %in% c('Shortening in Neurites','Lengthening in Neurites')]
dat.ag <- dat[, .N, by=list(Significance, Measure)]

source('scripts/colors.R')
colors <- c('Shortening in Neurites'=neurite_down_color,
            'Lengthening in Neurites'=neurite_color)

g1 <- ggplot(data=dat.ag, aes(x=Measure, y=N, fill=Significance)) +
    geom_bar(stat='identity', position='stack') +
    theme_few() +
    labs(x='', y='Number of Genes') +
    scale_fill_manual(values=colors) +
    theme(legend.position="bottom",legend.direction="vertical")


ggsave(g1, file='APA/APA-3methods-bar.png', width=16, height=8, units='cm')
ggsave(g1, file='APA/APA-3methods-bar.pdf', width=12, height=8, units='cm')

dat2 <- dat2[, .N, by=list(Significance,Measure)]
dat2 <- dat2[Significance !='NS']
colors2 <- c('IPA activation in Neurites'=neurite_down_color,
            'IPA suppression in Neurites'=neurite_color)

g2 <- ggplot(data=dat2, aes(x=Measure, y=N, fill=Significance)) +
    geom_bar(stat='identity', position='stack') +
    theme_few() +
    labs(x='', y='Number of Genes') +
    scale_fill_manual(values=colors2) +
    theme(legend.position="bottom",legend.direction="vertical")

ggsave(g2, file='APA/APA-IPA-bar.png', width=8, height=8, units='cm')
ggsave(g2, file='APA/APA-IPA-bar.pdf', width=8, height=8, units='cm')

