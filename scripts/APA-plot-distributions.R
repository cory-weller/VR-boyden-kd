#!/usr/bin/env Rscript

if (! require('ggplot2'))        { install.packages('ggplot2'); library(ggplot2) }
if (! require('ggthemes'))        { install.packages('ggthemes'); library(ggthemes) }
if (! require('data.table'))        { install.packages('data.table'); library(data.table) }
if (! require('cowplot'))        { install.packages('cowplot'); library(cowplot) }

facet_order <- c(
'REPAC\nTranscript length change',
"APAlyzer\nTranscript length change",
"APAlyzer\nIntron Polyadenylation"
)

# APAlyzer 3UTR
dat1 <- fread('APA/APAlyzer-3UTR.csv')
dat1 <- dat1[, .SD, .SDcols=c('gene_symbol','RED','p_adj','Significance')]
dat1[, 'Measure' := facet_order[2]]
dat1[, 'p_adj' := NULL]
setnames(dat1, c('site','Value','Significance', 'Measure'))
dat1[Significance == 'Shortened in Neurites', Significance := 'Decreased in Neurites']
dat1[Significance == 'Lengthened in Neurites', Significance := 'Increased in Neurites']

# APAlyzer IPA
dat2 <- fread('APA/APAlyzer-IPA.csv')
dat2 <- dat2[, .SD, .SDcols=c('PASid','RED','p_adj','Significance')]
dat2[, 'Measure' := facet_order[3]]
dat2[, p_adj := NULL]
setnames(dat2, c('site','Value','Significance', 'Measure'))
dat2[Significance == 'IPA activation in Neurites', Significance := 'Increased in Neurites']
dat2[Significance == 'IPA suppression in Neurites', Significance := 'Decreased in Neurites']


# REPAC
dat3 <- fread('APA/REPAC-gene-level-significance.csv')
setnames(dat3, 'gene_name', 'site')
dat3 <- dat3[, .SD, .SDcols=c('site','cFC','Significance')]
dat3[, 'Measure' := facet_order[1]]
setnames(dat3, c('site','Value','Significance', 'Measure'))

dat3[Significance == 'Shortened in Neurites', Significance := 'Decreased in Neurites']
dat3[Significance == 'Lengthened in Neurites', Significance := 'Increased in Neurites']

dat <- rbindlist(list(dat1, dat2, dat3))
dat[, Measure := factor(Measure, levels=rev(facet_order))]

dat[, Significance := factor(Significance, levels=rev(c('NS','Decreased in Neurites','Increased in Neurites')))]

dat.ag <- dat[, .N, by=list(Significance, Measure)]

g <- ggplot(data=dat.ag, aes(x=Measure, y=N, fill=Significance)) +
    geom_bar(stat='identity', position='stack') +
    theme_few() +
    labs(x='', y='Number of Genes') +
    coord_flip()


ggsave(g, file='APA/APA-stacked-bar.png', width=24, height=8, units='cm')