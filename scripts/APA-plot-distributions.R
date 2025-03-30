#!/usr/bin/env Rscript

if (! require('ggplot2'))        { install.packages('ggplot2'); library(ggplot2) }
if (! require('ggthemes'))        { install.packages('ggthemes'); library(ggthemes) }
if (! require('data.table'))        { install.packages('data.table'); library(data.table) }
if (! require('cowplot'))        { install.packages('cowplot'); library(cowplot) }

facet_order <- c(
"APAlyzer 3'UTR Polyadenylation\nRelative Expression Difference",
"APAlyzer Intron Polyadenylation\nRelative Expression Difference",
'REPAC\nCompositional Fold Change'
)

# APAlyzer 3UTR
dat1 <- fread('APA/APAlyzer-3UTR.csv')
dat1 <- dat1[, .SD, .SDcols=c('gene_symbol','RED','p_adj','Significance')]
dat1[, 'Measure' := facet_order[1]]
setnames(dat1, c('site','Value','P','Significance', 'Measure'))

# APAlyzer IPA
dat2 <- fread('APA/APAlyzer-IPA.csv')
dat2 <- dat2[, .SD, .SDcols=c('PASid','RED','p_adj','Significance')]
dat2[, 'Measure' := facet_order[2]]
setnames(dat2, c('site','Value','P','Significance', 'Measure'))



# REPAC
dat3 <- fread('APA/REPAC-results.csv')
dat3 <- dat3[, .SD, .SDcols=c('Site','cFC','p.adj','Significance')]
dat3[, 'Measure' := facet_order[3]]
setnames(dat3, c('site','Value','P','Significance', 'Measure'))



dat <- rbindlist(list(dat1, dat2, dat3))
dat[, Significance := factor(Significance, levels=c('Whole cell','NS','Neurite'))]
dat[, Measure := factor(Measure, levels=facet_order)]

g <- ggplot(data=dat, aes(x=Significance, y=Value)) +
    facet_grid(~Measure) +
    geom_boxplot() +
    theme_few() +
    geom_hline(yintercept=0, linetype='dashed', alpha=0.5)

ggsave(g, file='APA/APA-boxplot-combined.png', width=24, height=16, units='cm')