#!/usr/bin/env Rscript

#!/usr/bin/env Rscript

# R 4.3.2

library(data.table)
library(ggplot2)
library(ggthemes)
library(cowplot)
library(org.Hs.eg.db)
library(clusterProfiler)
library(viridis)
library(ggthemes)

####################################################################################################
enrich_pvalue <- 0.05

getGO <- function(entrez_list, all_gene_list, enrich_pvalue=0.05) {
    desired_cols <- c('ID','Description','GeneRatio','BgRatio','pvalue','p.adjust','qvalue','geneID','Count','group')
    d1 <- enrichGO(gene          = entrez_list,
                                    universe      = all_gene_list,
                                    OrgDb         = org.Hs.eg.db,
                                    ont           = "CC",
                                    pAdjustMethod = "BH",
                                    pvalueCutoff  = enrich_pvalue,
                                    qvalueCutoff  = enrich_pvalue,
                                    readable      = TRUE)
    if(is.null(d1)) {
        d1 <- data.table(matrix(ncol=10,  nrow=0))
        setnames(d1, desired_cols)
    } else {
        d1 <- as.data.table(d1@result)
    }
    d1[, 'group' := 'go_cc']

    d2 <- enrichGO(gene          = entrez_list,
                                    universe      = all_gene_list,
                                    OrgDb         = org.Hs.eg.db,
                                    ont           = "BP",
                                    pAdjustMethod = "BH",
                                    pvalueCutoff  = enrich_pvalue,
                                    qvalueCutoff  = enrich_pvalue,
                                    readable      = TRUE)
    if(is.null(d2)) {
        d2 <- data.table(matrix(ncol=10,  nrow=0))
        setnames(d2, desired_cols)
    } else {
        d2 <- as.data.table(d2@result)
    }
    d2[, 'group' := 'go_bp']

    d3 <- enrichGO(gene          = entrez_list,
                                    universe      = all_gene_list,
                                    OrgDb         = org.Hs.eg.db,
                                    ont           = "MF",
                                    pAdjustMethod = "BH",
                                    pvalueCutoff  = enrich_pvalue,
                                    qvalueCutoff  = enrich_pvalue,
                                    readable      = TRUE)
    if(is.null(d3)) {
        d3 <- data.table(matrix(ncol=10,  nrow=0))
        setnames(d3, desired_cols)
    } else {
        d3 <- as.data.table(d3@result)
    }
    d3[, 'group' := 'go_mf']

    d4 <- enrichKEGG(gene         = entrez_list,
                                    organism     = 'hsa',
                                    universe      = all_gene_list,
                                    pAdjustMethod = "BH",
                                    qvalueCutoff  = enrich_pvalue,
                                    pvalueCutoff = enrich_pvalue)
    if(is.null(d4)) {
        d4 <- data.table(matrix(ncol=10,  nrow=0))
        setnames(d4, desired_cols)
    } else {
        d4 <- as.data.table(d4@result)
    }
    d4[, 'group' := 'kegg']
    d4 <- d4[, .SD, .SDcols=desired_cols]
    rbindlist(list(d1,d2,d3,d4))
}

plot_enrich <- function(dat, filestem, lims=NULL, x.lbl=NULL, bestN=15) {
        # Concatenate into one table
        dat <-dat[qvalue <= 0.05]
        dat[group=='go_mf', 'Domain' := 'GO Molecular Function']
        dat[group=='go_cc', 'Domain' := 'GO Cellular Component']
        dat[group=='go_bp', 'Domain' := 'GO Biological Process']
        dat[group=='kegg', 'Domain' := 'KEGG Pathway']

        dat.hold <- copy(dat)

        for (grp in c('go_mf','go_cc','go_bp','kegg')) {
        dat <- copy(dat.hold)
        dat <- dat[group==grp]
        # continue with next if no signif
        if (nrow(dat) == 0) { next }

        # Create qsort column to sort within 'counts' ties
        dat[, 'qsort' := -1*qvalue]

        # Set order for plotting
        setkey(dat, 'Domain', Count, qsort)
        #dat[, Description := stringr::str_wrap(Description, width=35)]
        desc_order <- dat$Description
        dat[Count > 0, 'hj' := 1]
        dat[Count > 0, Description := paste0(Description, ' ')]
        dat[Count < 0, 'hj' := 0]
        dat[Count < 0, Description := paste0(' ', Description)]
        dat[, Description := factor(Description, levels=desc_order)]
        #dat[, rnk := frank(-Count, ties.method='random'), by=list(Count<0)]
        #dat <- dat[rnk <= bestN]


        g <- ggplot(dat, aes(x=Count, y=Description, fill=qvalue)) + geom_bar(stat='identity') +
                theme_few() +
                scale_fill_viridis(limits=c(0,0.05), direction=-1) +
                geom_vline(xintercept=0) +
                labs(x=x.lbl) +
                labs(y='GO') +
                #ggplot2::facet_grid(`GO domain` ~ ., drop=T, scales='free', space='free', switch='y') +
                theme(strip.text.y.left = element_text(angle = 0)) +
                geom_text(data=dat, aes(x=0, label=Description, hjust=hj)) +
                theme(axis.title.y=element_blank(),
                axis.text.y=element_blank(),
                axis.ticks.y=element_blank()) +
                scale_x_continuous(expand = expansion(mult = 0.1)) +
                scale_x_continuous(breaks=scales::pretty_breaks(), limits=lims) +
                facet_grid(.~Domain)

        height_cm <- 5 + 0.5 * nrow(dat)
        outfile <- paste0(filestem, '-', grp, '.pdf')
        ggsave(g, file=outfile, width=40, height=height_cm, units='cm', limitsize = FALSE)
    }
}


symbol_to_entrez <- function(symbols) {
    entrez <- mapIds(org.Hs.eg.db, keys = symbols, column = "ENTREZID",keytype="SYMBOL")
    return(unique(unlist(entrez)))
}

ensembl_to_entrez <- function(symbols) {
    entrez <- mapIds(org.Hs.eg.db, keys = symbols, column = "ENTREZID",keytype="ENSEMBL")
    return(unique(unlist(entrez)))
}



# Get sets of genes

## Load in gene up/down tables

files <- list.files('../RNA_DE/', pattern='*.csv', full.names=T)
files <- unlist(lapply(files, normalizePath))


file <- '/vf/users/CARD_ARDIS/users/wellerca/VR-boyden-kd/RNA_DE/hnRNPA1-Neurite.csv'
{
    bn <- basename(file)
    bn <- gsub('.csv','',bn)
    bn <- strsplit(bn, split='-')[[1]]
    experiment <- 'hnRNPA1'
    fraction <- 'Neurite'
    plot_grp <- paste0(experiment, '-', fraction)

    dat <- fread(file)
    
    up <- ensembl_to_entrez(dat[Significance == 'Up']$Geneid)
    down <- ensembl_to_entrez(dat[Significance == 'Down']$Geneid)
    background <- ensembl_to_entrez(dat$Geneid)
    
    up <- getGO(entrez_list = up, all_gene_list=background)
    down <- getGO(entrez_list = down, all_gene_list=background)
    down[, Count := -1 * Count]
    combined <- rbindlist(list(up, down))
    
    plot_enrich(combined, filestem=plot_grp, lims=c(-15,15), x.lbl='Enrichment in hnRNPA-1 knockdown, relative to non-targeting control')

}





file <- '/vf/users/CARD_ARDIS/users/wellerca/VR-boyden-kd/RNA_DE/hnRNPA1-Soma.csv'
{

    experiment <- 'hnRNPA1'
    fraction <- 'Soma'
    plot_grp <- paste0(experiment, '-', fraction)

    dat <- fread(file)
    
    up <- ensembl_to_entrez(dat[Significance == 'Up']$Geneid)
    down <- ensembl_to_entrez(dat[Significance == 'Down']$Geneid)
    background <- ensembl_to_entrez(dat$Geneid)
    
    up <- getGO(entrez_list = up, all_gene_list=background)
    down <- getGO(entrez_list = down, all_gene_list=background)
    down[, Count := -1 * Count]
    combined <- rbindlist(list(up, down))
    
    plot_enrich(combined, filestem=plot_grp, lims=c(-10,10), x.lbl='Enrichment in hnRNPA-1 knockdown, relative to non-targeting control')
}



file <- '/vf/users/CARD_ARDIS/users/wellerca/VR-boyden-kd/RNA_DE/TDP-43-Neurite.csv'
{

    experiment <- 'TDP-43'
    fraction <- 'Neurite'
    plot_grp <- paste0(experiment, '-', fraction)

    dat <- fread(file)
    
    up <- ensembl_to_entrez(dat[Significance == 'Up']$Geneid)
    down <- ensembl_to_entrez(dat[Significance == 'Down']$Geneid)
    background <- ensembl_to_entrez(dat$Geneid)
    
    up <- getGO(entrez_list = up, all_gene_list=background)
    down <- getGO(entrez_list = down, all_gene_list=background)
    down[, Count := -1 * Count]
    combined <- rbindlist(list(up, down))
    
    plot_enrich(combined, filestem=plot_grp, lims=c(-30,30), x.lbl='Enrichment in TDP-43 knockdown, relative to non-targeting control')

}





file <- '/vf/users/CARD_ARDIS/users/wellerca/VR-boyden-kd/RNA_DE/TDP-43-Soma.csv'
{

    experiment <- 'TDP-43'
    fraction <- 'Soma'
    plot_grp <- paste0(experiment, '-', fraction)

    dat <- fread(file)
    
    up <- ensembl_to_entrez(dat[Significance == 'Up']$Geneid)
    down <- ensembl_to_entrez(dat[Significance == 'Down']$Geneid)
    background <- ensembl_to_entrez(dat$Geneid)
    
    up <- getGO(entrez_list = up, all_gene_list=background)
    down <- getGO(entrez_list = down, all_gene_list=background)
    down[, Count := -1 * Count]
    combined <- rbindlist(list(up, down))
    
    plot_enrich(combined, filestem=plot_grp, lims=c(-30,30), x.lbl='Enrichment in TDP-43 knockdown, relative to non-targeting control')

}




file <- '/vf/users/CARD_ARDIS/users/wellerca/VR-boyden-kd/RNA_DE/WT-Soma-Neurite.csv'
{

    experiment <- 'WT'
    fraction <- 'Neurite-v-Soma'
    plot_grp <- paste0(experiment, '-', fraction)

    dat <- fread(file)
    
    up <- ensembl_to_entrez(dat[Significance == 'Up']$Geneid)
    down <- ensembl_to_entrez(dat[Significance == 'Down']$Geneid)
    background <- ensembl_to_entrez(dat$Geneid)
    
    up <- getGO(entrez_list = up, all_gene_list=background)
    down <- getGO(entrez_list = down, all_gene_list=background)
    down[, Count := -1 * Count]
    combined <- rbindlist(list(up, down))
    
    plot_enrich(combined, filestem=plot_grp, lims=c(-200,200), x.lbl='Enrichment in Neurites, relative to Whole cell')

}




quit()

