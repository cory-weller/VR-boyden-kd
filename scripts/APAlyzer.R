#!/usr/bin/env Rscript

if (! require('BiocManager'))   { install.packages("BiocManager") }
if (! require('APAlyzer'))      { BiocManager::install("APAlyzer"); library(APAlyzer) }
if (! require('repmis'))        { install.packages('repmis'); library(repmis) }
if (! require('ggplot2'))        { install.packages('ggplot2'); library(ggplot2) }
if (! require('ggthemes'))        { install.packages('ggthemes'); library(ggthemes) }


#2 Load Reference Genome, check if hg38 
URL="https://github.com/RJWANGbioinfo/PAS_reference_RData/blob/master/"
file ="hg38_REF.RData"
source_data(paste0(URL,file,"?raw=True"))

#3 Build 3' and Intron references
PASREF <- REF4PAS(refUTRraw_hg38,dfIPA_hg38,dfLE_hg38)
UTRdbraw <- PASREF$UTRdbraw
dfIPA <- PASREF$dfIPA
dfLE <- PASREF$dfLE



library(data.table)
bam_dir <- 'RNA_alignment'

# Create data frame of files
files <- data.table(path=unlist(lapply(list.files(bam_dir, pattern='soma-.*.bam|neurite-.*.bam', recursive=T, full.names=T), normalizePath)))
files[, SampleID := tstrsplit(basename(path), split='\\.')[1]]
files[, Condition := tstrsplit(SampleID, split='-')[1]]
files[, Condition := factor(Condition, levels=c('soma','neurite'))]
setkey(files, Condition)

# Set up named vector of files for APAlyzer
flsall <- files$path
names(flsall) <- files$SampleID

if (! file.exists('APA/APAlyzer-UTR.RDS')) {
# Calculate relative expression of 3'UTR APA
# Strandtype here is none because other options are forward or reverse sequencing only. Because it is paired end, we use NONE here. 
    UTR_APA_OUT <- PASEXP_3UTR(UTRdbraw, flsall, Strandtype = "NONE")
    saveRDS(UTR_APA_OUT, 'APA/APAlyzer-UTR.RDS')
} else {
    UTR_APA_OUT <- readRDS('APA/APAlyzer-UTR.RDS')
}
#6 save 3'UTR as a RDS file for local loading and further anlysis in R


if (! file.exists('APA/APAlyzer-IPA.RDS')) {

# Calculate relative expression of IPA
IPA_OUT=PASEXP_IPA(dfIPA, dfLE, flsall, nts=1, SeqType = "ThreeMostPairEnd", Strandtype= "NONE")
saveRDS(IPA_OUT, "APA/APAlyzer-IPA.RDS")

} else {
    IPA_OUT <- readRDS('APA/APAlyzer-IPA.RDS')
}

sampleTable <- files[, .SD, .SDcols=c('SampleID','Condition')]
setnames(sampleTable, c('samplename', 'condition'))
sampleTable <- as.data.frame(sampleTable)


# Significantly regulated APA in 3’UTRs

APA_3UTR <- APAdiff(sampleTable,
                    UTR_APA_OUT,
                    conKET='soma',
                    trtKEY='neurite',
                    PAS='3UTR',
                    CUTreads=5)

APA_3UTR <- as.data.table(APA_3UTR)
APA_3UTR[APAreg == 'DN', Significance := 'Whole cell']
APA_3UTR[APAreg == 'UP', Significance := 'Neurite']
APA_3UTR[APAreg == 'NC', Significance := 'NS']

writeLines(as.character(unique(APA_3UTR$gene_symbol)), con='APA/APAlyzer_3utr_set.txt')
writeLines(as.character(unique(APA_3UTR[Significance=='Neurite']$gene_symbol)), con='APA/APAlyzer_3utr_neurite.txt')
writeLines(as.character(unique(APA_3UTR[Significance=='Whole cell']$gene_symbol)), con='APA/APAlyzer_3utr_soma.txt')

source('scripts/colors.R')

g.3utr.volcano <- ggplot(data=APA_3UTR, aes(x=RED, y=-1*log10(p_adj), color=Significance)) + 
        geom_point() + 
        scale_color_manual(values=c('NS'=not_signif_color, 'Whole cell'=soma_color,'Neurite'=neurite_color)) +
        labs(x='Relative Expression Delta (RED)', y='-log10(P)', title="APAlyzer Significantly Regulated 3'UTRs") +
        theme_few() +
        geom_hline(yintercept=-1*log10(0.215), linetype='dashed', alpha=0.5)

ggsave(g.3utr.volcano , file='APA/APAlyzeer-3UTR-volcano.png', width=15, height=15, units='cm')
ggsave(g.3utr.volcano , file='APA/APAlyzeer-3UTR-volcano.pdf', width=15, height=15, units='cm')

# APAVolcano(APA_3UTR, PAS='3UTR', Pcol = "pvalue", plot_title='3UTR APA')
# png(file="APAlyzerResults/3UTRVolcanoplot.png")

# APABox(APA_3UTR, xlab = "APAreg", ylab = "RED", plot_title = NULL)
# png(file="APAlyzerResults/APABoxPlot.png")


# Significantly regulated APA in Intron
APA_IPA <- APAdiff(sampleTable,
                        IPA_OUT, 
                        conKET='soma',
                        trtKEY='neurite',
                        PAS='IPA',
                        CUTreads=5)
                        

APA_IPA <- as.data.table(APA_IPA)
APA_IPA[APAreg == 'DN', Significance := 'Whole cell']
APA_IPA[APAreg == 'UP', Significance := 'Neurite']
APA_IPA[APAreg == 'NC', Significance := 'NS']


writeLines(as.character(unique(APA_IPA$gene_symbol)), con='APA/APAlyzer_ipa_set.txt')
writeLines(as.character(unique(APA_IPA[Significance=='Neurite']$gene_symbol)), con='APA/APAlyzer_ipa_neurite.txt')
writeLines(as.character(unique(APA_IPA[Significance=='Whole cell']$gene_symbol)), con='APA/APAlyzer_ipa_soma.txt')

g.ipa.volcano <- ggplot(data=APA_IPA, aes(x=RED, y=-1*log10(p_adj), color=Significance)) + 
        geom_point() + 
        scale_color_manual(values=c('NS'=not_signif_color, 'Whole cell'=soma_color,'Neurite'=neurite_color)) +
        labs(x='Relative Expression Delta (RED)', y='-log10(P)', title="APAlyzer Significantly Regulated IPAs") +
        theme_few() +
        geom_hline(yintercept=-1*log10(0.15), linetype='dashed', alpha=0.5)

ggsave(g.ipa.volcano , file='APA/APAlyzeer-IPA-volcano.png', width=15, height=15, units='cm')
ggsave(g.ipa.volcano , file='APA/APAlyzeer-IPA-volcano.pdf', width=15, height=15, units='cm')