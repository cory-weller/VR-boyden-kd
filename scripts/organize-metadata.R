#!/usr/bin/env Rscript

library(data.table)

# hnRNPA1
dat <- fread('metadata/hnrnpa1.txt')
dat[, 'Experiment' := 'hnRNPA1']
dat <- dat[, .SD, .SDcols=c('Experiment','Sample ID','Long Name','Short Name')]
setnames(dat, 'Long Name', 'Long_Name')
setnames(dat, 'Short Name', 'Short_Name')
setnames(dat, 'Sample ID', 'Sample_ID')
dat[, c('sgRNA','Fraction','Replicate') := tstrsplit(Long_Name, split=' ')]
dat[, Replicate := as.numeric(Replicate)]
dat[sgRNA=='sg1126', Contrast := 'nt_control']
dat[sgRNA=='sg1128', Contrast := 'hnRNPA1_knockdown']
hnRNPA1 <- copy(dat)


HNRNPA1.files <- data.table('filename'=list.files('data/hnRNPA1-kd', pattern='*.fastq.gz'))
HNRNPA1.files[, 'filestem' := tstrsplit(filename, split='_')[1]]
HNRNPA1.files[, 'filename' := NULL]
HNRNPA1.files <- unique(HNRNPA1.files)
HNRNPA1.files[, 'Sample_ID' := tstrsplit(filestem, split='-')[3]]
HNRNPA1.files[, 'Sample_ID' := as.numeric(Sample_ID)]
hnRNPA1 <- merge(hnRNPA1, HNRNPA1.files, by='Sample_ID')
kept_order <- colnames(hnRNPA1)

# TDP-43
dat <- fread('metadata/tdp43.txt')
dat[, 'Experiment' := 'TDP-43']
dat <- dat[, .SD, .SDcols=c('Experiment','Sample ID','Long name','Short name')]
setnames(dat, 'Long name', 'Long_Name')
setnames(dat, 'Short name', 'Short_Name')
setnames(dat, 'Sample ID', 'Sample_ID')
dat[, c('sgRNA','Fraction','Replicate') := tstrsplit(Long_Name, split=' ')]
dat[, Replicate := as.numeric(Replicate)]
dat[sgRNA=='sg100', Contrast := 'nt_control']
dat[sgRNA=='sg200', Contrast := 'TDP43_knockdown']
tdp43 <- copy(dat)

TDP43.files <- data.table('filename'=list.files('data/TDP43-kd', pattern='*.fastq.gz'))
TDP43.files[, 'filestem' := tstrsplit(filename, split='_')[1]]
TDP43.files[, 'filename' := NULL]
TDP43.files <- unique(TDP43.files)
TDP43.files[, 'Sample_ID' := tstrsplit(filestem, split='-')[3]]
TDP43.files[, 'Sample_ID' := as.numeric(Sample_ID)]
tdp43 <- merge(tdp43, TDP43.files, by='Sample_ID')



# WT
dat <- CJ('Experiment'='WT', 'Sample_ID'=4:7, 'Fraction'=c('Neurite','Soma'))
dat[, Replicate := Sample_ID]
dat[, 'Long_Name' := paste0('WT ', Fraction, ' ', Replicate)]
dat[, sgRNA := NA]
dat[, Short_Name := NA]
dat[, Contrast := Fraction]
wt <- copy(dat)
wt[, filestem := paste0(tolower(Fraction), '-', Sample_ID)]
setcolorder(wt, kept_order)
setkey(wt, filestem)


# FUS
dat <- fread('metadata/fus.txt')
dat[, 'Experiment' := 'FUS']
dat <- dat[, .SD, .SDcols=c('Experiment','Sample ID','Long Name','Short Name')]
setnames(dat, 'Long Name', 'Long_Name')
setnames(dat, 'Short Name', 'Short_Name')
setnames(dat, 'Sample ID', 'Sample_ID')
dat[, c('sgRNA','Fraction','Replicate') := tstrsplit(Long_Name, split=' ')]
dat[, Replicate := as.numeric(Replicate)]
dat[sgRNA=='sg100', Contrast := 'nt_control']
dat[sgRNA=='sg1152', Contrast := 'FUS_knockdown']
dat[, Sample_ID := Sample_ID - 24]
fus <- copy(dat)


FUS.files <- data.table('filename'=list.files('data/FUS-kd', pattern='*.fastq.gz'))
FUS.files[, 'filestem' := tstrsplit(filename, split='_')[1]]
FUS.files[, 'filename' := NULL]
FUS.files <- unique(FUS.files)
FUS.files[, 'Sample_ID' := tstrsplit(filestem, split='-')[3]]
FUS.files[, 'Sample_ID' := as.numeric(Sample_ID)]
fus <- merge(fus, FUS.files, by='Sample_ID')


metadata <- rbindlist(list(wt, hnRNPA1, tdp43, fus))
fwrite(metadata, file='metadata/rnaseq-metadata.csv', quote=F, row.names=F, col.names=T, sep=',', na='NA')