# 20 samples collected across India with the help of BAIF, analysed/sequenced at NCCS
# Correlation analysis
# Microbe network analysis

library(dada2)
library(fastqcr)
# install.packages("fastqcr")
library(phyloseq)

# setwd('C:\\Users\\Admin\\Desktop\\NIWASm\\NEW')
path = 'C:\\Users\\Admin\\Desktop\\NIWASm\\NEW'
# list.files()
packageVersion('dada2')

# Forward and reverse fastq filenames have format: SAMPLENAME_R1_001.fastq and SAMPLENAME_R2_001.fastq
fnFs <- sort(list.files(path, pattern="_R1_001.fastq.gz", full.names = TRUE))
fnRs <- sort(list.files(path, pattern="_R2_001.fastq.gz", full.names = TRUE))
sample.names <- sapply(strsplit(basename(fnFs), "_"), `[`, 1)
sample.names
plotQualityProfile(fnRs[1:2])

# Place filtered files in filtered/ subdirectory
filtFs <- file.path(path, "filtered", paste0(sample.names, "_F_filt.fastq.gz"))
filtRs <- file.path(path, "filtered", paste0(sample.names, "_R_filt.fastq.gz"))
names(filtFs) <- sample.names
names(filtRs) <- sample.names


out <- filterAndTrim(fnFs, filtFs, fnRs, filtRs, truncLen=c(160,160),
                     maxN=0, maxEE=c(2,2), truncQ=2, rm.phix=TRUE,
                     compress=TRUE, multithread=FALSE) # On Windows set multithread=FALSE
head(out)


errF <- learnErrors(filtFs, multithread=TRUE)
errR <- learnErrors(filtRs, multithread=TRUE)
plotErrors(errF, nominalQ=TRUE)
plotErrors(errR, nominalQ=TRUE)

dadaFs <- dada(filtFs, err=errF, multithread=TRUE)
dadaRs <- dada(filtRs, err=errR, multithread=TRUE)
dadaFs[[1]]


mergers <- mergePairs(dadaFs, filtFs, dadaRs, filtRs, verbose=TRUE)
# Inspect the merger data.frame from the first sample
head(mergers[[1]])
mergers[[1]]

seqtab <- makeSequenceTable(mergers)
dim(seqtab)

# Inspect distribution of sequence lengths
table(nchar(getSequences(seqtab)))


seqtab.nochim <- removeBimeraDenovo(seqtab, method="consensus", multithread=TRUE, verbose=TRUE)
dim(seqtab.nochim)
sum(seqtab.nochim)/sum(seqtab)


getN <- function(x){
  sum(getUniques(x))
}
track <- cbind(out, sapply(dadaFs, getN), sapply(dadaRs, getN), sapply(mergers, getN), rowSums(seqtab.nochim))
# If processing a single sample, remove the sapply calls: e.g. replace sapply(dadaFs, getN) with getN(dadaFs)
colnames(track) <- c("input", "filtered", "denoisedF", "denoisedR", "merged", "nonchim")
rownames(track) <- sample.names
head(track)

taxa <- assignTaxonomy(seqtab.nochim, "C:/Users/Admin/Desktop/NIWASm/silva/silva_nr99_v138.1_train_set.fa.gz", multithread=TRUE)

taxa.print <- taxa # Removing sequence row names for display only
rownames(taxa.print) <- NULL
head(taxa.print)

saveRDS(taxa,'taxa.rds')
write.csv(track,'track.csv')
saveRDS(seqtab.nochim,'seqtab.nochim.rds')

# readRDS('taxa.rds')
# readRDS('seqtab.nochim.rds')




############# Make Phyloseq table ################
library(ggplot2)
metadata <- read.table("C:/Users/Admin/Desktop/NIWASm/sample_metadata.csv", sep='\t',header=T,row.names = 1,check.names=F)
View(metadata)
sample_data(ps) <- metadata
  
otu.tab = otu_table(object = seqtab.nochim, taxa_are_rows = F)
tax.tab = tax_table(object = taxa)

ps = phyloseq(otu.tab, tax.tab, metadata)
richness = plot_richness(ps,x='Group', measures=c("Shannon", "Simpson"),color = 'Group')
richness + geom_point(size=3)
box = richness + geom_boxplot()

ps.prop <- transform_sample_counts(ps, function(otu) otu/sum(otu))
ord.nmds.bray <- ordinate(ps.prop, method="NMDS", distance="bray")
plot_ordination(ps.prop, ord.nmds.bray, color='Group', title="Bray NMDS") + geom_point(size = 3.5)


save.image(file='niwasm1.RData')
load.image('niwasm1.RData')
######### Perform Rarefaction ###########







