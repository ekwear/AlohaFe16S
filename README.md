# AlohaFe16S
16S rRNA gene amplicon analysis code for project AlohaFe

1) ASV generation in DADA2, using R through server

library(dada2); packageVersion("dada2")
setwd("~/AlohaFE")
path <- "~/AlohaFE"
forward_reads <- sort(list.files(path,pattern="_R1_001.fastq",full.names=TRUE))
reverse_reads <- sort(list.files(path,pattern="_R2_001.fastq",full.names=TRUE))
sample.names <- sapply(strsplit(basename(forward_reads), "_"), `[`, 1)
pdf("~/AlohaFE/read1_quality.pdf")
plotQualityProfile(forward_reads[1:6])
dev.off()
pdf("~/AlohaFE/read2_quality.pdf")
plotQualityProfile(reverse_reads[1:6])
dev.off()
filtFs <- file.path(path, "filtered", paste0(sample.names, "_F_filt.fastq.gz"))
filtRs <- file.path(path, "filtered", paste0(sample.names, "_R_filt.fastq.gz"))
names(filtFs) <- sample.names
names(filtRs) <- sample.names
filtered_out <- filterAndTrim(forward_reads, filtFs, reverse_reads, filtRs, truncLen=c(245,160),maxN=0, maxEE=c(2,5), truncQ=8, rm.phix=TRUE, compress=TRUE,multithread=TRUE)
write.table(filtered_out, "filtered_out.txt")
errF <- learnErrors(filtFs, nbases=2e+08, multithread=TRUE)
errR <- learnErrors(filtRs, nbases=2e+08, multithread=TRUE)
pdf("~/AlohaFE/errF.pdf")
plotErrors(errF, nominalQ=TRUE)
dev.off()
derep_forward <- derepFastq(filtFs, verbose=TRUE)
names(derep_forward) <- sample.names
derep_reverse <- derepFastq(filtRs, verbose=TRUE)
names(derep_reverse) <- sample.names
dada_forward <- dada(derep_forward, err=errF, pool=TRUE)
dada_reverse <- dada(derep_reverse, err=errR, pool=TRUE)
merged_amplicons <- mergePairs(dada_forward, derep_forward, dada_reverse, derep_reverse, verbose=TRUE)
head(merged_amplicons)
seqtab <- makeSequenceTable(merged_amplicons)
class(seqtab)
dim(seqtab)
seqtab.nochim <- removeBimeraDenovo(seqtab, verbose=T)
sum(seqtab.nochim)/sum(seqtab)
getN <- function(x) sum(getUniques(x))
track <- cbind(filtered_out, sapply(dada_forward, getN), sapply(dada_reverse, getN), sapply(merged_amplicons, getN), rowSums(seqtab.nochim))
colnames(track) <- c("input", "filtered", "denoisedF", "denoisedR", "merged", "nonchim")
rownames(track) <- sample.names
head(track)
write.table(track, "track_seqnumbers.tsv",sep="\t")
taxa <- assignTaxonomy(seqtab.nochim, "~/AlohaFE/silva_nr_v138_train_set.fa.gz", tryRC=T) 
taxa.plus <- addSpecies(taxa,"~/AlohaFE/silva_species_assignment_v138.fa.gz", verbose=TRUE) 
unname(taxa.plus)
asv_headers <- vector(dim(seqtab.nochim)[2], mode="character")
for (i in 1:dim(seqtab.nochim)[2]) {
 asv_headers[i] <- paste(">ASV", i, sep="_")
 }
asv_seqs <- colnames(seqtab.nochim)
asv_fasta <- c(rbind(asv_headers, asv_seqs))
write.table(asv_fasta, "AlohaFE_ASV_April20.fasta", sep="\t",quote=F, col.names=NA)
asv_tab <- t(seqtab.nochim)
row.names(asv_tab) <- sub(">", "", asv_headers)
write.table(asv_tab,"AlohaFE_asv_counts_April20.tsv",sep="\t", quote=F, col.names=NA) 
asv_tax <- taxa.plus
row.names(asv_tax) <- sub(">", "", asv_headers)
write.table(asv_tax,"AlohaFE_ASVs_tax_April20.tsv", sep="\t", quote=F, col.names=NA)
