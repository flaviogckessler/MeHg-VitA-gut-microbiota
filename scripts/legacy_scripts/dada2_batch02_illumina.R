if (!requireNamespace("BiocManager", quietly=TRUE))
  install.packages("BiocManager")
BiocManager::install("DECIPHER")

library("dada2")
library(biomformat)
library(phyloseq)
library(Biostrings)
library(ggplot2)

######Dada2#######

sra_bioproject <- "PRJNA1489663"
raw_fastq_dir <- Sys.getenv("RAW_FASTQ_DIR", unset = "")

if (!nzchar(raw_fastq_dir) || !dir.exists(raw_fastq_dir)) {
  stop(
    paste0(
      "Raw FASTQ files are not distributed in this repository. ",
      "After BioProject ", sra_bioproject,
      " is released, download the reads from NCBI SRA and set RAW_FASTQ_DIR ",
      "to the directory containing the paired FASTQ files."
    )
  )
}

path <- normalizePath(raw_fastq_dir, mustWork = TRUE)
quality_output_dir <- Sys.getenv(
  "DADA2_QC_OUTPUT_DIR",
  unset = file.path(path, "dada2_qc")
)
dir.create(quality_output_dir, recursive = TRUE, showWarnings = FALSE)

fnFs <- sort(list.files(path, pattern="_R1_001.fastq", full.names = TRUE))
fnRs <- sort(list.files(path, pattern="_R2_001.fastq", full.names = TRUE))

sample.names <- sapply(strsplit(basename(fnFs), "_"), `[`, 1)
sample.names <- sapply(strsplit(basename(fnRs), "_"), `[`, 1)

l.images = list(1:9, 10:19, 20:29, 30:39, 40, 49, 50)
i = 0
for (l in l.images) {
  
  print(sample.names[l])
  # p.f <- plotQualityProfile(fnFs[l])
  p.r <- plotQualityProfile(fnRs[l])
  
  i = i+1
  
  # tiff(filename = paste0("output/", get("i"), "_plotquality_forward.tif"), width = 10, height = 8, units = "in", res = 300)
  # print(p.f)
  # dev.off()

  tiff(filename = file.path(quality_output_dir, paste0(get("i"), "_plotquality_reverse.tif")), width = 10, height = 8, units = "in", res = 300)
  print(p.r)
  dev.off()

}

filtered_dir <- file.path(path, "filtered")
dir.create(filtered_dir, recursive = TRUE, showWarnings = FALSE)
filtFs <- file.path(filtered_dir, paste0(sample.names, "_F_filt.fastq.gz"))
filtRs <- file.path(filtered_dir, paste0(sample.names, "_R_filt.fastq.gz"))

names(filtFs) <- sample.names
names(filtRs) <- sample.names

out <- filterAndTrim(fnFs, filtFs, fnRs, filtRs, truncLen=c(149, 70),
                     maxN=0, maxEE=2, truncQ=2, rm.phix=TRUE,
                     compress=TRUE, trimLeft= c(20, 0),multithread=TRUE)

errF <- learnErrors(filtFs, multithread=TRUE)
errR <- learnErrors(filtRs, multithread=TRUE)

plotErrors(errF, nominalQ=TRUE)
plotErrors(errR, nominalQ=TRUE)

dadaFs <- dada(filtFs, err=errF, pool="pseudo", multithread=TRUE)
dadaRs <- dada(filtRs, err=errR, pool="pseudo", multithread=TRUE)

mergers <- mergePairs(dadaFs, filtFs, dadaRs, filtRs, verbose=TRUE)

seqtab <- makeSequenceTable(mergers)
dim(seqtab)

seqtab.nochim <- removeBimeraDenovo(seqtab, multithread=TRUE, method="consensus", verbose=TRUE)
dim(seqtab.nochim)

getN <- function(x) sum(getUniques(x))
track <- cbind(out, sapply(dadaFs, getN), sapply(dadaRs, getN), sapply(mergers, getN), rowSums(seqtab.nochim))
# If processing a single sample, remove the sapply calls: e.g. replace sapply(dadaFs, getN) with getN(dadaFs)
colnames(track) <- c("input", "filtered", "denoisedF", "denoisedR", "merged", "nonchim")
rownames(track) <- sample.names
head(track)
View(as.data.frame(track))

taxa <- assignTaxonomy(seqtab.nochim, "silva_nr99_v138.1_train_set.fa.gz", multithread = TRUE)

  taxa <- addSpecies(taxa, "silva_species_assignment_v138.1.fa.gz") #verificar se vale a pena adicionar espécie.

######Phyloseq#######

theme_set(theme_bw())

ps <- phyloseq(otu_table(seqtab.nochim, taxa_are_rows=FALSE), tax_table(taxa))

mapfile <- import_qiime_sample_data("~/fernanda/microbiota/metadados_microbiota_test.csv") #deve estar delimitado por espaços

dna <- Biostrings::DNAStringSet(taxa_names(ps))
names(dna) <- taxa_names(ps)

ps <- phyloseq(otu_table(seqtab.nochim, taxa_are_rows=FALSE),
               sample_data(mapfile),
               tax_table(taxa),
               refseq(dna))

ps.m <- ps
taxa_names(ps) <- paste0("ASV", seq(ntaxa(ps)))

save.image(file = "dada2_com_primer.RData")

ps <- subset_samples(ps, SampleID != "102")
ps <- subset_samples(ps, SampleID != "29")

#Exporting data (asv == sequence)

tax<-as(tax_table(ps),"matrix")
tax_cols <- colnames(tax)
tax<-as.data.frame(tax)
tax$taxonomy<-do.call(paste, c(tax[tax_cols], sep=";"))
for(co in tax_cols) tax[co]<-NULL
write.table(tax, "output_com_primer/tax.txt", quote=FALSE, col.names=FALSE, sep="\t")

otu<-t(as(otu_table(ps),"matrix"))
otu_biom<-make_biom(data=otu)
write_biom(otu_biom,"output_com_primer/otu_biom.biom")

write.table(t(otu_table(ps)), "output_com_primer/seqtab.txt", sep="\t", row.names=TRUE, col.names=NA, quote=FALSE)

uniquesToFasta(seqtab.nochim, fout='output_com_primer/rep-seqs.fna', ids=colnames(seqtab.nochim))

alpha <- estimate_richness(ps, measures=c("Observed", "Simpson", "InvSimpson", "Shannon", "Chao1"))
write.csv2(alpha, "output_com_primer/alpha_diversity.csv", quote=FALSE)

write.table(sample_data(ps),"output_com_primer/sample-metadata.txt", sep="\t", row.names=FALSE, col.names=TRUE, quote=FALSE)

save.image("retry_dada2.RData")

#Exporting data (asv == asv_number)

tax<-as(tax_table(ps.m),"matrix")
tax_cols <- colnames(tax)
tax<-as.data.frame(tax)
tax$taxonomy<-do.call(paste, c(tax[tax_cols], sep="\t"))
for(co in tax_cols) tax[co]<-NULL
write.table(tax, "output_com_primer/tax_m.txt", quote=FALSE, col.names=FALSE, sep="\t")

uniquesToFasta(seqtab.nochim, fout='output_com_primer/rep-seqs.fna', ids=colnames(seqtab.nochim))

top20 <- names(sort(taxa_sums(ps), decreasing=TRUE))[1:20]
ps.top20 <- transform_sample_counts(ps, function(OTU) OTU/sum(OTU))
ps.top20 <- prune_taxa(top20, ps.top20)
plot_bar(ps.top20, x="GRUPO", fill="Family") + facet_wrap(~SEXO ~LINHAGEM, scales="free_x")
