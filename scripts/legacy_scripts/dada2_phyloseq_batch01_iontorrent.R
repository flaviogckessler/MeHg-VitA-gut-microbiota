#https://benjjneb.github.io/dada2/tutorial.html

#install.packages("BiocManager")
#BiocManager::install("dada2")
#BiocManager::install("biomformat")
#BiocManager::install("phyloseq")
#BiocManager::install("Biostrings")
#####################################
library("dada2")
library("biomformat")
library("phyloseq")
library("Biostrings")
######Dada2#######

path <- file.path("data", "external", "fastq", "batch01")

seqs <- sort(list.files(path, full.names = TRUE))

sample.names <- sapply(strsplit(basename(seqs), ".IonXpress"), `[`, 1)
sample.names
#plotQualityProfile(seqs[1:2])

#Criar diretório para filtragem
filtseqs <- file.path(path, "filtered", paste0(sample.names, "_filt.fastq.gz"))

names(filtseqs) <- sample.names

#Filtragem e trimerização
#TruncLen=275 (255 do tamanho V4 + 20 primer)
#trimLeft = 20 (tamanho do primer)
out <- filterAndTrim(seqs, filtseqs, truncLen=c(275),
                     maxN=0, maxEE=2, truncQ=2, rm.phix=TRUE,
                     compress=TRUE, trimLeft= c(20), multithread=TRUE)
head(out)
#Frequencia de erro
errseqs <- learnErrors(filtseqs, multithread=TRUE)

plotErrors(errseqs, nominalQ=TRUE)

#Para IonTorrent se recomenda os seguintes parametros HOMOPOLYMER_GAP_PENALTY=-1, BAND_SIZE=32
#https://benjjneb.github.io/dada2/faq.html#can-i-use-dada2-with-my-454-or-ion-torrent-data
dadaPseudos <- dada(filtseqs, err=errseqs, pool="pseudo", HOMOPOLYMER_GAP_PENALTY=-1, BAND_SIZE=32, multithread = TRUE)

seqtab <- makeSequenceTable(dadaPseudos)
dim(seqtab)

seqtab.nochim <- removeBimeraDenovo(seqtab, method="consensus", verbose=TRUE, multithread=TRUE)
dim(seqtab.nochim)
rm(seqtab)

taxa <- assignTaxonomy(seqtab.nochim, "silva_nr99_v138.1_train_set.fa.gz", multithread=TRUE)

taxa <- addSpecies(taxa, "silva_species_assignment_v138.1.fa.gz") #verificar se vale a pena adicionar espécie.

######Phyloseq#######

library("ggplot2")
theme_set(theme_bw())

ps <- phyloseq(otu_table(seqtab.nochim, taxa_are_rows=FALSE), tax_table(taxa))

meta_data <- read.csv("Metadata_mehg.csv", sep=";", encoding = 'Latin-1',header = TRUE,row.names = 1)
View(meta_data)


#meta_data_dariane <- read.csv("mapping_file.CSV", sep=";", encoding = 'Latin-1')
#View(meta_data_dariane)
sampledata = sample_data(data.frame(meta_data,
                                    stringsAsFactors = FALSE))
#sampledata = sample_data(data.frame(meta_data,
#                                    row.names=sample_names(ps),
#                                    stringsAsFactors = FALSE))
#sampledata = sample_data(data.frame(meta_data,
#                                    row.names=meta_data$sample_id,
#                                    stringsAsFactors = FALSE))

dna <- Biostrings::DNAStringSet(taxa_names(ps))
names(dna) <- taxa_names(ps)

ps <- phyloseq(otu_table(seqtab.nochim, taxa_are_rows=FALSE), 
                sampledata, 
                tax_table(taxa),
                refseq(dna))
View(ps)

ps.m <- ps
taxa_names(ps.m) <- paste0("ASV", seq(ntaxa(ps.m)))

#Exporting data (asv == sequence)

tax<-as(tax_table(ps),"matrix")
tax_cols <- colnames(tax)
tax<-as.data.frame(tax)
tax$taxonomy<-do.call(paste, c(tax[tax_cols], sep=";"))
for(co in tax_cols) tax[co]<-NULL
write.table(tax, "output/tax.txt", quote=FALSE, col.names=FALSE, sep="\t")

otu<-t(as(otu_table(ps),"matrix"))
otu_biom<-make_biom(data=otu)
write_biom(otu_biom,"output/otu_biom.biom")
rm(otu_biom)

write.table(t(otu_table(ps)), "output/seqtab.txt", sep="\t", row.names=TRUE, col.names=NA, quote=FALSE)

uniquesToFasta(seqtab.nochim, fout='output/rep-seqs.fna', ids=colnames(seqtab.nochim))

alpha <- estimate_richness(ps, measures=c("Observed", "Simpson", "InvSimpson", "Shannon", "Chao1"))
View(alpha)

write.csv2(alpha, "output/alpha_diversity.csv", quote=FALSE)

write.table(sample_data(ps),"output/sample-metadata.txt", sep="\t", row.names=FALSE, col.names=TRUE, quote=FALSE)

#Exporting data (asv == asv_number)

tax<-as(tax_table(ps.m),"matrix")
tax_cols <- colnames(tax)
tax<-as.data.frame(tax)
tax$taxonomy<-do.call(paste, c(tax[tax_cols], sep="\t"))
for(co in tax_cols) tax[co]<-NULL
write.table(tax, "output/tax_m.txt", quote=FALSE, col.names=FALSE, sep="\t")

uniquesToFasta(seqtab.nochim, fout='output/rep-seqs.fna', ids=colnames(seqtab.nochim))
#################################################################################
#################################################################################
for(i in 1:length(tax$taxonomy)) {
  buf = strsplit(tax$taxonomy[i],"\t")
  buf = paste("k__",buf[[1]][1],";p__",buf[[1]][2],";c__",buf[[1]][3],";o__",buf[[1]][4],";f__",buf[[1]][5],";g__",buf[[1]][6],";s__",buf[[1]][7],sep = "")
  #print(buf)
  tax$taxonomy[i] = buf
}
otu_d <- as.data.frame(otu)

library(tibble)
teste = add_column(otu_d, tax, .after = 0)
rm(otu_d)
View(teste)

names(teste)[names(teste) == 'taxonomy'] <- 'lineage'

write.csv(teste, file.path("output", "batch01_phyloseq_export.csv"))

library("ape")
random_tree = rtree(ntaxa(ps), rooted=TRUE, tip.label=taxa_names(ps))
plot(random_tree)
