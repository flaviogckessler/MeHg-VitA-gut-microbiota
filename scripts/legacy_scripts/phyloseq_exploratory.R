#Source: https://ycl6.github.io/16S-Demo/index.html

# Historical workstation path removed for privacy; run from the project root.
setwd(".")

library("data.table")
library("phyloseq")
library("DESeq2")
library("ggplot2")
library("ggbeeswarm")
library("ggrepel")
library("vegan")
library("tidyr")
library("dplyr")
library("gridExtra")
library("GUniFrac")
library("dada2")
source("taxa_summary.R", local = TRUE)


if(!dir.exists("outfiles")) { dir.create("outfiles") }
if(!dir.exists("images")) { dir.create("images") }

gc()
load("dada2.RData")

ps0 = subset_taxa(ps, Kingdom == "Bacteria" & !is.na(Phylum) & !is.na(Class) & 
                    Family != "Mitochondrias")

tdt = data.table(setDT(as.data.frame(tax_table(ps0))), 
                 TotalCounts = taxa_sums(ps0), SV = taxa_names(ps0))

tdt

ggplot(tdt, aes(TotalCounts)) + geom_histogram(bins = 50) + theme_bw() + 
  ggtitle("Histogram of Total Counts")

tdt[(TotalCounts == 1), .N]     # singletons
tdt[(TotalCounts == 2), .N]     # doubletons

# taxa cumulative sum
taxcumsum = tdt[, .N, by = TotalCounts]
setkey(taxcumsum, TotalCounts)
taxcumsum[, CumSum := cumsum(N)]

# Plot the cumulative sum of ASVs against the total counts
pCumSum = ggplot(taxcumsum, aes(TotalCounts, CumSum)) + geom_point() + theme_bw() + 
  xlab("Filtering Threshold") + ylab("ASV Filtered")

png("images/FilterTaxa-taxa-abundance.png", width = 8, height = 8, units = "in", res = 300)
gridExtra::grid.arrange(pCumSum, pCumSum + xlim(0, 500), 
                        pCumSum + xlim(0, 100), pCumSum + xlim(0, 50), nrow = 2, 
                        top = "ASVs that would be filtered vs. minimum taxa counts threshold")
invisible(dev.off())

mdt = fast_melt(ps0)

mdt

prevdt = mdt[, list(Prevalence = sum(count > 0), TotalCounts = sum(count), 
                    MaxCounts = max(count)), by = TaxaID]

prevdt

ggplot(prevdt, aes(Prevalence)) + geom_histogram(bins = 50) + theme_bw() +
  ggtitle("Histogram of Taxa Prevalence")

prevdt[(Prevalence == 1), .N]   # singletons
prevdt[(Prevalence == 2), .N]   # doubletons
ggplot(prevdt, aes(MaxCounts)) + geom_histogram(bins = 100) + xlim(0, 500) + theme_bw() + 
  ggtitle("Histogram of Maximum TotalCounts")

table(prevdt$MaxCounts)[1:50]
prevdt[(MaxCounts == 1)]    # singletons

# taxa cumulative sum
prevcumsum = prevdt[, .N, by = Prevalence]
setkey(prevcumsum, Prevalence)
prevcumsum[, CumSum := cumsum(N)]

# Plot the cumulative sum of ASVs against the prevalence
pPrevCumSum = ggplot(prevcumsum, aes(Prevalence, CumSum)) + geom_point(size = 2, alpha = 0.5) + 
  theme_bw() + xlab("Filtering Threshold") + ylab("ASVs Filtered") + 
  ggtitle("ASVs that would be filtered vs. minimum sample count threshold")

png("images/FilterTaxa-taxa-prevalence.png", width = 8, height = 8, units = "in", res = 300)
pPrevCumSum
invisible(dev.off())

png("images/FilterTaxa-Prevalence-TotalCounts.png", width = 8, height = 8, units = "in", res = 300)
ggplot(prevdt, aes(Prevalence, TotalCounts)) + geom_point(size = 2, alpha = 0.5) + 
  scale_y_log10() + theme_bw() + xlab("Prevalence [No. Samples]") + ylab("TotalCounts [Taxa]")
invisible(dev.off())

addPhylum = unique(copy(mdt[, list(TaxaID, Phylum)]))

# Join by TaxaID
setkey(prevdt, TaxaID)
setkey(addPhylum, TaxaID)
prevdt = addPhylum[prevdt]
setkey(prevdt, Phylum)

png("images/FilterTaxa-Prevalence-TotalCounts-Phylum.png", width = 8, height = 8, 
    units = "in", res = 300)
ggplot(prevdt, aes(Prevalence, TotalCounts, color = Phylum)) + 
  geom_point(size = 1, alpha = 0.5) + scale_y_log10() + theme_bw() + 
  facet_wrap(~Phylum, nrow = 4) + theme(legend.position="none") + 
  xlab("Prevalence [No. Samples]") + ylab("Total Abundance")
invisible(dev.off())

#Perform taxa filtering

prevalenceThreshold = 6
abundanceThreshold = 2
maxThreshold = 1

keepTaxa = prevdt[(Prevalence > prevalenceThreshold & 
                     TotalCounts > abundanceThreshold & 
                     MaxCounts > maxThreshold), TaxaID]

ps1 = prune_taxa(keepTaxa, ps0)

phy_tree(ps1) = ape::root(phy_tree(ps1), sample(taxa_names(ps1), 1), resolve.root = TRUE)

ps1

#Create FASTA file
dada2::uniquesToFasta(df[rownames(df) %in% taxa_names(ps1),], "outfiles/expr.asv.fasta", 
                      ids = df[rownames(df) %in% taxa_names(ps1),]$id)

#Create BIOM file
biomformat::write_biom(biomformat::make_biom(data = t(as.matrix(otu_table(ps1)))), 
                       "outfiles/expr.biom")

#Create ASV and Taxonomy tables
write.table(as.data.table(otu_table(ps1), keep.rownames = T), file = "outfiles/expr.otu_table.txt", 
            sep = "\t", quote = F, row.names = F, col.names = T)
write.table(as.data.table(tax_table(ps1), keep.rownames = T), file = "outfiles/expr.tax_table.txt", 
            sep = "\t", quote = F, row.names = F, col.names = T)

#Create raw abundance tables
write.table(ps1 %>% psmelt() %>% arrange(OTU) %>% rename(ASV = OTU) %>% 
              select(ASV, Kingdom, Phylum, Class, Order, Family, Genus, Sample, Abundance) %>% 
              spread(Sample, Abundance),
            file = "outfiles/expr.abundance.all.txt", sep = "\t", 
            quote = F, row.names = F, col.names = T)
write.table(ps1 %>% tax_glom(taxrank = "Phylum") %>% psmelt() %>% 
              select(Phylum, Sample, Abundance) %>% spread(Sample, Abundance),
            file = "outfiles/expr.abundance.abphy.txt", sep = "\t", 
            quote = F, row.names = F, col.names = T)
write.table(ps1 %>% tax_glom(taxrank = "Class") %>% psmelt() %>% 
              select(Class, Sample, Abundance) %>% spread(Sample, Abundance),
            file = "outfiles/expr.abundance.abcls.txt", sep = "\t", 
            quote = F, row.names = F, col.names = T)
write.table(ps1 %>% tax_glom(taxrank = "Family") %>% psmelt() %>%
              select(Family, Sample, Abundance) %>% spread(Sample, Abundance),
            file = "outfiles/expr.abundance.abfam.txt", sep = "\t", 
            quote = F, row.names = F, col.names = T)
write.table(ps1 %>% tax_glom(taxrank = "Genus") %>% psmelt() %>%
              select(Genus, Sample, Abundance) %>% spread(Sample, Abundance),
            file = "outfiles/expr.abundance.abgen.txt", sep = "\t", 
            quote = F, row.names = F, col.names = T)
  
#Create relative abundance tables
write.table(ps1 %>% transform_sample_counts(function(x) {x/sum(x)} ) %>% psmelt() %>% 
              arrange(OTU) %>% rename(ASV = OTU) %>% 
              select(ASV, Kingdom, Phylum, Class, Order, Family, Genus, Sample, Abundance) %>%
              spread(Sample, Abundance), 
            file = "outfiles/expr.relative_abundance.all.txt", 
            sep = "\t", quote = F, row.names = F, col.names = T)
write.table(ps1 %>% tax_glom(taxrank = "Phylum") %>% 
              transform_sample_counts(function(x) {x/sum(x)} ) %>% psmelt() %>%
              select(Phylum, Sample, Abundance) %>% spread(Sample, Abundance), 
            file = "outfiles/expr.relative_abundance.abphy.txt", 
            sep = "\t", quote = F, row.names = F, col.names = T)
write.table(ps1 %>% tax_glom(taxrank = "Class") %>%
              transform_sample_counts(function(x) {x/sum(x)} ) %>% psmelt() %>% 
              select(Class, Sample, Abundance) %>% spread(Sample, Abundance), 
            file = "outfiles/expr.relative_abundance.abcls.txt", 
            sep = "\t", quote = F, row.names = F, col.names = T)
write.table(ps1 %>% tax_glom(taxrank = "Family") %>%
              transform_sample_counts(function(x) {x/sum(x)} ) %>% psmelt() %>% 
              select(Family, Sample, Abundance) %>% spread(Sample, Abundance), 
            file = "outfiles/expr.relative_abundance.abfam.txt", 
            sep = "\t", quote = F, row.names = F, col.names = T)
write.table(ps1 %>% tax_glom(taxrank = "Genus") %>% 
              transform_sample_counts(function(x) {x/sum(x)} ) %>% psmelt() %>% 
              select(Genus, Sample, Abundance) %>% spread(Sample, Abundance), 
            file = "outfiles/expr.relative_abundance.abgen.txt", 
            sep = "\t", quote = F, row.names = F, col.names = T)

#################################################################################
# CHANGE METADATA NAMES IN ps1
library(tidyverse)
ps1@sam_data

ps1@sam_data$treatment <- sub('Control','CN',ps1@sam_data$treatment)
ps1@sam_data$treatment <- sub('MeHG','MeHg',ps1@sam_data$treatment)
ps1@sam_data$treatment <- sub('CoExpo','VitA+MeHg',ps1@sam_data$treatment)

ps1@sam_data
#################################################################################
# Transform to proportions (relative abundances)
ps1.rp = transform_sample_counts(ps1, function(OTU) OTU/sum(OTU))

# Top N taxa
N = 1000
topN = names(sort(taxa_sums(ps1), decreasing=TRUE))[1:N]
ps1.topN = transform_sample_counts(ps1, function(OTU) OTU/sum(OTU))
ps1.topN = prune_taxa(topN, ps1.topN)

ps1.topN



#At Phylum level
# All taxa
group_level <- c("CN", "VitA", "MeHg", "VitA+MeHg")
sex_level <- c("Female","Male")
phyla_level <- c('Firmicutes','Bacteroidota','Proteobacteria','Actinobacteriota','Campylobacterota','Desulfobacterota','Verrucomicrobiota','Patescibacteria')

phylum_ps <- tax_glom(ps1.rp, "Phylum", NArm = T)
phylum_ps_F1 = subset_samples(phylum_ps, generation == "F1")


ptaxa1 = plot_bar(phylum_ps_F1, x = "Sample", fill = "Phylum", 
                  title = paste(ntaxa(ps1.rp), "Taxa (Phylum)")) +
  geom_bar(stat = "identity", size = 0.1,width = 1) +
  facet_wrap(~ sex+treatment, scales = "free_x", nrow = 1) +
  #guides(fill = guide_legend(ncol = 1)) +
  scale_fill_brewer(palette = "Paired") +
  scale_y_continuous(labels = scales::percent) +
  theme(axis.text.x = element_text(angle = 90, size = 6),
        axis.text.y = element_text(color = "black"),
        strip.text = element_text(face = "bold"),
        strip.background = element_blank()) + 
  xlab("SampleID") + ylab("Relative abundance")
  
ptaxa1$data$treatment <- factor(ptaxa1$data$treatment, levels = group_level)
ptaxa1$data$sex <- factor(ptaxa1$data$sex, levels = sex_level)
ptaxa1$data$Phylum <- factor(ptaxa1$data$Phylum, levels = phyla_level)
ptaxa1

###########################################
tiff("images/Abundacias_phylum_F1.tiff", width = 10, height = 8, units = "in", res = 300)
ptaxa1
invisible(dev.off())
##########################################
# Top 200 taxa
ptaxa2 = plot_bar(ps1.topN, x = "Sample", fill = "Phylum", 
                  title = paste("Top",N, "Taxa (Phylum)")) + 
  geom_bar(stat = "identity", size = 0.1, color = "black") + 
  facet_wrap(~generation + treatment, scales = "free_x", nrow = 1) + guides(fill = guide_legend(ncol = 1)) + 
  scale_y_continuous(breaks = seq(0, 1, 0.1)) + 
  theme(legend.title = element_text(size = 8), legend.text = element_text(size = 7), 
        axis.text.x = element_text(size = 5, angle = 90, vjust = 0.5, hjust = 1), 
        axis.text.y = element_text(size = 12)) + xlab("Sample") + ylab("Relative Abundance")

png("images/new_plot_bar_phylum.png", width = 10, height = 8, units = "in", res = 300)
ptaxa1
invisible(dev.off())

#At Class level
class_ps <- tax_glom(ps1.rp, "Class", NArm = T)
ptaxa3 = plot_bar(class_ps, x = "Sample", fill = "Class", 
                  title = paste(ntaxa(ps1.rp), "Taxa (Class)")) +
  geom_bar(stat = "identity", size = 0.1, color = "black") +
  facet_wrap(~generation + treatment, scales = "free_x", nrow = 1) + guides(fill = guide_legend(ncol = 1)) +
  scale_y_continuous(breaks = seq(0, 1, 0.1)) +
  theme(legend.title = element_text(size = 8), legend.text = element_text(size = 7),
        axis.text.x = element_text(size = 5, angle = 90, vjust = 0.5, hjust = 1),
        axis.text.y = element_text(size = 12)) + xlab("Sample") + ylab("Relative abundance")
ptaxa3$data$treatment  <- factor(ptaxa3$data$treatment , levels = group_level)

ptaxa4 = plot_bar(ps1.rp, x = "Sample", fill = "Class", 
                  title = paste("Top",N, "Taxa (Class)")) +
  geom_bar(stat = "identity", size = 0.1, color = "black") +
  facet_wrap(~generation + treatment, scales = "free_x", nrow = 1) + guides(fill = guide_legend(ncol = 1)) +
  scale_y_continuous(breaks = seq(0, 1, 0.1)) +
  theme(legend.title = element_text(size = 8), legend.text = element_text(size = 7),
        axis.text.x = element_text(size = 5, angle = 90, vjust = 0.5, hjust = 1),
        axis.text.y = element_text(size = 12)) + xlab("Sample") + ylab("Relative Abundance")

png("images/plot_bar_class.png", width = 10, height = 8, units = "in", res = 300)
ptaxa3
invisible(dev.off())

#At Family level
family_ps <- tax_glom(ps1.rp, "Family", NArm = T)

ptaxa5 = plot_bar(family_ps, x = "Sample", fill = "Family", 
                  title = paste(ntaxa(ps1.rp), "Taxa (Family)")) +
  geom_bar(stat = "identity", size = 0.1, color = "black") +
  facet_wrap(~generation + treatment, scales = "free_x", nrow = 1) + guides(fill = guide_legend(ncol = 1)) +
  scale_y_continuous(breaks = seq(0, 1, 0.1)) + 
  theme(legend.title = element_text(size = 8), legend.text = element_text(size = 7),
        axis.text.x = element_text(size = 5, angle = 90, vjust = 0.5, hjust = 1),
        axis.text.y = element_text(size = 12)) + xlab("Sample") + ylab("Relative Abundance")
ptaxa5$data$treatment  <- factor(ptaxa5$data$treatment , levels = group_level)

png("images/plot_bar_family.png", width = 10, height = 8, units = "in", res = 300)
ptaxa5
invisible(dev.off())

#At Genus level
genus_ps <- tax_glom(ps1.rp, "Genus", NArm = T)
ptaxa6 = plot_bar(genus_ps, x = "Sample", fill = "Genus", 
                  title = paste(ntaxa(ps1.rp), "Taxa (Genus)")) +
  geom_bar(stat = "identity", size = 0.1, color = "black") +
  facet_wrap(~generation + treatment, scales = "free_x", nrow = 1) + guides(fill = guide_legend(ncol = 1)) +
  scale_y_continuous(breaks = seq(0, 1, 0.1)) + 
  theme(legend.title = element_text(size = 8), legend.text = element_text(size = 7),
        axis.text.x = element_text(size = 5, angle = 90, vjust = 0.5, hjust = 1),
        axis.text.y = element_text(size = 12)) + xlab("Sample") + ylab("Relative Abundance")
ptaxa6$data$treatment  <- factor(ptaxa6$data$treatment , levels = group_level)

png("images/plot_bar_genus.png", width = 10, height = 8, units = "in", res = 300)
ptaxa6
invisible(dev.off())


#Beta diversity
unifracs = GUniFrac::GUniFrac(otu_table(ps1), phy_tree(ps1), alpha = c(0, 0.5, 1))$unifracs

# Calculate geometric means prior to estimate size factors
gm_mean = function(x, na.rm = TRUE){
  exp(sum(log(x[x > 0]), na.rm = na.rm) / length(x))
}

# Perform variance stabilizing transformation
ps1.vsd = subset_samples(ps1, generation == "F1")

dds = phyloseq_to_deseq2(ps1, ~ treatment)

geoMeans = apply(counts(dds), 1, gm_mean)
dds = estimateSizeFactors(dds, geoMeans = geoMeans)
dds = estimateDispersions(dds)

abund.vsd = getVarianceStabilizedData(dds)
abund.vsd[abund.vsd < 0] = 0    # set negative values to 0
otu_table(ps1.vsd) = otu_table(abund.vsd, taxa_are_rows = TRUE)

# Create distance objects
dist_un = as.dist(unifracs[, , "d_UW"])       # Unweighted UniFrac
attr(dist_un, "method") = "Unweighted UniFrac"

dist_wu = as.dist(unifracs[, , "d_1"])        # Weighted UniFrac
attr(dist_wu, "method") = "Weighted UniFrac"

# dist_vu = as.dist(unifracs[, , "d_VAW"])      # Variance-adjusted-weighted UniFrac
# attr(dist_vu, "method") = "Variance-adjusted-weighted UniFrac"

dist_gu = as.dist(unifracs[, , "d_0.5"])      # GUniFrac with alpha 0.5
attr(dist_gu, "method") = "GUniFrac with alpha 0.5"

dist_bc = phyloseq::distance(ps1.vsd, "bray") # Bray-Curtis
dist_dc = phyloseq::distance(ps1.vsd, "dpcoa") # DPCoA
dist_eu = phyloseq::distance(ps1.vsd, "euclidean") # euclidean distance


# Perform ordination
ord_un = ordinate(ps1, method = "PCoA", distance = dist_un)
ord_wu = ordinate(ps1, method = "PCoA", distance = dist_wu)
# ord_vu = ordinate(ps1, method = "PCoA", distance = dist_vu)
ord_gu = ordinate(ps1, method = "PCoA", distance = dist_gu)
set.seed(123)
ord_bc = ordinate(ps1, method = "NMDS", distance = dist_bc)
ord_dc = ordinate(ps1, method = "PCoA", distance = dist_eu)

tiff("images/f1_plot_ordination.tiff", width = 6, height = 5, units = "in", res = 300)
plot_ordination(ps1, ord_un, type = "samples", color = "treatment", shape = "sex", 
                title = "PCoA on unweighted UniFrac") + 
  theme_bw() + coord_fixed(ratio = 1) + 
  stat_ellipse(aes(group = sex), type = "t", size = 0.5, linetype = "dashed")

plot_ordination(ps1, ord_wu, type = "samples", color = "treatment", shape = "sex", 
                title = "PCoA on weighted UniFrac") + 
  theme_bw() + coord_fixed(ratio = 1) + 
  stat_ellipse(aes(group = treatment), type = "t", size = 0.5, linetype = "dashed")

plot_ordination(ps1, ord_gu, type = "samples", color = "treatment", shape = "sex", 
                title = "PCoA on GUniFrac with alpha 0.5") + 
  theme_bw() + coord_fixed(ratio = 1) + 
  stat_ellipse(aes(group = treatment), type = "t", size = 0.5, linetype = "dashed")

plot_ordination(ps1, ord_bc, type = "samples", color = "treatment", shape = "sex", 
                title = "NMDS on Bray-Curtis distance (VST)") + 
  theme_bw() + coord_fixed(ratio = 1) + 
  stat_ellipse(aes(group = treatment), type = "t", size = 0.5, linetype = "dashed")

invisible(dev.off())

#Multivariate analyses
#PERMANOVA
set.seed(123)
adonis_un <- adonis2(dist_un ~ Stress + Treatment, data.frame(sample_data(ps1.vsd)), permutations = 1000)

set.seed(123)
adonis_wu <- adonis2(dist_wu ~ Stress + Treatment, data.frame(sample_data(ps1.vsd)), permutations = 1000)

# set.seed(123)
# adonis_vu <- adonis(dist_vu ~ Sex_Stress_Treatment, data.frame(sample_data(ps1)), permutations = 1000)

set.seed(123)
adonis_gu <- adonis2(dist_gu ~ Stress + Treatment, data.frame(sample_data(ps1.vsd)), permutations = 1000)

set.seed(123)
adonis_bc <- adonis2(dist_bc ~ Stress + Treatment, data.frame(sample_data(ps1.vsd)), permutations = 1000)

#Multivariate dispersion
disp1 = betadisper(dist_un, group = data.frame(sample_data(ps1.vsd))$final_group)
set.seed(123)
permutest(disp1, permutations = 1000)

disp2 = betadisper(dist_wu, group = data.frame(sample_data(ps1.vsd))$final_group) 
set.seed(123)
permutest(disp2, permutations = 1000)

# disp3 = betadisper(dist_vu, group = data.frame(sample_data(ps1))$Sex_Stress_Treatment)
# set.seed(123)
# permutest(disp3, permutations = 1000)

disp4 = betadisper(dist_gu, group = data.frame(sample_data(ps1.vsd))$final_group)
set.seed(123)
permutest(disp4, permutations = 1000)

disp5 = betadisper(dist_bc, group = data.frame(sample_data(ps1.vsd))$final_group)
set.seed(123)
permutest(disp5, permutations = 1000)

pdf("images/female_plot_dispersion.pdf", width = 6, height = 5, pointsize = 12)
plot(disp1, hull = FALSE, ellipse = TRUE)
plot(disp2, hull = FALSE, ellipse = TRUE)
# plot(disp3, hull = FALSE, ellipse = TRUE)
plot(disp4, hull = FALSE, ellipse = TRUE)
plot(disp5, hull = FALSE, ellipse = TRUE)
invisible(dev.off())

save.image(file = "phyloseq.RData")


# Select alpha-diversity measures
divIdx = c("Observed", "Chao1", "Shannon", "Simpson")

png("images/plot_richness.png", width = 6, height = 5, units = "in", res = 300)
plot_richness(ps1, x = "Sex_Stress_Treatment", measures = divIdx, color = "Sex_Stress_Treatment", nrow = 1) + 
  geom_point(size = 0.8) + theme_bw() + 
  theme(legend.position = "none", 
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) + 
  labs(x = "Group", y = "Alpha Diversity Measure")
invisible(dev.off())

ad = estimate_richness(ps1, measures = divIdx)
rownames(ad) = lapply(rownames(ad), gsub, pattern="X", replace="") #remove "x" in front of rownames
ad = merge(data.frame(sample_data(ps1)), ad, by = "row.names")
ad = ad %>% select(Sample.id, Sex_Stress_Treatment, all_of(divIdx)) %>% 
  gather(key = "alpha", value = "measure", -c(Sample.id, Sex_Stress_Treatment))

png("images/plot_richness_boxplot.png", width = 6, height = 5, units = "in", res = 300)
ggplot(ad, aes(Sex_Stress_Treatment, measure, color = Sex_Stress_Treatment)) + 
  geom_boxplot(outlier.shape = NA, size = 0.8, width = 0.8) + 
  geom_quasirandom(size = 0.8, color = "black") + theme_bw() + 
  facet_wrap(~ alpha, scales = "free_y", nrow = 1) + 
  theme(legend.position = "none", 
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) + 
  labs(x = "Group", y = "Alpha Diversity Measure")
invisible(dev.off())

