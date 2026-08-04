#Source: https://ycl6.github.io/16S-Demo/index.html
#BiocManager::install("ALDEx2")

#Running picrust2 on console
#recomended installation
# mamba create -n picrust2 -c bioconda -c conda-forge picrust2=2.5.0
# it is recommended to install mamba into conda's base environment

#Use terminal to run the following commands
conda activate picrust2
picrust2_pipeline.py --study_fasta outfiles/expr.asv.fasta --input outfiles/expr.biom \
--output picrust2_out_stratified --processes 4 --stratified --remove_intermediate --verbose

#Working on R
library("data.table")
library("phyloseq")
library("ALDEx2")
library("tidyverse")
load("lefse.RData")

if(!dir.exists("outfiles")) { dir.create("outfiles") }
if(!dir.exists("images")) { dir.create("images") }

#Importing files from picrust2
picrust2 = "picrust2_out_stratified"
# list.files(picrust2, recursive = TRUE)

p2_EC = paste0(picrust2, "/EC_metagenome_out/pred_metagenome_unstrat.tsv.gz")
p2_KO = paste0(picrust2, "/KO_metagenome_out/pred_metagenome_unstrat.tsv.gz")
p2_PW = paste0(picrust2, "/pathways_out/path_abun_unstrat.tsv.gz")

mapfile = "/home/luiz/miniconda3/envs/picrust2/lib/python3.8/site-packages/picrust2/default_files/description_mapfiles"
# list.files(mapfile, recursive = TRUE)
mapfile_EC = paste0(mapfile, "/ec_level4_info.tsv.gz")
mapfile_KO = paste0(mapfile, "/ko_info.tsv.gz")
mapfile_PW = paste0(mapfile, "/metacyc_pathways_info.txt.gz")


#Preparing input
mapEC = as.data.frame(fread(mapfile_EC, header = FALSE))
colnames(mapEC) = c("function","description")
mapKO = as.data.frame(fread(mapfile_KO, header = FALSE, sep = "\t"))
colnames(mapKO) = c("function","description")
mapPW = as.data.frame(fread(mapfile_PW, header = FALSE))
colnames(mapPW) = c("pathway","description")

p2EC = as.data.frame(fread(p2_EC))
rownames(p2EC) = p2EC$"function"
p2EC = as.matrix(p2EC[,-1])
p2EC = round(p2EC)

p2KO = as.data.frame(fread(p2_KO))
rownames(p2KO) = p2KO$"function"
p2KO = as.matrix(p2KO[,-1])
p2KO = round(p2KO)

p2PW = as.data.frame(fread(p2_PW))
rownames(p2PW) = p2PW$"pathway"
p2PW = as.matrix(p2PW[,-1])
p2PW = round(p2PW)

#Subset control vs depressive
p2EC1 = p2EC[,sample_names(ps1a)]
p2KO1 = p2KO[,sample_names(ps1a)]
p2PW1 = p2PW[,sample_names(ps1a)]

#Subset control vs resilient
p2EC2 = p2EC[,sample_names(ps1b)]
p2KO2 = p2KO[,sample_names(ps1b)]
p2PW2 = p2PW[,sample_names(ps1b)]

#Subset depressive vs resilient
p2EC3 = p2EC[,sample_names(ps1c)]
p2KO3 = p2KO[,sample_names(ps1c)]
p2PW3 = p2PW[,sample_names(ps1c)]

#SUBSET-1
#EC
set.seed(12345)
system.time({
  aldex2_EC1 = aldex(p2EC1, sample_data(ps1a)$Description, mc.samples = 500, test = "t", 
                     effect = TRUE, denom = "iqlr", verbose = TRUE)
})

# KO
set.seed(12345)
system.time({
  aldex2_KO1 = aldex(p2KO1, sample_data(ps1a)$Description, mc.samples = 500, test = "t", 
                     effect = TRUE, denom = "iqlr", verbose = TRUE)
})

# Pathway
set.seed(12345)
system.time({
  aldex2_PW1 = aldex(p2PW1, sample_data(ps1a)$Description, mc.samples = 500, test = "t", 
                     effect = TRUE, denom = "iqlr", verbose = TRUE)
})

#SUBSET2
# EC
set.seed(12345)
system.time({
  aldex2_EC2 = aldex(p2EC2, sample_data(ps1b)$Description, mc.samples = 500, test = "t", 
                     effect = TRUE, denom = "iqlr", verbose = TRUE)
})

# KO
set.seed(12345)
system.time({
  aldex2_KO2 = aldex(p2KO2, sample_data(ps1b)$Description, mc.samples = 500, test = "t", 
                     effect = TRUE, denom = "iqlr", verbose = TRUE)
})

# Pathway
set.seed(12345)
system.time({
  aldex2_PW2 = aldex(p2PW2, sample_data(ps1b)$Description, mc.samples = 500, test = "t", 
                     effect = TRUE, denom = "iqlr", verbose = TRUE)
})

#SUBSET3
# EC
set.seed(12345)
system.time({
  aldex2_EC3 = aldex(p2EC3, sample_data(ps1c)$Description, mc.samples = 500, test = "t", 
                     effect = TRUE, denom = "iqlr", verbose = TRUE)
})

# KO
set.seed(12345)
system.time({
  aldex2_KO3 = aldex(p2KO3, sample_data(ps1c)$Description, mc.samples = 500, test = "t", 
                     effect = TRUE, denom = "iqlr", verbose = TRUE)
})

# Pathway
set.seed(12345)
system.time({
  aldex2_PW3 = aldex(p2PW3, sample_data(ps1c)$Description, mc.samples = 500, test = "t", 
                     effect = TRUE, denom = "iqlr", verbose = TRUE)
})

#Check estimated effect size
png("images/ALDEx2_picrust2_effect_1.png", width = 6, height = 6, units = "in", res = 300)
par(mfrow = c(2,2))
hist(aldex2_EC1$effect, breaks = 20, xlab = "effect size", col = "yellow", main = "EC")
hist(aldex2_KO1$effect, breaks = 20, xlab = "effect size", col = "yellow", main = "KO")
hist(aldex2_PW1$effect, breaks = 20, xlab = "effect size", col = "yellow", main = "Pathway")
invisible(dev.off())

png("images/ALDEx2_picrust2_effect_2.png", width = 6, height = 6, units = "in", res = 300)
par(mfrow = c(2,2))
hist(aldex2_EC2$effect, breaks = 20, xlab = "effect size", col = "yellow", main = "EC")
hist(aldex2_KO2$effect, breaks = 20, xlab = "effect size", col = "yellow", main = "KO")
hist(aldex2_PW2$effect, breaks = 20, xlab = "effect size", col = "yellow", main = "Pathway")
invisible(dev.off())

png("images/ALDEx2_picrust2_effect_3.png", width = 6, height = 6, units = "in", res = 300)
par(mfrow = c(2,2))
hist(aldex2_EC3$effect, breaks = 20, xlab = "effect size", col = "yellow", main = "EC")
hist(aldex2_KO3$effect, breaks = 20, xlab = "effect size", col = "yellow", main = "KO")
hist(aldex2_PW3$effect, breaks = 20, xlab = "effect size", col = "yellow", main = "Pathway")
invisible(dev.off())

#Create MW and MA plots
#SUBSET1
png("images/ALDEx2_picrust2_MW_MA_1.png", width = 6, height = 8, units = "in", res = 300)
par(mfrow = c(3,2))
aldex.plot(aldex2_EC1, type = "MW", test = "wilcox", all.cex = 0.4, rare.cex = 0.4, 
           called.cex = 0.6, cutoff = 0.05, xlab = "Dispersion", ylab = "Difference")
title(main = "(EC) MW Plot")

aldex.plot(aldex2_EC1, type = "MA", test = "wilcox", all.cex = 0.4, rare.cex = 0.4, 
           called.cex = 0.6, cutoff = 0.05, xlab = "Log-ratio abundance", ylab = "Difference")
title(main = "(EC) MA Plot")

aldex.plot(aldex2_KO1, type = "MW", test = "wilcox", all.cex = 0.4, rare.cex = 0.4, 
           called.cex = 0.6, cutoff = 0.05, xlab = "Dispersion", ylab = "Difference")
title(main = "(KO) MW Plot")

aldex.plot(aldex2_KO1, type = "MA", test = "wilcox", all.cex = 0.4, rare.cex = 0.4, 
           called.cex = 0.6, cutoff = 0.05, xlab = "Relative abundance", ylab = "Difference")
title(main = "(KO) MA Plot")

aldex.plot(aldex2_PW1, type = "MW", test = "wilcox", all.cex = 0.4, rare.cex = 0.4, 
           called.cex = 0.6, cutoff = 0.05, xlab = "Dispersion", ylab = "Difference")
title(main = "(PW) MW Plot")

aldex.plot(aldex2_PW1, type = "MA", test = "wilcox", all.cex = 0.4, rare.cex = 0.4, 
           called.cex = 0.6, cutoff = 0.05, xlab = "Relative abundance", ylab = "Difference")
title(main = "(PW) MA Plot")
invisible(dev.off())

#SUBSET2
png("images/ALDEx2_picrust2_MW_MA_2.png", width = 6, height = 8, units = "in", res = 300)
par(mfrow = c(3,2))
aldex.plot(aldex2_EC2, type = "MW", test = "wilcox", all.cex = 0.4, rare.cex = 0.4, 
           called.cex = 0.6, cutoff = 0.05, xlab = "Dispersion", ylab = "Difference")
title(main = "(EC) MW Plot")

aldex.plot(aldex2_EC2, type = "MA", test = "wilcox", all.cex = 0.4, rare.cex = 0.4, 
           called.cex = 0.6, cutoff = 0.05, xlab = "Relative abundance", ylab = "Difference")
title(main = "(EC) MA Plot")

aldex.plot(aldex2_KO2, type = "MW", test = "wilcox", all.cex = 0.4, rare.cex = 0.4, 
           called.cex = 0.6, cutoff = 0.05, xlab = "Dispersion", ylab = "Difference")
title(main = "(KO) MW Plot")

aldex.plot(aldex2_KO2, type = "MA", test = "wilcox", all.cex = 0.4, rare.cex = 0.4, 
           called.cex = 0.6, cutoff = 0.05, xlab = "Relative abundance", ylab = "Difference")
title(main = "(KO) MA Plot")

aldex.plot(aldex2_PW2, type = "MW", test = "wilcox", all.cex = 0.4, rare.cex = 0.4, 
           called.cex = 0.6, cutoff = 0.05, xlab = "Dispersion", ylab = "Difference")
title(main = "(PW) MW Plot")

aldex.plot(aldex2_PW2, type = "MA", test = "wilcox", all.cex = 0.4, rare.cex = 0.4, 
           called.cex = 0.6, cutoff = 0.05, xlab = "Relative abundance", ylab = "Difference")
title(main = "(PW) MA Plot")
invisible(dev.off())

#SUBSET3
png("images/ALDEx2_picrust2_MW_MA_3.png", width = 6, height = 8, units = "in", res = 300)
par(mfrow = c(3,2))
aldex.plot(aldex2_EC3, type = "MW", test = "wilcox", all.cex = 0.4, rare.cex = 0.4, 
           called.cex = 0.6, cutoff = 0.05, xlab = "Dispersion", ylab = "Difference")
title(main = "(EC) MW Plot")

aldex.plot(aldex2_EC3, type = "MA", test = "wilcox", all.cex = 0.4, rare.cex = 0.4, 
           called.cex = 0.6, cutoff = 0.05, xlab = "Relative abundance", ylab = "Difference")
title(main = "(EC) MA Plot")

aldex.plot(aldex2_KO3, type = "MW", test = "wilcox", all.cex = 0.4, rare.cex = 0.4, 
           called.cex = 0.6, cutoff = 0.05, xlab = "Dispersion", ylab = "Difference")
title(main = "(KO) MW Plot")

aldex.plot(aldex2_KO3, type = "MA", test = "wilcox", all.cex = 0.4, rare.cex = 0.4, 
           called.cex = 0.6, cutoff = 0.05, xlab = "Relative abundance", ylab = "Difference")
title(main = "(KO) MA Plot")

aldex.plot(aldex2_PW3, type = "MW", test = "wilcox", all.cex = 0.4, rare.cex = 0.4, 
           called.cex = 0.6, cutoff = 0.05, xlab = "Dispersion", ylab = "Difference")
title(main = "(PW) MW Plot")

aldex.plot(aldex2_PW3, type = "MA", test = "wilcox", all.cex = 0.4, rare.cex = 0.4, 
           called.cex = 0.6, cutoff = 0.05, xlab = "Relative abundance", ylab = "Difference")
title(main = "(PW) MA Plot")
invisible(dev.off())

#Relationship between effect, difference, and P values

#SUBSET1
png("images/ALDEx2_picrust2_P_adjP_1.png", width = 6, height = 8, units = "in", res = 300)
par(mfrow = c(3,2))
plot(aldex2_EC1$effect, aldex2_EC1$we.ep, log = "y", cex = 0.4, col = "blue", pch = 19, 
     xlab = "Effect size", ylab = "P value", main = "(EC) Effect size plot")
points(aldex2_EC1$effect, aldex2_EC1$we.eBH, cex = 0.5, col = "red", pch = 19)
abline(h = 0.05, lty = 2, col = "grey")
legend("bottom", legend = c("P value", "BH-adjusted"), pch = 19, col = c("blue", "red"))

plot(aldex2_EC1$diff.btw, aldex2_EC1$we.ep, log = "y", cex = 0.4, col = "blue", pch = 19, 
     xlab = "Difference", ylab = "P value", main = "(EC) Volcano plot")
points(aldex2_EC1$diff.btw, aldex2_EC1$we.eBH, cex = 0.5, col = "red", pch = 19)
abline(h = 0.05, lty = 2, col = "grey")

plot(aldex2_KO1$effect, aldex2_KO1$we.ep, log = "y", cex = 0.4, col = "blue", pch = 19, 
     xlab = "Effect size", ylab = "P value", main = "(KO) Effect size plot")
points(aldex2_KO1$effect, aldex2_KO1$we.eBH, cex = 0.5, col = "red", pch = 19)
abline(h = 0.05, lty = 2, col = "grey")
legend("bottom", legend = c("P value", "BH-adjusted"), pch = 19, col = c("blue", "red"))

plot(aldex2_KO1$diff.btw, aldex2_KO1$we.ep, log = "y", cex = 0.4, col = "blue", pch = 19, 
     xlab = "Difference", ylab = "P value", main = "(KO) Volcano plot")
points(aldex2_KO1$diff.btw, aldex2_KO1$we.eBH, cex = 0.5, col = "red", pch = 19)
abline(h = 0.05, lty = 2, col = "grey")

plot(aldex2_PW1$effect, aldex2_PW1$we.ep, log = "y", cex = 0.4, col = "blue", pch = 19, 
     xlab = "Effect size", ylab = "P value", main = "(PW) Effect size plot")
points(aldex2_PW1$effect, aldex2_PW1$we.eBH, cex = 0.5, col = "red", pch = 19)
abline(h = 0.05, lty = 2, col = "grey")
legend("bottom", legend = c("P value", "BH-adjusted"), pch = 19, col = c("blue", "red"))

plot(aldex2_PW1$diff.btw, aldex2_PW1$we.ep, log = "y", cex = 0.4, col = "blue", pch = 19, 
     xlab = "Difference", ylab = "P value", main = "(PW) Volcano plot")
points(aldex2_PW1$diff.btw, aldex2_PW1$we.eBH, cex = 0.5, col = "red", pch = 19)
abline(h = 0.05, lty = 2, col = "grey")
invisible(dev.off())

#SUBSET2
png("images/ALDEx2_picrust2_P_adjP_2.png", width = 6, height = 8, units = "in", res = 300)
par(mfrow = c(3,2))
plot(aldex2_EC2$effect, aldex2_EC2$we.ep, log = "y", cex = 0.4, col = "blue", pch = 19, 
     xlab = "Effect size", ylab = "P value", main = "(EC) Effect size plot")
points(aldex2_EC2$effect, aldex2_EC2$we.eBH, cex = 0.5, col = "red", pch = 19)
abline(h = 0.05, lty = 2, col = "grey")
legend("bottom", legend = c("P value", "BH-adjusted"), pch = 19, col = c("blue", "red"))

plot(aldex2_EC2$diff.btw, aldex2_EC2$we.ep, log = "y", cex = 0.4, col = "blue", pch = 19, 
     xlab = "Difference", ylab = "P value", main = "(EC) Volcano plot")
points(aldex2_EC2$diff.btw, aldex2_EC2$we.eBH, cex = 0.5, col = "red", pch = 19)
abline(h = 0.05, lty = 2, col = "grey")

plot(aldex2_KO2$effect, aldex2_KO2$we.ep, log = "y", cex = 0.4, col = "blue", pch = 19, 
     xlab = "Effect size", ylab = "P value", main = "(KO) Effect size plot")
points(aldex2_KO2$effect, aldex2_KO2$we.eBH, cex = 0.5, col = "red", pch = 19)
abline(h = 0.05, lty = 2, col = "grey")
legend("bottom", legend = c("P value", "BH-adjusted"), pch = 19, col = c("blue", "red"))

plot(aldex2_KO2$diff.btw, aldex2_KO2$we.ep, log = "y", cex = 0.4, col = "blue", pch = 19, 
     xlab = "Difference", ylab = "P value", main = "(KO) Volcano plot")
points(aldex2_KO2$diff.btw, aldex2_KO2$we.eBH, cex = 0.5, col = "red", pch = 19)
abline(h = 0.05, lty = 2, col = "grey")

plot(aldex2_PW2$effect, aldex2_PW2$we.ep, log = "y", cex = 0.4, col = "blue", pch = 19, 
     xlab = "Effect size", ylab = "P value", main = "(PW) Effect size plot")
points(aldex2_PW2$effect, aldex2_PW2$we.eBH, cex = 0.5, col = "red", pch = 19)
abline(h = 0.05, lty = 2, col = "grey")
legend("bottom", legend = c("P value", "BH-adjusted"), pch = 19, col = c("blue", "red"))

plot(aldex2_PW2$diff.btw, aldex2_PW2$we.ep, log = "y", cex = 0.4, col = "blue", pch = 19, 
     xlab = "Difference", ylab = "P value", main = "(PW) Volcano plot")
points(aldex2_PW2$diff.btw, aldex2_PW2$we.eBH, cex = 0.5, col = "red", pch = 19)
abline(h = 0.05, lty = 2, col = "grey")
invisible(dev.off())

#SUBSET3
png("images/ALDEx2_picrust2_P_adjP_3.png", width = 6, height = 8, units = "in", res = 300)
par(mfrow = c(3,2))
plot(aldex2_EC3$effect, aldex2_EC3$we.ep, log = "y", cex = 0.4, col = "blue", pch = 19, 
     xlab = "Effect size", ylab = "P value", main = "(EC) Effect size plot")
points(aldex2_EC3$effect, aldex2_EC3$we.eBH, cex = 0.5, col = "red", pch = 19)
abline(h = 0.05, lty = 2, col = "grey")
legend("bottom", legend = c("P value", "BH-adjusted"), pch = 19, col = c("blue", "red"))

plot(aldex2_EC3$diff.btw, aldex2_EC3$we.ep, log = "y", cex = 0.4, col = "blue", pch = 19, 
     xlab = "Difference", ylab = "P value", main = "(EC) Volcano plot")
points(aldex2_EC3$diff.btw, aldex2_EC3$we.eBH, cex = 0.5, col = "red", pch = 19)
abline(h = 0.05, lty = 2, col = "grey")

plot(aldex2_KO3$effect, aldex2_KO3$we.ep, log = "y", cex = 0.4, col = "blue", pch = 19, 
     xlab = "Effect size", ylab = "P value", main = "(KO) Effect size plot")
points(aldex2_KO3$effect, aldex2_KO3$we.eBH, cex = 0.5, col = "red", pch = 19)
abline(h = 0.05, lty = 2, col = "grey")
legend("bottom", legend = c("P value", "BH-adjusted"), pch = 19, col = c("blue", "red"))

plot(aldex2_KO3$diff.btw, aldex2_KO3$we.ep, log = "y", cex = 0.4, col = "blue", pch = 19, 
     xlab = "Difference", ylab = "P value", main = "(KO) Volcano plot")
points(aldex2_KO3$diff.btw, aldex2_KO3$we.eBH, cex = 0.5, col = "red", pch = 19)
abline(h = 0.05, lty = 2, col = "grey")

plot(aldex2_PW3$effect, aldex2_PW3$we.ep, log = "y", cex = 0.4, col = "blue", pch = 19, 
     xlab = "Effect size", ylab = "P value", main = "(PW) Effect size plot")
points(aldex2_PW3$effect, aldex2_PW3$we.eBH, cex = 0.5, col = "red", pch = 19)
abline(h = 0.05, lty = 2, col = "grey")
legend("bottom", legend = c("P value", "BH-adjusted"), pch = 19, col = c("blue", "red"))

plot(aldex2_PW3$diff.btw, aldex2_PW3$we.ep, log = "y", cex = 0.4, col = "blue", pch = 19, 
     xlab = "Difference", ylab = "P value", main = "(PW) Volcano plot")
points(aldex2_PW3$diff.btw, aldex2_PW3$we.eBH, cex = 0.5, col = "red", pch = 19)
abline(h = 0.05, lty = 2, col = "grey")
invisible(dev.off())

#Merge with map file data
#SUBSET1
df_EC1 = aldex2_EC1 %>% tibble::rownames_to_column(var = "EC") %>% 
  inner_join(mapEC, by = c("EC" = "function")) %>% arrange(EC)

df_KO1 = aldex2_KO1 %>% tibble::rownames_to_column(var = "KO") %>% 
  inner_join(mapKO, by = c("KO" = "function")) %>% arrange(KO)

df_PW1 = aldex2_PW1 %>% tibble::rownames_to_column(var = "Pathway") %>% 
  inner_join(mapPW, by = c("Pathway" = "pathway")) %>% arrange(Pathway)

#SUBSET2
df_EC2 = aldex2_EC2 %>% tibble::rownames_to_column(var = "EC") %>%
  inner_join(mapEC, by = c("EC" = "function")) %>% arrange(EC)

df_KO2 = aldex2_KO2 %>% tibble::rownames_to_column(var = "KO") %>%
  inner_join(mapKO, by = c("KO" = "function")) %>% arrange(KO)

df_PW2 = aldex2_PW2 %>% tibble::rownames_to_column(var = "Pathway") %>%
  inner_join(mapPW, by = c("Pathway" = "pathway")) %>% arrange(Pathway)

#SUBSET3
df_EC3 = aldex2_EC3 %>% tibble::rownames_to_column(var = "EC") %>%
  inner_join(mapEC, by = c("EC" = "function")) %>% arrange(EC)

df_KO3 = aldex2_KO3 %>% tibble::rownames_to_column(var = "KO") %>%
  inner_join(mapKO, by = c("KO" = "function")) %>% arrange(KO)

df_PW3 = aldex2_PW3 %>% tibble::rownames_to_column(var = "Pathway") %>%
  inner_join(mapPW, by = c("Pathway" = "pathway")) %>% arrange(Pathway)

#Output to file
#SUBSET1
write.table(df_EC1, file = "outfiles/ALDEx2_picrust2_EC_results_1.txt", sep = "\t", quote = F, 
            row.names = F, col.names = T)
write.table(df_KO1, file = "outfiles/ALDEx2_picrust2_KO_results_1.txt", sep = "\t", quote = F, 
            row.names = F, col.names = T)
write.table(df_PW1, file = "outfiles/ALDEx2_picrust2_Pathway_results_1.txt", sep = "\t", quote = F, 
            row.names = F, col.names = T)

#SUBSET2
write.table(df_EC2, file = "outfiles/ALDEx2_picrust2_EC_results_2.txt", sep = "\t", quote = F, 
            row.names = F, col.names = T)
write.table(df_KO2, file = "outfiles/ALDEx2_picrust2_KO_results_2.txt", sep = "\t", quote = F, 
            row.names = F, col.names = T)
write.table(df_PW2, file = "outfiles/ALDEx2_picrust2_Pathway_results_2.txt", sep = "\t", quote = F, 
            row.names = F, col.names = T)

#SUBSET3
write.table(df_EC3, file = "outfiles/ALDEx2_picrust2_EC_results_3.txt", sep = "\t", quote = F, 
            row.names = F, col.names = T)
write.table(df_KO3, file = "outfiles/ALDEx2_picrust2_KO_results_3.txt", sep = "\t", quote = F, 
            row.names = F, col.names = T)
write.table(df_PW3, file = "outfiles/ALDEx2_picrust2_Pathway_results_3.txt", sep = "\t", quote = F, 
            row.names = F, col.names = T)

save.image(file = "picrust2.RData")
