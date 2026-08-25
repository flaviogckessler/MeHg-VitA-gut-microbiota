library(readr)
library(dplyr)
library(tidyverse)
library(ggplot2)
library(rstatix)
library(ggpubr)

#install.packages("metacoder")
library("metacoder")
#install.packages("agricolae")
library(agricolae)
library(reshape2)
######################################################################################################
p.adjust_bonferroni <- function (p, method = p.adjust.methods, n) 
{
  method <- match.arg(method)
  if (method == "fdr") 
    method <- "BH"
  nm <- names(p)
  p <- as.numeric(p)
  p0 <- setNames(p, nm)
  if (all(nna <- !is.na(p))) 
    nna <- TRUE
  else p <- p[nna]
  lp <- length(p)
  #stopifnot(n >= lp)
  if (n <= 1) 
    return(p0)
  if (n == 2 && method == "hommel") 
    method <- "hochberg"
  p0[nna] <- switch(method,
                    bonferroni = pmin(1, n * p),
                    none = p)
  return(p0)
}
#######
###########################################################################################################
group_level <- c("CN", "VitA", "MeHg", "VitA+MeHg")
sex_level <- c("Female","Male")
phyla_level <- c('Firmicutes','Bacteroidota','Proteobacteria','Actinobacteriota','Campylobacterota','Desulfobacterota','Verrucomicrobiota','Patescibacteria')
###########################################################################################################
metadata_1year <- read.csv(file.path("data", "metadata", "legacy", "metadata.csv"))
rownames(metadata_1year) <- metadata_1year$sampleID

abund_relative_1year <- read_delim("outfiles/expr.relative_abundance.all.txt", delim = "\t")
rownames(abund_relative_1year) <- abund_relative_1year$ASV
###############################################
alpha_f0 <- subset(alpha1,Geracao == "F0")
#####################################################################################
#####################################################################################
abund_relative_1year$lineage <- NA
for(i in 1:length(abund_relative_1year$ASV)) {
  buf = paste("k__",abund_relative_1year$Kingdom[i],";p__",abund_relative_1year$Phylum[i],";c__",abund_relative_1year$Class[i],
              ";o__",abund_relative_1year$Order[i],";f__",abund_relative_1year$Family[i],";g__",abund_relative_1year$Genus[i],sep = "")
  abund_relative_1year$lineage[i] <- buf
}

#write.csv(teste,"output/RelativeAbundance_1year_lineage.csv")
#####################################################################################
#####################################################################################
# https://grunwaldlab.github.io/metacoder_documentation/example.html
l <- metadata_1year$sampleID
l <- append(c("ASV","lineage"),l)
mydata_l_1year <- abund_relative_1year[,l]
mydata_l_1year

obj <- parse_tax_data(mydata_l_1year, class_cols = "lineage", class_sep = ";",           #AQUI mydata
                      class_key = c(tax_rank = "taxon_rank", tax_name = "taxon_name"),
                      class_regex = "^(.+)__(.+)$")

obj
View(obj$data$tax_data)
View(obj$data$class_data)

# Until now, we have per ASV abundance, not taxon abundance. The following function calculate taxons abundances
obj$data$tax_abund <- calc_taxon_abund(obj, "tax_data", cols = as.character(metadata_1year$sampleID))    
obj

# Calculate group differences
obj$data$tax_occ <- calc_n_samples(obj, "tax_abund", groups = metadata_1year$treatment)                                      #AQUI
print(obj)
#Comparing any number of treatments/groups
obj$data$diff_table <- compare_groups(obj, data = "tax_abund",
                                      cols = as.character(metadata_1year$sampleID), # What columns of sample data to use #AQUI
                                      groups = metadata_1year$treatment) # What category each sample is assigned to          #AQUI
obj

# Correção de comparações múltiplas
obj <- mutate_obs(obj, "diff_table",
                  fdr = p.adjust(wilcox_p_value, method = "fdr"))
obj
obj$data$diff_table$fdr <- p.adjust_bonferroni(obj$data$diff_table$wilcox_p_value,"bonferroni",n=3)
obj <- mutate_obs(obj, "cleaned_names", gsub(taxon_names, pattern = "\\[|\\]", replacement = ""))
View(obj$data$diff_table)

# Zera p valores ajustados não significativos
obj$data$diff_table$log2_median_ratio[obj$data$diff_table$fdr > 0.05] <- 0
head(obj$data$diff_table,10)

#Plotagem arvore
obj %>%
  metacoder::filter_taxa(taxon_ranks == "g", supertaxa = TRUE, reassign_obs = FALSE) %>%
  metacoder::filter_taxa(grepl(cleaned_names, pattern = "^[a-zA-Z]+$"), reassign_obs = FALSE) %>%
  heat_tree_matrix(data = "diff_table",
                   node_label = cleaned_names,
                   node_size = n_obs, # number of OTUs
                   node_color = log2_median_ratio, # difference between groups
                   node_color_trans = "linear",
                   node_color_interval = c(-3, 3), # symmetric interval
                   edge_color_interval = c(-3, 3), # symmetric interval
                   node_color_range = diverging_palette(), # diverging colors
                   node_size_axis_label = "OTU count",
                   node_color_axis_label = "Log 2 ratio of median counts",
                   layout = "da", initial_layout = "re",
                   key_size = 0.67,
                   seed = 2,
                   output_file = "images/Differential_Tree_1year.pdf")
##########################################################################################################
##########################################################################################################
x_meta <- obj[["data"]][["class_data"]]
x_meta <- dplyr::select(obj[["data"]][["class_data"]], -input_index)
x_meta <- na.omit(x_meta)
x_meta <- unique(x_meta)
#x_meta[x_meta$tax_name == "NA"] <- NA

x_meta[x_meta == "NA"] <- NA
x_meta <- na.omit(x_meta)
#View(x_meta)
x <- obj[["data"]][["tax_abund"]]
#y <- obj[["data"]][["diff_table"]]
abund_taxons <- x_meta %>% left_join(x, by = "taxon_id")
#abund_taxons = merge(x,x_meta,by="taxon_id")
rm(x, x_meta)
write.csv(abund_taxons, file = "images/AbundanciasTaxons_1year.csv")
#####################################################################################
# ************************************************************************************* #
abund_taxons_offs <- read_delim(file.path("data", "metadata", "legacy", "AbundanciasTaxons_Todos.csv"), delim = ",")
metadata_offs <- read_delim(file.path("data", "metadata", "legacy", "Metadata_mehg.csv"), delim = ";")

x <- abund_taxons_offs %>% filter(tax_rank != "s")
x <- dplyr::select(x, -taxon_id, -tax_rank, -tax_name,-Num)
x<-as.data.frame(x)

x<-subset(x, regex_match!="NA")

row.names(x)<-x$regex_match
x <- dplyr::select(x, -regex_match)
x <- x[rowSums(x) > 0, , drop=FALSE]
rowSums(x)[order(rowSums(x),decreasing = TRUE)]

tx <- t(x)
tx <- as.data.frame(tx)
tx$SampleID <- row.names(tx)
#taxons = merge(tx,metadata,by="SampleID")
#tx$SampleID <- as.numeric(tx$SampleID)
taxons <- metadata_offs %>% left_join(tx, by = "SampleID")
rownames(taxons) <- taxons$SampleID
rm(tx,x)
#rm(abund_taxons)
###############################################################################################################
taxons <- dplyr::select(taxons, -treatment,-sex)
taxons <- subset(taxons,SampleID != "10F")
taxons <- subset(taxons,SampleID != "22F")
taxons <- taxons %>% mutate(treatment = factor(Grupo, levels = group_level))
taxons <- taxons %>% mutate(sex = factor(Sexo, levels = sex_level))

###############################################################################################################
###############################################################################################################
# Plotagem Abundancia gênero entre geração/treatment ##############################################################
#####################################################################
taxons_1year_F1 <- subset(taxons, generation == "F1")
taxons_1year_F2 <- subset(taxons, generation == "F2")

taxons <- taxons_1year_F2
#####################################################################
# TESTE homocedasticity
teste_shapiro <-taxons %>% 
  subset(Geracao == "F1") %>%
  mutate(ratio_log = p__Proteobacteria) %>%
  #group_by(sex) %>%
  #group_by(sex, treatment) %>%
  mutate(N_Samples = n()) %>%
  nest() %>%
  mutate(Shapiro = map(data, ~shapiro.test(.x$ratio_log)))

library(broom)
teste_shapiro <- teste_shapiro %>%
  mutate(glance_shapiro = Shapiro %>% map(glance)) %>%
  unnest(glance_shapiro)

teste_shapiro
########################################
#https://www.datanovia.com/en/lessons/kruskal-wallis-test-in-r/
#####################################################################
# NONPARAMETRIC
#####################################################################
residue <- taxons %>% 
  subset(Geracao == "F1") %>%
  wilcox_test(p__Proteobacteria ~ Sexo)
#group_by(sex) %>%
#kruskal_test(Firmicutes ~ treatment)
residue

#Effect size
taxons %>%
  subset(Geracao == "F1") %>%
  #group_by(sex) %>%
  kruskal_effsize(p__Firmicutes ~ Grupo)

stat.test <- taxons %>%
  subset(Geracao == "F1") %>%
  #group_by(sex) %>% 
  dunn_test(p__Firmicutes ~ Grupo, p.adjust.method = "fdr") 
stat.test
stat.test <- stat.test[!stat.test$p.adj.signif=="ns",]
stat.test

# Box plots with p-values
stat.test <- stat.test %>% add_y_position()
#text = paste(get_test_label(residue, detailed = TRUE,type = "text")[[1]],
#             "                                                                                                   ",
#             get_test_label(residue, detailed = TRUE,type = "text")[[2]])
text <- paste("Sex comparison: ",get_test_label(residue, detailed = TRUE,type = "text"))

g <- ggboxplot(taxons %>% subset(Geracao == "F1"), x = "Grupo", y = "p__Firmicutes",palette = c("#98F5FF","#B3EE3A","#CD2626","#4D4D4D"),
               outlier.shape = NA,
               fill = "Grupo"
               #, facet.by = "Sexo"
               ) +
  stat_pvalue_manual(stat.test, label = "p.adj.signif", tip.length = 0.01) +
  labs(title = "Firmicutes",
       subtitle = text,
       caption = get_pwc_label(stat.test))+
  theme(axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        #legend.title = element_blank(),
        legend.direction='vertical',
        legend.position = "right",
        plot.subtitle = element_text(hjust = 1),
        panel.background = element_rect(fill = "gray93",colour = "grey", linetype = "solid"),
        panel.grid.major = element_line(size = 0.5, linetype = 'solid',colour = "white"), )+
  xlab("") +
  ylab("Firmicutes") +
  guides(fill=guide_legend(title="Group"))
g
tiff("images/abundancias_filo/Firmicutes_nonparametric_F2_1year.tiff",width = 6,height = 5,units = "in", res=150)
g
dev.off()
##########################################################################################################
#####################################################################
# PARAMETRIC
#####################################################################
residue <- taxons %>% 
  t_test(Firmicutes ~ sex)
residue

stat.test <- taxons %>%
  #group_by(sex) %>% 
  tukey_hsd(Firmicutes ~ treatment) 
stat.test
stat.test <- stat.test[!stat.test$p.adj.signif=="ns",]
stat.test

# Box plots with p-values
stat.test <- stat.test %>% add_y_position()
text <- paste("Sex comparison: ",get_test_label(residue, detailed = TRUE,type = "text"))

g <- ggboxplot(taxons, x = "treatment", y = "Firmicutes",palette = c("#98F5FF","#B3EE3A","#CD2626","#4D4D4D"),
               outlier.shape = NA,
               fill = "treatment", facet.by = "sex") +
  stat_pvalue_manual(stat.test, label = "p.adj.signif", tip.length = 0.01) +
  labs(title = "Firmicutes",
       subtitle = text,
       caption = get_pwc_label(stat.test))+
  theme(axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        #legend.title = element_blank(),
        legend.direction='vertical',
        legend.position = "right",
        plot.subtitle = element_text(hjust = 1),
        panel.background = element_rect(fill = "gray93",colour = "grey", linetype = "solid"),
        panel.grid.major = element_line(size = 0.5, linetype = 'solid',colour = "white"), )+
  xlab("") +
  ylab("Firmicutes") +
  guides(fill=guide_legend(title="Group"))

tiff("images/abundancias_filo/Firmicutes_parametric_F2_1year.tiff",width = 6,height = 5,units = "in", res=150)
g
dev.off()
###############################################################################################################
###############################################################################################################
# for loop para grafico de abundancio por nivel de taxon escolhido
#
##########################################################################################################
##########################################################################################################
head(taxons_1year_F1[6:13])
taxons_1year_F1_F <- subset(taxons_1year_F1, sex == "Female")
taxons_1year_F1_M <- subset(taxons_1year_F1, sex == "Male")
taxons_1year_F2_F <- subset(taxons_1year_F2, sex == "Female")
taxons_1year_F2_M <- subset(taxons_1year_F2, sex == "Male")

#teste <- lapply(taxons_1year_F1_F[6:183],function(x) aov(x ~ taxons_1year_F1_F$treatment, data = taxons_1year_F1_F))
#teste <- lapply(taxons_1year_F1_M[6:183],function(x) aov(x ~ taxons_1year_F1_M$treatment, data = taxons_1year_F1_M))
#teste <- lapply(taxons_1year_F2_F[6:183],function(x) aov(x ~ taxons_1year_F2_F$treatment, data = taxons_1year_F2_F))
#teste <- lapply(taxons_1year_F2_M[6:183],function(x) aov(x ~ taxons_1year_F2_M$treatment, data = taxons_1year_F2_M))

####
head(taxons[1:10])

taxons_dams <- subset(taxons, Geracao_Sexo == "F0_Dams")
taxons_F1Female <- subset(taxons, Geracao_Sexo == "F1_Female")
taxons_F1Male <- subset(taxons, Geracao_Sexo == "F1_Male")

#teste <- lapply(taxons_F1Male[7:153],function(x) aov(x ~ taxons_F1Male$treatment, data = taxons_F1Male))
#teste <- lapply(taxons_F1Female[7:153],function(x) aov(x ~ taxons_F1Female$treatment, data = taxons_F1Female))
#teste <- lapply(taxons_dams[7:153],function(x) aov(x ~ taxons_dams$treatment, data = taxons_dams))

##########################################################################################################
# Data Assignment

taxon_ref <- taxons_F1Male
names(taxon_ref) <- make.names(names(taxon_ref))
teste <- lapply(taxon_ref[7:260],function(x) aov(x ~ taxon_ref$Grupo, data = taxon_ref))
buf_name <- "PND21_F1_Male_"
##########################################################################################################

buf <- NA
for(i in 1:length(teste)){
  #print(summary(teste[[i]])[[1]][["Pr(>F)"]][1])
  buf <- append(buf,summary(teste[[i]])[[1]][["Pr(>F)"]][1])
  if(summary(teste[[i]])[[1]][["Pr(>F)"]][1] == 'NaN'){
    next
  }else{
    if(summary(teste[[i]])[[1]][["Pr(>F)"]][1] <= 0.05){
      print(names(teste[i]))
      print(summary(teste[[i]])[[1]][["Pr(>F)"]][1])
    }
  }
}
buf <- na.omit(buf)
#buf[order(buf)]
######################################
# PARAMETRIC
for(i in 1:length(teste)){
  if(summary(teste[[i]])[[1]][["Pr(>F)"]][1] == 'NaN'){
    next
  }else{
    if(summary(teste[[i]])[[1]][["Pr(>F)"]][1] <= 0.05){
      print(names(teste[i]))
      print(summary(teste[[i]])[[1]][["Pr(>F)"]][1])
      #tukey_result<-HSD.test(teste[[i]],"taxon_ref$Grupo", group= TRUE)                  #AQUI
      #group_data <- tukey_result$groups[order(rownames(tukey_result$groups)),]
      #print(group_data)
      
      #group_data$Grupo <- rownames(group_data)
      #group_data <- group_data %>% mutate(Grupo = factor(Grupo, levels = group_level))
      #colnames(group_data) <- c("taxon","groups","Grupo")
      #####
      tukey_result <- tukey_hsd(teste[[i]])
      abcd <- paste(names(teste[i]),"~ Grupo")
      #tukey_result <- do.call("dunn_test",list(formula=formula(abcd),data=quote(taxon_ref),p.adjust.method="fdr"))
      tukey_result <- do.call("add_y_position",list(test = quote(tukey_result),formula=formula(abcd),data=quote(taxon_ref)))
      tukey_result <- tukey_result[!tukey_result$p.adj.signif=="ns",]
      
      buf <- paste("Anova, p= ",formatC(summary(teste[[i]])[[1]][["Pr(>F)"]][1], format = "e", digits = 2))
      #####
      g<-ggboxplot(taxon_ref,
                     x = "Grupo", y = names(teste[i]),
                     fill = "Grupo", 
                     outlier.shape = NA,
                     palette = c("#98F5FF","#B3EE3A","#CD2626","#4D4D4D")
                     #ylim = c(0, 45)
        )+
          theme(axis.text.x=element_blank(),
                axis.ticks.x=element_blank(),
                legend.title = element_blank(),
                legend.direction='vertical',
                legend.position = "right",
                panel.background = element_rect(fill = "gray93",colour = "grey", linetype = "solid"),
                panel.grid.major = element_line(size = 0.5, linetype = 'solid',colour = "white"), )+
          ggtitle(names(teste[i]))+
          xlab("Groups") +
          ylab("Relative Abundance") + 
          stat_pvalue_manual(tukey_result, label = "p.adj.signif", tip.length = 0.01) +
          labs(title = names(teste[i]),
               subtitle = buf,
               caption = get_pwc_label(tukey_result))
        ggsave(filename= paste("images/Abundancias/",buf_name,names(teste[i]),".tiff"),            #AQUI
               plot = g, device = "tiff",
               width = 6,height = 5,units = "in", dpi=100)
    }
  }
}
i
g
######################################
# NONPARAMETRIC
teste <- lapply(taxon_ref[7:260],function(x) kruskal_test(x ~ taxon_ref$Grupo, data = taxon_ref))
teste$p__Firmicutes$p

buf <- NA
for(i in 1:length(teste)){
  #print(summary(teste[[i]])[[1]][["Pr(>F)"]][1])
  buf <- append(buf,teste[[i]]$p)
  if(teste[[i]]$p == 'NaN'){
    next
  }else{
    if(teste[[i]]$p <= 0.05){
      print(names(teste[i]))
      print(teste[[i]]$p)
    }
  }
}
buf <- na.omit(buf)
#buf[order(buf)]


for(i in 1:length(teste)){
  if(teste[[i]]$p == 'NaN'){
    next
  }else{
    if(teste[[i]]$p <= 0.05){
      print(names(teste[i]))
      print(teste[[i]]$p)
      #####
      abcd <- paste(names(teste[i]),"~ Grupo")
      dunn_result <- do.call("dunn_test",list(formula=formula(abcd),data=quote(taxon_ref),p.adjust.method="fdr"))
      dunn_result <- do.call("add_y_position",list(test = quote(dunn_result),formula=formula(abcd),data=quote(taxon_ref)))
      dunn_result <- dunn_result[!dunn_result$p.adj.signif=="ns",]
      #####
      g<-ggboxplot(taxon_ref,
                   x = "Grupo", y = names(teste[i]),
                   fill = "Grupo", 
                   outlier.shape = NA,
                   palette = c("#98F5FF","#B3EE3A","#CD2626","#4D4D4D")
                   #ylim = c(0, 45)
      )+
        theme(axis.text.x=element_blank(),
              axis.ticks.x=element_blank(),
              legend.title = element_blank(),
              legend.direction='vertical',
              legend.position = "right",
              panel.background = element_rect(fill = "gray93",colour = "grey", linetype = "solid"),
              panel.grid.major = element_line(size = 0.5, linetype = 'solid',colour = "white"), )+
        ggtitle(names(teste[i]))+
        xlab("Groups") +
        ylab("Relative Abundance") + 
        stat_pvalue_manual(dunn_result, label = "p.adj.signif", tip.length = 0.01) +
        labs(title = names(teste[i]),
             subtitle = get_test_label(teste[[i]], detailed = TRUE,type = "text"),
             caption = get_pwc_label(dunn_result))
      ggsave(filename= paste("images/Abundancias/",buf_name,names(teste[i]),".tiff"),            #AQUI
             plot = g, device = "tiff",
             width = 6,height = 5,units = "in", dpi=100)
    }
  }
}

##########################################################################################################
#####################################################################
# Plot abundancies by Phylum without sex facets
#####################################################################
#####################################################################
MajorPhyla <- c("p__Firmicutes","p__Bacteroidota","p__Proteobacteria","p__Actinobacteriota")
append(c("SampleID","Grupo","Sexo","Geracao"),MajorPhyla)

#taxons_majorphyla <- taxons %>% select(append(c("SampleID","Grupo","Sexo","Geracao"),MajorPhyla))
taxons_majorphyla <- taxons_1year_F2 %>% select(append(c("sampleID","treatment","sex","generation"),MajorPhyla))
taxons_majorphyla <- taxons_majorphyla %>% subset(Geracao == "F2")
#abcd <- taxons_majorphyla$SampleID
abcd <- taxons_majorphyla$sampleID
#taxons_majorphyla <- taxons_majorphyla %>% select(-Sexo, -Geracao, -SampleID)
taxons_majorphyla <- taxons_majorphyla %>% select(-sex, -generation, -sampleID)
rownames(taxons_majorphyla) <- abcd
taxons_majorphyla

#taxons_majorphyla <- melt(taxons_majorphyla, id.vars = "Grupo")
taxons_majorphyla <- melt(taxons_majorphyla, id.vars = "treatment")
head(taxons_majorphyla)
#taxons_majorphyla <- taxons_majorphyla %>% mutate(Grupo = factor(Grupo, levels = group_level))
taxons_majorphyla <- taxons_majorphyla %>% mutate(treatment = factor(treatment, levels = group_level))
taxons_majorphyla$variable <- str_replace_all(taxons_majorphyla$variable,"p__","")
head(taxons_majorphyla)
taxons_majorphyla <- taxons_majorphyla %>% mutate(variable = factor(variable, levels = phyla_level))

# TESTE homocedasticity
teste_shapiro <-taxons_1year_F2 %>% 
  #subset(Geracao == "F0") %>%
  #mutate(ratio_log = p__Firmicutes) %>%
  #mutate(ratio_log = p__Bacteroidota) %>%
  mutate(ratio_log = p__Proteobacteria) %>%
  #mutate(ratio_log = p__Actinobacteriota) %>%
  #group_by(sex) %>%
  #group_by(sex, treatment) %>%
  mutate(N_Samples = n()) %>%
  nest() %>%
  mutate(Shapiro = map(data, ~shapiro.test(.x$ratio_log)))

library(broom)
teste_shapiro <- teste_shapiro %>%
  mutate(glance_shapiro = Shapiro %>% map(glance)) %>%
  unnest(glance_shapiro)

#teste_shapiro
#shapiro <-NA
shapiro <- append(shapiro,teste_shapiro$p.value)
shapiro

shapiro <- na.omit(shapiro)

########################################
#https://www.datanovia.com/en/lessons/kruskal-wallis-test-in-r/
#####################################################################
# NONPARAMETRIC
#####################################################################
residue <- taxons_1year_F2 %>% 
  #subset(Geracao == "F1") %>%
  #wilcox_test(p__Firmicutes ~ sex)
  #wilcox_test(p__Bacteroidota ~ sex)
  #wilcox_test(p__Proteobacteria ~ sex)
  wilcox_test(p__Actinobacteriota ~ sex)
residue

#Effect size
taxons %>% subset(Geracao == "F1") %>%
  #group_by(sex) %>%
  kruskal_effsize(p__Firmicutes ~ Sexo)

stat.test <- taxons_majorphyla %>%
  #subset(variable == "Actinobacteriota") %>%
  group_by(variable) %>%
  dunn_test(value ~ treatment, p.adjust.method = "fdr") 
  #tukey_hsd(value ~ treatment, p.adjust.method = "fdr") 
View(stat.test)
stat.test <- stat.test[!stat.test$p.adj.signif=="ns",]
stat.test

# Box plots with p-values
stat.test <- stat.test %>% add_y_position()

g <- ggplot(taxons_majorphyla, aes(x = variable, y = value),
               outlier.shape = NA) +
  geom_boxplot(aes(fill=treatment),outliers = TRUE) +
  scale_fill_manual(values = c("#98F5FF","#B3EE3A","#CD2626","#4D4D4D")) +
  #stat_pvalue_manual(stat.test, label = "p.adj.signif", tip.length = 0.01) +
  labs(title = "F2 - 1 year",
       subtitle = "Sex comparisons for all phyla, Wilcoxon test p > 0.05",
       caption = "Shapiro-Wilk test: all phyla, p < 0.05"
       )+
  xlab("Phylum") +
  ylab("Relative Abundance")+
  theme(axis.text.x=element_blank(),axis.ticks.x=element_blank()) +
  facet_wrap(~ variable, scales = "free") +
  guides(fill=guide_legend(title="Group"))

g


tiff("images/abundancias_filo/MajorPhyla_nonparametric_F2_1year.tiff",width = 6,height = 5,units = "in", res=300)
g
dev.off()
