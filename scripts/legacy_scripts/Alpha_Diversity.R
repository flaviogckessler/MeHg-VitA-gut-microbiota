library(readr)
library(dplyr)
library(tidyverse)
library(ggplot2)
library(rstatix)
library(ggpubr)
###########################################################################################################
group_level <- c("CN", "VitA", "MeHg", "VitA+MeHg")
sex_level <- c("Female","Male")
phyla_level <- c('Firmicutes','Bacteroidota','Proteobacteria','Actinobacteriota','Campylobacterota','Desulfobacterota','Verrucomicrobiota','Patescibacteria')
###########################################################################################################
# F0 and F1 PND21
metadata <- read.csv("C:/Users/flavi/OneDrive/Documentos/R_Microbiota/Alberto_Hg/Data_metadata/Metadata_mehg.csv")
rownames(metadata) <- metadata$SampleID

alpha1<- read_csv("C:/Users/flavi/OneDrive/Documentos/R_Microbiota/Alberto_Hg/Data_metadata/Alpha_filtered_Metadata_F0F1.csv")
rownames(alpha1) <- alpha1$SampleID

# 1 year
metadata <- read.csv("C:/Users/flavi/OneDrive/Documentos/R_Microbiota/Alberto_Hg/Data_metadata/metadata.csv")
rownames(metadata) <- metadata$sampleID

alpha<- read_csv("C:/Users/flavi/OneDrive/Documentos/R_Microbiota/Alberto_Hg/Data_metadata/alpha_diversity.csv")
rownames(alpha) <- alpha$sampleID
###############################################
#alpha_f0 <- subset(alpha1,Geracao == "F0")
alpha_f1 <- subset(alpha,generation == "F1")
alpha_f2 <- subset(alpha,generation == "F2")

alpha <- alpha_f1
rownames(alpha) <- alpha$sampleID

#alpha <- alpha %>%  mutate(Grupo = factor(Grupo, levels = group_level))
#alpha <- alpha %>%  mutate(Sexo = factor(Sexo, levels = sex_level))
alpha <- alpha %>%  mutate(treatment = factor(treatment, levels = group_level))
alpha <- alpha %>%  mutate(sex = factor(sex, levels = sex_level))
colnames(alpha)

alpha <- as.data.frame(alpha)
#rownames(alpha) <- alpha$SampleID
rownames(alpha) <- alpha$sampleID
##################################################################################
##################################################################################
teste_shapiro <-alpha %>% 
  mutate(ratio_log = InvSimpson) %>%
  #group_by(Sexo) %>%
  #group_by(Sexo, Grupo) %>%
  mutate(N_Samples = n()) %>%
  nest() %>%
  mutate(Shapiro = map(data, ~shapiro.test(.x$ratio_log)))

library(broom)
teste_shapiro <- teste_shapiro %>%
  mutate(glance_shapiro = Shapiro %>% map(glance)) %>%
  unnest(glance_shapiro)

teste_shapiro

#####################################################################
#https://www.datanovia.com/en/lessons/kruskal-wallis-test-in-r/
#####################################################################
# NONPARAMETRIC
#####################################################################
residue <- alpha %>% 
  wilcox_test(InvSimpson ~ sex)
  #group_by(Sexo) %>%
  #kruskal_test(InvSimpson ~ Grupo)
residue

#Effect size
alpha %>% 
  group_by(sex) %>%
  kruskal_effsize(InvSimpson ~ treatment)

stat.test <- alpha %>%
  group_by(sex) %>%
  dunn_test(InvSimpson ~ treatment, p.adjust.method = "fdr") 
stat.test
stat.test <- stat.test[!stat.test$p.adj.signif=="ns",]
stat.test

# Box plots with p-values
stat.test <- stat.test %>% add_y_position()
#text = paste(get_test_label(residue, detailed = TRUE,type = "text")[[1]],
#             "                                                                                                   ",
#             get_test_label(residue, detailed = TRUE,type = "text")[[2]])
text <- paste("Sex comparison: ",get_test_label(residue, detailed = TRUE,type = "text"))

g <- ggboxplot(alpha, x = "treatment", y = "InvSimpson",palette = c("#98F5FF","#B3EE3A","#CD2626","#4D4D4D"),
               outlier.shape = NA,
               fill = "treatment", facet.by = "sex") +
  stat_pvalue_manual(stat.test, label = "p.adj.signif", tip.length = 0.01) +
  labs(title = "InvSimpson",
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
  ylab("InvSimpson") +
  guides(fill=guide_legend(title="Group"))
g

tiff("images/InvSimpson_nonparametric_F1_DPN21.tiff",width = 8,height = 5,units = "in", res=150)
g
dev.off()
##########################################################################################################
#####################################################################
# PARAMETRIC
#####################################################################
residue <- alpha %>% 
  t_test(InvSimpson ~ sex)
residue

stat.test <- alpha %>%
  group_by(sex) %>% 
  tukey_hsd(InvSimpson ~ treatment) 
stat.test
stat.test <- stat.test[!stat.test$p.adj.signif=="ns",]
stat.test

# Box plots with p-values
stat.test <- stat.test %>% add_y_position()
text <- paste("Sex comparison: ",get_test_label(residue, detailed = TRUE,type = "text"))

g <- ggboxplot(alpha, x = "treatment", y = "InvSimpson",palette = c("#98F5FF","#B3EE3A","#CD2626","#4D4D4D"),
               outlier.shape = NA,
               fill = "treatment", facet.by = "sex") +
  stat_pvalue_manual(stat.test, label = "p.adj.signif", tip.length = 0.01) +
  labs(title = "InvSimpson",
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
  ylab("InvSimpson") +
  guides(fill=guide_legend(title="Group"))
g

tiff("images/InvSimpson_parametric_F1_1year.tiff",width = 8,height = 5,units = "in", res=150)
g
dev.off()