#GERANDO GRÁFICO DE BARRAS MICROBIOTA

#Vídeo tutorial: https://www.youtube.com/watch?v=siIoupAnILk
#Repositório: https://github.com/surh/scip_barplot
#No repositório, os dados são encontrados na pasta data/rhizo
#No repositório, os script utilizado no vídeo tem o seguinte nome: video_code.Rmd
#no mesmo repositório, o arquivo "extended_example.Rmd" contém o mesmo script tutorial mas com explicação completa

  
library(tidyverse)
library("dplyr")
getwd()
setwd("D:/a.UFRGS/BOLSA/bioq.meHg/MeHg_Epigenetica_Microbiota/Microbiota")

Tab <- read_tsv("seqtab.txt",
                col_types = cols(OTU = col_character(),
                                 .default = col_number()))
Tab
rowSums(Tab[2:55])
colSums(Tab[2:55])[order(colSums(Tab[2:55]))]

Tab <- select(Tab, -c("10F","22F"))
colSums(Tab[2:53])[order(colSums(Tab[2:53]))]
1697

Tab <- filter(Tab, rowSums(Tab[2:53]) > 2)
Tab
#rowSums(Tab[2:109])
#########################################################################################################
dat <- Tab %>%
  pivot_longer(-OTU, names_to = "SampleID", values_to = "count")
dat
#########################################################################################################
Tax <- read_tsv("Data_Metadata/tax.txt",
                col_types = cols(.default = "character"))
Tax
#########################################################################################################
dat <- dat %>%
  left_join(Tax, by = "OTU")
dat

dat <- dat %>% select(-c("Kingdom"))
dat

#dat <- filter(dat, !is.na(dat$Phylum))
#dat
#########################################################################################################
#Meta <- read_tsv("Data_Metadata/sample-metadata.txt",
#                 col_types = cols(.default = col_character()))
Meta <- read_delim("Data_Metadata/Metadata_mehg.csv",delim=';',
                 col_types = cols(.default = col_character()))
Meta

#########################################################################################################
dat <- dat %>%
  left_join(Meta, by = "SampleID")
dat
#########################################################################################################
dat %>%
  ggplot(aes(x = SampleID, y = count)) +
  facet_grid(scales = "free_x", space = "free_x") +
  geom_bar(aes(fill = Phylum), stat = "identity", position = "fill", width = 1)

dat %>%
  ggplot(aes(x = SampleID, y = count)) +
  facet_grid(~ Geracao_Sexo + Grupo, scales = "free_x", space = "free_x") +
  geom_bar(aes(fill = Phylum), stat = "identity", position = "fill", width = 1)

dat %>%
  ggplot(aes(x = SampleID, y = count)) +
  facet_grid(~ Geracao_Sexo + Grupo, scales = "free_x", space = "free_x") +
  geom_bar(aes(fill = Phylum), stat = "identity", position = "fill", width = 1) +
  scale_fill_brewer(palette = "Paired")
#########################################################################################################
group_level <- c("CN","VitA","MeHg","VitA+MeHg")
geracao_sexo_level <- c("F0_Dams","F1_Female","F1_Male")

#######################################
#Determinação dos filos mais abundantes
mean_freqs <- dat %>%
  group_by(SampleID, Phylum) %>%
  summarise(phyl_count = sum(count),
            .groups = 'drop') %>%
  group_by(SampleID) %>%
  summarise(phyl_freq = phyl_count / sum(phyl_count),
            Phylum = Phylum,
            .groups = 'drop') %>%
  group_by(Phylum) %>%
  summarise(mean = mean(phyl_freq),
            .groups = 'drop') %>%
  arrange(desc(mean))
print(mean_freqs)
phyla_level <- mean_freqs$Phylum
phyla_level <- phyla_level[!is.na(phyla_level)]
phyla_level <- append(phyla_level, NA)
phyla_level
#############################################################################
dat <- dat %>%
  mutate(Phylum = factor(Phylum, levels = phyla_level))

dat <- dat %>%
  mutate(Grupo = factor(Grupo, levels = group_level))

dat <- dat %>%
  mutate(Geracao_Sexo = factor(Geracao_Sexo, levels = geracao_sexo_level))


#dat_chow <- subset(dat, dat$LINHAGEM == "CHOW")
#dat_chow <- select(dat_chow, -c("Species"))
#sum(dat_chow$count)
#dat_chow <- filter(dat_chow, !is.na(dat_chow$Phylum))
#sum(dat_chow$count)
#i <- 25042-79

#dat_chow <- subset(dat_chow, count > 0)
#########################################################################################################
plot <- dat %>%
          ggplot(aes(x = SampleID, y = count)) +
          facet_grid(~ Geracao_Sexo + Grupo, scales = "free_x", space = "free_x") +
          geom_bar(aes(fill = Phylum), stat = "identity", position = "fill", width = 1) +
          scale_fill_brewer(palette = "Paired") +
          scale_y_continuous(name = "Relative abundance",
                             labels = scales::percent) +
          theme(axis.text.x = element_text(angle = 90, size = 6),
                axis.text.y = element_text(color = "black"),
                strip.text = element_text(face = "bold"),
                strip.background = element_blank())
plot

tiff("Figures/Abundancias_Filo_Mestrado.tiff", width = 12,height = 7 , units = "in", res=300)
plot
dev.off()

#ggsave("Figures/Abundancias_Filo_FernandaSantos.png",plot,width = 1920, height = 991, device = "png", dpi = 300, units = "px", limitsize = FALSE)
#ggsave("Figures/Abundancias_Filo_FernandaSantos2.pdf",plot, device = "pdf",width = 12, height = 7, dpi = 300, units = "in")
