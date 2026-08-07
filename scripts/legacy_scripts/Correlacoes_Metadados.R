library(readr)
library(dplyr)
#################
group_level <- c("CN","VitA","MeHg","VitA+MeHg")
geracao_sexo_level <- c("F0_Dams","F1_Female","F1_Male")
#########################################################################################################
metadata<-read_delim("Data_Metadata/Metadata_mehg.csv", delim = ";")
rownames(metadata) <- metadata$SampleID
row.names(metadata)
metadata <- metadata[!(row.names(metadata) %in% c("10F","22F")),]
rownames(metadata) <- metadata$SampleID

metadata <- metadata %>%  mutate(Grupo = factor(Grupo, levels = group_level))
metadata <- metadata %>%  mutate(Geracao_Sexo = factor(Geracao_Sexo, levels = geracao_sexo_level))
rownames(metadata) <- metadata$SampleID
#########################################################################################################
meta_hg_hplc <- read_delim("Data_Metadata/Metadata_THg_HPLC.csv", delim = ",",col_types = cols(SampleID = col_character()))
rownames(meta_hg_hplc) <- meta_hg_hplc$SampleID
row.names(meta_hg_hplc)
meta_hg_hplc <- meta_hg_hplc[!(row.names(meta_hg_hplc) %in% c("10F","22F")),]
rownames(meta_hg_hplc) <- meta_hg_hplc$SampleID

meta_hg_hplc <- meta_hg_hplc %>%  mutate(Grupo = factor(Grupo, levels = group_level))
meta_hg_hplc <- meta_hg_hplc %>%  mutate(Geracao_Sexo = factor(Geracao_Sexo, levels = geracao_sexo_level))
rownames(meta_hg_hplc) <- meta_hg_hplc$SampleID
#########################################################################################################

abund_taxons <- read_delim("Figures/AbundanciasTaxons_Todos_EditTexto.csv", delim = ",")

x <- abund_taxons %>% filter(tax_rank != "s")
x <- x %>% filter(tax_rank != "k")
#x <- abund_taxons
x <- dplyr::select(x, -taxon_id, -tax_rank, -tax_name)
x<-as.data.frame(x)

x<-subset(x, regex_match!="NA")

row.names(x)<-x$regex_match
x <- dplyr::select(x, -regex_match,-Num)
x <- x[rowSums(x) > 0, , drop=FALSE]

tx <- t(x)
tx <- as.data.frame(tx)
tx$SampleID <- row.names(tx)
#taxons = merge(tx,metadata,by="SampleID")
#tx$SampleID <- as.numeric(tx$SampleID)
taxons_hg_hplc <- meta_hg_hplc %>% left_join(tx, by = "SampleID")
row.names(taxons_hg_hplc) <- taxons_hg_hplc$SampleID
rm(tx,x)

#########################################################################################################
###### CORRELAÇÃO #########################################################################
#########################################################################################################
SampleID <- taxons_hg_hplc$SampleID

# Amostras devem estar nas linhas e dados nas colunas, caso estejam ao contrário será necessário transpor a matriz
#t_otumat <-t(otumat)

#Análise de correlação entre aspectos dos individuos com os taxons, com p valor ajustado
#Análise faz correlações de uma matriz, metadados, com outra matriz, abundância dos taxons. Ambas matrizes compartilham mesma identificações de indivíduos como nome das linhas

library(ggbiplot)
library(ggplot2)
#install.packages("BiocManager")
#library(BiocManager)
BiocManager::install("microbiome")
library(microbiome)
library("knitr")

x <- as.matrix(select(taxons_hg_hplc, -colnames(taxons_hg_hplc)[1:21]))
y <- as.matrix(select(taxons_hg_hplc, -colnames(taxons_hg_hplc)[1:6], -colnames(taxons_hg_hplc)[22:length(colnames(taxons_hg_hplc))]))
rownames(x) <- taxons_hg_hplc$SampleID
rownames(y) <- taxons_hg_hplc$SampleID


correlation.table <- associate(y,x, method = "spearman", mode = "table", 
                               p.adj.threshold = 0.05, n.signif = 1,filter.self.correlations = TRUE)
correlation.table <- filter(correlation.table, p.adj <= 0.05)

#write.csv(correlation.table, file= "C:/Users/flavi/OneDrive/Documentos/R_Microbiota/Pedrotti/correlacoes_metadados-metadados_Pedrotti.csv")
#correlation.table <- filter(correlation.table, Correlation > 0.3)
kable(head(correlation.table))

p <- heat(correlation.table, "X1", "X2", 
          fill = "Correlation", 
          star = "p.adj", 
          p.adj.threshold = 0.13)
p + theme(text=element_text(size=10), 
          axis.text.x = element_text(angle = 90, hjust = 1), 
          legend.key.size = unit(1.3, "cm"))
p

#https://www.youtube.com/watch?v=eXczd0ewVgE
install.packages("ggstatsplot")
library(ggstatsplot)
g<- ggscatterstats(data = taxons_hg_hplc,
                   x = Hg,y = g__Dorea,
                   type = "nonparametric") # or "np" or "r"

tiff("Figures/Dorea_hg_correlation_nonparametric.tiff", units="in", width=5.5, height=4, res=150)
g
dev.off()

#install.packages("ggpubr")
library(ggpubr)

###### perason correlation for parametric
###### spearman correlation for nonparametric

g <- ggscatter(taxons_hg_hplc, x = "Hg", y = "g__Dorea",
          add = "reg.line",                         # Add regression line
          conf.int = TRUE,                          # Add confidence interval
          color = "Sexo", 
          palette = "npg",           # Color by groups "cyl"
          #shape = "Sexo",                          # Change point shape by groups "cyl"
          #facet.by = "Sexo",
          cor.method = "spearman",xlab = "THg, µg/g", ylab = "Dorea sp.") +
  theme(axis.title.y = element_text(face = 'italic'))+
  stat_cor(aes(color = Sexo),label.x = c(0.15,1.6),label.y = c(0.0042,0.0075),
           show.legend = FALSE, method = "spearman")           # Add correlation coefficient

ggpar(g,legend = "right", legend.title = "Sex")

tiff("Figures/Dorea_SEX_hg_correlation_nonparametric.tiff", units="in", width=5.5, height=4, res=150)
ggpar(g,legend = "right", legend.title = "Sex")
dev.off()

################################################################################################################
