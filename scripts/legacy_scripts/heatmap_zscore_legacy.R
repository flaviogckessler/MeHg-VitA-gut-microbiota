#install.packages("gplots")  #Instalar gplots#
#install.packages("rlist")
library(gplots) #rodar a cada inicializa??o do Rstudio#
library(readr)
library(dplyr)
library(rstatix)
library(rlist)

###########################################################################################################
group_level <- c("CN", "VitA", "MeHg", "VitA+MeHg")
sex_level <- c("Female","Male")
phyla_level <- c('Firmicutes','Bacteroidota','Proteobacteria','Actinobacteriota','Campylobacterota','Desulfobacterota','Verrucomicrobiota','Patescibacteria')
###########################################################################################################
#importar dataset, pode ser atrav?s do ?cone ou usando read.csv#
taxons_offs <- read_delim("images/AbundanciasTaxons_Todos.csv", delim = ",")
taxons_1year <- read_delim("images/AbundanciasTaxons_1year.csv", delim = ",")

# Metadata
metadata_offs <- read_delim(file.path("data", "metadata", "legacy", "Metadata_mehg.csv"), delim = ";")
metadata_1year <- read_delim(file.path("data", "metadata", "legacy", "metadata.csv"), delim = ",")
rownames(metadata_offs) <- metadata_offs$SampleID
rownames(metadata_1year) <- metadata_1year$sampleID
###########################################################################################################
x <- taxons_offs %>% filter(tax_rank != "s")
x <- dplyr::select(x, -taxon_id, -tax_rank, -tax_name,-Num)
x<-as.data.frame(x)

x<-subset(x, regex_match!="NA")

row.names(x)<-x$regex_match
x <- dplyr::select(x, -regex_match)
x <- x[rowSums(x) > 0, , drop=FALSE]
#rowSums(x)[order(rowSums(x),decreasing = TRUE)]

tx <- t(x)
tx <- as.data.frame(tx)
tx$SampleID <- row.names(tx)
taxons_offs <- metadata_offs %>% left_join(tx, by = "SampleID")

taxons_offs <- dplyr::select(taxons_offs, -Geracao_Sexo,-Grupo_Geracao_Sexo,-k__Bacteria)
taxons_offs <- subset(taxons_offs,SampleID != "10F")
taxons_offs <- subset(taxons_offs,SampleID != "22F")
taxons_offs <- taxons_offs %>% mutate(Grupo = factor(Grupo, levels = group_level))
taxons_offs <- taxons_offs %>% mutate(Sexo = factor(Sexo, levels = sex_level))
rownames(taxons_offs) <- taxons_offs$SampleID
###
x <- taxons_1year %>% filter(tax_rank != "s")
x <- dplyr::select(x, -taxon_id, -tax_rank, -tax_name,-Num)
x<-as.data.frame(x)

x<-subset(x, regex_match!="NA")

row.names(x)<-x$regex_match
x <- dplyr::select(x, -regex_match)
x <- x[rowSums(x) > 0, , drop=FALSE]
#rowSums(x)[order(rowSums(x),decreasing = TRUE)]

tx <- t(x)
tx <- as.data.frame(tx)
tx$sampleID <- row.names(tx)
taxons_1year <- metadata_1year %>% left_join(tx, by = "sampleID")

taxons_1year <- dplyr::select(taxons_1year,-obs,-k__Bacteria)
taxons_1year <- taxons_1year %>% mutate(treatment = factor(treatment, levels = group_level))
taxons_1year <- taxons_1year %>% mutate(sex = factor(sex, levels = sex_level))
rownames(taxons_1year) <- taxons_1year$sampleID
rm(x,tx)
###########################################################################################################
# Subset generations
taxons_dams <- subset(taxons_offs,Geracao == "F0") %>% select(-Sexo,-Geracao)
taxons_offs_F1 <- subset(taxons_offs,Geracao == "F1")

taxons_1year_F1 <- subset(taxons_1year,generation == "F1")
taxons_1year_F2 <- subset(taxons_1year,generation == "F2")
###########################################################################################################
#taxons_ref <- taxons_dams
#taxons_ref <- taxons_offs_F1 %>% subset(Sexo == "Female")
#taxons_ref <- taxons_offs_F1 %>% subset(Sexo == "Male")
#taxons_ref <- taxons_1year_F1 %>% subset(sex == "Female")
#taxons_ref <- taxons_1year_F1 %>% subset(sex == "Male")
#taxons_ref <- taxons_1year_F2 %>% subset(sex == "Female")
taxons_ref <- taxons_1year_F2 %>% subset(sex == "Male")
#names(taxons_ref) <- make.names(names(taxons_ref))
#teste <- lapply(taxons_ref[5:length(taxons_ref)],function(x) kruskal_test(x ~ taxons_ref$Grupo, data = taxons_ref))
teste <- lapply(taxons_ref[5:length(taxons_ref)],function(x) kruskal_test(x ~ taxons_ref$treatment, data = taxons_ref))

buf <- NA
for(i in 1:length(teste)){
  if(teste[[i]]$p == 'NaN'){
    next
  }else{
    if(teste[[i]]$p <= 0.05){
      buf <- append(buf,names(teste[i]))
      print(names(teste[i]))
      print(teste[[i]]$p)
    }
  }
}
buf <- na.omit(buf)
buf

#Signif_dams <- buf
#Signif_offs_F1Female <- buf
#Signif_offs_F1Male <- buf
#Signif_1year_F1Female <- buf
#Signif_1year_F1Male <- buf
#Signif_1year_F2Female <- buf
Signif_1year_F2Male <- buf

rm(taxons_ref)
###########################################################################################################
# Order significant taxons 
Signif_dams <- Signif_dams[order(Signif_dams)]
Signif_dams <- list.append(Signif_dams[grep('p__',Signif_dams)],
                                    Signif_dams[grep('c__',Signif_dams)],
                                    Signif_dams[grep('o__',Signif_dams)],
                                    Signif_dams[grep('f__',Signif_dams)],
                                    Signif_dams[grep('g__',Signif_dams)])

Signif_offs_F1Female <- Signif_offs_F1Female[order(Signif_offs_F1Female)]
Signif_offs_F1Female <- list.append(Signif_offs_F1Female[grep('p__',Signif_offs_F1Female)],
                                     Signif_offs_F1Female[grep('c__',Signif_offs_F1Female)],
                                     Signif_offs_F1Female[grep('o__',Signif_offs_F1Female)],
                                     Signif_offs_F1Female[grep('f__',Signif_offs_F1Female)],
                                     Signif_offs_F1Female[grep('g__',Signif_offs_F1Female)])

Signif_offs_F1Male <- Signif_offs_F1Male[order(Signif_offs_F1Male)]
Signif_offs_F1Male <- list.append(Signif_offs_F1Male[grep('p__',Signif_offs_F1Male)],
                                   Signif_offs_F1Male[grep('c__',Signif_offs_F1Male)],
                                   Signif_offs_F1Male[grep('o__',Signif_offs_F1Male)],
                                   Signif_offs_F1Male[grep('f__',Signif_offs_F1Male)],
                                   Signif_offs_F1Male[grep('g__',Signif_offs_F1Male)])


Signif_1year_F1Female <- Signif_1year_F1Female[order(Signif_1year_F1Female)]
Signif_1year_F1Female <- list.append(Signif_1year_F1Female[grep('p__',Signif_1year_F1Female)],
                                     Signif_1year_F1Female[grep('c__',Signif_1year_F1Female)],
                                     Signif_1year_F1Female[grep('o__',Signif_1year_F1Female)],
                                     Signif_1year_F1Female[grep('f__',Signif_1year_F1Female)],
                                     Signif_1year_F1Female[grep('g__',Signif_1year_F1Female)])

Signif_1year_F1Male <- Signif_1year_F1Male[order(Signif_1year_F1Male)]
Signif_1year_F1Male <- list.append(Signif_1year_F1Male[grep('p__',Signif_1year_F1Male)],
                                     Signif_1year_F1Male[grep('c__',Signif_1year_F1Male)],
                                     Signif_1year_F1Male[grep('o__',Signif_1year_F1Male)],
                                     Signif_1year_F1Male[grep('f__',Signif_1year_F1Male)],
                                     Signif_1year_F1Male[grep('g__',Signif_1year_F1Male)])

Signif_1year_F2Female <- Signif_1year_F2Female[order(Signif_1year_F2Female)]
Signif_1year_F2Female <- list.append(Signif_1year_F2Female[grep('p__',Signif_1year_F2Female)],
                                   Signif_1year_F2Female[grep('c__',Signif_1year_F2Female)],
                                   Signif_1year_F2Female[grep('o__',Signif_1year_F2Female)],
                                   Signif_1year_F2Female[grep('f__',Signif_1year_F2Female)],
                                   Signif_1year_F2Female[grep('g__',Signif_1year_F2Female)])

Signif_1year_F2Male <- Signif_1year_F2Male[order(Signif_1year_F2Male)]
Signif_1year_F2Male <- list.append(Signif_1year_F2Male[grep('p__',Signif_1year_F2Male)],
                                     Signif_1year_F2Male[grep('c__',Signif_1year_F2Male)],
                                     Signif_1year_F2Male[grep('o__',Signif_1year_F2Male)],
                                     Signif_1year_F2Male[grep('f__',Signif_1year_F2Male)],
                                     Signif_1year_F2Male[grep('g__',Signif_1year_F2Male)])
###########################################################################################################

#Signif <- data.frame(tax = Signif_offs_F1Female,CN = NA,VitA=NA, MeHg=NA, "VitA+MeHg"= NA,med=NA,sd=NA)
#Signif <- data.frame(tax = Signif_offs_F1Male,CN = NA,VitA=NA, MeHg=NA, "VitA+MeHg"= NA,med=NA,sd=NA)
#Signif <- data.frame(tax = Signif_dams,CN = NA,VitA=NA, MeHg=NA, "VitA+MeHg"= NA,med=NA,sd=NA)
#Signif <- data.frame(tax = Signif_1year_F1Female,CN = NA,VitA=NA, MeHg=NA, "VitA+MeHg"= NA,med=NA,sd=NA)
#Signif <- data.frame(tax = Signif_1year_F1Male,CN = NA,VitA=NA, MeHg=NA, "VitA+MeHg"= NA,med=NA,sd=NA)
#Signif <- data.frame(tax = Signif_1year_F2Female,CN = NA,VitA=NA, MeHg=NA, "VitA+MeHg"= NA,med=NA,sd=NA)
Signif <- data.frame(tax = Signif_1year_F2Male,CN = NA,VitA=NA, MeHg=NA, "VitA+MeHg"= NA,med=NA,sd=NA)
colnames(Signif)[5] <-"VitA+MeHg"
rownames(Signif) <- Signif$tax

for (i in 1:length(Signif$tax)) {
  buf <- taxons_offs_F1 %>% subset(Sexo == "Male"&Grupo =="CN") %>% select(Signif$tax[i])
  Signif[i,2] <-buf[[Signif$tax[i]]] %>% mean()

  buf <- taxons_offs_F1 %>% subset(Sexo == "Male"&Grupo =="VitA") %>% select(Signif$tax[i])
  Signif[i,3] <-buf[[Signif$tax[i]]] %>% mean()
  
  buf <- taxons_offs_F1 %>% subset(Sexo == "Male"&Grupo =="MeHg") %>% select(Signif$tax[i])
  Signif[i,4] <-buf[[Signif$tax[i]]] %>% mean()
  
  buf <- taxons_offs_F1 %>% subset(Sexo == "Male"&Grupo =="VitA+MeHg") %>% select(Signif$tax[i])
  Signif[i,5] <-buf[[Signif$tax[i]]] %>% mean()
  
  Signif[i,6] <- unlist(Signif[i,2:5]) %>% mean()
  Signif[i,7] <- unlist(Signif[i,2:5]) %>% sd()
  }
for (i in 1:length(Signif$tax)) {
  buf <- taxons_dams %>% subset(Grupo =="CN") %>% select(Signif$tax[i])
  Signif[i,2] <-buf[[Signif$tax[i]]] %>% mean()
  
  buf <- taxons_dams %>% subset(Grupo =="VitA") %>% select(Signif$tax[i])
  Signif[i,3] <-buf[[Signif$tax[i]]] %>% mean()
  
  buf <- taxons_dams %>% subset(Grupo =="MeHg") %>% select(Signif$tax[i])
  Signif[i,4] <-buf[[Signif$tax[i]]] %>% mean()
  
  buf <- taxons_dams %>% subset(Grupo =="VitA+MeHg") %>% select(Signif$tax[i])
  Signif[i,5] <-buf[[Signif$tax[i]]] %>% mean()
  
  Signif[i,6] <- unlist(Signif[i,2:5]) %>% mean()
  Signif[i,7] <- unlist(Signif[i,2:5]) %>% sd()
}
for (i in 1:length(Signif$tax)) {
  buf <- taxons_1year_F2 %>% subset(sex == "Male"&treatment =="CN") %>% select(Signif$tax[i])
  Signif[i,2] <-buf[[Signif$tax[i]]] %>% mean()
  
  buf <- taxons_1year_F2 %>% subset(sex == "Male"&treatment =="VitA") %>% select(Signif$tax[i])
  Signif[i,3] <-buf[[Signif$tax[i]]] %>% mean()
  
  buf <- taxons_1year_F2 %>% subset(sex == "Male"&treatment =="MeHg") %>% select(Signif$tax[i])
  Signif[i,4] <-buf[[Signif$tax[i]]] %>% mean()
  
  buf <- taxons_1year_F2 %>% subset(sex == "Male"&treatment =="VitA+MeHg") %>% select(Signif$tax[i])
  Signif[i,5] <-buf[[Signif$tax[i]]] %>% mean()
  
  Signif[i,6] <- unlist(Signif[i,2:5]) %>% mean()
  Signif[i,7] <- unlist(Signif[i,2:5]) %>% sd()
}

#Signif_offs_F1Female <- Signif
#Signif_offs_F1Male <- Signif
#Signif_dams <- Signif
#Signif_1year_F1Female <- Signif
#Signif_1year_F1Male <- Signif
#Signif_1year_F2Female <- Signif
Signif_1year_F2Male <- Signif
###########################################################################################################
# create Z-Score dataframes
ZScore_list <- list(Signif_dams = Signif_dams,
                    Signif_offs_F1Female = Signif_offs_F1Female,
                    Signif_offs_F1Male = Signif_offs_F1Male,
                    Signif_1year_F1Female = Signif_1year_F1Female,
                    Signif_1year_F1Male = Signif_1year_F1Male,
                    Signif_1year_F2Female = Signif_1year_F2Female,
                    Signif_1year_F2Male = Signif_1year_F2Male)

for (i in 1:length(ZScore_list)) {
  buf <- ZScore_list[[i]]
  for (j in 1:length(buf$tax)) {
    buf[j,]$CN <- (ZScore_list[[i]][j,]$CN - ZScore_list[[i]][j,]$med)/ZScore_list[[i]][j,]$sd
    buf[j,]$VitA <- (ZScore_list[[i]][j,]$VitA - ZScore_list[[i]][j,]$med)/ZScore_list[[i]][j,]$sd
    buf[j,]$MeHg <- (ZScore_list[[i]][j,]$MeHg - ZScore_list[[i]][j,]$med)/ZScore_list[[i]][j,]$sd
    buf[j,]$`VitA+MeHg` <- (ZScore_list[[i]][j,]$`VitA+MeHg` - ZScore_list[[i]][j,]$med)/ZScore_list[[i]][j,]$sd
  }
  ZScore_list[[i]] <- buf %>% select(-tax,-med,-sd)
}

rm(i,j,buf,Signif,teste,teste_shapiro,vec)
###########################################################################################################
#x <- ZScore_list[["Signif_dams"]]
#x <- rbind(ZScore_list[["Signif_offs_F1Female"]],ZScore_list[["Signif_offs_F1Male"]])
#x <- rbind(ZScore_list[["Signif_1year_F1Female"]],ZScore_list[["Signif_1year_F1Male"]])
x <- rbind(ZScore_list[["Signif_1year_F2Female"]],ZScore_list[["Signif_1year_F2Male"]])
x

y <- data.matrix(x) #transformar o objeto numa matriz#
y
#aqui pode transformar diretamente a planilha numa matriz, sem transform?-la num objeto antes#
heatmap.2(y) #pedir heatmap da matriz# #erro de margens muito grandes ? comum#

yb <-colorRampPalette(c("blue", "white", "red"))
heatmap.2(y,key= TRUE, col=yb, margins = c(10,25), Rowv = FALSE, Colv = FALSE,
          dendrogram = NULL, trace = 'none', denscol = NULL, cexRow = 1.7, cexCol = 1.7,
          key.title = "Z-Score",density.info = 'none',key.xlab = NA, keysize = 1)

heatmap.2(y,key= TRUE, col=yb, margins = c(10,25), Rowv = FALSE, Colv = FALSE,
          dendrogram = NULL, trace = 'none', denscol = NULL, cexRow = 1.7, cexCol = 1.7,
          key.title = "Z-Score",density.info = 'none',key.xlab = NA, keysize = 1,
          rowsep = c(29),sepwidth = c(0.05,0.15))

##
tiff('images/ZScore/1year_F2_Heatmap_Zscore.tiff', width = 10,height = 10, units = "in", res=200)
heatmap.2(y,key= TRUE, col=yb, margins = c(10,25), Rowv = FALSE, Colv = FALSE,
          dendrogram = NULL, trace = 'none', denscol = NULL, cexRow = 1.7, cexCol = 1.7,
          key.title = "Z-Score",density.info = 'none',key.xlab = NA, keysize = 1,
          rowsep = c(8),sepwidth = c(0.5,0.15)
          )
dev.off()
