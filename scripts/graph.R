library("phyloseq")
library("ggplot2")
library("broom")
library("ggpubr")
library("agricolae")
library("rstatix")
library("multcompView")
library("dplyr")
library("tidyr")
library("purrr")
setwd("~/Desktop/thiaguinho")
load("lefse.RData")
if(!dir.exists("significant/")) { dir.create("significant/") }
mapfile <- import_qiime_sample_data("qiime/manifest.csv") #deve estar delimitado por espaços
sample_names(mapfile) <- mapfile$Sample.id
sample_data(ps1) <- mapfile
taxa_names(ps.m) <- paste0("ASV", seq(ntaxa(ps.m)))
ps.t <- ps1
ps.t.m <- subset_samples(ps1, Sex == "Male")
ps.t.f <- subset_samples(ps1, Sex == "Female")

sample_data(ps.t)$final_group <- gsub("-", "", sample_data(ps.t)$final_group)
taxonomic_levels <- names(as.data.frame(tax_table(ps.t.f)))
taxonomic_levels <- taxonomic_levels[taxonomic_levels != "Kingdom"]

#versão automatizada
for (level in taxonomic_levels) {
  print(level)
  glom <- tax_glom(ps.t , get("level"), NArm = TRUE)
  print(glom)
  glom <- transform_sample_counts(glom, function(x) x / sum(x))
  
  glommelt <- psmelt(glom)[, -c(11:12)]
  taxon_level <- tail(names(glommelt), 1)
  
  if(!dir.exists(paste0("significant/", get("taxon_level")))) { dir.create(paste0("significant/", get("taxon_level"))) }
  
  group_level <- c("M.I_O.S", "M.I_O.3", "M.M.S_O.S", "M.M.S_O.3", "F.I_O.S", "F.I_O.3", "F.M.S_O.S", "F.M.S_O.3")
  
  anova <- glommelt %>%
    nest(data = -get("taxon_level")) %>%
    mutate(test = map(.x=data, ~aov(Abundance~final_group, data=.x) %>% tidy)) %>%
    unnest(test) %>%
    mutate(p.adjust = p.adjust(p.value, method="BH")) %>%
    filter(p.adjust < 0.05)
  
  tukey <- list()
  
  for (taxon in anova[[get("taxon_level")]]) {
    
    temp <- HSD.test(aov(Abundance~final_group, glommelt[glommelt[[get("taxon_level")]] %in% taxon, ]), "final_group", group = TRUE)
    temp$groups$final_group <- rownames(temp$groups)
    temp$groups <- temp$groups %>% arrange(final_group, levels = group_level)
    
    if (length(unique(temp$groups$groups))!=1) {
      
      tukey[[taxon]] <- as.data.frame(temp$groups)
      
    }
    
  }
  
  for (taxon in names(tukey))  {
    
    print(taxon)
    print(max((glommelt %>% filter(., get(taxon_level) %in% taxon))$Abundance)*1.05)
    print(tukey[[taxon]]$groups)
    
    p <- glommelt %>% filter(., get(taxon_level) %in% taxon) %>%
      ggplot(data = ., aes(x = final_group, y = Abundance)) +
      geom_boxplot(outlier.shape  = NA,
                   aes(color = final_group)) +
      labs(color=NULL)+
      theme(axis.text.x = element_blank(), axis.ticks.x = element_blank()) +
      # geom_jitter(aes(color = final_group), height = 0, width = .2) +
      labs(x = "", y = "Relative abundance\n") +
      ggtitle(label = get("taxon"), subtitle = paste("Anova, p-value =", formatC((anova %>% filter(., get(taxon_level) %in% taxon))$p.adjust, format = "e", digits = 2))) +
      theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 15),
            plot.subtitle = element_text(hjust = 0.5)) +
      geom_text(data = data.frame(),
                aes(x = tukey[[taxon]]$final_group,
                    y = (max((glommelt %>% filter(., get(taxon_level) %in% taxon))$Abundance)*1.05),
                    label = tukey[[taxon]]$groups),
                col = 'black',size = 3.5,inherit.aes = FALSE )
    
    p$data$final_group <- factor(p$data$final_group, levels = group_level)
    
    tiff(paste0("significant/", get("taxon_level"), "/unstratified_",  get("taxon_level"), "_", get("taxon") ,".tif"), width = 9, height = 7, units = "in", res = 450)
    print(p)
    dev.off()
    
  }
  
}

###versão manual
glom <- tax_glom(ps.t, "Kingdom", NArm = TRUE)
glom <- transform_sample_counts(glom, function(x) x / sum(x))


glommelt <- psmelt(glom)[, -c(10:11)]
taxon_level <- tail(names(glommelt), 1)

if(!dir.exists(paste0("significant/", get("taxon_level")))) { dir.create(paste0("significant/", get("taxon_level"))) }


group_level <- c("WISTARNAIVE21","WISTARNAIVE30","WISTARNAIVE40",
                 "WISTARSHAM21","WISTARSHAM30","WISTARSHAM40",
                 "WISTARMALF21","WISTARMALF30","WISTARMALF40",
                 "WARNAIVE21","WARNAIVE30","WARNAIVE40",
                 "WARSHAM21","WARSHAM30","WARSHAM40",
                 "WARMALF21","WARMALF30","WARMALF40")

lmresults <- glommelt %>%
  nest(data = -get("taxon_level")) %>%
  mutate(test = map(.x=data, ~lm(Abundance~final_group, data=.x) %>% tidy)) %>%
  unnest(test) %>%
  mutate(p.adjust = p.adjust(p.value, method="BH")) %>%
  filter(p.adjust < 0.05)

tukey <- list()
for (taxon in anova[[get("taxon_level")]]) {
  
  temp <- HSD.test(lm(Abundance~final_group, glommelt[glommelt[[get("taxon_level")]] %in% taxon, ]), "final_group", group = TRUE)
  temp$groups$final_group <- rownames(temp$groups)
  temp$groups <- temp$groups %>% arrange(final_group, levels = group_level)
  if (length(unique(temp$groups$groups))!=1) {
    
    tukey[[taxon]] <- as.data.frame(temp$groups)
    
  }
  
}


for (taxon in names(tukey)) {
  print(taxon)
  print(max((glommelt %>% filter(., get(taxon_level) %in% taxon))$Abundance)*1.05)
  print(tukey[[taxon]]$groups)
  # print(glommelt %>% filter(., Phylum %in% "taxon"))
  
  p <- glommelt %>% filter(., get(taxon_level) %in% taxon) %>%
    ggplot(data = ., aes(x = final_group, y = Abundance)) +
    geom_boxplot(outlier.shape  = NA,
                 aes(color = final_group)) +
    theme(axis.text.x = element_blank(), axis.ticks.x = element_blank()) +
    # geom_jitter(aes(color = final_group), height = 0, width = .2) +
    labs(x = "", y = "Relative abundance\n") +
    ggtitle(label = get("taxon"), subtitle = paste("Anova, p-value =", formatC((anova %>% filter(., get(taxon_level) %in% taxon))$p.adjust, format = "e", digits = 2))) +
    theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 17),
          plot.subtitle = element_text(hjust = 0.5)) +
    geom_text(data = data.frame(),
              aes(x = tukey[[taxon]]$final_group,
                  y = (max((glommelt %>% filter(., get(taxon_level) %in% taxon))$Abundance)*1.05),
                  label = tukey[[taxon]]$groups),
              col = 'black',size = 4,inherit.aes = FALSE )
  # stat_compare_means(method = "anova", label.x = 17,label.y = max((glommelt %>% filter(., get(taxon_level) %in% taxon))$Abundance)*1)
  p
  p$data$final_group <- factor(p$data$final_group, levels = group_level)
  print(p)
  
  tiff(paste0("significant/", get("taxon_level"), "/",  get("taxon_level"), "_", get("taxon") ,".tif"), width = 10, height = 8, units = "in", res = 600)
  print(p)
  dev.off()
}

# 
# images
# 
# geom_text(data = data.frame(),
#           aes(x = group_data$LGI,
#               y = max(taxons[,i+6])*1.05,
#               label = group_data$groups,
#               color = LGI),
#           #position = position_dodge(width = 1),vjust = -7,
#           col = 'black',size = 3,inherit.aes = FALSE )
# 
# p
# 
# 
# p$data$final_group <- factor(p$data$final_group, levels = group_level)
# tiff("Phylum_boxplot_fernanda_santos.tiff", height = 15, width = 25, units = "in", res = 300)
# p
# dev.off()
# 
# # taxonlst <- unique(glommelt$Phylum)

# res <- list()
# for (taxon in anova$Phylum) {
#   
#   temp <- glommelt %>%
#     filter(Phylum == taxon)  %>%
#     tukey_hsd(Abundance~final_group) %>%
#     filter(p.adj.signif != "ns")
#   
#   res[[taxon]] <- p
# 
# }
