library(readr)
library(stringr)

# Read and process KO data
liness <- read.delim2("./KO_source_kegg")
lines <- liness$X.D

# Initialize variables
level_A <- level_B <- level_C <- NA
data <- list()

# Process hierarchy
for (line in lines) {
  line <- trimws(line)
  
  if (startsWith(line, "A")) {
    content <- trimws(substr(line, 2, nchar(line))) 
    level_A <- substring(content, 2)
  } else if (startsWith(line, "B")) {
    content <- trimws(substr(line, 2, nchar(line))) 
    level_B <- substring(content, 2)
  } else if (startsWith(line, "C")) {
    content <- trimws(substr(line, 2, nchar(line))) 
    level_C <- substring(content, 2)
  } else if (startsWith(line, "D")) {
    content <- trimws(substr(line, 2, nchar(line))) 
    level_D <- substring(content, 2)
    data <- append(data, list(c(level_A, level_B, level_C, level_D)))
  }
}

# Convert to data frame
data <- do.call(rbind, data)
colnames(data) <- c("Level1", "Level2", "Level3", "KO")
data2 <- as.data.frame(data, stringsAsFactors = FALSE)

# Flattened data frame
print(data2)
dim(data2) #  61655     4

# Extract EC, Ko, Kegg access number
data3 <- data2
data3$KEGG <- paste0("ko", sub(".*(\\d{5})\\]$", "\\1", data3$Level3))
data3$EC <- paste0("EC:", str_extract(data3$KO, "(?<=EC:)[^]]+"))
data3$KO_num <- paste0("K", substr(x = data3$KO , start = 1, stop = 6))
data3[is.na(data3)] <- ""

library(data.table)
fwrite(data3, "./Kegg_KO_17112025_ff.csv")

kegg_pathways <- data3 %>% dplyr::select(Level1,Level2,Level3,KEGG)
kegg_pathways <- unique(kegg_pathways)
rownames(kegg_pathways) <- kegg_pathways$KEGG

rm(liness,lines,level_A,level_B,level_C,data,data2,data3)
