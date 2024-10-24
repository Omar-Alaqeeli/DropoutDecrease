library(xtable)
 
library(scales) 
 
library("RColorBrewer")

library("gplots")

library("RColorBrewer")


library(lattice)

library(grid)

library(ComplexHeatmap)

library(colorRamp2)

library("gplots")

library("ggplot2")

library("heatmap3")

library(purrr)

library(imager)

require(ggplot2)

library(fields)

library(hrbrthemes)

library(geomtextpath)

library(plyr)

library(dplyr)

library(arm)

library(geodiv)


suppressPackageStartupMessages(library(SummarizedExperiment))
suppressPackageStartupMessages(library(MultiAssayExperiment))

files <- list.files("~/Desktop/Datasets", pattern = ".", all.files = FALSE, recursive = TRUE, full.names = TRUE)

zero_increment <- c()
where = 0

for (k in files) {
  k
  dataset <- readRDS(k)
  dataset_file <- basename(k) #get full dataset name
  dataset_id <- tools::file_path_sans_ext(dataset_file) #get full dataset name without extension
  where = where + 1
  dataset <- updateObject(dataset)
  experiments(dataset)
  (dataset <- experiments(dataset)[["gene"]])
  dataset <- assays(dataset)[["TPM"]]
  dim(dataset)
  
i = 1
j = 1
zero_count = 0
entry_count = 0

zero_increment <- array(, dim = c((nrow(dataset) * ncol(dataset))))

for (i in 1:nrow(dataset)) {
  for (j in 1:ncol(dataset)) {
    entry_count = entry_count + 1
    if(dataset[i,j] == "0") {
      zero_count = zero_count + 1
    }
    zero_increment[entry_count] = zero_count
    print(zero_count)
    print(paste0("-------------------------------> ", where))
  }
}
assign(paste(dataset_id), zero_increment) 
}
length(zero_increment)


write.table(dataset_id, paste0("dataset_id.txt"), quote = FALSE, col.names = FALSE, row.names = FALSE)
d <- read.delim("Zeros/dataset_id.txt")


png(file = paste0("dataset_id.png"))
#plot(1:length(dataset_id), dataset_id, type = "l", main = dataset_id, xlab="Number of gene abundance entries", ylab="The total number of zeros", col = "red")
plot(peaks)
dev.off()


df <- data.frame(id = 1:length(dataset_id),
                 zeros = c(dataset_id))
valley = 0
count = 0
valleys <- data.frame(zeroValley = c(),
                      position = c())

for (i in df$id) {
  if(valley == df$zeros[i+1])
    count = count + 1
  else
    valley = df$zeros[i+1]
    count = 0
    valleys <- rbind(valleys, list(count, i))
   
}

png(file = paste0("dataset_id.png"))
plot(1:length(EMTAB2805), EMTAB2805, type = "l", main = "EMTAB2805", xlab="Number of gene abundance entries", ylab="The total number of zeros", col = "red")
dev.off()

