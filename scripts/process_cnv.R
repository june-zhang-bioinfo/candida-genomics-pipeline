#!/usr/bin/env Rscript
options(repos = "https://cloud.r-project.org/")

# 1. PARSE COMMAND-LINE ARGUMENTS
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 3) {
  stop("❌ Error: Missing required arguments.\nUsage: Rscript process_cnv.R <depth_file> <gff_file> <genes_csv>")
}

samtools_depth_file <- args[1]
gff_file            <- args[2]
genes_csv_file      <- args[3]

print(paste("📦 Processing sample depth metrics for:", samtools_depth_file))
print(paste("🧬 Utilizing structural annotation file:", gff_file))
print(paste("📊 Utilizing targeted candidate gene list:", genes_csv_file))

# 2. LOAD DEPENDENCIES CLEANLY
library(dplyr)
library(data.table) 
library(rtracklayer)

# Load resources dynamically from parameters
genes_of_interest <- read.csv(genes_csv_file)

# 3. DEFINE FEATURE LOOKUP FUNCTION
samtools_depth_mean_coverage <- function(samtools_depth_file, prokka_gff_file) {
  cluster <- read.table(samtools_depth_file)
  track <- import(prokka_gff_file)
  
  # Create a dataframe of sequence name, start, end coordinates, and gene name
  cluster_ranges <- data.frame(seqnames(track), start(track), end(track), track$Name)
  
  # Unite gff file with depth information
  all_range_means <- numeric()
  for (i in 1:nrow(cluster_ranges)) {
    range_of_interest <- cluster_ranges$start.track[i] : cluster_ranges$end.track[i]
    subset1 <- subset(cluster, subset = cluster_ranges$seqnames.track[i] == cluster$V1 & cluster$V2 %in% range_of_interest)
    range_mean <- mean(subset1$V3)
    all_range_means <- c(all_range_means, range_mean)
  }
  
  cluster_ranges <- cbind(cluster_ranges, all_range_means)
  return(cluster_ranges)
}

# 4. EXECUTE ANALYSIS LOGIC
cluster_ranges <- samtools_depth_mean_coverage(samtools_depth_file, gff_file)
cluster_ranges <- merge(cluster_ranges, genes_of_interest, by.x = "track.Name", by.y = "Cparap.ID", all.x = TRUE)

# Chromosomal normalization framework
df <- cluster_ranges %>%
  group_by(seqnames.track.) %>%
  mutate(
    avg_read_depth = mean(all_range_means, na.rm = TRUE), 
    depth_ratio = all_range_means / avg_read_depth
  ) %>%
  ungroup() %>%
  arrange(desc(depth_ratio)) %>%
  distinct(track.Name, .keep_all = TRUE) %>%
  filter(!is.na(track.Name))

# 5. WRITE OUTPUT ADJACENT TO THE DEPTH FILE
output_csv <- paste0(samtools_depth_file, ".csv")
write.csv(df, output_csv, row.names = FALSE)
print(paste("✅ Profile matching sequence complete. Generated:", output_csv))