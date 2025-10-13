#!/usr/bin/env Rscript

# --- Load libraries ---
args <- commandArgs(trailingOnly = F)  
scriptPath <- normalizePath(dirname(sub("^--file=", "", args[grep("^--file=", args)])))
source(paste(scriptPath, "/R_functions.R", sep=""))

add_dependency("ggplot2")
add_dependency("dplyr") 
add_dependency("gghighlight")
add_dependency("optparse")
add_dependency("stringr")
add_dependency("viridis")
add_dependency("data.table")

# --- Parse options ---
option_list = list(
  make_option(c("-d", "--dataframe"), default=NULL, help="Dataframe file name"),
  make_option(c("-g", "--highlight"), default=NULL, help="Highlight this GO term (optional)"),
  make_option(c("-o", "--output"), default=NULL, help="Output PDF"),
  make_option(c("-i", "--height"), default=8, help="Height of pdf (inches)"),
  make_option(c("-w", "--width"), default=12, help="Width of pdf (inches)")
)

opt_parser = OptionParser(usage = "Usage: %prog -d <dataframe>", option_list=option_list)
opt = parse_args(opt_parser)

if (is.null(opt$dataframe)) {
  print_help(opt_parser)
  stop("Dataframe must be supplied (input file)", call.=FALSE)
}

# --- Read and clean data ---
if (!file.exists(opt$dataframe)) stop("Error: dataframe does not exist", call.=FALSE)

library(data.table)
#data <- read.table(opt$dataframe, sep='\t', header=TRUE, stringsAsFactors = FALSE)
data <- fread(opt$dataframe, sep = "\t", header = TRUE, data.table = FALSE)

# Column names
#colnames(data) <- c("contig","feature","desc","strand","upstream_distance","downstream_distance","start_pos","uplog10","downlog10")

# Replace NA and zero values
data$desc[data$desc == "" | is.na(data$desc)] <- "none"
data$upstream_distance[data$upstream_distance == 0] <- 0.5
data$downstream_distance[data$downstream_distance == 0] <- 0.5

# Recalculate log10 values
data$uplog10 <- log10(data$upstream_distance)
data$downlog10 <- log10(data$downstream_distance)

# Check how many secreted genes there are (use from secreted plot step)
secreted_count <- 4372912

# Downsample non-secreted to match secreted count
set.seed(42)  # For reproducibility
if (nrow(data) > secreted_count) {
  data <- data[sample(nrow(data), secreted_count), ]
}

# Output file
output <- ifelse(is.null(opt$output),
                 paste(opt$dataframe, "-densityplot_viridis.pdf", sep=""),
                 opt$output)

# Compute medians
#median_up <- median(data$uplog10, na.rm=TRUE)
#median_down <- median(data$downlog10, na.rm=TRUE)

# medians for all genes =
# perl -ane 'push @a, $F[7] if $. > 1; END { @a = sort { $a <=> $b } @a; print $a[@a/2], "\n" }' All_eukaryotes_distance_file_with_unique_contigs_and_gene_ids2_no_dup_lines.tab
median_up <- 2.99694924849538
# perl -ane 'push @a, $F[8] if $. > 1; END { @a = sort { $a <=> $b } @a; print $a[@a/2], "\n" }' All_eukaryotes_distance_file_with_unique_contigs_and_gene_ids2_no_dup_lines.tab
median_down <- 3.09131515969722


# Create plot
pdf(output, height=opt$height, width=opt$width)

plot <- ggplot(data, aes(x=downlog10, y=uplog10)) +
  geom_bin2d(bins = 60) +
  geom_density_2d(color = "black", alpha = 0.3, bins = 10) +  # Density contours
  geom_vline(xintercept = median_down, linetype = "dashed", linewidth = 1, color = "dodgerblue4") +
  geom_hline(yintercept = median_up, linetype = "dashed", linewidth = 1, color = "dodgerblue4") +
  scale_fill_gradientn(
    colours = viridis::viridis(10),
    trans = "log10",
    breaks = c(1, 10, 100, 1000, 10000),
    labels = c("1", "10", "100", "1K", "10K+"),
    limits = c(1, 10000),  # <<< THIS fixes the gradient range
    oob = scales::squish,   # Prevents values >10K or <1 from erroring
    name = "Gene Count"
  ) +
  labs(
    x = "3′ log10 FIR (bp)",
    y = "5′ log10 FIR (bp)",
    title = "Flanking Intergenic Regions (FIR) Density Plot"
  ) +
  coord_cartesian(xlim = c(1, 7), ylim = c(1, 7)) +
  theme_minimal(base_size = 14) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "right"
  )

print(plot)
dev.off()