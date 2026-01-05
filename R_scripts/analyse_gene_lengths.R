#!/usr/bin/env Rscript

suppressMessages(library(optparse))
suppressMessages(library(ggplot2))
suppressMessages(library(readr))  # Use read_tsv for safer import
suppressMessages(library(dplyr))

option_list <- list(
  make_option(c("-f", "--file"), type="character", help="Input file: all_gene_fir_lengths.tsv"),
  make_option(c("-p", "--plot"), action="store_true", default=FALSE, help="Generate scatter plots [optional]")
)

opt_parser <- OptionParser(option_list=option_list)
opt <- parse_args(opt_parser)

if (is.null(opt$file)) {
  stop("Please provide input file using -f")
}

# ---- Robust file reading with readr ----
cat("Reading file:", opt$file, "\n")

df <- tryCatch({
  read_tsv(opt$file, col_names=c(
    "Assembly", "Contig", "GeneID", "GeneLength", "UpstreamDistance", "DownstreamDistance", "Secreted"
  ), col_types = cols(
    Assembly = col_character(),
    Contig = col_character(),
    GeneID = col_character(),
    GeneLength = col_double(),
    UpstreamDistance = col_double(),
    DownstreamDistance = col_double(),
    Secreted = col_integer()
  ))
}, error = function(e) {
  stop("Failed to read TSV file: ", e$message)
})

cat("✅ Read", nrow(df), "rows with", ncol(df), "columns\n")

# Check for NA rows in critical fields
na_rows <- df %>% filter(is.na(GeneLength) | is.na(UpstreamDistance) | is.na(DownstreamDistance))
if (nrow(na_rows) > 0) {
  cat("⚠️ Skipping", nrow(na_rows), "rows with missing numeric values\n")
  df <- df %>% filter(!is.na(GeneLength) & !is.na(UpstreamDistance) & !is.na(DownstreamDistance))
}

# Convert Secreted to factor
df$Secreted <- as.factor(df$Secreted)

# ---- Summary statistics ----
summary_stats <- df %>%
  group_by(Secreted) %>%
  summarise(
    Mean = mean(GeneLength),
    Median = median(GeneLength),
    SD = sd(GeneLength),
    Count = n()
  )

print("📊 Gene length summary by 'Secreted':")
print(summary_stats)

# ---- Welch's t-test ----
cat("\n📈 Welch's t-test comparing gene lengths:\n")
print(t.test(GeneLength ~ Secreted, data=df))

# ---- Correlation ----
cat("\n🔗 Correlation: Gene length vs 5′ upstream distance\n")
print(cor.test(df$GeneLength, df$UpstreamDistance, method="pearson"))

cat("\n🔗 Correlation: Gene length vs 3′ downstream distance\n")
print(cor.test(df$GeneLength, df$DownstreamDistance, method="pearson"))

# ---- Optional plots ----
if (opt$plot) {
  g1 <- ggplot(df, aes(x=UpstreamDistance, y=GeneLength)) +
    geom_point(alpha=0.2) +
    theme_minimal() +
    labs(title="Gene Length vs 5′ FIR", x="5′ Intergenic Distance", y="Gene Length")

  g2 <- ggplot(df, aes(x=DownstreamDistance, y=GeneLength)) +
    geom_point(alpha=0.2) +
    theme_minimal() +
    labs(title="Gene Length vs 3′ FIR", x="3′ Intergenic Distance", y="Gene Length")

  ggsave("GeneLength_vs_UpstreamFIR.png", g1, width=6, height=4, dpi=300)
  ggsave("GeneLength_vs_DownstreamFIR.png", g2, width=6, height=4, dpi=300)
}
