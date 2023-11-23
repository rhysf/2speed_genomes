#!/usr/bin/Rscript
args <- commandArgs(trailingOnly = F)  
scriptPath <- normalizePath(dirname(sub("^--file=", "", args[grep("^--file=", args)])))
R_functions <- paste(scriptPath, "/R_functions.R", sep="")
source(R_functions)
add_dependency("dplyr") 
add_dependency("gghighlight")
add_dependency("optparse")
add_dependency("hash")

### r.farrer@exeter.ac.uk 
# modified from Theresa Wacker script: 
# https://bitbucket.org/Theresa_42/wackeretal_2022_bsal_2speedgenome/src/master/HypergeometricTest_DensityPlots_and_Histo_inverse_cmdline_looped_cor.R

# Opening commands
option_list = list(
    make_option(c("-d", "--dataframe"), default=NULL, help="Dataframe file name")
);
opt_parser = OptionParser(usage = "Usage: %prog -d <dataframe>
Optional: -o Output [opt$dataframe-Enrichment_tests.tab]

Notes: Dataframe should look like this:
contig  feature desc    strand  upstream_distance       downstream_distance     uplog10 downlog10
scaffold0001    BSAL2_000002  cat1,cat2,cat3          -       1271    221     3.10414555055401        2.34439227368511
scaffold0001    BSAL2_000003  cat4            +       1069    221     3.02897770520878        2.34439227368511", option_list=option_list);
opt = parse_args(opt_parser);
if(is.null(opt$dataframe)) {
  print_help(opt_parser)
  stop("Dataframe must be supplied (input file)", call.=FALSE)
}
options(warn = 0)
if(!file.exists(opt$dataframe)) { stop("Error: dataframe does not appear to be valid", call.=FALSE) }
data = read.table(opt$dataframe, com='', sep='\t', quote = "", header=T)

# outfile
name = opt$dataframe
outfile <- paste(opt$dataframe, "-Enrichment_tests.tab", sep="")
if(!is.null(opt$output)) { outfile <- opt$output }

# outfile2 significant hits (with multiple correction)
outfile2 <- paste(outfile, "-significant-hyper-BH.tab", sep="")
outfile3 <- paste(outfile, "-significant-chi-BH.tab", sep="")

#calculate median for 5' and for 3' genes
median_up = median(data$uplog10, na.rm=TRUE)
median_down= median(data$downlog10, na.rm=TRUE)

# IMPORTANT: quadrants as follows: 
# ------------------------------
#|  Q1          |           Q2 |
#|              |              |
#---------median uplog---------|
#|  Q4          |           Q3 |
#|              |              |
#|--------------|--------------|
#               |<- median downlog

### Quadrants for all genes

# Q1
Q1_temp = data%>%filter(uplog10 > median_up)
Q1 = Q1_temp %>% filter(downlog10 < median_down)
# Q2
Q2 = Q1_temp %>% filter(downlog10 > median_down)
# Q3
Q3_temp = data %>% filter(downlog10 > median_down)
Q3 = Q3_temp %>% filter(uplog10 < median_up)
# Q4
Q4_temp = data %>% filter(uplog10 < median_up)
Q4 = Q4_temp %>% filter(downlog10 < median_down)

# k = sample size (The number of genes in quadrant)
k1 = nrow(Q1)
k2 = nrow(Q2)
k3 = nrow(Q3)
k4 = nrow(Q4)

# N = total number of genes, n the number of non-successes (all non-M36s)
N = nrow(data)

# Split categories
all_categories <- strsplit(data$desc, ",")
unique_categories <- unique(unlist(all_categories))

# Save everything to memory in order to to do padjust and print in a clear way!
GO_sum <- hash()
GO_everything_else_sum <- hash()
GO_UL_sum <- hash()
GO_UR_sum <- hash()
GO_LR_sum <- hash()
GO_LL_sum <- hash()
GO_hyper_p_UL <- hash()
GO_hyper_p_UR <- hash()
GO_hyper_p_LR <- hash()
GO_hyper_p_LL <- hash()
GO_chi_p_UL <- hash()
GO_chi_p_UR <- hash()
GO_chi_p_LR <- hash()
GO_chi_p_LL <- hash()

# Loop over categories
for (category in unique_categories) {
	#print(paste("GO term:", category))
	table = data %>% select(uplog10,downlog10,desc,feature)
  
	# Save lines with the matching category
	matching_rows <- data %>% filter(grepl(category, data$desc))

	# m = the number of possible successes in the test, aka for instance all M36s
	m = nrow(matching_rows)
	GO_sum[[category]] <- m

	# n = the number of non-successes 
	n = (N-m)
	GO_everything_else_sum[[category]] <- n

	# x = number of marked elements in the selection
	# Q1
	x1_temp  = matching_rows %>% filter(uplog10 > median_up)
	x1_temp2 = x1_temp %>% filter(downlog10 < median_down)
	x1 = nrow(x1_temp2)
	GO_UL_sum[[category]] <- x1

	# Q2
	x2_temp = x1_temp %>% filter(downlog10 > median_down)
	x2 = nrow(x2_temp)
	GO_UR_sum[[category]] <- x2

	# Q3
	x3_temp  = matching_rows %>% filter(downlog10 > median_down)
	x3_temp2 = x3_temp %>% filter(uplog10 < median_up)
	x3 = nrow(x3_temp2)
	GO_LR_sum[[category]] <- x3

	# Q4
	x4_temp  = matching_rows %>% filter(uplog10 < median_up)
	x4_temp2 = x4_temp %>% filter(downlog10 < median_down)
	x4 = nrow(x4_temp2)
	GO_LL_sum[[category]] <- x4

	# Hypergeometric Test and Chisquared Test 
	x_values = c(x1, x2, x3, x4)
	k_values = c(k1, k2, k3, k4)
	hypergeometrics = c('na','na','na','na')
	chisquared = c('na','na','na','na')
	for (i in 1:4){
		xval = x_values[i]
		kval = k_values[i]

		# Hypergeometric Test
		hypergeometrics[i] = phyper(xval-1,m,n,kval,lower.tail=FALSE)

		# chisquared
		noselection = kval-xval
		none_x = m-xval
		nothing = N-kval-none_x
		row1 = cbind(xval,none_x)
		row2 = cbind(noselection,nothing)
		table = rbind(row1, row2)
		names = c("quadrant","none-quadrant")
		colnames(table) = names
		rows = c("successes","none_successes")
		rownames(table) = rows
		chisq = chisq.test(table)
		chisquared[i] = chisq$p.value
	}

	# save
	GO_hyper_p_UL[[category]] <- hypergeometrics[1]
	GO_hyper_p_UR[[category]] <- hypergeometrics[2]
	GO_hyper_p_LR[[category]] <- hypergeometrics[3]
	GO_hyper_p_LL[[category]] <- hypergeometrics[4]
	GO_chi_p_UL[[category]] <- chisquared[1]
	GO_chi_p_UR[[category]] <- chisquared[2]
	GO_chi_p_LR[[category]] <- chisquared[3]
	GO_chi_p_LL[[category]] <- chisquared[4]
}

# convert p values in hash to vector
GO_categories_vector <- as.vector(keys(GO_hyper_p_UL))
GO_hyper_p_UL_vector <- as.vector(values(GO_hyper_p_UL))
GO_hyper_p_UR_vector <- as.vector(values(GO_hyper_p_UR))
GO_hyper_p_LR_vector <- as.vector(values(GO_hyper_p_LR))
GO_hyper_p_LL_vector <- as.vector(values(GO_hyper_p_LL))
GO_chi_p_UL_vector <- as.vector(values(GO_chi_p_UL))
GO_chi_p_UR_vector <- as.vector(values(GO_chi_p_UR))
GO_chi_p_LR_vector <- as.vector(values(GO_chi_p_LR))
GO_chi_p_LL_vector <- as.vector(values(GO_chi_p_LL))

# p.adjust
GO_hyper_p_UL_adj <- p.adjust(GO_hyper_p_UL_vector, method="BH")
GO_hyper_p_UR_adj <- p.adjust(GO_hyper_p_UR_vector, method="BH")
GO_hyper_p_LR_adj <- p.adjust(GO_hyper_p_LR_vector, method="BH")
GO_hyper_p_LL_adj <- p.adjust(GO_hyper_p_LL_vector, method="BH")
GO_chi_p_UL_adj <- p.adjust(GO_chi_p_UL_vector, method="BH")
GO_chi_p_UR_adj <- p.adjust(GO_chi_p_UR_vector, method="BH")
GO_chi_p_LR_adj <- p.adjust(GO_chi_p_LR_vector, method="BH")
GO_chi_p_LL_adj <- p.adjust(GO_chi_p_LL_vector, method="BH")

# Print

# save the numbers k per quadrant & write to file
print_blank_line = paste("", sep = '')
print_header = paste("#ALL genes Q_UL:", k1, "ALL genes Q_UR:", k2, "ALL genes Q_LR:", k3, "ALL genes Q_LL:", k4, sep = '\t')
print_line = paste("---------------", sep = '')
cat(print_header, file = outfile, append = TRUE, sep="\n")

# Loop over categories
i=0
#for (category in unique_categories) {
for (category in GO_categories_vector) {
	i=i+1

	# Counts of genes with GO-term and genes without that GO-term
	print_sum = paste(category, " sum (m):\t", GO_sum[[category]], sep = '')
	print_sum2 = paste("not-", category, " (n):\t", GO_everything_else_sum[[category]], sep = '')
	cat(print_line, file = outfile, append = TRUE, sep="\n")
	cat(print_sum, file = outfile, append = TRUE, sep="\n")
	cat(print_sum2, file = outfile, append = TRUE, sep="\n")
	cat(print_blank_line, file = outfile, append = TRUE, sep="\n")

	# Genes with GO term in each quadrant
	print_GO_per_quadrant1 = paste(category, " genes Q_UL:\t", GO_UL_sum[[category]], sep = '')
	print_GO_per_quadrant2 = paste(category, " genes Q_UR:\t", GO_UR_sum[[category]], sep = '')
	print_GO_per_quadrant3 = paste(category, " genes Q_LR:\t", GO_LR_sum[[category]], sep = '')
	print_GO_per_quadrant4 = paste(category, " genes Q_LL:\t", GO_LL_sum[[category]], sep = '')
	cat(print_GO_per_quadrant1, file = outfile, append = TRUE, sep="\n")
	cat(print_GO_per_quadrant2, file = outfile, append = TRUE, sep="\n")
	cat(print_GO_per_quadrant3, file = outfile, append = TRUE, sep="\n")
	cat(print_GO_per_quadrant4, file = outfile, append = TRUE, sep="\n")
	cat(print_blank_line, file = outfile, append = TRUE, sep="\n")

	# Print results table for hyper
	print_hyper1 = paste(category, " p(hypergeometric) Q_UL:\t", GO_hyper_p_UL[[category]], sep = '')
	print_hyper2 = paste(category, " p(hypergeometric) Q_UR:\t", GO_hyper_p_UR[[category]], sep = '')
	print_hyper3 = paste(category, " p(hypergeometric) Q_LR:\t", GO_hyper_p_LR[[category]], sep = '')
	print_hyper4 = paste(category, " p(hypergeometric) Q_LL:\t", GO_hyper_p_LL[[category]], sep = '')
	cat(print_hyper1, file = outfile, append = TRUE, sep="\n")
	cat(print_hyper2, file = outfile, append = TRUE, sep="\n")
	cat(print_hyper3, file = outfile, append = TRUE, sep="\n")
	cat(print_hyper4, file = outfile, append = TRUE, sep="\n")
	cat(print_blank_line, file = outfile, append = TRUE, sep="\n")

	# Print results table for chi2
	print_chi1 = paste(category, " p(chi-squared) Q_UL:\t", GO_chi_p_UL[[category]], sep = '')
	print_chi2 = paste(category, " p(chi-squared) Q_UR:\t", GO_chi_p_UR[[category]], sep = '')
	print_chi3 = paste(category, " p(chi-squared) Q_LR:\t", GO_chi_p_LR[[category]], sep = '')
	print_chi4 = paste(category, " p(chi-squared) Q_LL:\t", GO_chi_p_LL[[category]], sep = '')
	cat(print_chi1, file = outfile, append = TRUE, sep="\n")
	cat(print_chi2, file = outfile, append = TRUE, sep="\n")
	cat(print_chi3, file = outfile, append = TRUE, sep="\n")
	cat(print_chi4, file = outfile, append = TRUE, sep="\n")
	cat(print_blank_line, file = outfile, append = TRUE, sep="\n")

	# Multiple correction for hyper
	print_hyper_adj1 = paste(category, " p(hypergeometric) Q_UL (Benjamini, Hochberg):\t", GO_hyper_p_UL_adj[i], sep = '')
	print_hyper_adj2 = paste(category, " p(hypergeometric) Q_UR (Benjamini, Hochberg):\t", GO_hyper_p_UR_adj[i], sep = '')
	print_hyper_adj3 = paste(category, " p(hypergeometric) Q_LR (Benjamini, Hochberg):\t", GO_hyper_p_LR_adj[i], sep = '')
	print_hyper_adj4 = paste(category, " p(hypergeometric) Q_LL (Benjamini, Hochberg):\t", GO_hyper_p_LL_adj[i], sep = '')
	cat(print_hyper_adj1, file = outfile, append = TRUE, sep="\n")
	cat(print_hyper_adj2, file = outfile, append = TRUE, sep="\n")
	cat(print_hyper_adj3, file = outfile, append = TRUE, sep="\n")
	cat(print_hyper_adj4, file = outfile, append = TRUE, sep="\n")
	cat(print_blank_line, file = outfile, append = TRUE, sep="\n")

	# Multiple correction for chi2
	print_chi_adj1 = paste(category, " p(chi-squared) Q_UL (Benjamini, Hochberg):\t", GO_chi_p_UL_adj[i], sep = '')
	print_chi_adj2 = paste(category, " p(chi-squared) Q_UR (Benjamini, Hochberg):\t", GO_chi_p_UR_adj[i], sep = '')
	print_chi_adj3 = paste(category, " p(chi-squared) Q_LR (Benjamini, Hochberg):\t", GO_chi_p_LR_adj[i], sep = '')
	print_chi_adj4 = paste(category, " p(chi-squared) Q_LL (Benjamini, Hochberg):\t", GO_chi_p_LL_adj[i], sep = '')
	cat(print_chi_adj1, file = outfile, append = TRUE, sep="\n")
	cat(print_chi_adj2, file = outfile, append = TRUE, sep="\n")
	cat(print_chi_adj3, file = outfile, append = TRUE, sep="\n")
	cat(print_chi_adj4, file = outfile, append = TRUE, sep="\n")
	cat(print_blank_line, file = outfile, append = TRUE, sep="\n")	

	# significant hyper hits
	if(!is.na(GO_hyper_p_UL_adj[i]) && GO_hyper_p_UL_adj[i] < 0.01) {
		cat(print_hyper_adj1, file = outfile2, append = TRUE, sep="\n")
	}
	if(!is.na(GO_hyper_p_UR_adj[i]) && GO_hyper_p_UR_adj[i] < 0.01) {
		cat(print_hyper_adj2, file = outfile2, append = TRUE, sep="\n")
	}
	if(!is.na(GO_hyper_p_LR_adj[i]) && GO_hyper_p_LR_adj[i] < 0.01) {
		cat(print_hyper_adj3, file = outfile2, append = TRUE, sep="\n")
	}
	if(!is.na(GO_hyper_p_LL_adj[i]) && GO_hyper_p_LL_adj[i] < 0.01) {
		cat(print_hyper_adj4, file = outfile2, append = TRUE, sep="\n")
	}

	# significant chi2 hits
	if(!is.na(GO_chi_p_UL_adj[i]) && GO_chi_p_UL_adj[i] < 0.01) {
		cat(print_chi_adj1, file = outfile3, append = TRUE, sep="\n")
	}
	if(!is.na(GO_chi_p_UR_adj[i]) && GO_chi_p_UR_adj[i] < 0.01) {
		cat(print_chi_adj2, file = outfile3, append = TRUE, sep="\n")
	}
	if(!is.na(GO_chi_p_LR_adj[i]) && GO_chi_p_LR_adj[i] < 0.01) {
		cat(print_chi_adj3, file = outfile3, append = TRUE, sep="\n")
	}
	if(!is.na(GO_chi_p_LL_adj[i]) && GO_chi_p_LL_adj[i] < 0.01) {
		cat(print_chi_adj4, file = outfile3, append = TRUE, sep="\n")
	}
}

# Check for warnings and print them
if (length(warnings()) > 0) {
  cat("Warnings:\n")
  warnings()
}
