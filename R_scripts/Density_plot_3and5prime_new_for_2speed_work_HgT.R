#!/usr/bin/Rscript
args <- commandArgs(trailingOnly = F)  
scriptPath <- normalizePath(dirname(sub("^--file=", "", args[grep("^--file=", args)])))
R_functions <- paste(scriptPath, "/R_functions.R", sep="")
source(R_functions)
add_dependency("dplyr") 
add_dependency("gghighlight")
add_dependency("optparse")

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
scaffold0001    BSAL2_000002            -       1271    221     3.10414555055401        2.34439227368511
scaffold0001    BSAL2_000003            +       1069    221     3.02897770520878        2.34439227368511", option_list=option_list);
opt = parse_args(opt_parser);
if(is.null(opt$dataframe)) {
  print_help(opt_parser)
  stop("Dataframe must be supplied (input file)", call.=FALSE)
}
if(!file.exists(opt$dataframe)) { stop("Error: dataframe does not appear to be valid", call.=FALSE) }
data = read.table(opt$dataframe, com='', sep='\t', quote = "", header=T)

# outfile
name = opt$dataframe
outfile <- paste(opt$dataframe, "-Enrichment_tests.tab", sep="")
if(!is.null(opt$output)) { outfile <- opt$output }

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

# Print gene IDs for quadrants as lists in files 
#outfileQ1=paste(outfile,name, "_Q1_gene_IDs.txt", sep="")
#outfileQ2=paste(outfile,name, "_Q2_gene_IDs.txt", sep="")
#outfileQ3=paste(outfile,name, "_Q3_gene_IDs.txt", sep="")
#outfileQ4=paste(outfile,name, "_Q4_gene_IDs.txt", sep="")
#write.table(as.character(Q1$feature), outfileQ1, row.names = FALSE, sep="\t", col.names = FALSE)
#write.table(as.character(Q2$feature), outfileQ2, row.names = FALSE, sep="\t", col.names = FALSE)
#write.table(as.character(Q3$feature), outfileQ3, row.names = FALSE, sep="\t", col.names = FALSE)
#write.table(as.character(Q4$feature), outfileQ4, row.names = FALSE, sep="\t", col.names = FALSE)

# save the numbers k per quadrant & write to file
#outfileQs=paste(outfile, "_numbers_in_quadrants.txt", sep="")
print0 = paste("", sep = '')
print1 = paste("ALL genes Q_UL:", k1, sep = '\t')
print2 = paste("ALL genes Q_UR:", k2, sep = '\t')
print3 = paste("ALL genes Q_LR:", k3, sep = '\t')
print4 = paste("ALL genes Q_LL:", k4, sep = '\t')
cat(print1, file = outfile, append = TRUE, sep="\n")
cat(print2, file = outfile, append = TRUE, sep="\n")
cat(print3, file = outfile, append = TRUE, sep="\n")
cat(print4, file = outfile, append = TRUE, sep="\n")
cat(print0, file = outfile, append = TRUE, sep="\n")

# N = total number of genes, n the number of non-successes (all non-M36s)
N = nrow(data)

#vectorize
#categories = as.vector(categories)#only works properly as a vector

#loop over categories (secreted etc.)
categories = c("cat1")
for (var in categories){
  table = data %>% select(uplog10,downlog10,desc,feature)
  
  #filter the categories
  #subset=data%>%filter(data$UQ(var)==!!(var))
  subset = data %>% filter(table$desc == 'cat1,')

  # m = the number of possible successes in the test, aka for instance all M36s
  m = nrow(subset)
  # n = the number of non-successes 
  n = (N-m)

  print1 = paste(var, " sum (m):\t", m, sep = '')
  print2 = paste("not-", var, " (n):\t", n, sep = '')
  cat(print1, file = outfile, append = TRUE, sep="\n")
  cat(print2, file = outfile, append = TRUE, sep="\n")
  cat(print0, file = outfile, append = TRUE, sep="\n")
  
  # x = number of marked elements in the selection
  # Q1
  x1_temp  = subset %>% filter(uplog10 > median_up)
  x1_temp2 = x1_temp %>% filter(downlog10 < median_down)
  x1 = nrow(x1_temp2)

  # Q2
  x2_temp = x1_temp %>% filter(downlog10 > median_down)
  x2 = nrow(x2_temp)
  
  # Q3
  x3_temp  = subset %>% filter(downlog10 > median_down)
  x3_temp2 = x3_temp %>% filter(uplog10 < median_up)
  x3 = nrow(x3_temp2)
  
  # Q4
  x4_temp  = subset %>% filter(uplog10 < median_up)
  x4_temp2 = x4_temp %>% filter(downlog10 < median_down)
  x4 = nrow(x4_temp2)

  # write to file
  print1 = paste(var, " genes Q_UL:\t", x1, sep = '')
  print2 = paste(var, " genes Q_UR:\t", x2, sep = '')
  print3 = paste(var, " genes Q_LR:\t", x3, sep = '')
  print4 = paste(var, " genes Q_LL:\t", x4, sep = '')
  cat(print1, file = outfile, append = TRUE, sep="\n")
  cat(print2, file = outfile, append = TRUE, sep="\n")
  cat(print3, file = outfile, append = TRUE, sep="\n")
  cat(print4, file = outfile, append = TRUE, sep="\n")
  cat(print0, file = outfile, append = TRUE, sep="\n")

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

  # Generate results table for hyper and chi
  print1 = paste(var, " p(hypergeometric) Q_UL:\t", hypergeometrics[1], sep = '')
  print2 = paste(var, " p(hypergeometric) Q_UR:\t", hypergeometrics[2], sep = '')
  print3 = paste(var, " p(hypergeometric) Q_LR:\t", hypergeometrics[3], sep = '')
  print4 = paste(var, " p(hypergeometric) Q_LL:\t", hypergeometrics[4], sep = '')
  cat(print1, file = outfile, append = TRUE, sep="\n")
  cat(print2, file = outfile, append = TRUE, sep="\n")
  cat(print3, file = outfile, append = TRUE, sep="\n")
  cat(print4, file = outfile, append = TRUE, sep="\n")
  cat(print0, file = outfile, append = TRUE, sep="\n")

  print1 = paste(var, " p(chi-squared) Q_UL:\t", chisquared[1], sep = '')
  print2 = paste(var, " p(chi-squared) Q_UR:\t", chisquared[2], sep = '')
  print3 = paste(var, " p(chi-squared) Q_LR:\t", chisquared[3], sep = '')
  print4 = paste(var, " p(chi-squared) Q_LL:\t", chisquared[4], sep = '')
  cat(print1, file = outfile, append = TRUE, sep="\n")
  cat(print2, file = outfile, append = TRUE, sep="\n")
  cat(print3, file = outfile, append = TRUE, sep="\n")
  cat(print4, file = outfile, append = TRUE, sep="\n")
}