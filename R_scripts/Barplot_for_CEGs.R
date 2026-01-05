#!/usr/bin/Rscript

#source("/home/rhys/Code/bioinformatics_r/R_functions.R")

# CEGannoteP.R
# Makes barplots of CEG analyses from CEGannoteP.pl. Multiple gene sets may be compared.
# CEGs are widely conserved single copy eukaryotic orthologs.
# These provide a sampling allowing the coverage and fragmentation
# of a genome to be assessed. 
# The input for CEGannoteP.R (which comes from CEGannoteP.pl) looks like this:

#ngenes	ncegs	nuniqcegs	ngenie	below50	n50to60	n60to70	n70to80	n80to90	n90to100	nmissing
#Annc_alge_insect_USDA_JJB_V2_PASA_UPDATES_1_FIXED_COMBO_1_UTRTRIM_1	3598	209	138	120	4	2	12	18	35	67	110
#Annc_alge_PRA109_V2_PASA_UPDATES_1_MOD3_1_UTRTRIM_1	3232	156	81	44	27	2	8	15	10	19	167

##### NEW
#ncegs   below70 n70to80 n80to90 n90to100        nmissing
#GCA_013398475.1 216     22      11      18      187     10
#GCA_020499815.1 244     2       5       15      224     2
#GCA_019434415.1 235     11      6       18      211     2
#GCA_022606255.1 244     3       4       13      227     1


# There is one less column name than columns (the gene set column is not named)
# This tells R to use the gene set names to name table rows.

# Read this table from disk:

# Opening commands
args<-commandArgs(trailingOnly = TRUE)
if(length(args) != 1) { 
	cat("Usage: Rscript script_name.R <CEGSummary.txt>\n") 
	stop("Notes: If empty pdf, may need to run on cluster\n")
}
CEGSum.tbl <- read.table(args[1], com='', sep='\t', header=T)
#CEGSum.tbl <- read.table("CEGSummary.txt",header = TRUE)

# Output
output <- paste(args[1], ".pdf", sep="")

# Some kind of processing of data
cat("Plotting CEGomes\n")
cat(output)
cat("\n")
trans <- t(CEGSum.tbl)


# Can i make columns of maybe 50?
Grp <- 0:(length(trans) - 1) %/% 50
names(Grp) <- trans

Grp[0]


#below70 <-  0
#for ( i in 1:dim(trans)[2] ) {
#   below70[i] <- sum(trans[6:8,i])
#}
#print (below70)
#trans <- rbind(trans,below70)

# Number of lines
HEIGHT<-nrow(CEGSum.tbl)

# Specify output
pdf(output, onefile=TRUE, width=8, height=HEIGHT)
#pdf(output, onefile=TRUE, width=12, height=12)
#pdf(output)

# margin (bottom, left, top and right margins )
#par(mai=c(0.5,1.5,1.5,1.5))
#par(xpd=TRUE)

# Rainbow plot of CEG content of input genomes
#par(mai=c(4,1,0.25,1.5))
legend <- c("missing","below 70%","70-80%","80-90%","90-100%")
colors=c("white","lightcoral","lightgreen","lightblue","mediumpurple1")
revcolors <- rev(colors)
#rows2plot <- c(11:6,12)
#rows2plot <- c(6,2,3,4,5)
rows2plot <- c(5,4,3,2,6)

#rows2plot
#trans[rows2plot,]
#Grp1[rows2plot,]

barplot(trans[rows2plot,],
	las = 2,
	xlab="Number CEGs",
	col=rev(colors),
	cex.names=c(0.8),axes=FALSE,
	font=3,
        horiz = TRUE)

#usr <- par("usr")
#par(usr=c(usr[1:2], 0, 250))

# axis 3 = above
axis(3,at=c(0,50,100,150,200,250), outer=FALSE)
legend("topright",
       legend=legend,
       ncol=1,
       fill=colors,
       cex=.5,
       title="CEG coverage",
       inset = c(0,0))
dev.off()
#proc.time()
