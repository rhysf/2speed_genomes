#!/usr/bin/Rscript
args <- commandArgs(trailingOnly = F)  
scriptPath <- normalizePath(dirname(sub("^--file=", "", args[grep("^--file=", args)])))
R_functions <- paste(scriptPath, "/R_functions.R", sep="")
source(R_functions)
add_dependency("ggplot2")
#add_dependency("gridExtra") # plot placement
#add_dependency("gtable") # same
#add_dependency("grid") # ggplot rendering 
#add_dependency("reshape2") # data manipulation
add_dependency("optparse")

# Opening commands
option_list = list(
		   make_option(c("-d", "--dataframe"), default=NULL, help="Dataframe file name"),
		   make_option(c("-q", "--ymin"), default=0, help="Ymin (take into account stacked min value)"),
		   make_option(c("-y", "--ymax"), default=1, help="Ymax (take into account stacked max value)"),
		   make_option(c("-i", "--height"), default=5, help="Height (inches)"),
		   make_option(c("-w", "--width"), default=8, help="Width (inches)"),
		   make_option(c("-r", "--orientation"), default="v", help="Orientation of plot"),
		   make_option(c("-o", "--output"), default=NULL, help="Output PDF")
		   );
opt_parser = OptionParser(usage = "Usage: %prog -d <dataframe>
Optional: -r Orientation (v=verticle, h=horizontal) [v]
          -q Ymin [0]
          -y Ymax [1]
          -i Height (inches) [5]
          -w Width (inches) [8]
          -o Output [opt$dataframe-boxplot.pdf]\n
Notes: Dataframe should look like this:
Quadrant	Value
QLL	0.340495766099381
QLL	0.0970764268128472
QLL	0.842806500852561", option_list=option_list);
opt = parse_args(opt_parser);
if(is.null(opt$dataframe)) {
	print_help(opt_parser)
	stop("Dataframe must be supplied (input file)", call.=FALSE)
}
if(!file.exists(opt$dataframe)) { stop("Error: dataframe does not appear to be valid", call.=FALSE) }
data = read.table(opt$dataframe, com='', sep='\t', header=T)

# Output
output <- paste(opt$dataframe, "-boxplot.pdf", sep="")
if(!is.null(opt$output)) { output <- opt$output }
pdf(output, paper="special", height=opt$height, width=opt$width)

# http://www.sthda.com/english/wiki/ggplot2-box-plot-quick-start-guide-r-software-and-data-visualization

# Plot 
# Change  automatically color by groups
bp <- ggplot(data, aes(x=Quadrant, y=Value, fill=Quadrant)) + 
  geom_boxplot()+
  labs(title="HgT p-values for Eukaryotic secreted genes",x="Quadrant", y = "p-value") + #+ scale_y_continuous(trans='log10')
  scale_y_log10(limits = c(2.22e-300, 1), breaks = c(1e-300, 1e-200, 1e-100, 1))

# Discrete colors
Plot <- bp + 
    scale_fill_brewer(palette="Dark2") + 
    theme_minimal() +
    geom_hline(yintercept=2.13e-6,linetype=2, color="#A9A9A9")


Plot

dev.off()
