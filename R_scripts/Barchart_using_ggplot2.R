#!/usr/bin/Rscript
args <- commandArgs(trailingOnly = F)  
scriptPath <- normalizePath(dirname(sub("^--file=", "", args[grep("^--file=", args)])))
R_functions <- paste(scriptPath, "/R_functions.R", sep="")
source(R_functions)
add_dependency("ggplot2")
add_dependency("gridExtra") # plot placement
add_dependency("gtable") # same
add_dependency("grid") # ggplot rendering 
add_dependency("reshape2") # data manipulation
add_dependency("optparse")

# Opening commands
option_list = list(
		   make_option(c("-d", "--dataframe"), default=NULL, help="Dataframe file name"),
		   make_option(c("-q", "--ymin"), default=0, help="Ymin (take into account stacked min value)"),
		   make_option(c("-y", "--ymax"), default=20, help="Ymax (take into account stacked max value)"),
		   make_option(c("-i", "--height"), default=5, help="Height (inches)"),
		   make_option(c("-w", "--width"), default=8, help="Width (inches)"),
		   make_option(c("-r", "--orientation"), default="v", help="Orientation of plot"),
		   make_option(c("-o", "--output"), default=NULL, help="Output PDF")
		   );
opt_parser = OptionParser(usage = "Usage: %prog -d <dataframe>
Optional: -r Orientation (v=verticle, h=horizontal) [v]
          -q Ymin [0]
          -y Ymax [20]
          -i Height (inches) [5]
          -w Width (inches) [8]
	  -o Output [opt$dataframe-barchart-stacked.pdf]\n
Notes: Dataframe should look like this:
File	kingdom	phylum	Species_name	secreted_genes_pc
GCA_020617735.1_ASM2061773v1	Fungi	Chytridiomycota	Batrachochytrium salamandrivorans	20.54304173
GCA_018290095.1_Pv_5.2	Metazoa	Arthropoda	Polypedilum vanderplanki	20.3209954
GCA_002217885.1_Fsol_1.0	N/A	Bacillariophyta	Fistulifera solaris	19.7080292
GCA_021223935.1_ASM2122393v1	Metazoa	Arthropoda	Drosophila pseudotakahashii	19.20489297", option_list=option_list);
opt = parse_args(opt_parser);
if(is.null(opt$dataframe)) {
	print_help(opt_parser)
	stop("Dataframe must be supplied (input file)", call.=FALSE)
}
if(!file.exists(opt$dataframe)) { stop("Error: dataframe does not appear to be valid", call.=FALSE) }
data = read.table(opt$dataframe, com='', sep='\t', header=T)

# Output
output <- paste(opt$dataframe, "-barchart-stacked.pdf", sep="")
if(!is.null(opt$output)) { output <- opt$output }
pdf(output, paper="special", height=opt$height, width=opt$width)

# Plot 
#(limits = order) , limits=data$secreted_genes_pc
tmpPlot <- ggplot(data, 
	   aes(x=reorder(go_term_and_desc, Enriched_in_this_many_genomes),
	       y=Enriched_in_this_many_genomes
	      ), 
	   color=kingdom) +
	   scale_y_continuous("Genomes with GO Terms significantly enriched in QUR",limits=c(opt$ymin,opt$ymax)) +
	   scale_x_discrete("GO Terms and Description") +
	   geom_bar(stat="identity", position="stack") 

# adding breaks line below, will give control of axis tick marks and numbers e.g. from 0 to 6, spacing every 1. 6 could be replaced with whatever the limit is to get all.
#scale_y_continuous("Total length (Mb)",limits=c(0,17),breaks=seq(0, 6, 1)) +

#Plot <- tmpPlot + theme(axis.text.x=element_text(angle = 45, hjust = 1),
Plot <- tmpPlot + theme(axis.text.x=element_text(angle = 90, hjust = 1),
	   axis.text.y=element_text(angle = 0, hjust=1))

if(opt$orientation == 'h') {
	Plot <- Plot + coord_flip()
}

Plot

dev.off()
