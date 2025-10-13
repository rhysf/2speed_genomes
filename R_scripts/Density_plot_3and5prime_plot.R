#!/usr/bin/Rscript
args <- commandArgs(trailingOnly = F)  
scriptPath <- normalizePath(dirname(sub("^--file=", "", args[grep("^--file=", args)])))
R_functions <- paste(scriptPath, "/R_functions.R", sep="")
source(R_functions)
add_dependency("ggplot2")
add_dependency("dplyr") 
add_dependency("gghighlight")
add_dependency("optparse")

### r.farrer@exeter.ac.uk 

# Opening commands
option_list = list(
    make_option(c("-d", "--dataframe"), default=NULL, help="Dataframe file name"),
    make_option(c("-o", "--output"), default=NULL, help="Output PDF"),
    make_option(c("-l", "--limit"), default=30, help="Max number of genes per square"),
    make_option(c("-i", "--height"), default=8, help="Height of pdf (inches)"),
    make_option(c("-w", "--width"), default=12, help="Width of pdf (inches)")
);
opt_parser = OptionParser(usage = "Usage: %prog -d <dataframe>
Optional: -o Output [opt$dataframe-densityplot_recolored.pdf]
          -i Height of pdf [8]
          -w Width of pdf [12]

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
data = read.table(opt$dataframe, com='', sep='\t', header=T)

output <- paste(opt$dataframe, "-densityplot_recolored.pdf", sep="")
if(!is.null(opt$output)) { output <- opt$output }
pdf(output, paper="special", height=opt$height, width=opt$width)

#calculate median for 5' and for 3' and then filter data for above and below (median because skewed data)
median_up = median(data$uplog10,na.rm=TRUE)
median_low= median(data$downlog10, na.rm=TRUE)

#Density Plots-----------------------------------------------------------------------------------------
# calculate selection k and number of genes -------------------------------------------------------------------------

# IMPORTANT: quadrants as follows: 
# ------------------------------
#|  Q1          |           Q2 |
#|              |              |
#---------median uplog---------|
#|  Q4          |           Q3 |
#|              |              |
#|--------------|--------------|
#               |<- median downlog
#Q1 is synonymous to q_ul, Q2 to q_ur, Q3 to q_lr and Q4 to q_ll

FIRS=data

#select for log plot and rename
# See the !!; they mean “hey R, remember the expression I stored recently? Now take it, and ‘unquote’ it, that is, just run it!”. The double exclamation mark is just syntactic sugar for that phrase.
#table = FIRS %>% select(uplog10,downlog10,!!(var),feature)
table = FIRS %>% select(uplog10,downlog10,desc,feature)

#plot with different bins
plot=ggplot(table, aes(x=downlog10, y=uplog10)) +
  labs(x="3' log10 FIR (bp)", y="5' log10 FIR (bp)") +
  #bin for heatmap
  geom_bin2d(bins = 60) +
  #color scheme
  scale_fill_gradient2(low="dodgerblue3",mid="darkorange2", high="red3",midpoint=15, limits=c(0,opt$limit)) +
  #format text
  theme(axis.text.x = element_text( color="black",size=16)) +
  theme(axis.text.y = element_text( color="black",size=16)) +
  theme(axis.title=element_text(size=16)) +
  labs(fill="Number of x")
  
# add quadrants
plot = plot + geom_vline(data = FIRS, aes(xintercept = median_low), linetype="dashed", linewidth=1, color = "dodgerblue4")
plot = plot + geom_hline(data = FIRS, aes(yintercept = median_up), linetype="dashed", linewidth=1, color = "dodgerblue4")
  
# Higlight genes of interest 
#plot2 = plot + gghighlight(FIRS$UQ(var) == !!(var))
plot2 = plot + gghighlight(table$desc == 'cat1,')

# save each plot
ggsave(output, height=opt$height, width=opt$width, units='in', dpi=600, pointsize=16)
