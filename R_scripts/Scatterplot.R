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
add_dependency("ggrepel") # text labels


### r.farrer@exeter.ac.uk

# Opening commands
option_list = list(
		   make_option(c("-d", "--dataframe"), default=NULL, help="Dataframe file name"),
		   make_option(c("-a", "--xvar"), default=NULL, help="x variable"),
		   make_option(c("-b", "--yvar"), default=NULL, help="y variable"),
		   make_option(c("-c", "--title"), default="title", help="Title"),
		   make_option(c("-t", "--type"), default="b", help="Scatterplot type (b=based, g=grouped by variable)"),
		   make_option(c("-s", "--size_point"), default=2, help="Size of points"),
		   make_option(c("-g", "--groupvar"), default=NULL, help="if type=g, group by what variable/header"),
		   make_option(c("-l", "--labelsvar"), default=NULL, help="if type=g, add labels from a variable/header"),
		   make_option(c("-r", "--regression"), default="n", help="add a line of best fit (n/y)"),
		   make_option(c("-z", "--r2"), default="y", help="add r2 value to regression line (n/y)"),
		   make_option(c("-i", "--height"), default=5, help="Height (inches)"),
		   make_option(c("-w", "--width"), default=8, help="Width (inches)"),
		   make_option(c("-o", "--output"), default=NULL, help="Output PDF")
		   );
opt_parser = OptionParser(usage = "Usage: %prog -d <dataframe> -a <x variable> -b <y variable>
Optional: -t Type (b=based, g=grouped by variable) [b]
          -s Size of point [2]
          -g Group variable (one of the columns) []
          -l Labels variable (one of the columns) if grouped
          -r Add regression line (y/n) [n]
          -z Add r2 value to regression line (y/n) [y]
          -i Height (inches) [5]
          -w Width (inches) [8]
          -c Title [title]
          -o Output [opt$scatterplot1.pdf]

Notes: Dataframe should look like mtcars e.g.:
                  wt  mpg cyl qsec
 Mazda RX4      2.62 21.0   6 16.5
 Mazda RX4 Wag  2.88 21.0   6 17.0
 Datsun 710     2.32 22.8   4 18.6
 Hornet 4 Drive 3.21 21.4   6 19.4", option_list=option_list);
opt = parse_args(opt_parser);
if(is.null(opt$dataframe)) {
	print_help(opt_parser)
	stop("Dataframe must be supplied (input file)", call.=FALSE)
}
options(warn=1)
if(!file.exists(opt$dataframe)) { stop("Error: dataframe does not appear to be valid", call.=FALSE) }
data = read.table(opt$dataframe, com='', sep='\t', header=T)
#attach(data)
#names(data)

# Get the x and y values wanted
xvar1 <- opt$xvar
yvar1 <- opt$yvar

# Get the data for x and y-values
xvar2 <- data[[xvar1]]
yvar2 <- data[[yvar1]]

# Correlation and regression stats:
# Remove missing values
keep <- complete.cases(xvar2, yvar2)
x <- xvar2[keep]
y <- yvar2[keep]
n <- length(x)

cat("\n Pearson correlation:\n")
ct <- cor.test(x, y, method = "pearson")
print(ct)

cat("\n Linear model (y ~ x):\n")
fit <- lm(y ~ x)
s <- summary(fit)
cat("n =", n, "\n")
cat("R-squared =", s$r.squared, "\n")
cat("Slope p-value =", s$coefficients[2, 4], "\n")

# Output
output <- paste(opt$dataframe, "-scatterplot1.pdf", sep="")
if(!is.null(opt$output)) { output <- opt$output }
pdf(output, paper="special", height=opt$height, width=opt$width)

# https://www.datanovia.com/en/lessons/ggplot-scatter-plot/

# size of points
point_size <- opt$size_point

# Initiate a ggplot
b <- ggplot(data, aes(x = xvar2, y = yvar2)) +
	 labs(title = opt$title, x = xvar1, y = yvar1) +
	 theme(plot.title = element_text(hjust = 0.5)) +
	 geom_point(size = point_size)

# Add regression line
if(opt$regression == 'y') {
	b <- b + geom_smooth(method = "lm", se = FALSE, color = "black")

	if(opt$r2 == 'y') {
		add_dependency("ggpubr") # for r2 value

		b <- b + stat_cor(
			aes(label = after_stat(rr.label)), # R² value
				label.x.npc = "center", 
				label.y.npc = "center",
				size = 5
            )
	}
}

# basic plot
if(opt$type == 'b') {
	# Basic scatter plot
	
	# no labels
	if(is.null(opt$labelsvar)) {
		Plot <- b + geom_point()
	} else {
		varlab <- opt$labelsvar
		varlab2 <- data[[varlab]]

		Plot <- b + geom_point() +
		geom_label_repel(aes(label = varlab2,  color = varnew2), size = 3)
	}
	# Change color, shape and size
	# b + geom_point(color = "#00AFBB", size = 2, shape = 23)

	# plot
	#b
}

# Grouped
if(opt$type == 'g') {
	# grouping factor
	varnew <- opt$groupvar
	varnew2 <- data[[varnew]]
	varnew2 <- as.factor(varnew2)
	
	var_count <- length(unique(varnew2))
	colors <- c("red", "green", "blue", "purple", "orange", "yellow", "pink", "brown")

	# Subset colors to match number of unique strings
	color_subset <- colors[1:var_count]

	# no labels
	if(is.null(opt$labelsvar)) {

		print("Making grouped scatterplot without labels...")

		# Make all points the same size (required scale_shape_identity)
		Plot <- b + geom_point(aes(shape = 19, color = varnew2)) +
			scale_color_manual(values = color_subset) +
			#scale_color_manual(values = c("#00AFBB", "#E7B800")) +
			scale_shape_identity() +
			labs(color = "Key")
	}
	# labels
	else {

		print("Making grouped scatterplot with labels...")

		varlab <- opt$labelsvar
		varlab2 <- data[[varlab]]


		# Change point shapes and colours by grouping factor
		#scale_color_manual(values = c("#00AFBB", "#E7B800", "#FC4E07"))
		# without rectangle = geom_text_repel(aes(label = varlab2,  color = varnew2), size = 3)+
		Plot <- b + geom_point(aes(shape = varnew2, color = varnew2)) +
		geom_label_repel(aes(label = varlab2,  color = varnew2), size = 3)+
		scale_color_manual(values = c("#00AFBB", "#E7B800"))
	}

}

# Plot
Plot

dev.off()
