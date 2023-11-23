### r.farrer@exeter.ac.uk

# Check, install, or add libraries
#lib_location = '/home/unix/rfarrer/R/x86_64-unknown-linux-gnu-library/3.1'
#lib_location = '/Library/Frameworks/R.framework/Resources/library'

find_os <- function() {
	lib_location = ""
	OS <- Sys.info()['sysname']
	if(OS == "Darwin") {
		lib_location = '/Library/Frameworks/R.framework/Resources/library'
	}
	else if(OS == "Linux") {
		#lib_location = '/usr/lib/R/library'
		lib_location = '~/R/x86_64-pc-linux-gnu-library/4.1'
	}
	else {
		stop("ERROR: R_functions.R -> find_os: OS cannot be determined", call.=FALSE)
	}
	return(lib_location)
	#lib_location = .libPaths()
	#return(lib_location[1])
}

add_dependency <- function(library_wanted) {
	if (!require(library_wanted, character.only = TRUE)) {
		message("installing ", library_wanted, " ...")

		lib_location <- find_os()
		install.packages(library_wanted, lib_location, repos = 'http://cran.us.r-project.org', dependencies = TRUE)
		library(library_wanted, character.only = TRUE)
	}
}

add_dependency_Bioconductor <- function(library_wanted) {
	if (!require(library_wanted, character.only = TRUE)) {

		# Get BiocManager if needed
		if (!requireNamespace("BiocManager", quietly = TRUE)) {
		    install.packages("BiocManager")
		}

		message("installing with BiocManager: ", library_wanted, " ...")
		BiocManager::install(library_wanted)
	}
}

# Calculate and collect features of a y-axis in relation to another y-axis
CalcFudgeAxis = function(y1, y2=y1) {

	# x gets mapped to range of ylim2
	Cast2To1 = function(x) ((ylim1[2]-ylim1[1])/(ylim2[2]-ylim2[1])*x) 
	ylim1 <- c(min(y1),max(y1))
	ylim2 <- c(min(y2),max(y2))    
	yf <- Cast2To1(y2)
	labelsyf <- pretty(y2)  

	return(list(
		yf=yf,
		labels=labelsyf,
		breaks=Cast2To1(labelsyf)
	))
}

# Plot second y-axis
# http://stackoverflow.com/questions/3099219/plot-with-2-y-axes-one-y-axis-on-the-left-and-another-y-axis-on-the-right
# based on: https://rpubs.com/kohske/dual_axis_in_ggplot2
PlotWithFudgeAxis = function( plot1, FudgeAxis) {
	plot2 <- plot1 + with(FudgeAxis, scale_y_continuous( breaks=breaks, labels=labels))

	#extract gtable
	g1<-ggplot_gtable(ggplot_build(plot1))
	g2<-ggplot_gtable(ggplot_build(plot2))

	#overlap the panel of the 2nd plot on that of the 1st plot
	pp<-c(subset(g1$layout, name=="panel", se=t:r))
	g<-gtable_add_grob(g1, g2$grobs[[which(g2$layout$name=="panel")]], pp$t, pp$l, pp$b,pp$l)

	ia <- which(g2$layout$name == "axis-l")
	ga <- g2$grobs[[ia]]
	ax <- ga$children[[2]]
	ax$widths <- rev(ax$widths)
	ax$grobs <- rev(ax$grobs)
	ax$grobs[[1]]$x <- ax$grobs[[1]]$x - unit(1, "npc") + unit(0.15, "cm")
	g <- gtable_add_cols(g, g2$widths[g2$layout[ia, ]$l], length(g$widths) - 1)
	g <- gtable_add_grob(g, ax, pp$t, length(g$widths) - 1, pp$b)

	#grid.draw(g)
	return(g)
}
