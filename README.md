## 2speed_genomes

Intergenic distances and gene function have been associated with genome organisation, most strikingly in instances of “two-speed genomes”, a feature described in pathogens such as Batrachochytrium salamandrivorans (1) and Phytophthora infestans (2). Here, I develop a set of scripts useful for testing and visualising intergenic distances for any given genome assembly.

1.  T. Wacker, et al, Two-speed genome evolution drives pathogenicity in fungal pathogens of animals. Proc Natl Acad Sci U S A 120, e2212633120 (2023).
2.  S. Raffaele et al, Genome evolution following host jumps in the Irish potato famine pathogen lineage. Science 330, 1540–1543 (2010).

## Dependencies

* perl dependencies available from cpan/cpanm:
```bash
Getopt::Std
File::Basename
Data::Dumper
Getopt::Std
FindBin
Math::Matrix
```

* R dependencies (will try and download on first run):
```bash
ggplot2
gridExtra
gtable
grid
reshape2
optparse
dplyr
gghighlight
hash
stringr
```

## Approach to visualising signatures of 2-speed genomes

Genes are classified into quadrants that are defined as:

  Q1 (upper left) | Q2 (upper right)
  
  Q4 (lower left) | Q3 (lower right)

  
  1 = < median log10(downstream length) && > median log10(upstream length)

  2 = > median log10(downstream length) && > median log10(upstream length)

  3 = > median log10(downstream length) && < median log10(upstream length)

  4 = < median log10(downstream length) && < median log10(upstream length)

## Accessory scripts (Perl)

* Bootstrap_for_random_genes_chosen_for_HgT.pl
* GFF_and_genome_FASTA_to_gene_FASTA.pl
* HgT_p_values_to_R_dataframe.pl
* Join_all_stats_bootstraps.pl
* Markov_results_GFF_and_SigP_file.pl
* Markov_results_GFF_and_SigP_file_output_genes.pl
* RepeatMasker_output_parse.pl
* Repeat_folder_and_taxonomy_dataset_and_genome_length_to_dataframe3.pl
* Repeat_folder_and_taxonomy_dataset_to_dataframe3_vs_sigp4.pl

## Visualising scripts (R)

* Barchart_using_ggplot2.R - Generic barchart plotter
* Barchart_using_ggplot2_second_axis.R - Generic barchart plotter with 2nd axis
* Boxplot_using_ggplot2.R - Generic boxplot plotter
* Density_plot_3and5prime_HgT.R - Performs hypergeometric tests, comparing the numbers of genes in each of the quadrants.
* Density_plot_3and5prime_plot.R - Plots density plot of genes (either all, or a selected gene)
* Density_plot_3and5prime_plot_no_header-with-contours.R - Plots density plot of genes (either all, or a selected gene) with contours. No headers expected in the input file.
* Density_plot_3and5prime_plot_no_header.R - Plots density plot of genes (either all, or a selected gene). No headers expected in the input file.

## Steps for running the markov chain approach to identifying consecutive genes 

1. GFF-upstream-downstream-distances.pl (reads in GFF's and outputs tab files describing intergenic distances for every gene)
2. Intergenic_distances_to_median_quadrants.pl (outputs tab files of gene locations and quadrant information) 
3. Markov_chain_from_quadrant_file.pl (identifies significant numbers of consecutive genes of a given quadrant) - Cumulative_probability_of_consecutive_genes is the probability of finding that number of consecutive genes in that given quadrant.
