# 2speed_genomes
Code for identifying signatures of 2 speed genomes

# Steps

1. GFF-upstream-downstream-distances.pl (outputs tab file of intergenic distances)
2. Intergenic_distances_to_median_quadrants.pl (outputs tab of gene locations and quadrant) 
3. Markov_chain_from_quadrant_file.pl (identifies significant numbers of consecutive genes of a given quadrant)

Quadrants are defined as:

Q1 (upper left) | Q2 (upper right)
Q4 (lower left) | Q3 (lower right)

1 = < median log10(downstream length) && > median log10(upstream length)
2 = > median log10(downstream length) && > median log10(upstream length)
3 = > median log10(downstream length) && < median log10(upstream length)
4 = < median log10(downstream length) && < median log10(upstream length)
