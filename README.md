# 2speed_genomes
Code for identifying signatures of 2 speed genomes

# Steps

1. GFF-upstream-downstream-distances.pl (outputs tab file of intergenic distances)
2. Intergenic_distances_to_median_quadrants.pl (outputs tab of gene locations and quadrant) 
3. Markov_chain_from_quadrant_file.pl (identifies significant numbers of consecutive genes of a given quadrant)

Quadrants are defined as such (ASCII plot made by TW):
  \--------------------------------

  \|   Q1          |           Q2 |

  \|               |              |
  \----------median uplog---------|
  \|   Q4          |           Q3 |
  \|               |              |
  \|---------------|--------------|
  \                |<- median downlog
