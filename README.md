# Bionim
[![test](https://github.com/Unaimend/bionim/actions/workflows/main_test.yml/badge.svg?branch=main)](https://github.com/Unaimend/bionim/actions/workflows/main_test.yml)

This package tries to provide a lot of the most useful data structures and algorithms need in the different subfields of bioinformatics

The apis for most of the algorithms we be quite unstable until version 1.0.0 is reached, also this package will grrow quite a lot in the next few months

# Instalation
Requires Nim to be installed on your system. See https://nim-lang.org/

Installation with the nimble package manager is recommended:

`nimble install bionim`


# Tests
To run test just run `nimble test` in the bionim directory



Currently implemented alogrithms and data structures
# Alignment
Algorithms and data structures usefull in the alignment of genomic and protein based sequences.
 ## Algorithms
 * [Needleman-Wunsch](https://en.wikipedia.org/wiki/Needleman-Wunsch_algorithm "Wikipedia page of the Needleman-Wunsch algorithm") Basic algorithm used for global alignment
 * [Smith-Waterman](https://en.wikipedia.org/wiki/Smith-Waterman_algorithm) Basic algorithm used for local alignment.

# Phylogenetics
 * Fitch's algorithm for constructing phylogenetic trees from aligned sequences.

# Assembly
 * DeBruijn Graph building

# Credits
In this section are all packages listed on which bionim depends and which were not created by me.
 * [PhylogeNi](https://github.com/kerrycobb/PhylogeNi) A library for working with phylogenetic trees
 
# Contributions
If you want to contribute or think that a specific algorithm or data structure should be included in this package just create an issue.
Also I am always happy if someone suggest code improvements and stuff like that.
# Todo
## Parsing
 * FASTA and FASTQ parsing
## Alignment
 * Calculation of all possible optimal alignment for Needleman-Wunsch
 * Usage of similiarity matrices in Needleman-Wunsch
     * Customizable gap penalties for Needleman-Wunsch
## Minimum skew problem
