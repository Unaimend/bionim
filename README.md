# Bionim
[![test](https://github.com/Unaimend/bionim/actions/workflows/main_test.yml/badge.svg?branch=main)](https://github.com/Unaimend/bionim/actions/workflows/main_test.yml)

This package tries to provide a lot of the most useful data structures and algorithms needed in the different subfields of bioinformatics

The apis for most of the algorithms we be quite unstable until version 1.0.0 is reached, also this package will grrow quite a lot in the next few months

# Instalation
Requires Nim to be installed on your system. See https://nim-lang.org/

Installation with the nimble package manager is recommended:

`nimble install bionim`


# Tests
To run test just run `nimble test` in the bionim directory

# File Types
 * FASTA Alignments, those are fasta files with added constraint that all sequences must be of the same length
 * General FASTA Files
 * Reading FastQ files

Currently implemented algorithms and data structures
# Alignment
Algorithms and data structures usefull in the alignment of genomic and protein based sequences.
 ## Algorithms
 * [Needleman-Wunsch](https://en.wikipedia.org/wiki/Needleman-Wunsch_algorithm "Wikipedia page of the Needleman-Wunsch algorithm") Basic algorithm used for global alignment
 * [Smith-Waterman](https://en.wikipedia.org/wiki/Smith-Waterman_algorithm) Basic algorithm used for local alignment.

# Phylogenetics
 * Fitch's algorithm for constructing phylogenetic trees from aligned sequences.

# Assembly
 * DeBruijn Graph building and rendering
 * Paired-End DeBruijn grap building and rendering

# Documentation
All algorithms and data structures are documented in the projects wiki.
Documentation for `PhylogeNi` and `bio_seq` can be found on their respective git hub readme's

# Credits
In this section are all packages listed on which bionim depends and which were not created by me.
 * [PhylogeNi](https://github.com/kerrycobb/PhylogeNi) A library for working with phylogenetic trees
 * [bio_seq](https://github.com/kerrycobb/nim-dna-seq) A library for working with bioinformatic sequence types
 
# Contributions
If you want to contribute or think that a specific algorithm or data structure should be included in this package just create an issue.
Also I am always happy if someone suggest code improvements and stuff like that.

# Contact 
If you have any questions or just want to contact me regarding this package. Please write to
Unaimend+bionim@gmail.com
# Todo
## Parsing
 * FASTQ parsing
## Alignment
 * Calculation of all possible optimal alignment for Needleman-Wunsch
 * Usage of similiarity matrices in Needleman-Wunsch
     * Customizable gap penalties for Needleman-Wunsch
## Minimum skew problem
