# Package

version       = "0.0.4"
author        = "Thomas Dost"
description   = "A collection of useful algorithms and data structures for bioinformatics"
license       = "MIT"
srcDir        = "src"
binDir        = "bin"
installExt    = @["nim"]


# Dependencies

requires "nim >= 1.4.4", "phylogeni >= 0.0.2", "bio_seq >= 0.0.5"
