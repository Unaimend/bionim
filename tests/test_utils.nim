import unittest
import bionim/utils
import tables


proc testCount()=
  suite "Test counting proc":
    test "test":
      let c: CountTable[string] = countKmers("ATGACTAATGTACCTGAATGT", 4)
  
      
        


testCount()
