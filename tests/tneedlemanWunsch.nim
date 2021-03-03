import std/unittest
import bionim/needlemanWunsch
import bionim/utils

proc needlemanWunschMain() =
  suite "Main test":
    setup:
      var seq1 = "WHAT"
      var seq2 = "WHY"

      let indel_penal: int8 = -1
      let mismatch: int8 = -1
      let match: int8 = 1

      var grid = needlemanWunsch(seq1, seq2, indel_penal, match, mismatch)
    test "Build grid":
      var test: Matrix 
      new(test)
      test[] = newSeq[seq[int]](seq2.len+1)
      test[0] = @[0,-1,-2,-3, -4]
      test[1] = @[-1,1,0,-1, -2]
      test[2] = @[-2,0,2,1, 0]
      test[3] = @[-3,-1,1,1, 0]
      check(grid[] == test[])
    

    test "Align normale":
      let alignment = calculateAlignment(grid, seq1, seq2, indel_penal, match, mismatch)
      check(alignment[0] == "WHAT")
      check(alignment[1] == "WH-Y")

    test "Build grid first single":
      seq1 = "W"
      seq2 = "WHY"
      grid = needlemanWunsch(seq1, seq2, indel_penal, match, mismatch)
      var test: Matrix 
      new(test)
      test[] = newSeq[seq[int]](seq2.len+1)
      test[0] = @[0, -1]
      test[1] = @[-1, 1]
      test[2] = @[-2, 0]
      test[3] = @[-3, -1]
      check(grid[] == test[])


needlemanWunschMain()






