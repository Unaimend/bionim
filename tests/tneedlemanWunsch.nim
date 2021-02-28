import std/unittest
import bionim/needlemanWunsch


proc needlemanWunschMain() =
  suite "Main test":
    setup:
      let seq1 = "WHAT"
      let seq2 = "WHY"

      let indel_penal: int8 = -1
      let gap_penal: int8 = -1
      let match: int8 = 1

      let grid = needlemanWunschGlobal(seq1, seq2, gap_penal, match, indel_penal)
    test "Build grid":
      var test: Matrix 
      new(test)
      test[] = newSeq[seq[int]](seq2.len+1)
      test[0] = @[0,-1,-2,-3, -4]
      test[1] = @[-1,1,0,-1, -2]
      test[2] = @[-2,0,2,1, 0]
      test[3] = @[-3,-1,1,1, 0]
      printGrid(test, "WHAT", "WHY")
      check(grid[] == test[])
    test "Align":
      let alignment = calculateAlignment(grid, seq1, seq2, gap_penal, match, indel_penal)
      check(alignment[0] == "WHAT")
      check(alignment[1] == "WH-Y")


needlemanWunschMain()






