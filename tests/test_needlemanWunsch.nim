import std/unittest
import bionim/needlemanWunsch
import bionim/utils
import bio_seq/io/fasta
import bio_seq/iupac_uint8

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
    

    test "Align":
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


    test "Build grid second single":
      seq1 = "WHY"
      seq2 = "W"
      grid = needlemanWunsch(seq1, seq2, indel_penal, match, mismatch)
      var test: Matrix 
      new(test)
      test[] = newSeq[seq[int]](seq2.len+1)
      test[0] = @[0, -1, -2, -3]
      test[1] = @[-1, 1, 0, -1]
      check(grid[] == test[])

    test "Build grid first empty":
      seq1 = ""
      seq2 = "WHY"
      grid = needlemanWunsch(seq1, seq2, indel_penal, match, mismatch)
      var test: Matrix 
      new(test)
      test[] = newSeq[seq[int]](seq2.len+1)
      test[0] = @[0 ]
      test[1] = @[-1]
      test[2] = @[-2]
      test[3] = @[-3 ]
      check(grid[] == test[])


    test "Build grid second empty":
      seq1 = "WHY"
      seq2 = ""
      grid = needlemanWunsch(seq1, seq2, indel_penal, match, mismatch)
      var test: Matrix 
      new(test)
      test[] = newSeq[seq[int]](seq2.len+1)
      test[0] = @[0, -1, -2, -3 ]
      check(grid[] == test[])

      #TODO Add more alignment tests
      #Test the NeedlemanWunschOptions type

proc needlemanWunschFile()=
  suite "Main test":
    setup:
      let indel_penal: int8 = -1
      let mismatch: int8 = -1
      let match: int8 = 1
    
    test "needleman wunsch from file":
      # this test is only here to check that the integration from bio_seq
      # and this algorithms works
      var seqs = parseFastaFile("tests/files/nm.fasta")

      var seq1 = seqs.seqs[0].sequence
      var seq2 = seqs.seqs[1].sequence

      
      var grid = needlemanWunsch(seq1, seq2, indel_penal, match, mismatch)   
      check(true)

needlemanWunschMain()
needlemanWunschFile()






