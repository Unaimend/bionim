## Implements the Smiah-Waterman algorithm used for local alignment
import utils
#TODO Find a way to do this better
const DEBUG = false

type SmithWatermanConfig* = ref object
  ## Type which reprensents all need information for the Smith-Waterman algorithm
  sequence1*: string
    ## First sequence used in the alignment
  sequence2*: string
    ## Second sequence used in the alignment
  indel_penal*: int8
    ## The penaltie which for insertions and deletions aka a gap
  match*: int8
    ## The reward when a match happens
  mismatch*: int8
    ## The penaltie for a mismatch

proc smithWaterman*(sequence1: string, sequence2: string, indel_penal: int8, match: int8, mismatch: int8 ): Matrix =
  ## This procedure builds and returns the 2D-Arrays used in the later steps of the algorithm. 
  let seq1_len = sequence1.len+1
  let seq2_len = sequence2.len+1

  #Create 2D-Array for the algorithm
  var grid: Matrix
  new(grid)
  grid[] = newSeq[seq[int]](seq2_len)
  for y in 0 ..< seq2_len:
    #Add sequences of length sequence1 for each character of sequence2
    grid[y] = newSeq[int](seq1_len)
    #initilize the 0 colum
    grid[y][0] = 0
    for x in 0 ..< seq1_len:
      if y==0:
        #initialize the 0 row
        grid[y][x] = 0
      elif x > 0 and y > 0:
        #calculate score values
        let s = if sequence1[x-1] == sequence2[y-1]: match else: mismatch
        grid[y][x] = max(grid[y][x-1] + indel_penal, 
                        max(grid[y-1][x] + indel_penal, 
                         max(grid[y-1][x-1]+s, 0)))
  result = grid

proc smithWaterman(options: SmithWatermanConfig): Matrix=
  smithWaterman(options.sequence1, options.sequence2, options.indel_penal, options.match, options.mismatch)




when isMainModule:
  let a = smithWaterman("GGCTCAATCA", "ACCTAAGG", -2, 2, -1)
  printGrid(a, "GGCTCAATCA", "ACCTAAGG")
