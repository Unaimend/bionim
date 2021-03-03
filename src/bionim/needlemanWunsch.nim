##Module description

import algorithm
import utils
const DEBUG = true
## Type alias for working with 2D-Matrix
#TODO Should probably renamed to Matrix2D

type NeedlemanWunschConfig* = ref object
  ## Type which reprensents all need information for the Needleman-Wunsch algorithm
  sequence1*: string
    ## First sequence used in the alignment
  sequence2*: string
    ## Second sequence used in the alignment
  indel_penal*: int8
    ## The penaltie which for insertions and deletions
  match*: int8
    ## The reward when a match happens
  mismatch*: int8
    ## The penaltie for a mismatch

proc needlemanWunsch*(sequence1: string, sequence2: string, indel_penal: int8, match: int8, mismatch: int8): Matrix =
  ## This procedure build and returns the 2D-Arrays used in the algorithm. 
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
    grid[y][0] = indel_penal*y
    for x in 0 ..< seq1_len:
      if y==0:
        #initialize the 0 row
        grid[y][x] = indel_penal*x
      elif x > 0 and y > 0:
        #calculate score values
        let s = if sequence1[x-1] == sequence2[y-1]: match else:mismatch 
        grid[y][x] = max(grid[y][x-1] + indel_penal, 
                     max(grid[y-1][x] + indel_penal, 
                         grid[y-1][x-1]+s))

  when DEBUG:
    echo "Length of first sequence ", seq1_len, " Length of first array index ", grid[].len
    echo "Length of second sequence ", seq2_len, " Length of second array index ", grid[0].len
  result = grid

proc needlemanWunsch*(options: NeedlemanWunschConfig): Matrix=
  # Helper function for use with the NeedlemanWunschConfig type
  needlemanWunsch(options.sequence1, options.sequence2, options.indel_penal, options.match, options.mismatch)


proc calculateAlignment*(grid: Matrix, sequence1: string, sequence2: string, indel_penal: int8, match: int8, mismatch: int8): (string, string) = ## Calculates the actual allignment, returs a Tuple of (string, string) in which t[0] represents the first aligned sequence and t[1] the second aligned sequence.
  var x = sequence1.len
  var y = sequence2.len

  var alignA: string
  var alignB: string
  #TODO Add comment why i do this weird stuff
  proc lazy(): int8 =
    if sequence1[x-1] == sequence2[y-1]: match else: mismatch
     
  while x > 0 or y > 0:
    var current = grid[y][x]
    #find out if there was a match or an indel in sequence1 or an indel in sequence2
    #Check if we can still go to the upper left, and if we came from the upperleft
    if x > 0 and y > 0 and current == grid[y-1][x-1] + lazy():
        #Match
        alignA.add(sequence1[x-1])
        alignB.add(sequence2[y-1])
        x = x - 1 
        y = y - 1
    elif y > 0 and current == grid[y-1][x] + indel_penal:
      alignA.add("-")
      alignB.add(sequence2[y-1])
      y = y - 1
    elif x > 0 and current == grid[y][x-1] + indel_penal:
      alignA.add(sequence1[x-1])
      alignB.add("-")
      x = x - 1
    else:
      echo "THIS SHOULD NEVER HAPPEN"
  #Since we got those sequences by backtracking we have to reverse them
  alignA.reverse
  alignB.reverse
  result = (alignA, alignB)

proc calculateAlignment*(grid: Matrix, options: NeedlemanWunschConfig) : (string, string) =
  # Helper function for use with the NeedlemanWunschConfig type
  calculateAlignment(grid, options.sequence1, options.sequence2, options.indel_penal, options.match, options.mismatch) 

when isMainModule:
  let a = needlemanWunsch("G", "WHY", -1, 1, -1)
  printGrid(a, "G", "WHY")
  let b = calculateAlignment(a, "G", "WHY", -1,1,-1)
  echo b[0]
  echo b[1]

