import utils
const DEBUG = true

type Matrix* = ref seq[seq[int]]

type SmithWatermanConfig* = ref object
  sequence1*: string
  sequence2*: string
  indel_penal*: int8
  match*: int8
  gap_penal*: int8

proc smithWaterman*(sequence1: string, sequence2: string, gap_penal: int8, match: int8, mismatch: int8 ): Matrix =
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
        echo "TEST"
        let s = if sequence1[x-1] == sequence2[y-1]: match else: mismatch
        grid[y][x] = max(grid[y][x-1] + gap_penal, 
                        max(grid[y-1][x] + gap_penal, 
                         max(grid[y-1][x-1]+s, 0)))
                         

  when DEBUG:
    echo "Length of first sequence ", seq1_len, " Length of first array index ", grid[].len
    echo "Length of second sequence ", seq2_len, " Length of second array index ", grid[0].len
  result = grid



when isMainModule:
  let a = smithWaterman("GGCTCAATCA", "ACCTAAGG", -2, 2, -1)
  printGrid(a, "GGCTCAATCA", "ACCTAAGG")
