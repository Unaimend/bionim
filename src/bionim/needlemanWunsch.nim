import strformat
import algorithm
const DEBUG = true

type Matrix* = ref seq[seq[int]]

type NeedlemanWunschConfig* = ref object
  sequence1*: string
  sequence2*: string
  match*: int8
  gap_penal*: int8
  indel_penal*: int8

proc needlemanWunsch*(sequence1: string, sequence2: string, gap_penal: int8, match: int8, indel_penal: int8): Matrix =
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
    grid[y][0] = gap_penal*y
    for x in 0 ..< seq1_len:
      if y==0:
        #initialize the 0 row
        grid[y][x] = gap_penal*x
      elif x > 0 and y > 0:
        #calculate score values
        let s = if sequence1[x-1] == sequence2[y-1]: match else: indel_penal
        grid[y][x] = max(grid[y][x-1] + gap_penal, 
                     max(grid[y-1][x] + gap_penal, 
                         grid[y-1][x-1]+s))

  when DEBUG:
    echo "Length of first sequence ", seq1_len, " Length of first array index ", grid[].len
    echo "Length of second sequence ", seq2_len, " Length of second array index ", grid[0].len
  result = grid

proc needlemanWunsch*(options: NeedlemanWunschConfig): Matrix=
  needlemanWunsch(options.sequence1, options.sequence2, options.gap_penal, options.match, options.indel_penal)


proc calculateAlignment*(grid: Matrix, sequence1: string, sequence2: string, gap_penal: int8, match: int8, indel_penal: int8): (string, string) = 
  var x = sequence1.len
  var y = sequence2.len

  var alignA: string
  var alignB: string
  while x > 0 or y > 0:
    var current = grid[y][x]
    #find out if there was a match or an indel in sequence1 or an indel in sequence2
    let s = if sequence1[x-1] == sequence2[y-1]: match else: indel_penal
    #Check if we can still go to the upper left, and if we came from the upperleft
    if x > 0 and y > 0 and current == grid[y-1][x-1] + s:
      #Match
      alignA.add(sequence1[x-1])
      alignB.add(sequence2[y-1])
      x = x - 1 
      y = y - 1
    elif y > 0 and current == grid[y-1][x] + gap_penal:
      alignA.add("-")
      alignB.add(sequence2[y-1])
      y = y - 1
    elif x > 0 and current == grid[y][x-1] + gap_penal:
      alignA.add(sequence1[x-1])
      alignB.add("-")
      x = x - 1
  alignA.reverse
  alignB.reverse
  result = (alignA, alignB)

proc calculateAlignment*(grid: Matrix, options: NeedlemanWunschConfig) : (string, string) =
  calculateAlignment(grid, options.sequence1, options.sequence2, options.gap_penal, options.match, options.indel_penal) 

proc printGrid*(grid: Matrix, sequence1: string, sequence2: string): void = 
  stdout.write "    -"
  for i in sequence1:
    let s: string = $i
    #TODO Make this not static
    stdout.write fmt"{s:>4}"
  echo ""
  let altSeq = "-" & sequence2
  for i in 0 ..< grid[].len:
    stdout.write altSeq[i]
    for j in 0 ..< grid[i].len:
      let s: string = $grid[i][j]
      stdout.write fmt"{s:>4}"
    echo ""

when isMainModule:
  let a = needlemanWunsch("WHAT", "WHY", -1, 1, -1)
  let b = calculateAlignment(a, "WHAT", "WHY", -1,1,-1)
  echo b[0]
  echo b[1]

  printGrid(a, "WHAT", "WHY")
