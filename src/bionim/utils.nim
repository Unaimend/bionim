import strformat

type Matrix* = ref seq[seq[int]]

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
