import std/strformat
import std/[tables]


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

iterator kmers*(sequence: string, kmerLength: int): string =
  for i in 0 ..< sequence.len - (kmerLength - 1):
    yield sequence[i .. i + kmerLength - 1]


proc countKmers*(sequence: string, kmerLength: int): CountTable[string]=
  var table: CountTable[string]  
  #make this an iterator
  for i in kmers(sequence, kmerLength):
    table.inc(i)
  table

proc countNucleotides*(sequence: string): CountTableRef[char]=
  var table:  CountTableRef[char] = newCountTable(sequence)
  table

proc prefix*(s: string): string {.inline.} = s[0 .. s.len-2]
proc suffix*(s: string) : string {.inline.} = s[1 .. s.len-1]

