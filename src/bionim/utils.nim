import std/strformat
import std/[tables]

# get all kmers of length kmerLength in sequence
iterator kmers*(sequence: string, kmerLength: int): string =
  for i in 0 ..< sequence.len - (kmerLength - 1):
    yield sequence[i .. i + kmerLength - 1]

# Calculates all kmers(of length kmerLength) of a string and their corresponding amount
proc countKmers*(sequence: string, kmerLength: Positive): CountTable[string]=
  var table: CountTable[string]  
  #make this an iterator
  for i in kmers(sequence, kmerLength):
    table.inc(i)
  table

proc prefix*(s: string): string {.inline.} = s[0 .. s.len-2]
proc suffix*(s: string) : string {.inline.} = s[1 .. s.len-1]


