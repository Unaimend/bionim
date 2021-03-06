{.experimental: "views".}
{.experimental: "strictFuncs".}
import std/strformat
import std/[tables, hashes]


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

type ID = int64
type
  StrView = object
    start, finish: int
    base: ptr string

type DeBruijnGraph* = ref object
  source*: string
  kmers*: seq[StrView]
  edgesOut*: seq[seq[ID]]
  edgesIn*: seq[seq[ID]]
  kMerMap*: Table[StrView, ID]

proc toStrView(str: var string, start, finish: int): StrView =
  StrView(base: addr str, start: start, finish: finish)

#@Discord haxscramper, helped me fixing the view stuff i.e implemente the StrView type since I could get it to work
#with the nim buil in toOpenArray

iterator intRange(start, finish: int): int =
  if start <= finish:
    for idx in start .. finish:
      yield idx

  else:
    for idx in countdown(start, finish):
      yield idx

iterator items(view: StrView): char =
 for ch in intRange(max(0, view.start), max(0, view.finish)):
   yield view.base[][ch]

proc `$`(view: StrView): string =
  for ch in view:
    result.add ch
proc len(v1: StrView): int = abs(v1.finish - v1.start)
proc `[]`(v: StrView, idx: int): char = v.base[][v.start + idx]


proc `==`*(v1, v2: StrView): bool =
  if v1.len != v2.len:
    return false

  for i in 0 ..< len(v1):
    if v1[i] != v2[i]:
      return false

  return true

proc hash*(view: StrView): Hash =
  var h: Hash = 0
  for ch in view:
    h = h !& hash(ch)

  result = !$h

proc findOrCreateKmer(graph: var DeBruijnGraph, kMer: StrView): ID=
  var index: ID
  #TODO tab.mgetOrPut(key, @[]).add(val) can i do this in one lookup
  if not (kMer in graph.kMerMap):
    #echo "CALLED '", kMer, "'"
    
    #for x,y in graph.kMerMap:
      #echo "  ", x.start..x.finish, ": '", x, "' -> ", graph.kMerMap[x], " hash=", hash(x)

    #echo "adding element with hash", hash(kMer)
    graph.kMerMap[kMer] = graph.kMerMap.len 
    index = graph.kMermap.len - 1
    #for x,y in graph.kMerMap:
      #echo "  ", kMer.start..kMer.finish, ": '", kMer, "' -> ", graph.kMerMap[x], " hash=", hash(kMer)

    graph.edgesOut.add(@[])
    graph.edgesIn.add(@[])
    graph.kmers.add(kMer)
  else:
    index = graph.kMerMap[kMer]
  #echo "Index ",  index
  index


proc build*(sourceString: string, kmerLength: int): DeBruijnGraph =
  var t: DebruijnGraph
  new(t)
  t = DeBruijnGraph(source: sourceString,
                        kmers: @[], edgesOut: @[],
                        edgesIn: @[],
                        kMerMap: initTable[StrView, ID]())
  var length = (t.source.len() - (kmerLength))
  for i in 0 ..< length:
    let kmerL = t.source.toStrView(i, i+(kmerLength-1) )
    let kmerR = t.source.toStrView(i+1, (i+1)+(kmerLength-1) )
    #TODO In debug histrogramm of k-mers
    #echo "WTF", sourceString[kmerL.start..kmerL.finish], "...", sourceString[kmerR.start..kmerR.finish]
    let nodeL: ID = t.findOrCreateKmer(kmerL)
    let nodeR: ID = t.findOrCreateKmer(kmerR)

    #echo "ID",  nodeL, " and ", nodeR
    t.edgesOut[nodeL].add(nodeR)
    t.edgesIn[nodeR].add(nodeL)


  #[echo "LEN", t.kMerMap.len
  for x,y in t.kMerMap:
    echo  x, "' -> ", t.kMerMap[x] #, " hash=", hash(x)
  
  for x in t.kmers:
    echo x.base[x.start..x.finish]
  
  for k in t.kmers:
    echo k, "->", t.edgesOut[t.kmerMap[k]]
  ]#
  t

proc toDot*(g: DeBruijnGraph, filename: string)=
  let f = open(filename & ".dot", fmWrite)
  f.writeLine("digraph {")

  for s in 0 ..< g.edgesOut.len:
    for e in g.edgesOut[s]:
      f.writeLine(g.kmers[s], "->", g.kmers[e])
  f.writeLine("}")
  f.close()
