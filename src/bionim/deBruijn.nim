{.experimental: "views".}
{.experimental: "strictFuncs".}

#@Discord haxscramper, helped me fixing the view stuff i.e implemente the StrView type since I could get it to work
#with the nim buil in toOpenArray

import std/[tables, hashes]

## Typedef for ID's
## Each node in the DeBruijn-Graph has a unique id of this type
## At the moment its just an increment integer
type ID = int64
## String View type
## This are basically slices which still reference the original
## string instead of making copies
type
  StrView* = object
    start, finish: int
    base: ptr string


## Creates a string view from a give src string
## You have to make sure that the given src does not change is address, while a string view exists
proc toStrView*(str: var string, start, finish: int): StrView =
  if (start < 0):
    raise newException(Exception, "Your starting index must be >= 0")

  if (finish >= str.len):
    raise newException(Exception, "your finish cannot be greater then the length of the underlying string")
  StrView(base: addr str, start: start, finish: finish)

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

proc `$`*(view: StrView): string =
  for ch in view:
    result.add ch

proc len*(v1: StrView): int = abs(v1.finish - v1.start)+1
proc `[]`*(v: StrView, idx: int): char = v.base[][v.start + idx] 

proc `==`*(v1, v2: StrView): bool =
  if v1.len != v2.len:
    return false
  #Check content 
  for i in 0 ..< len(v1):
    if v1[i] != v2[i]:
      return false

  return true


proc `==`*(v1: string ,v2: StrView): bool =
  if v1.len != v2.len:
    return false
  #Check content 
  for i in 0 ..< len(v1):
    if v1[i] != v2[i]:
      return false

  return true


proc `==`*(v1: StrView,v2: string): bool =
  if v1.len != v2.len:
    return false
  #Check content 
  for i in 0 ..< len(v1):
    if v1[i] != v2[i]:
      return false

  return true

proc hash*(view: StrView): Hash =
  var h: Hash = 0
  for ch in view:
    h = h !& hash(ch)

  result = !$h

##DeBruijn Graph Type
type DeBruijnGraph* = ref object
  ## The source string, this has to stay in scope since all StrViews 
  ## reference this string
  source*: string
  ## Collection of all available strings view, to get a specifig
  ## string view just do kmers[ID]
  kmers*: seq[StrView]
  ## Collection of all available outgoing edges for a specific node
  edgesOut*: seq[seq[ID]]
  ## Collection of all available incoming edges for a specific node
  edgesIn*: seq[seq[ID]]
  ## Maps StrViews to ID, this enables you to get the ID for a specific kmer
  kMerMap*: Table[StrView, ID]



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
