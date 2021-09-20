{.experimental: "views".}
{.experimental: "strictFuncs".}

#@Discord haxscramper, helped me fixing the view stuff i.e implemente the StrView type since I could get it to work
#with the nim buil in toOpenArray

import std/[tables, hashes]
import utils
import sequtils
import bio_seq
import sets

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


template dbg(data: untyped) = 
  when not defined(release):
    echo data

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

proc `==`*(v1, v2: seq[StrView]): bool =
  for i in countup(0, len(v1)):
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
  source*: seq[string]
  ## Collection of all available strings view, to get a specifig
  ## string view just do kmers[ID]
  kmers*: seq[StrView]
  ## Collection of all available outgoing edges for a specific node
  edgesOut*: seq[seq[ID]]
  ## Collection of all available incoming edges for a specific node
  edgesIn*: seq[seq[ID]]
  ## Maps StrViews to ID, this enables you to get the ID for a specific kmer
  kMerMap*: Table[StrView, ID]


type WeightedDeBruijnGraph* = ref object
  ## The source string, this has to stay in scope since all StrViews 
  ## reference this string
  source*: seq[string]
  ## Collection of all available strings view, to get a specifig
  ## string view just do kmers[ID]
  kmers*: seq[StrView]
  ## Collection of all available outgoing edges for a specific node
  edgesOut*: Table[ID, Table[ID, int]]
  ## Collection of all available incoming edges for a specific node
  edgesIn*: Table[ID, Table[ID, int]]
  ## Maps StrViews to ID, this enables you to get the ID for a specific kmer
  kMerMap*: Table[StrView, ID]

type TempDebruijnGraph* = ref object
  ## The source string, this has to stay in scope since all StrViews 
  ## reference this string
  source*: seq[string]
  ## Collection of all available strings view, to get a specifig
  ## string view just do kmers[ID]
  contigs*: Table[string, ID]
  ## Collection of all available outgoing edges for a specific node
  edgesOut*: Table[ID, Table[ID, int]]
  ## Collection of all available incoming edges for a specific node
  edgesIn*: Table[ID, Table[ID, int]]
  ## Maps StrViews to ID, this enables you to get the ID for a specific kmer
  kMerMap*: Table[string, ID]


proc suffix(str: StrView): StrView=
  StrView(base: str.base, start: str.start+1, finish: str.finish)

proc prefix(str: StrView): StrView=
  StrView(base: str.base, start: str.start, finish: str.finish-1)


proc findOrCreateKmer(graph: var DeBruijnGraph, kMer: StrView): ID=
  var index: ID
  var pre = prefix(kMer)
  var suff = suffix(kMer)
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


proc findOrCreateKmerW(graph: var WeightedDeBruijnGraph, kMer: StrView): ID=
  var index: ID
  var pre = prefix(kMer)
  var suff = suffix(kMer)
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

    graph.edgesOut[index] = initTable[ID, int]()
    graph.edgesIn[index] = initTable[ID, int]()
    graph.kmers.add(kMer)
  else:
    index = graph.kMerMap[kMer]
  #echo "Index ",  index
  index


proc build*(sourceString: string, kmerLength: int): DeBruijnGraph =
  var t: DebruijnGraph
  new(t)
  t = DeBruijnGraph(source: @[sourceString],
                        kmers: @[], edgesOut: @[],
                        edgesIn: @[],
                        kMerMap: initTable[StrView, ID]())
  var length = (t.source.len() - (kmerLength))
  for i in 0 ..< length:
    let kmerL = t.source[0].toStrView(i, i+(kmerLength-1) )
    let kmerR = t.source[0].toStrView(i+1, (i+1)+(kmerLength-1) )
    #TODO In debug histrogramm of k-mers
    #echo "WTF", sourceString[kmerL.start..kmerL.finish], "...", sourceString[kmerR.start..kmerR.finish]
    let nodeL: ID = t.findOrCreateKmer(kmerL)
    let nodeR: ID = t.findOrCreateKmer(kmerR)

    #echo "ID",  nodeL, " and ", nodeR
    t.edgesOut[nodeL].add(nodeR)
    t.edgesIn[nodeR].add(nodeL)
  t

proc build*(reads: var seq[string], kmerLength: int): DeBruijnGraph =
  var t: DebruijnGraph
  new(t)
  t = DeBruijnGraph(source: reads,
                        kmers: @[], edgesOut: @[],
                        edgesIn: @[],
                        kMerMap: initTable[StrView, ID]())
  for read in reads.mitems:
    var length = (read.len() - (kmerLength-1))
    for i in 0 ..< length:
      let kmer = read.toStrView(i, i+(kmerLength-1) )
      let suff = kmer.suffix
      let pre = kmer.prefix
     # echo kmer
     # echo pre 
     # echo suff

      let nodeL: ID = t.findOrCreateKmer(pre)
      let nodeR: ID = t.findOrCreateKmer(suff)
     # echo "-----"

      #echo "ID",  nodeL, " and ", nodeR
      t.edgesOut[nodeL].add(nodeR)
      t.edgesIn[nodeR].add(nodeL)
  t

proc buildWeighted*(reads: var seq[string], kmerLength: int): WeightedDeBruijnGraph =
  var weighted: WeightedDeBruijnGraph 
  new(weighted)
  weighted = WeightedDeBruijnGraph( 
    source: reads,
    kmers: @[],
    edgesOut: initTable[ID, Table[ID, int]](), 
    edgesIn: initTable[ID, Table[ID, int]](), 
    kMerMap: initTable[StrView, ID]())
  
  for read in reads.mitems:
    var length = (read.len() - (kmerLength-1))
    for i in 0 ..< length:
      let kmer = read.toStrView(i, i+(kmerLength-1) )
      let suff = kmer.suffix
      let pre = kmer.prefix
     # echo kmer
     # echo pre 
     # echo suff

      let nodeL: ID = weighted.findOrCreateKmerW(pre)
      let nodeR: ID = weighted.findOrCreateKmerW(suff)

      #t.edgesOut[nodeL].add(nodeR)
      #t.edgesIn[nodeR].add(nodeL)
      #
      
      # If node L is already in our List then 
      if nodeL in weighted.edgesOut:
        #We still have to find out if there exists an entry for NodeR
        if nodeR in weighted.edgesOut[nodeL]:
          # If there is such an entry we just icrement it because we just got a kmer for ,at least, the second time
          weighted.edgesOut[nodeL][nodeR] += 1
        else:
          # Else nodeR is new and because we encountered it for the first time we set its weight to 1
          weighted.edgesOut[nodeL][nodeR]  = 1
      else:
        # We first have to init the Table for the nodeL entry
        weighted.edgesOut[nodeL] = initTable[ID, int]()
        #and then still do the check for nodeR
        if nodeR in weighted.edgesOut[nodeL]:
          # If there is such an entry we just icrement it because we just got a kmer for ,at least, the second time
          weighted.edgesOut[nodeL][nodeR] += 1
        else:
          # Else nodeR is new and because we encountered it for the first time we set its weight to 1
          weighted.edgesOut[nodeL][nodeR]  = 1

      # If nodeR is already in our List then 
      if nodeR in weighted.edgesIn:
        #We still have to find out if there exists an entry for NodeR
        if nodeL in weighted.edgesIn[nodeR]:
          # If there is such an entry we just icrement it because we just got a kmer for ,at least, the second time
          weighted.edgesIn[nodeR][nodeL] += 1
        else:
          # Else nodeR is new and because we encountered it for the first time we set its weight to 1
          weighted.edgesIn[nodeR][nodeL] = 1
      else:
        # We first have to init the Table for the nodeR entry
        weighted.edgesIn[nodeR] = initTable[ID, int]()
        #and then still do the check for nodeR
        if nodeL in weighted.edgesIn[nodeR]:
          # If there is such an entry we just icrement it because we just got a kmer for ,at least, the second time
          weighted.edgesIn[nodeR][nodeL] += 1
        else:
          # Else nodeR is new and because we encountered it for the first time we set its weight to 1
          weighted.edgesIn[nodeR][nodeL] = 1

  #echo weighted.edgesI
  #echo weighted.edgesOut

     # echo "-----"
  weighted

proc toDot*(g: TempDeBruijnGraph, filename: string)=
  let f = open(filename & ".dot", fmWrite)
  f.writeLine("digraph {")
  for contig,id in g.contigs:
    f.writeline(contig, " ", "[label=\"" & $contig & "," &  $id & "\"]")
  f.writeLine("}")
  f.close()

proc contigs*(graph: WeightedDeBruijnGraph): seq[string]=
  var t : TempDebruijnGraph
  new(t)
  var contigs: seq[string]
  var contig_counter: int64 = 0
  # Stores the end k-mer and the coresspondig contig id
  # USAGE:
  # To restore connection information you input a specific k-mer you then get all the ids to which the contigs coreesponding to
  # that `end k-mer` shoukd connect
  # Keep in mind that those are the IDs from the WeightedDBG
  # The second truple member corresponds to the edge weight
  var outgoing_kmers_map: Table[string, seq[(ID, int)]] 
  var incoming_kmers_map: Table[string, seq[ID]]
  
  t = TempDebruijnGraph( 
    source: graph.source,
    contigs: initTable[string, ID](),
    edgesOut: initTable[ID, Table[ID, int]](), 
    edgesIn: initTable[ID, Table[ID, int]](), 
    kMerMap: initTable[string, ID]())
  var marked =  newSeq[bool]()
  marked.setLen(len(graph.kMerMap))

  for str, id in graph.kMerMap:
    #TODO check that all entries are false by default
    if marked[id] == true:
      continue
    marked[id] = true
    # DIE COND FUCKT UP
    if( len(graph.edgesOut[id]) == 1 and len(graph.edgesIn[id]) == 1):
      # BECAUSE THIS ID SHOULD NEVER AGAIN OCCUR IT SHOULD BESAFE TO USE
      var newId: ID = id
      # We are in a path which can both ways
      dbg (id,":", $str, "     MIDDLE PATH")
      ## We now go forward until we reach a "bad" note 

      ##################### GOING FORWARD STEP ##########################################
      var outgoing_edges = graph.edgesOut[id]
      var next_node: ID
      var forward_contig: string = ""
      next_node = id
      # TODO I CANT ACCES NODES BY INDEX IN A TABLE
      # FOR LOOP HACK, but I know there is only one
      #for curr_id,weight in outgoing_edges:
      #   next_node = curr_id
      dbg ("NODE ID FORWARD PATH START", next_node, "  ",graph.kmers[next_node])
      while true:
        #TODO CHECK THAT; EX T3 is wrong with this inplace
        #if graph.edgesIn[next_node].len > 1:
          #break
        #DO STUFF 
        
        forward_contig.add($graph.kmers[next_node])
        marked[next_node] = true
        dbg ("CURRENT NODE  ", graph.kmers[next_node])
        
        if graph.edgesOut[next_node].len > 1:
          # This results in bubble sources inluded in contigs
          dbg ("FORWARD PATH BUUBBLES UP", next_node ,"  ", graph.kmers[next_node])
          dbg(" FORWARD CONNECTION ")
          dbg(graph.edgesOut[next_node])
          if not outgoing_kmers_map.contains($graph.kmers[next_node]):
            outgoing_kmers_map[$graph.kmers[next_node]] = @[]
          for k,v in graph.edgesOut[next_node]:
            outgoing_kmers_map[$graph.kmers[next_node]].add((k,v))
          dbg(" MAP")
          dbg(outgoing_kmers_map)
          break
        if graph.edgesOut[next_node].len == 0:
          dbg ("END OF FORWARD PATH", next_node ,"  ", graph.kmers[next_node])
          break
      
        # AND GO ON
        outgoing_edges = graph.edgesOut[next_node]
        var old_id = next_node
        for curr_id,weight in outgoing_edges:
          next_node = curr_id
          dbg ("next node ID ", next_node, " next node KMER: ",  graph.kmers[next_node])
        # WE ARE STILL ON A PATH
        if graph.edgesIn[next_node].len > 1:
          dbg ("NOT ALLOWED TO INCLUDE", next_node ,"  ", graph.kmers[next_node])
          ## IF WE GET "KILLED" HERE WE HAVE TO RETRIEVE THE CONNECTION INFO
          dbg(" FORWARD CONNECTION ")
          dbg(graph.edgesOut[old_id])
          if not outgoing_kmers_map.contains($graph.kmers[old_id]):
            outgoing_kmers_map[$graph.kmers[old_id]] = @[]
          for k,v in graph.edgesOut[old_id]:
            outgoing_kmers_map[$graph.kmers[old_id]].add((k,v))
          dbg(" MAP")
          dbg(outgoing_kmers_map)
          break
        dbg ("----------------")
      dbg ("FOWARD CONTIG ", forward_contig)
      ##################### GOING BACK STEP ##########################################
      ##
      var incoming_edges = graph.edgesIn[id]
      var previous_node: ID
      var backwards_contig: string = ""


      # TODO I CANT ACCES NODES BY INDEX IN A TABLE
      # FOR LOOP HACK, but I know there is only one
      for curr_id,weight in incoming_edges:
         previous_node = curr_id

      dbg ("NODE ID BEFORE BACKWARD STEP", previous_node, "  ",graph.kmers[previous_node])
      ## TODO When the previous node has more then one outgoing edge I would ignore that
      while true:
        if graph.edgesOut[previous_node].len > 1:
          dbg ("END OF BACKWARDSPATH, because SPLIT", next_node ,"  ", graph.kmers[next_node])
          break
        # DO STUFF 
        backwards_contig = ($graph.kmers[previous_node]) & backwards_contig
        marked[previous_node] = true
        dbg ("CURRENT NODE  ", graph.kmers[previous_node])

        if graph.edgesIn[previous_node].len >= 1:
          dbg ("SINK OF BUBBLE UP", next_node ,"  ", graph.kmers[next_node])
          break

        if graph.edgesIn[previous_node].len == 0:
          dbg ("END OF BACKWARDSPATH", next_node ,"  ", graph.kmers[next_node])
          break
      
        # AND GO ON
        incoming_edges = graph.edgesIn[previous_node]
        for curr_id,weight in incoming_edges:
          previous_node = curr_id
          dbg ("next node ID ", previous_node, " next node KMER: ",  graph.kmers[previous_node])
      ##
      ##
      ##
      dbg ("BACKWARD CONTIG ", backwards_contig)

      dbg ("FULL CONTIG: ", backwards_contig & forward_contig)
      #var s: Sequence = newSequence($contig_counter,  backwards_contig & forward_contig)
      contigs.add(backwards_contig & forward_contig)
      contig_counter += 1
      t.contigs[backwards_contig & forward_contig] = newId
      dbg ("------ PATH END  ---------")
      #TODO DO I NEED BACKWARDS/FORWARD PATH DISTINCTION
    elif( len(graph.edgesOut[id]) <= 1 and len(graph.edgesIn[id]) == 1):
      # We are in a path which can only go backwards
      dbg (id,":", $str, "   BACKWARDPATH")
      continue 
    elif ( len(graph.edgesOut[id]) ==  1 and len(graph.edgesIn[id]) <= 1):
      continue
    else:
      dbg (id,":", $str, "   BAD NODE")
      let l = t.kMerMap.len 
      #Safe "bad" kmer with old id
      t.kMerMap[$str] = id
      # Safe old edges which will later be used in the translation step
      t.edgesOut[l] = graph.edgesOut[id]
      t.edgesIn[l] = graph.edgesIn[id] 
  dbg ("-----------------CONTIGS----------")
  echo (contigs)
  t.toDot("test3")
  #contigs.writeFastaFile("contigs.fasta")
  return contigs

proc toDot*(g: DeBruijnGraph, filename: string)=
  let f = open(filename & ".dot", fmWrite)
  f.writeLine("digraph {")

  for s in 0 ..< g.edgesOut.len:
    for e in g.edgesOut[s]:
      f.writeLine(g.kmers[s], "->", g.kmers[e])
  f.writeLine("}")
  f.close()

proc toDot*(g: WeightedDeBruijnGraph, filename: string)=
  let f = open(filename & ".dot", fmWrite)
  f.writeLine("digraph {")
  for kmer in g.kmers:
    var k = g.kMerMap[kmer]
    f.writeline(kmer, " ", "[label=\"" & ($kmer) & "," &  $k & "\"]")

  f.writeline("")
  for src_id, edges in g.edgesOut:
    #echo src_id, " ", edges
    for dest_id, weight in edges:
      #echo dest_id
      f.writeLine(g.kmers[src_id], "->", g.kmers[dest_id], " ", "[label=\"" & $weight & "\"]")
  f.writeLine("}")
  f.close()



##DeBruijn Graph Type
type PairedDeBruijnGraph* = ref object
  ## The source string, this has to stay in scope since all StrViews 
  ## reference this string
  source*: seq[tuple[read1: string, read2: string]]
  ## Collection of all available strings view, to get a specifig
  ## string view just do kmers[ID]
  kmers*: seq[tuple[read1: StrView, read2: StrView]]
  ## Collection of all available outgoing edges for a specific node
  edgesOut*: seq[seq[ID]]
  ## Collection of all available incoming edges for a specific node
  edgesIn*: seq[seq[ID]]
  ## Maps StrViews to ID, this enables you to get the ID for a specific kmer
  kMerMap*: Table[tuple[read1: StrView, read2: StrView], ID]

proc findOrCreateKmer(graph: var PairedDeBruijnGraph, kMer: (StrView, StrView)): ID=
  var index: ID
  #TODO tab.mgetOrPut(key, @[]).add(val) can i do this in one lookup
  if not (kMer in graph.kMerMap):
    graph.kMerMap[kMer] = graph.kMerMap.len 
    index = graph.kMermap.len - 1

    graph.edgesOut.add(@[])
    graph.edgesIn.add(@[])
    graph.kmers.add(kMer)
  else:
    index = graph.kMerMap[kMer]
  index



proc build*(sourceString: seq[(string, string)], kmerLength: int): PairedDeBruijnGraph =
  var t: PairedDebruijnGraph
  new(t)
  t = PairedDeBruijnGraph(source: sourceString,
                        kmers: @[], edgesOut: @[],
                        edgesIn: @[],
                        kMerMap: initTable[(StrView, StrView), ID]())
  #For every paired end in our file                     
  for paired_reads in t.source.mitems: 
    #We assume that both reads have the same size
    #Maybe not true in real word. also we assume that the gap is always the same size
    var length = (paired_reads.read1.len() - (kmerLength-1))
    #For every read create all kmers
    for i in 0 ..< length:
      #Those two should always overlapt
      #KmerL
      let r1prefix = paired_reads.read1.toStrView(i, i+(kmerLength-2) )
      let r2prefix = paired_reads.read2.toStrView(i, i+(kmerLength-2) )
      #KmerR
      echo i+1+(kmerLength-1)
      let r1suffix = paired_reads.read1.toStrView(i+1, (i+1)+(kmerLength-2) ) 
      let r2suffix = paired_reads.read2.toStrView(i+1, (i+1)+(kmerLength-2) ) 
      let kmerL = (r1prefix, r2prefix)
      let kmerR = (r1suffix, r2suffix)
      let nodeL: ID = t.findOrCreateKmer(kmerL)
      let nodeR: ID = t.findOrCreateKmer(kmerR)
      echo "ID",  $nodeL, " and ", $nodeR
      t.edgesOut[nodeL].add(nodeR)
      t.edgesIn[nodeR].add(nodeL)
    result = t
    

proc prefix(pair :(string,string)):(string, string)=
  (pair[0].prefix, pair[1].prefix)

proc suffix(pair :(string,string)):(string, string)=
  (pair[0].suffix, pair[1].suffix)

proc toReadPairs(src: string, kmerLength, gap: int): seq[(string, string)]=
  var length = (src.len() - (kmerLength+gap))
  for i in 0 ..< length-(kmerLength - 1):
    var x: string = src[i .. i+(kmerLength-1)]
    var y: string = src[(i+kmerLength+gap) .. ((i-1)+gap+(2*kmerLength))]
    result.add((x,y))
  result

#TODO Change method to kmers instead of edges for iteration, this would fail to draw a node with degree 0
proc toDot*(g: PairedDeBruijnGraph, filename: string)=
  let f = open(filename & ".dot", fmWrite)
  f.writeLine("digraph {")
  for s in 0 ..< g.edgesOut.len:
    for e in g.edgesOut[s]:
      let node1 = $g.kmers[s].read1 & "\n" & $g.kmers[s].read2
      f.writeLine( s, " " , "[label=\"", node1, "\"]")

  for s in 0 ..< g.edgesOut.len:
    for e in g.edgesOut[s]:
      f.writeLine( s , "->", e,  ";")
  f.writeLine("}")
  f.close()

when isMainModule:
  #[ var b: PairedDebruijnGraph
  new(b)
  var seq1 = "TAATGCCATGGGATGTT"
  let pairs = toReadPairs(seq1, 3, 1)
  echo map(pairs, suffix)
  b = build(pairs, 3)
  b.toDot("paired")]#

  var b: DeBruijnGraph
  new(b)
  #var t: seq[string] = @["AAC", "ACB", "CBB", "BBA"]
  #var t: seq[string] = @["AAC", "ACB", "CBB", "BBA", "BBT"]
  #
  # EAXAMPLE T3
  #var t: seq[string] = @["AAC", "ACB", "CBB", "BBA", "BBT", "BTQ", "TQA"]
  #b = build(t, 3)
  var w: WeightedDeBruijnGraph 
  new(w)
  #var BA1 = @["AABCS", "ABCSS", "ABCSS", "ABDSS"] #interessant mit kmer=4
  var BA1 = @["TTAABCS", "TAABCSS", "AABCSSQ", "ABDSSQR"] #interessant mit kmer=4
  #BA1 = @["TTAABCS", "TAABCSS", "AABCSSQ", "ABDSSQR", "BDRTT"] #interessant mit kmer=4


  #let seq = parseFastQFile("reads/50.reads.fna.fq.gz")
  #w = buildWeighted(BA1, 3101) 
  #var str: seq[string] = @[]
  #for rec in seq.sequences:
  #  str.add($rec.nucs)
    
  w = buildWeighted(BA1, 3) 
  w.toDot("graph")
  discard w.contigs
  #discard b.toWeighted

