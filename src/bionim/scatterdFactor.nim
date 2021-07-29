## Implements the Smiah-Waterman algorithm used for local alignment
import utils
import algorithm



when isMainModule:
  var query = "anna"
  var seq = "banana"
  #var a = smithWaterman(query, seq, -1, 1, -2)
  printGrid(a, query, seq)
  #echo calculateAlignment(a, query, seq, -1, 1, -2)
  query = "anna"
  #seq = "banatnaannta"
  seq = "anttna"
 
  #a = smithWaterman(query, seq, -1, 1, -99)
  printGrid(a, query, seq)
  echo calculateAlignment(a, query, seq, -1, 1, -99)
