import unittest
import bionim/utils
import bionim/deBruijn
import bionim/utils
import tables


proc `==`(v1: seq[string],v2: seq[StrView]): bool =
  if v1.len != v2.len:
    return false
  #Check content 
  for i in 0 ..< len(v1):
    if v1[i] != v2[i]:
      return false

  return true


proc testDebruijn()=
  suite "Test DeBruijn":
    test "Test StrView":
      var s  = "TEST"
      var sView = s.toStrView(0, len(s)-1)
      check sView.len == s.len
      check s == sView
      check sView == s
      for i in 0 ..< s.len:
        check s[i] == sView[i]

      sView = s.toStrView(0, 2)
      check sView.len == s[0..2].len


      sView = s.toStrView(0, 1)
      check sView.len == s[0..1].len

      sView = s.toStrView(0, 0)
      check sView.len == s[0..0].len
      
      expect Exception:
        sView = s.toStrView(0, s.len)

      expect Exception:
        sView = s.toStrView(-1, 1)
    
    test "DeBruijn init":
      var g: DeBruijnGraph = build("ACTGTC", 3)
      check g.source == "ACTGTC"
      var s = "ACTGTC"
      var kmers: seq[string] = @[]
      for k in s.kmers(3):
        kmers.add(k)  
      check kmers == g.kmers
      
        


testDeBruijn()
