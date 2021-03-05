import deques
import tree

#[type
  Node* = ref object
    parent*: Node
    children*: seq[Node]
    score*: int
    data*: set[char]
    name*: string
    length*: float
]#
type Data = ref object
  data: set[char]
  score: int

proc buildTree(sequences: seq[string], site: Natural): Node[Data] =
  let n = sequences.len() 

  var q = initDeque[Node[Data]]()
 
 
  for i in 0 ..< n:
    var data: Data
    new(data)
    data[].score = 0
    data.data = {sequences[i][site]} 
    q.addLast(Node[Data](parent: nil, children: @[], data: data, name: sequences[i], length: 0))
  
  while q.len > 1:
    let left  = q.popFirst
    let right = q.popFirst
    var data: Data
    new(data)
    
    if left.data[].data * right.data[].data == {}:
      data[].score = left.data.score + right.data.score + 1
      data[].data = left.data.data + right.data.data
      q.addLast(Node[Data](parent: nil, children: @[left, right], data: data, name: "", length: 0))
    else:
      data[].score = left.data.score + right.data.score 
      data[].data = left.data.data * right.data.data
      q.addLast(Node[Data](parent: nil, children: @[left, right], data: data, name: "", length: 0))
  
  result = q.popFirst

      
  
when isMainModule:
  let sequences = @["AACGT", "AATGT", "AATTT", "ACATT", "ACTTT"]
  var q = buildTree(sequences, 2) 
  
  proc print(q: Node)=
    echo $q.data.data
    for c in q.children:
      print c

  print q
      

  

