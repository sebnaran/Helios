import pickle
from Functions import *
mesh = open("locrefs_quads_2.mesh","r")
RNodes        = []
ElementNodes = []
ElementEdges = []
EdgeNodes    = []

for i, lin  in enumerate(mesh):
    if i>8 and i<2867+8+1:
        Node = lin.split()
        RNodes.append([float(Node[0]),float(Node[1])])
    if i>2867+8+2 and i<2704+2867+8+3:
        linel    = lin.split()
        numnodes = int(linel[0])
        Element  = []
        for j in range(numnodes):
            Element.append(int(linel[j+1])-1)
        ElementNodes.append(Element)
    if i>2704+2867+8+3+1 and i<2*2704+2867+8+2+1+2:
        linel    = lin.split()
        numnodes = int(linel[0])
        Element  = []
        for j in range(numnodes):
            Element.append(int(linel[j+1])-1)
        ElementEdges.append(Element)
    if i>2*2704+2867+8+2+1+3 and i<2*2704+2867+8+2+2+3+5570:
        lined = lin.split()
        Edge  = [int(lined[0])-1,int(lined[1])-1]
        EdgeNodes.append(Edge)
mesh.close()
Nodes = []
for RNode in RNodes:
    Node = [2*RNode[0]-1,2*RNode[1]-1]
    Nodes.append(Node)
Orientations = []
for Element in ElementEdges:
    Ori = Orientation(Element,EdgeNodes,Nodes)
    Orientations.append(Ori)
i = 0
BoundaryNodes=[]
for Node in Nodes:
    x,y = Node[0],Node[1]
    if abs(x-1)<1E-5 or abs(x+1)<1E-5 or abs(y-1)<1E-5 or abs(y+1)<1E-5:
        BoundaryNodes.append(i)
    i = i+1
BottomToTop = []
LeftToRight = []
for i in BoundaryNodes:
    Node = Nodes[i]
    x,y  = Node[0],Node[1]
    if abs(x+1)<1E-5 and abs(y-1)>1E-5 and abs(y+1>1E-5):
        for j in BoundaryNodes:
            NNode = Nodes[j]
            xN,yN = NNode[0],NNode[1]
            if abs(xN-1)<1E-5 and abs(y-yN)<1E-5:
                LeftToRight.append([i,j])
    if abs(y+1)<1E-5 and abs(x-1)>1E-5 and abs(x+1>1E-5):
        for j in BoundaryNodes:
            NNode = Nodes[j]
            xN,yN = NNode[0],NNode[1]
            if abs(yN-1)<1E-5 and abs(xN-x)<1E-5:
                BottomToTop.append([i,j])
Corners = []
for i in BoundaryNodes:
    Node = Nodes[i]
    x,y  = Node[0],Node[1]
    if (abs(x-1)<1E-5 and abs(y-1)<1E-5) or (abs(x+1)<1E-5 and abs(y-1)<1E-5) or (abs(x-1)<1E-5 and abs(y+1)<1E-5) or (abs(x+1)<1E-5 and abs(y+1)<1E-5):
        Corners.append(i)
print(len(LeftToRight))
print(len(BottomToTop))
print(len(Corners))
with open('AMRmesh.txt', "wb") as fp:
    pickle.dump((Nodes,EdgeNodes,ElementEdges,BoundaryNodes,Orientations,BottomToTop,LeftToRight,Corners),fp)