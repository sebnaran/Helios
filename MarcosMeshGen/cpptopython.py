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
with open('AMRmesh.txt', "wb") as fp:
        pickle.dump((Nodes,EdgeNodes,ElementEdges,Orientations),fp)
#print(NodesPos)