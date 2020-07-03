import numpy as np
#Attributes:
#Nodes is a list of the coordinates of the nodes
#EdgeNodes are a list of the edges, each element of this list is a pair with the each component being the position of the node in Nodes
#ElementEdges is a list of the cells by edges, each member is itself a list of edges. The numbering follows its appearance in EdgeNodes
#BoundaryNodes is list of the positions in Nodes of the nodes along the boundary of the domain
#Each element in Ortientations corresponds to the element in the same spot in ElementEdges. A 1 is placed in the ordering of the edge
#accords with the divergence Theorem. A -1 is placed if this is not the case.
class HeliosMesh(object):
    def __init__(self,Nodes,EdgeNodes,ElementEdges,NumBoundaryNodes,Orientations):
        self.Nodes            = Nodes
        self.EdgeNodes        = EdgeNodes
        self.ElementEdges     = ElementEdges
        self.NumBoundaryNodes = NumBoundaryNodes
        self.BNodes           = [Nodes[i] for i in NumBoundaryNodes]
        self.Orientations  = Orientations


        self.MakeDictionaries()
        self.MakeNumIntNodes()
    #MakeDictionaries creates two lists NodestoCells and EdgestoCells. 
    #NodestoCells will, given the position of a node in Nodes, return a list of the cells that have such a node.
    #EdgestoCells will, likewise, return the list of cells that have each edge.
    def MakeDictionaries(self):
        self.NodestoCells = [[] for i in self.Nodes]
        self.EdgestoCells = [[] for j in self.EdgeNodes]

        for c in range(len(self.ElementEdges)):
            Cell = self.ElementEdges[c]
            for Edge in Cell:
                self.EdgestoCells[Edge].append(c)
                Node1 = self.EdgeNodes[Edge][0]
                Node2 = self.EdgeNodes[Edge][1]
                if c not in self.NodestoCells[Node1]:
                    self.NodestoCells[Node1].append(c)
                if c not in self.NodestoCells[Node2]:
                    self.NodestoCells[Node2].append(c)
    
    def MakeNumIntNodes(self):
        numnodes           = len(self.Nodes)
        AllNodes           = [i for i in range(numnodes)] 
        self.NumInternalNodes = np.setdiff1d(AllNodes,self.NumBoundaryNodes)