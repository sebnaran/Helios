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
        self.Orientations     = Orientations


        self.MakeDictionaries()
        self.ComputeMidponts() #This adds an array with the midpoints of every edge. It is in the same order as the edges
        self.ComputeBMidpoints() #Computes the boundary midpoints
        self.MakeNumIntNodes()   #Computes the internal midpoints of edges and internal nodes
                               
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
        numnodes                 = len(self.Nodes)
        AllNodes                 = [i for i in range(numnodes)] 
        self.NumInternalNodes    = np.setdiff1d(AllNodes,self.NumBoundaryNodes)

        numedges                 = len(self.EdgeNodes)
        AllEdges                 = [i for i in range(numedges)]
        self.NumInternalMidNodes = np.setdiff1d(AllEdges,self.NumBMidNodes)

    def ComputeMidponts(self):
        self.MidNodes = []
        for Edge in self.EdgeNodes:
            Node1, Node2 = self.Nodes[Edge[0]], self.Nodes[Edge[1]]
            x1,y1,x2,y2  = Node1[0],Node1[1],Node2[0],Node2[1]
            self.MidNodes.append([(x1+x2)/2,(y1+y2)/2])
    
    def ComputeBMidpoints(self):
        self.BMidNodes    = []
        self.NumBMidNodes = []
        i = 0
        for Node in self.MidNodes:
            x,y = Node[0],Node[1]
            if (abs(abs(x)-1) <1E-5 or abs(abs(y)-1)<1E-5):
                self.NumBMidNodes.append(i)
                self.BMidNodes.append(Node)
            i = i+1
        

    def StandardElement(self,Element,Ori):
    #This routine will reorient, if necessary, the edges of the element to agree with Gauss's theorem,
    #This is to say that the edges will be reoriented in such a way that the element will be traversed in the
    #Counterclockwise direction and rotation by pi/2 in the counterclockwise direction of the tangential vector
    #will result in an outward normal vector.
    #The last vertex,edge will be the first. This is in order to complete the loop.    
        N                = len(Element)
        OrientedEdges    = [0]*(N+1)
        OrientedVertices = [0]*(N+1)
        for i in range(N):
            if Ori[i] == 1:
                OrientedEdges[i] = self.EdgeNodes[Element[i]] #If they are "well-oriented" then do not alter them
            else:
                [v1,v2]          = self.EdgeNodes[Element[i]] #Otherwise reverse the order of their vertices
                OrientedEdges[i] = [v2,v1]

            OrientedVertices[i]  = self.Nodes[OrientedEdges[i][0]]
            OrientedEdges[N]     = OrientedEdges[0]
            OrientedVertices[N]  = OrientedVertices[0]
        return OrientedVertices,OrientedEdges

    def Centroid(self,Element,Ori):
        N  = len(Element)
        Cx = 0
        Cy = 0
        A  = 0
        Vertices,Edges = self.StandardElement(Element,Ori)
        for i in range(N):
            xi = Vertices[i][0]
            yi = Vertices[i][1]
            xiplusone = Vertices[i+1][0]
            yiplusone = Vertices[i+1][1]
            Cx        = Cx+(xi+xiplusone)*(xi*yiplusone-xiplusone*yi) #This formula is in Wikipedia
            Cy        = Cy+(yi+yiplusone)*(xi*yiplusone-xiplusone*yi)
            A         = A+xi*yiplusone-xiplusone*yi
        A = 0.5*A
        Cx = Cx/(6*A)
        Cy = Cy/(6*A)
        return Cx,Cy,A,Vertices,Edges

    def Area(self,Element,Ori):
        N              = len(Element)
        A              = 0
        Vertices,Edges = self.StandardElement(Element,Ori)
        for i in range(N):
            xi        = Vertices[i][0]
            yi        = Vertices[i][1]
            xiplusone = Vertices[i+1][0]
            yiplusone = Vertices[i+1][1]
            A         = A+xi*yiplusone-xiplusone*yi
        return 0.5*A,Vertices,Edges