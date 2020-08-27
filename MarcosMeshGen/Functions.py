from numpy import pi
import numpy as np
import math
#from sympy import Matrix
import pylab
#import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
#from scipy.interpolate import Rbf
import pickle
from scipy.sparse import csr_matrix
from scipy.sparse import lil_matrix
from scipy.sparse.linalg import spsolve
from EnergyClass import Energy
import matplotlib.pyplot as plt

#Data 
#Here the forcing diffusion coefficient, forcing terms, initial conditions, Dirichlet and Neumann conditions of the
#problem:
#Find E,B such that
#B_t=-curl E
#E=-nu curl B

C = 5

def Conductivity(x,y):
    #This is the diffusion coefficient
    #this function is scalar valued
    #Example 1
    #return 1
    #Energyplots 
    return 1

def EssentialBoundaryCond(x,y,t):
    #These are the essecial boundary conditions
    #example 1
    return -( 50*(math.exp(x)-math.exp(y))+math.cos(x*y)+math.sin(x*y) )*math.exp(-t)     
    #Energy Plots
    #return C*math.exp(C*t)*(50*math.exp(x)-50*math.exp(y)-math.sin(x*y)-math.cos(x*y))
    #Hartmann Flow 
    #return ( 2*math.sinh(0.5)-math.cosh(0.5) )/( 2*math.sinh(0.5) )

def InteriorCond(x,y):
    if abs(x-1)<10**(-5) or abs(x+1)<10**(-5) or abs(y-1)<10**(-5) or abs(y+1)<10**(-5):
        return 0
    else:
        return ( 2*math.sinh(0.5)-math.cosh(0.5) )/( 2*math.sinh(0.5) )
def InitialCond(x,y):
    #These are the initial condition on the magnetic field
    #r must be a 2 dimensional array
    #this function is vector valued
    #It must be divergence free
    #Example 1
    Bx = 50*math.exp(y)+x*math.sin(x*y)-x*math.cos(x*y)
    By = 50*math.exp(x)-y*math.sin(x*y)+y*math.cos(x*y)
    return Bx,By 
    #EnergyPlots
    #Bx = 50*math.exp(y)-x*math.sin(x*y)+x*math.cos(x*y)
    #By = 50*math.exp(x)+y*math.sin(x*y)-y*math.cos(x*y)
    #return Bx,By
    #Hartmann Flow
    #return (math.sinh(y)-2*y*math.sinh(0.5))/(2*math.sinh(0.5)),1
    #return 1,1
def ExactB(x,y,t):
    #This is the exact Magnetic field
    #Must be divergence free
    #Example 1
    Bx = ( 50*math.exp(y)+x*math.sin(x*y)-x*math.cos(x*y) )*math.exp(-t)
    By = ( 50*math.exp(x)-y*math.sin(x*y)+y*math.cos(x*y) )*math.exp(-t)
    return Bx,By#Solution 6 includes JxB
    #Example2
    #Bx,By = InitialCond(x,y)
    #return Bx*math.exp(-t),By*math.exp(-t)
    #EnergyPlots
    #Bx = (50*math.exp(y)-x*math.sin(x*y)+x*math.cos(x*y))*math.exp(C*t)
    #By = (50*math.exp(x)+y*math.sin(x*y)-y*math.cos(x*y))*math.exp(C*t)
    #return Bx,By
    #Hartmann Flow
    #return (math.sinh(y)-2*y*math.sinh(0.5))/(2*math.sinh(0.5)),1

def ExactE(x,y,t):
    #This is the exact Electric field
    #Example 1
    return -( 50*(math.exp(x)-math.exp(y))+math.cos(x*y)+math.sin(x*y) )*math.exp(-t)
    #EnergyPlots
    #return C*math.exp(C*t)*(50*math.exp(x)-50*math.exp(y)-math.sin(x*y)-math.cos(x*y))
    #Hartmann Flow
    #return ( 2*math.sinh(0.5)-math.cosh(0.5) )/( 2*math.sinh(0.5) )

def J(x,y):
    #Example 1
    Jx = ( (x**2+y**2-1)*(math.sin(x*y)+math.cos(x*y))-100*math.exp(x)+100*math.exp(y) )/( 2*(50*math.exp(x)-y*math.sin(x*y)+y*math.cos(x*y)) )
    Jy = -( (x**2+y**2-1)*(math.sin(x*y)+math.cos(x*y))-100*math.exp(x)+100*math.exp(y) )/( 2*(50*math.exp(y)+x*math.sin(x*y)-x*math.cos(x*y)) )
    return Jx,Jy 
    #Energyplots2, recall that due to how the code is written J = -u
    #Jx =  ( (-x**2-y**2-1)*(math.sin(x*y)+math.cos(x*y)) )/( 2*(50*math.exp(x)+y*math.sin(x*y)-y*math.cos(x*y)) )
    #Jy = -( (-x**2-y**2-1)*(math.sin(x*y)+math.cos(x*y)) )/( 2*(50*math.exp(y)-x*math.sin(x*y)+x*math.cos(x*y)) )
    #return C*Jx,C*Jy
    #Hartmann Flow
    #return -(math.cosh(0.5)-math.cosh(y))/(2*math.sinh(0.5)),0

def Poly1(x,y):
    return 1,0
def Poly2(x,y):
    return 0,1
def Poly3(x,y):
    return x,0
def Poly4(x,y):
    return y,0
def Poly5(x,y):
    return 0,x
def Poly6(x,y):
    return 0,y
def Poly(x,y):
    return x,y
















#MeshRetrievalFunctions


def GetMesh(file):
    #This function will, provided a text file in the format of meshes in Mathematica,
    #return the coordinates of the nodes, the Edge Nodes and the Nodes of each element
    UnprocessedMesh = open(file).read()
    FirstCut        = UnprocessedMesh.split('}}')
    Nodes           = FirstCut[0]+'}'
    EdgeNodes,Elements,trash = FirstCut[1].split(']}')
    for rep in (('{','['),('}',']'),('*^','*10**'),('Line[',''),('Polygon[',''),(']]',']')):
        Nodes     = Nodes.replace(rep[0],rep[1])
        EdgeNodes = EdgeNodes.replace(rep[0],rep[1])  #Replace the characters in the first position of the
        Elements  = Elements.replace(rep[0],rep[1])     #parenthesis with the second character
                                                    

    Nodes=Nodes+']'
    EdgeNodes = EdgeNodes+']'
    Elements  = Elements+']' #add the last ]
    Nodes     = eval(Nodes)
    EdgeNodes = eval(EdgeNodes)# turn string to list
    Elements  = eval(Elements)
    EdgeNodes = [(np.array(y)-1).tolist() for y in EdgeNodes]
    Elements  = [(np.array(y)-1).tolist() for y in Elements]
    #EdgeNodes=np.array(EdgeNodes)-1
    #Elements=np.array(Elements)-1 #subtract one from each position in the array(lists in mathematica begin with 1)
    return Nodes,EdgeNodes,Elements


#def EdgesElement(EdgeNodes,Elements):
    #This function will return the Edges of an element provided a list of the Nodes of the edges 
    #and the nodes of the elements
#    NumberElements=len(Elements)
#    NumberEdges=len(EdgeNodes)
#    ElementEdges=[0]*NumberElements
#    for i in range(NumberElements): #loop over elements
#        NumberNodesEdges=len(Elements[i]) #keep in mind there are the same number of Edges as there a
                                     #of nodes
#        ElementEdge=[0]*NumberNodesEdges    
#        Element=[0]*(NumberNodesEdges+1)
#        Element[0:NumberNodesEdges]=Elements[i] #pick a particular element
        
#        Element[NumberNodesEdges]=Element[0] # add the first vertex as the last vertex (so that edges become every pair)
#        for j in range(NumberNodesEdges): #run over every pair of vertices(edges)
#            Edge=[Element[j],Element[j+1]]
#            for ell in range(NumberEdges): #run over all edges
#                if Edge==EdgeNodes[ell] or [Edge[1],Edge[0]]==EdgeNodes[ell]: #if an the edge agrees, in any direction,
#                    ElementEdge[j]=ell #with one in the list of edges then we have
#        
#        ElementEdges[i]=ElementEdge
#    return ElementEdges                                                       #identified the edge

def EdgesElement(EdgeNodes,Elements):
    #This function will return the Edges of an element provided a list of the Nodes of the edges 
    #and the nodes of the elements
    NumberElements = len(Elements)
    NumberEdges = len(EdgeNodes)
    ElementEdges = [0]*NumberElements
    for i in range(NumberElements): #loop over elements
        NumberNodesEdges = len(Elements[i]) #keep in mind there are the same number of Edges as there a
                                     #of nodes
        ElementEdge = [0]*NumberNodesEdges    
        Element = [0]*(NumberNodesEdges+1)
        Element[0:NumberNodesEdges] = Elements[i] #pick a particular element
        
        Element[NumberNodesEdges] = Element[0] # add the first vertex as the last vertex (so that edges become every pair)
        for j in range(NumberNodesEdges): #run over every pair of vertices(edges)
            Edge = [Element[j],Element[j+1]]
            if Edge in EdgeNodes:
                ElementEdge[j] = EdgeNodes.index(Edge)
            else:
                Edge = [Element[j+1],Element[j]]
                ElementEdge[j] = EdgeNodes.index(Edge)
        ElementEdges[i]=ElementEdge
    return ElementEdges  


   
def Boundary(Nodes):
    #Given an array of Nodes this routine gives the nodes that lie on the boundary
    NumberNodes=len(Nodes)
    BoundaryNodes=[-1]*NumberNodes
    NumberBoundaryNodes=0
    for i in range(NumberNodes):
        Node=Nodes[i]
        if abs(Node[0]-1)<10**-10 or abs(Node[0]+1)<10**-10 or abs(Node[1]-1)<10**-10 or abs(Node[1]+1)<10**-10:
            BoundaryNodes[NumberBoundaryNodes]=i
            NumberBoundaryNodes=NumberBoundaryNodes+1
    return BoundaryNodes[0:NumberBoundaryNodes]

def Mesh(file):
    #Provided a file with a mesh in the language of mathematica this routine will return
    #four lists, Nodes, EdgeNodes,ElementEdges,BoundaryNodes
    Nodes,EdgeNodes,Elements=GetMesh(file) 
    ElementEdges=EdgesElement(EdgeNodes,Elements)
    BoundaryNodes=Boundary(Nodes)
    
    return Nodes,EdgeNodes,ElementEdges,BoundaryNodes

def FindVertecesEdges(Nodes,EdgeNodes):
    #This function, given a set of Edges, will return an array
    #the ith element of this array is the set of all edges that
    #have the ith vertex as an edpoint
    NumberNodes=len(Nodes)
    VertecesEdges=[[]]*NumberNodes
    i=0
    for Edge in EdgeNodes:
        v1=Edge[0]
        v2=Edge[1]
        VertecesEdges[v1]=list(set().union(VertecesEdges[v1],[i]))
        VertecesEdges[v2]=list(set().union(VertecesEdges[v2],[i]))
        i=i+1
    return VertecesEdges


def ProcessedMesh(Pfile):
    with open(Pfile, "rb") as fp:   # Unpickling
        N,E,EE,B,O = pickle.load(fp)
    return N,E,EE,B,O



















#AuxiliaryFunctions

def ElementCoordinates(Element,EdgeNodes,Nodes):
    #provided an element of a mesh this function returns its vertices 
    #as an array of the dimension of the number of edges on the element
    N=len(Element)
    Vertices=[0]*N #The number of vertices agrees with the number of edges
    Edge=Element[0]
    Vertices[0]=Nodes[EdgeNodes[Edge][0]]
    Vertices[1]=Nodes[EdgeNodes[Edge][1]] #The first two vertices are those in the first edge
    for i in range(1,N-1):
        Edge=EdgeNodes[Element[i]]
        v1=Nodes[Edge[0]]
        v2=Nodes[Edge[1]]    #new vertex added is the one on the i+1 edge that is not in the ith edge
        if Vertices[i-1]==v1:
            Vertices[i+1]=v2
        elif Vertices[i-1]==v2:
            Vertices[i+1]=v1
        elif Vertices[i]==v1:
            Vertices[i+1]=v2
        else:
            Vertices[i+1]=v1
    return Vertices

def Orientation(Element,EdgeNodes,Nodes):
    #Provided an element of a mesh this function returns a vector of dimension
    #as large as the number of edges with the orientation of the normal vector
    #1 if the orientation is outward and -1 if it is inward.
    #This algorithm assumes that the element is convex
    
    N=len(Element)
    #first we need to find a point inside the element to give indication of the orientation
    Vertices=ElementCoordinates(Element,EdgeNodes,Nodes)
    xinside=0
    yinside=0
    for i in range(N):
        xinside=xinside+Vertices[i][0]
        yinside=yinside+Vertices[i][1]
    xinside=xinside/N
    yinside=yinside/N   #since the element is convex any convex combination of the vertices lies inside
    
    #Now we move on to finding the orientation of the edges
    Ori=[0]*(N+1)
    
    for i in range(N):
        Edge=EdgeNodes[Element[i]]
        [x1,y1]=Nodes[Edge[0]]
        [x2,y2]=Nodes[Edge[1]]
        #This is the inner product between n, the vector t=(x2,y2)-(x1,y1) rotated pi/2 counterclockwise
        # and the vector u=(xinside,yinside)-(x1,y1) which should point towards the interior of the element
        #if the result is negative it means that the angle between t and u is larger than 90 which implies that
        #n points outside of the element
        sign=(y2-y1)*(xinside-x1)+(x1-x2)*(yinside-y1) 
        if sign<0:
            Ori[i]=1
        else: 
            Ori[i]=-1
    Ori[N]=Ori[0]
    return Ori

def StandardElement(Element,EdgeNodes,Nodes,Ori):
    #This routine will reorient, if necessary, the edges of the element to agree with stokes theorem,
    #This is to say that the edges will be reoriented in such a way that the element will be traversed in the
    #Counterclockwise direction and rotation by pi/2 in the counterclockwise direction of the tangential vector
    #will result in an outward normal vector.
    #The last vertex,edge will be the first. This is in order to complete the loop.
    N                = len(Element)
    OrientedEdges    = [0]*(N+1)
    OrientedVertices = [0]*(N+1)
    
    for i in range(N):
        if Ori[i]==1:
            OrientedEdges[i] = EdgeNodes[Element[i]] #If they are "well-oriented" then do not alter them
        else:
            [v1,v2]          = EdgeNodes[Element[i]] #Otherwise reverse the order of their vertices
            OrientedEdges[i] = [v2,v1]

        OrientedVertices[i]  = Nodes[OrientedEdges[i][0]]
    OrientedEdges[N]    = OrientedEdges[0]
    OrientedVertices[N] = OrientedVertices[0]
    return OrientedVertices,OrientedEdges


def Centroid(Element,EdgeNodes,Nodes,Ori):
    #This function, when provided with an element, will return its centroid or barycenter.
    N  = len(Element)
    Cx = 0
    Cy = 0
    A  = 0
    Vertices,Edges = StandardElement(Element,EdgeNodes,Nodes,Ori)
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

def InternalObjects(Boundary,Objects):
    #provided a set of geometrical objects, say vertices or edges, this routine returns those that
    #are internal. 
    N=len(Objects)
    Internal=np.sort(np.array(list(set(np.arange(N))-set(Boundary))))  
    NumberInternal=len(Internal)
    return Internal,NumberInternal
    


#Assembly


        
def LocprojE(Func,Element,EdgeNodes,Nodes):
    #This function will, provided a function a set of nodes and edges, compute the 
    #projection onto the space of edge-based functions. The direction of the unit normal
    #will be assumed to be the clockwise rotation of the tangential vector.
    
    N=len(Element)
    proj=np.zeros((N,1))
    j = 0
    for i in Element:
        x1=Nodes[EdgeNodes[i][0]][0]
        y1=Nodes[EdgeNodes[i][0]][1]
        x2=Nodes[EdgeNodes[i][1]][0]
        y2=Nodes[EdgeNodes[i][1]][1]
        lengthe=math.sqrt((x2-x1)**2+(y2-y1)**2)
        xmid=0.5*(x1+x2)
        ymid=0.5*(y1+y2)
        etimesnormal=[y2-y1,x1-x2]
        Fx,Fy=Func(xmid,ymid)
        proj[j]=(etimesnormal[0]*Fx+etimesnormal[1]*Fy)*lengthe**-1  #midpoint rule
        j = j+1
    return proj




def MFDMEWEMVWV(J,Basis,Element,EdgeNodes,Nodes,Ori):
    #This routine will compute the local mass matrix in the edge-based space E
    #Here we must ensure that the orientation of the elements is such that
    #We have an orientation for the edges that respects stoke's theorem
    n = len(Element)
    Dim = len(Basis)
    xP,yP,A,Vertices,Edges = Centroid(Element,EdgeNodes,Nodes,Ori)
    nu = Conductivity(xP,yP)
    NE = np.zeros((n,2))
    RE = np.zeros((n,2))
    NV = np.ones((n,1))*nu
    RV = np.zeros((n,1))
    
    for i in range(n):
        x1 = Vertices[i][0]
        y1 = Vertices[i][1]
        x2 = Vertices[i+1][0]
        y2 = Vertices[i+1][1]
        lengthEdge = math.sqrt((x2-x1)**2+(y2-y1)**2)
        #print('dif in x')
        #print(x2-x1)
        #print('dif in y')
        #print(y2-y1)
        #print('Length e')
        #print(lengthEdge)
        NE[i][0] = (y2-y1)*Ori[i]*lengthEdge**-1
        NE[i][1] = (x1-x2)*Ori[i]*lengthEdge**-1
        RE[i][0] = (0.5*(x1+x2)-xP)*Ori[i]*lengthEdge #These formulas are derived in the tex-document
        RE[i][1] = (0.5*(y1+y2)-yP)*Ori[i]*lengthEdge
    ME = LocalMassMatrix(NE,RE,n,A,1)
    #WE=LocalMassMatrix(RE,NE,n,A,1)
    
    
    
    x1n = Vertices[n-1][0] #first vertex of n-1th-edge
    y1n = Vertices[n-1][1]
    
    x2n = Vertices[0][0]
    y2n = Vertices[0][1] #second vertex of n-1th edge
    
    x11 = x2n #first vertex of first edge
    y11 = y2n
    
    x21 = Vertices[1][0]
    y21 = Vertices[1][1]  #second vertex of first edge
    
    omegan2 = (x2n-x1n)*((yP-y2n)+(2*yP-y1n-y2n))/6
    omega11 = (x21-x11)*((yP-y11)+(2*yP-y11-y21))/6
    RV[0] = omegan2+omega11
   
    for i in range(1,n):
        
        x1iminusone = Vertices[i-1][0] #first vertex of i-1th-edge
        y1iminusone = Vertices[i-1][1]
               
        x2iminusone = Vertices[i][0]    
        y2iminusone = Vertices[i][1] #second vertex of i-1th edge

        x1i = x2iminusone #first vertex of i+1 edge
        y1i = y2iminusone
    
        x2i = Vertices[i+1][0] #second vertex of i+1 edge
        y2i = Vertices[i+1][1]  
        
        omega2iminusone = (x2iminusone-x1iminusone)*((yP-y2iminusone)+\
                                                   (2*yP-y1iminusone-y2iminusone))/6
        omega1i = (x2i-x1i)*((yP-y1i)+(2*yP-y2i-y1i))/6   
        
        RV[i] = omega2iminusone+omega1i
        
    MV = LocalMassMatrix(NV,RV,n,A,nu)
    #WV=LocalMassMatrix(RV,NV,n,A,nu)
    
    NJ = np.zeros((Dim,n))
    
    
    for i in range(Dim):        
        NJ[i,:] = np.transpose( LocprojE(Basis[i],Element,EdgeNodes,Nodes) )
    
    
    NJ = np.transpose(NJ)
    #print(NJ)
    b = np.transpose(NJ).dot(ME)
    #print(b)
    #print(ME)
    #print(NJ)
    #print(np.transpose(NJ).dot(ME).dot(NJ))
    #print(np.linalg.inv( np.transpose(NJ).dot(ME).dot(NJ) ) )
    
    MJ = np.linalg.pinv( np.transpose(NJ).dot(ME).dot(NJ) )
    #print(MJ)
    MJ = MJ.dot(b)
    #print(MJ)



    
    PolyCoordinates = np.zeros((2*(len(Vertices)-1),Dim))
    JMatrix = np.zeros( (len(Vertices)-1,2*(len(Vertices)-1)) ) 
                       
    l = 0
    k = 0
    for Polynomial in Basis:
        
        
        for j in range(len(Vertices)-1):
            Vertex = Vertices[j]
            x = Vertex[0]
            y = Vertex[1]
            
            Px,Py = Polynomial(x,y)
            
            PolyCoordinates[2*j,l] = Px
            PolyCoordinates[2*j+1,l] = Py
            
            
    
            if k==0:
                Jx,Jy = J(x,y)
                JMatrix[j,2*j] = -Jy
                JMatrix[j,2*j+1] = Jx
                
            j = j+1
        k = 1        
        l = l+1
    
    MJ = JMatrix.dot(PolyCoordinates).dot(MJ)
    #MJ = np.transpose(MJ).dot(MV)
    MJ = MV.dot(MJ)
    #print(MJ)
    return ME,MV,MJ,Edges    
#The following three functions are necessary for the construction of the 
#mass matrix in the nodal space.
def m1(x,y,xP,yP):
        return 1

def m2(x,y,xP,yP):
        return x-xP

def m3(x,y,xP,yP):
        return y-yP

def P0(ElNodes,func,xP,yP):
        po = 0
        for node in ElNodes:
                x  = node[0]
                y  = node[1]
                po = po+func(x,y,xP,yP)
        return po/len(ElNodes)


def NewLocalMEWEMVWV(J,Basis,Element,EdgeNodes,Nodes,Ori):
    #This routine will compute the local mass matrix in the edge-based space E
    #Here we must ensure that the orientation of the elements is such that
    #We have an orientation for the edges that respects stoke's theorem
    n                      = len(Element)
    Dim                    = len(Basis)
    xP,yP,A,Vertices,Edges = Centroid(Element,EdgeNodes,Nodes,Ori)
    NE                     = np.zeros((n,2))
    RE                     = np.zeros((n,2))
    
    for i in range(n):
        x1         = Vertices[i][0]
        y1         = Vertices[i][1]
        x2         = Vertices[i+1][0]
        y2         = Vertices[i+1][1]
        lengthEdge = math.sqrt((x2-x1)**2+(y2-y1)**2)
        NE[i][0]   = (y2-y1)*Ori[i]*lengthEdge**-1
        NE[i][1]   = (x1-x2)*Ori[i]*lengthEdge**-1
        RE[i][0]   = (0.5*(x1+x2)-xP)*Ori[i]*lengthEdge #These formulas are derived in the tex-document
        RE[i][1]   = (0.5*(y1+y2)-yP)*Ori[i]*lengthEdge
    ME = LocalMassMatrix(NE,RE,n,A,1)
    #WE=LocalMassMatrix(RE,NE,n,A,1)
    
    #########################
    #Here we will construct the local nodal mass matrix

    OrVert,OrEdg = StandardElement(Element,EdgeNodes,Nodes,Ori)
    ElNodes      = OrVert[0:n]    
    ElEdges      = OrEdg[0:n]
    G            = np.zeros((3,3))

    G[0,0] = 1

    G[0,1] = P0(ElNodes,m2,xP,yP)
    G[1,1] = A
    
    G[0,2] = P0(ElNodes,m3,xP,yP)
    G[2,2] = A

    B      = np.ones((3,n))/n
    H      = np.zeros((3,3))
    H[0,0] = A
            
    for i in range(n):
        x1         = Vertices[i][0]
        y1         = Vertices[i][1]
        x2         = Vertices[i+1][0]
        y2         = Vertices[i+1][1]
        
        lengthedge  = math.sqrt((x2-x1)**2+(y2-y1)**2)
        taux        = (x2-x1)/lengthedge
        tauy        = (y2-y1)/lengthedge

        B[1,i]      =  0.5*lengthedge*tauy
        B[2,i]      = -0.5*lengthedge*taux
        h           = lengthedge/3
        nx          = tauy
        ny          = -taux
        costheta    = (x2-x1)/lengthedge
        sintheta    = (y2-y1)/lengthedge
        
        xot         = x1+h*costheta
        yot         = y1+h*sintheta
        
        xtt         = x1+2*h*costheta
        ytt         = y1+2*h*sintheta
        
        H[1,1]      = H[1,1] + h*nx*( m2(x1,y1,xP,yP)**3+\
                                 3*m2(xot,yot,xP,yP)**3+\
                                 3*m2(xtt,ytt,xP,yP)**3+\
                                   m2(x2,y2,xP,yP)**3 )/8  

        H[2,2]      = H[2,2] + h*nx*( m3(x1,y1,xP,yP)**3+\
                                 3*m3(xot,yot,xP,yP)**3+\
                                 3*m3(xtt,ytt,xP,yP)**3+\
                                   m3(x2,y2,xP,yP)**3 )/8       
   

        H[1,2]      = H[1,2] + 3*h*nx*( m3(x1,y1,xP,yP)*m2(x1,y1,xP,yP)**2+\
                                   3*m3(xot,yot,xP,yP)*m2(xot,yot,xP,yP)**2+\
                                   3*m3(xtt,ytt,xP,yP)*m2(xtt,ytt,xP,yP)**2+\
                                     m3(x2,y2,xP,yP)*m2(x2,y2,xP,yP)**2 )/16

        H[2,1]      = H[2,1] + 3*h*nx*( m3(x1,y1,xP,yP)*m2(x1,y1,xP,yP)**2+\
                                   3*m3(xot,yot,xP,yP)*m2(xot,yot,xP,yP)**2+\
                                   3*m3(xtt,ytt,xP,yP)*m2(xtt,ytt,xP,yP)**2+\
                                     m3(x2,y2,xP,yP)*m2(x2,y2,xP,yP)**2 )/16

  
       
    D      = np.ones((n,3))

    D[:,1] = [m2(x,y,xP,yP) for [x,y] in ElNodes]    
    D[:,2] = [m3(x,y,xP,yP) for [x,y] in ElNodes]

    Pistar = np.linalg.inv(G).dot(B)
    Pi     = D.dot(Pistar)
    Id     = np.identity(n)
    MV     = np.transpose(Pistar).dot(H.dot(Pistar))+A*np.transpose(Id-Pi).dot(Id-Pi)    
    MV     = MV*Conductivity(xP,yP)
    NJ = np.zeros((Dim,n))
    
    
    for i in range(Dim):        
        NJ[i,:] = np.transpose( LocprojE(Basis[i],Element,EdgeNodes,Nodes) )
    
    
    NJ = np.transpose(NJ)
    #print(NJ)
    b = np.transpose(NJ).dot(ME)
    #print(b)
    #print(ME)
    #print(NJ)
    #print(np.transpose(NJ).dot(ME).dot(NJ))
    #print(np.linalg.inv( np.transpose(NJ).dot(ME).dot(NJ) ) )
    
    MJ = np.linalg.pinv( np.transpose(NJ).dot(ME).dot(NJ) )
    #print(MJ)
    MJ = MJ.dot(b)
    #print(MJ)



    
    PolyCoordinates = np.zeros((2*(len(Vertices)-1),Dim))
    JMatrix = np.zeros( (len(Vertices)-1,2*(len(Vertices)-1)) ) 
                       
    l = 0
    k = 0
    for Polynomial in Basis:
        
        
        for j in range(len(Vertices)-1):
            Vertex = Vertices[j]
            x = Vertex[0]
            y = Vertex[1]
            
            Px,Py = Polynomial(x,y)
            
            PolyCoordinates[2*j,l] = Px
            PolyCoordinates[2*j+1,l] = Py

            if k==0:
                Jx,Jy = J(x,y)
                JMatrix[j,2*j] = -Jy
                JMatrix[j,2*j+1] = Jx
                
            j = j+1
        k = 1        
        l = l+1
    
    MJ = JMatrix.dot(PolyCoordinates).dot(MJ)
    MJ = MV.dot(MJ)
    return ME,MV,MJ,Edges

def LeastSquaresLocalMEWEMVWV(J,Basis,Element,EdgeNodes,Nodes,Ori):
    #This routine will compute the local mass matrix in the edge-based space E
    #Here we must ensure that the orientation of the elements is such that
    #We have an orientation for the edges that respects stoke's theorem
    n                      = len(Element)
    Dim                    = len(Basis)
    xP,yP,A,Vertices,Edges = Centroid(Element,EdgeNodes,Nodes,Ori)
    
    NE                     = np.zeros((n,2))
    RE                     = np.zeros((n,2))
    
    for i in range(n):
        x1         = Vertices[i][0]
        y1         = Vertices[i][1]
        x2         = Vertices[i+1][0]
        y2         = Vertices[i+1][1]
        lengthEdge = math.sqrt((x2-x1)**2+(y2-y1)**2)
        NE[i][0]   = (y2-y1)*Ori[i]*lengthEdge**-1
        NE[i][1]   = (x1-x2)*Ori[i]*lengthEdge**-1
        RE[i][0]   = (0.5*(x1+x2)-xP)*Ori[i]*lengthEdge #These formulas are derived in the tex-document
        RE[i][1]   = (0.5*(y1+y2)-yP)*Ori[i]*lengthEdge
    ME = LocalMassMatrix(NE,RE,n,A,1)
    #WE=LocalMassMatrix(RE,NE,n,A,1)
    
    #########################
    #Here we will construct the local nodal mass matrix

    OrVert,OrEdg = StandardElement(Element,EdgeNodes,Nodes,Ori)
    ElNodes      = OrVert[0:n]    
    ElEdges      = OrEdg[0:n]

    H      = np.zeros((3,3))
    H[0,0] = A
            
    for i in range(n):
        x1         = Vertices[i][0]
        y1         = Vertices[i][1]
        x2         = Vertices[i+1][0]
        y2         = Vertices[i+1][1]
        
        lengthedge  = math.sqrt((x2-x1)**2+(y2-y1)**2)
        taux        = (x2-x1)/lengthedge
        tauy        = (y2-y1)/lengthedge

        h           = lengthedge/3
        nx          = tauy
        ny          = -taux
        costheta    = (x2-x1)/lengthedge
        sintheta    = (y2-y1)/lengthedge
        
        xot         = x1+h*costheta
        yot         = y1+h*sintheta
        
        xtt         = x1+2*h*costheta
        ytt         = y1+2*h*sintheta
        
        H[1,1]      = H[1,1] + h*nx*( m2(x1,y1,xP,yP)**3+\
                                 3*m2(xot,yot,xP,yP)**3+\
                                 3*m2(xtt,ytt,xP,yP)**3+\
                                   m2(x2,y2,xP,yP)**3 )/8  

        H[2,2]      = H[2,2] + h*nx*( m3(x1,y1,xP,yP)**3+\
                                 3*m3(xot,yot,xP,yP)**3+\
                                 3*m3(xtt,ytt,xP,yP)**3+\
                                   m3(x2,y2,xP,yP)**3 )/8       
   

        H[1,2]      = H[1,2] + 3*h*nx*( m3(x1,y1,xP,yP)*m2(x1,y1,xP,yP)**2+\
                                   3*m3(xot,yot,xP,yP)*m2(xot,yot,xP,yP)**2+\
                                   3*m3(xtt,ytt,xP,yP)*m2(xtt,ytt,xP,yP)**2+\
                                     m3(x2,y2,xP,yP)*m2(x2,y2,xP,yP)**2 )/16

        H[2,1]      = H[2,1] + 3*h*nx*( m3(x1,y1,xP,yP)*m2(x1,y1,xP,yP)**2+\
                                   3*m3(xot,yot,xP,yP)*m2(xot,yot,xP,yP)**2+\
                                   3*m3(xtt,ytt,xP,yP)*m2(xtt,ytt,xP,yP)**2+\
                                     m3(x2,y2,xP,yP)*m2(x2,y2,xP,yP)**2 )/16

  
       
    D      = np.ones((n,3))

    D[:,1] = [m2(x,y,xP,yP) for [x,y] in ElNodes]    
    D[:,2] = [m3(x,y,xP,yP) for [x,y] in ElNodes]

    Pistar = np.linalg.inv(np.transpose(D).dot(D)).dot(np.transpose(D))
    Pi     = D.dot(Pistar)
    Id     = np.identity(n)
    MV     = np.transpose(Pistar).dot(H.dot(Pistar))+A*np.transpose(Id-Pi).dot(Id-Pi)    
   
    NJ = np.zeros((Dim,n))
    
    
    for i in range(Dim):        
        NJ[i,:] = np.transpose( LocprojE(Basis[i],Element,EdgeNodes,Nodes) )
    
    
    NJ = np.transpose(NJ)
    #print(NJ)
    b = np.transpose(NJ).dot(ME)
    #print(b)
    #print(ME)
    #print(NJ)
    #print(np.transpose(NJ).dot(ME).dot(NJ))
    #print(np.linalg.inv( np.transpose(NJ).dot(ME).dot(NJ) ) )
    
    MJ = np.linalg.pinv( np.transpose(NJ).dot(ME).dot(NJ) )
    #print(MJ)
    MJ = MJ.dot(b)
    #print(MJ)



    
    PolyCoordinates = np.zeros((2*(len(Vertices)-1),Dim))
    JMatrix = np.zeros( (len(Vertices)-1,2*(len(Vertices)-1)) ) 
                       
    l = 0
    k = 0
    for Polynomial in Basis:
        
        
        for j in range(len(Vertices)-1):
            Vertex = Vertices[j]
            x = Vertex[0]
            y = Vertex[1]
            
            Px,Py = Polynomial(x,y)
            
            PolyCoordinates[2*j,l] = Px
            PolyCoordinates[2*j+1,l] = Py
            
            
    
            if k==0:
                Jx,Jy = J(x,y)
                JMatrix[j,2*j] = -Jy
                JMatrix[j,2*j+1] = Jx
                
            j = j+1
        k = 1        
        l = l+1
    
    MJ = JMatrix.dot(PolyCoordinates).dot(MJ)
    MJ = MV.dot(MJ)
    return ME,MV,MJ,Edges

def psi3(x,y):
    return 1-x-y
def psi2(x,y):
    return y
def psi1(x,y):
    return x

def convexcoeffs(Element,EdgeNodes,Nodes,Ori):
    n                      = len(Element)
    OrVert,OrEdg           = StandardElement(Element,EdgeNodes,Nodes,Ori)
    xP,yP,A,Vertices,Edges = Centroid(Element,EdgeNodes,Nodes,Ori)

    ElNodes  = OrVert[0:n] 
    xs       = 0
    ys       = 0
    for node in ElNodes:
        xs = xs+node[0]
        ys = ys+node[1]
    
    xs = xs/n
    ys = ys/n
    Areas = [0]*n
    #2*Area = Ax(By - Cy) + Bx(Cy - Ay) + Cx(Ay - By)
    #A=starred B = kth C=k+1th
    xk         = ElNodes[n-1][0]
    yk         = ElNodes[n-1][1]
    xkp1       = ElNodes[0][0]
    ykp1       = ElNodes[0][1]
    Areas[n-1] = abs(xs*(yk-ykp1)+xk*(ykp1-ys)+xkp1*(ys-yk))/2
    for i in range(n-1):
        xk       = ElNodes[i][0]
        yk       = ElNodes[i][1]
        xkp1     = ElNodes[i+1][0]
        ykp1     = ElNodes[i+1][1]
        Areas[i] = abs(xs*(yk-ykp1)+xk*(ykp1-ys)+xkp1*(ys-yk))/2
    w    = [0]*n
    w[0] = ( (n+1)*(Areas[0]+Areas[n-1])+sum(Areas[1:n-1]) )/(3*n*A)
    
    for i in range(1,n):
        A1         = Areas[0:i-1]
        A2         = Areas[i+1:n]
        s          = sum(A1)
        p          = sum(A2)
        w[i] = ( (n+1)*(Areas[i]+Areas[i-1])+s+p )/(3*n*A)

    return w, Areas
def PieceWiseLocalMEWEMVWV(J,Basis,Element,EdgeNodes,Nodes,Ori):
    #This routine will compute the local mass matrix in the edge-based space E
    #Here we must ensure that the orientation of the elements is such that
    #We have an orientation for the edges that respects stoke's theorem
    n                      = len(Element)
    Dim                    = len(Basis)
    xP,yP,A,Vertices,Edges = Centroid(Element,EdgeNodes,Nodes,Ori)
    
    NE                     = np.zeros((n,2))
    RE                     = np.zeros((n,2))
    
    for i in range(n):
        x1         = Vertices[i][0]
        y1         = Vertices[i][1]
        x2         = Vertices[i+1][0]
        y2         = Vertices[i+1][1]
        lengthEdge = math.sqrt((x2-x1)**2+(y2-y1)**2)
        NE[i][0]   = (y2-y1)*Ori[i]*lengthEdge**-1
        NE[i][1]   = (x1-x2)*Ori[i]*lengthEdge**-1
        RE[i][0]   = (0.5*(x1+x2)-xP)*Ori[i]*lengthEdge #These formulas are derived in the tex-document
        RE[i][1]   = (0.5*(y1+y2)-yP)*Ori[i]*lengthEdge
    ME = LocalMassMatrix(NE,RE,n,A,1)
    #WE=LocalMassMatrix(RE,NE,n,A,1)
    
    #########################
    #Here we will construct the local nodal mass matrix

    OrVert,OrEdg = StandardElement(Element,EdgeNodes,Nodes,Ori)
    ElNodes      = OrVert[0:n]    
    #ElEdges      = OrEdg[0:n]

    H            = np.zeros((n+1,n+1))
    xstar        = xP
    ystar        = yP

    #We will map each triangle to the standard triangle with vertices at
    #(0,0), (1,0) and (0,1). We will always map v* to (0,0), 
    #the leftmost point (in the list) will be mapped to (1,0) while
    #the remaining point will be mapped to (0,1)
    #in the write-up this translates to

    #v <--> xstar,ystar <---> (0,0) <---> Psi3
    #w <--> xn,yn       <---> (1,0) <---> Psi1
    #z <--> xnp1,ynp1   <---> (0,1) <---> Psi2
    for i in range(n):
        noden   = OrVert[i]
        xn      = noden[0]
        yn      = noden[1]

        nodenp1 = OrVert[i+1]
        xnp1    = nodenp1[0]
        ynp1    = nodenp1[1]

        #Tn      = (1/6)*( 1/( (xn-xstar)*(ynp1-ystar)-(yn-ystar)*(xnp1-xstar) ) )
        Tn      = (1/6)*( (xn-xstar)*(ynp1-ystar)-(yn-ystar)*(xnp1-xstar)  )

        if i == n-1:
            #LeftmostNode
            H[i,i]   = (psi1(0,0.5)**2+psi1(0.5,0)**2+psi1(0.5,0.5)**2)*Tn + H[i,i] #Integrated against itself
            H[i,0]   = (psi1(0,0.5)*psi2(0,0.5)+psi1(0.5,0)*psi2(0.5,0)+psi1(0.5,0.5)*psi2(0.5,0.5))*Tn + H[i,0] #against the rightmost
            H[i,n]   = (psi1(0,0.5)*psi3(0,0.5)+psi1(0.5,0)*psi3(0.5,0)+psi1(0.5,0.5)*psi3(0.5,0.5))*Tn + H[i,n] #against v*

            #Rightmostnode   
            H[0,i]   = H[i,0] #Against the leftmost
            H[0,0]   = (psi2(0,0.5)**2+psi2(0.5,0)**2+psi2(0.5,0.5)**2)*Tn + H[0,0] #Agaist itself
            H[0,n]   = (psi2(0,0.5)*psi3(0,0.5)+psi2(0.5,0)*psi3(0.5,0)+psi2(0.5,0.5)*psi3(0.5,0.5))*Tn + H[0,n] #against v*

            #v* node
            H[n,i]   = H[i,n] #Against the leftmost
            H[n,0]   = H[0,n] #Against the rightmost
            H[n,n]   = (psi3(0,0.5)**2+psi3(0.5,0)**2+psi3(0.5,0.5)**2)*Tn + H[n,n] #Against itself
        else:
            #LeftmostNode
            H[i,i]   = (psi1(0,0.5)**2+psi1(0.5,0)**2+psi1(0.5,0.5)**2)*Tn + H[i,i] #Integrated agaist itself
            H[i,i+1] = (psi1(0,0.5)*psi2(0,0.5)+psi1(0.5,0)*psi2(0.5,0)+psi1(0.5,0.5)*psi2(0.5,0.5))*Tn + H[i,i+1] #against the rightmost
            H[i,n]   = (psi1(0,0.5)*psi3(0,0.5)+psi1(0.5,0)*psi3(0.5,0)+psi1(0.5,0.5)*psi3(0.5,0.5))*Tn + H[i,n] #against v*

            #Rightmostnode   
            H[i+1,i]   = H[i,i+1] #Against the leftmost
            H[i+1,i+1] = (psi2(0,0.5)**2+psi2(0.5,0)**2+psi2(0.5,0.5)**2)*Tn + H[i+1,i+1] #Agaist itself
            H[i+1,n]   = (psi2(0,0.5)*psi3(0,0.5)+psi2(0.5,0)*psi3(0.5,0)+psi2(0.5,0.5)*psi3(0.5,0.5))*Tn + H[i+1,n] #agaist v*

            #v* node
            H[n,i]   = H[i,n] #Against the leftmost
            H[n,i+1] = H[i+1,n] #Agaist the rightmost
            H[n,n]   = (psi3(0,0.5)**2+psi3(0.5,0)**2+psi3(0.5,0.5)**2)*Tn + H[n,n] #Agaist itself
    #I  = np.eye(n+1)
    #B            = np.eye(n+1)
    #for i in range(n+1):
    #    B[n,i] = 1/n
    w,Areas = convexcoeffs(Element,EdgeNodes,Nodes,Ori)
    B       = np.zeros((n+1,n))
    for i in range(n):
        B[i,i] = 1
        B[n,i] = w[i]
        #B[n,i] = 2*(i+1)/(n*(n+1))

    MV = np.transpose(B).dot(H.dot(B))#+A*np.transpose(I-B).dot(I-B)
    
    NJ = np.zeros((Dim,n))
    
    
    for i in range(Dim):        
        NJ[i,:] = np.transpose( LocprojE(Basis[i],Element,EdgeNodes,Nodes) )
    
    
    NJ = np.transpose(NJ)
    #print(NJ)
    b = np.transpose(NJ).dot(ME)
    #print(b)
    #print(ME)
    #print(NJ)
    #print(np.transpose(NJ).dot(ME).dot(NJ))
    #print(np.linalg.inv( np.transpose(NJ).dot(ME).dot(NJ) ) )
    
    MJ = np.linalg.pinv( np.transpose(NJ).dot(ME).dot(NJ) )
    #print(MJ)
    MJ = MJ.dot(b)
    #print(MJ)



    
    PolyCoordinates = np.zeros((2*(len(Vertices)-1),Dim))
    JMatrix = np.zeros( (len(Vertices)-1,2*(len(Vertices)-1)) ) 
                       
    l = 0
    k = 0
    for Polynomial in Basis:
        for j in range(len(Vertices)-1):
            Vertex = Vertices[j]
            x = Vertex[0]
            y = Vertex[1]
            
            Px,Py = Polynomial(x,y)
            
            PolyCoordinates[2*j,l]   = Px
            PolyCoordinates[2*j+1,l] = Py
            
            
    
            if k == 0:
                Jx,Jy = J(x,y)
                JMatrix[j,2*j] = -Jy
                JMatrix[j,2*j+1] = Jx
                
            j = j+1
        k = 1        
        l = l+1
    
    MJ = JMatrix.dot(PolyCoordinates).dot(MJ)
    MJ = MV.dot(MJ)
    return ME,MV,MJ,Edges

def MFDAssembly(J,Nodes,EdgeNodes,ElementEdges,Orientations):
    #This routine takes a mesh and assembles the global mass matrices and their inverses
    Basis          = [Poly1,Poly2,Poly]
    
    NumberElements = len(ElementEdges)
    NumberEdges    = len(EdgeNodes)
    NumberNodes    = len(Nodes)
    
    
    ME = lil_matrix((NumberEdges,NumberEdges))
    #ME = np.zeros((NumberEdges,NumberEdges))
    #WE=np.zeros((NumberEdges,NumberEdges))
    MV = lil_matrix((NumberNodes,NumberNodes))
    #MV = np.zeros((NumberNodes,NumberNodes))
    #WV=np.zeros((NumberNodes,NumberNodes))
    MJ = lil_matrix((NumberNodes,NumberEdges))
    #MJ = np.zeros((NumberNodes,NumberEdges))
    
    #loop over the elements
    k = 0
    for Element in ElementEdges: 
        #Compute the local mass and stiffness matrices
        #LocME,LocWE,LocMV,LocWV,Edges=LocalMEWEMVWV(Element,EdgeNodes,Nodes) 
        Ori = Orientations[k]
        k = k+1
        LocME,LocMV,LocMJ,Edges = MFDMEWEMVWV(J,Basis,Element,EdgeNodes,Nodes,Ori) 
        
        #The assembly of the edge-based functions is easier since Element
        #Contains the order in which the edges ought to be assembled
       
        for j in range(len(Element)):
            ME[Element[j],Element] = ME[Element[j],Element] + LocME[j]
            #WE[Element[j],Element] = WE[Element[j],Element] + LocWE[j]
        n=len(Edges)-1
        ElementVertices = [0]*n
        
        for i in range(n):
            ElementVertices[i] = Edges[i][0]
         
        for j in range(len(Element)):
            MV[ElementVertices[j],ElementVertices] = MV[ElementVertices[j],ElementVertices]+LocMV[j]
            MJ[ElementVertices[j],Element] = MJ[ElementVertices[j],Element]+LocMJ[j]
            #WV[ElementVertices[j],ElementVertices]=WV[ElementVertices[j],ElementVertices]+LocWV[j]
        
        #print(MJ)
        #i = 0
        #print(np.linalg.norm(LocMJ))
        #for Edge in Element:
        #    MJ[Edge,ElementVertices] = MJ[Edge,ElementVertices]+LocMJ[i]
        #    i = i+1
        #for j in range(len(Element)):
            #MJ[ElementVertices[j],Element] = MJ[ElementVertices[j],Element]+LocMJ[j]
    #return ME,WE,MV,WV
    
    ME = ME.tocsr()
    MV = MV.tocsr()
    MJ = MJ.tocsr()
    
    return ME,MV,MJ

def EAssembly(J,Nodes,EdgeNodes,ElementEdges,Orientations):
    #This routine takes a mesh and assembles the global mass matrices and their inverses
    Basis          = [Poly1,Poly2,Poly]
    
    NumberElements = len(ElementEdges)
    NumberEdges    = len(EdgeNodes)
    NumberNodes    = len(Nodes)
    
    
    ME = lil_matrix((NumberEdges,NumberEdges))
    #ME = np.zeros((NumberEdges,NumberEdges))
    #WE=np.zeros((NumberEdges,NumberEdges))
    MV = lil_matrix((NumberNodes,NumberNodes))
    #MV = np.zeros((NumberNodes,NumberNodes))
    #WV=np.zeros((NumberNodes,NumberNodes))
    MJ = lil_matrix((NumberNodes,NumberEdges))
    #MJ = np.zeros((NumberNodes,NumberEdges))
    
    #loop over the elements
    k = 0
    for Element in ElementEdges: 
        #Compute the local mass and stiffness matrices
        #LocME,LocWE,LocMV,LocWV,Edges=LocalMEWEMVWV(Element,EdgeNodes,Nodes) 
        Ori = Orientations[k]
        k = k+1
        LocME,LocMV,LocMJ,Edges = NewLocalMEWEMVWV(J,Basis,Element,EdgeNodes,Nodes,Ori) 
        
        #The assembly of the edge-based functions is easier since Element
        #Contains the order in which the edges ought to be assembled
       
        for j in range(len(Element)):
            ME[Element[j],Element] = ME[Element[j],Element] + LocME[j]
            #WE[Element[j],Element] = WE[Element[j],Element] + LocWE[j]
        n=len(Edges)-1
        ElementVertices = [0]*n
        
        for i in range(n):
            ElementVertices[i] = Edges[i][0]
         
        for j in range(len(Element)):
            MV[ElementVertices[j],ElementVertices] = MV[ElementVertices[j],ElementVertices]+LocMV[j]
            MJ[ElementVertices[j],Element] = MJ[ElementVertices[j],Element]+LocMJ[j]
            #WV[ElementVertices[j],ElementVertices]=WV[ElementVertices[j],ElementVertices]+LocWV[j]
        
        #print(MJ)
        #i = 0
        #print(np.linalg.norm(LocMJ))
        #for Edge in Element:
        #    MJ[Edge,ElementVertices] = MJ[Edge,ElementVertices]+LocMJ[i]
        #    i = i+1
        #for j in range(len(Element)):
            #MJ[ElementVertices[j],Element] = MJ[ElementVertices[j],Element]+LocMJ[j]
    #return ME,WE,MV,WV
    
    ME = ME.tocsr()
    MV = MV.tocsr()
    MJ = MJ.tocsr()
    
    return ME,MV,MJ

def LSAssembly(J,Nodes,EdgeNodes,ElementEdges,Orientations):
    #This routine takes a mesh and assembles the global mass matrices and their inverses
    
    Basis  = [Poly1,Poly2,Poly]
    NumberElements = len(ElementEdges)
    NumberEdges = len(EdgeNodes)
    NumberNodes = len(Nodes)
    
    
    ME = lil_matrix((NumberEdges,NumberEdges))
    #ME = np.zeros((NumberEdges,NumberEdges))
    #WE=np.zeros((NumberEdges,NumberEdges))
    MV = lil_matrix((NumberNodes,NumberNodes))
    #MV = np.zeros((NumberNodes,NumberNodes))
    #WV=np.zeros((NumberNodes,NumberNodes))
    MJ = lil_matrix((NumberNodes,NumberEdges))
    #MJ = np.zeros((NumberNodes,NumberEdges))
    
    #loop over the elements
    k = 0
    for Element in ElementEdges: 
        #Compute the local mass and stiffness matrices
        #LocME,LocWE,LocMV,LocWV,Edges=LocalMEWEMVWV(Element,EdgeNodes,Nodes) 
        Ori = Orientations[k]
        k = k+1
        LocME,LocMV,LocMJ,Edges = LeastSquaresLocalMEWEMVWV(J,Basis,Element,EdgeNodes,Nodes,Ori) 
        
        #The assembly of the edge-based functions is easier since Element
        #Contains the order in which the edges ought to be assembled
       
        for j in range(len(Element)):
            ME[Element[j],Element] = ME[Element[j],Element] + LocME[j]
            #WE[Element[j],Element] = WE[Element[j],Element] + LocWE[j]
        n=len(Edges)-1
        ElementVertices = [0]*n
        
        for i in range(n):
            ElementVertices[i] = Edges[i][0]
         
        for j in range(len(Element)):
            MV[ElementVertices[j],ElementVertices] = MV[ElementVertices[j],ElementVertices]+LocMV[j]
            MJ[ElementVertices[j],Element] = MJ[ElementVertices[j],Element]+LocMJ[j]
            #WV[ElementVertices[j],ElementVertices]=WV[ElementVertices[j],ElementVertices]+LocWV[j]
        
        #print(MJ)
        #i = 0
        #print(np.linalg.norm(LocMJ))
        #for Edge in Element:
        #    MJ[Edge,ElementVertices] = MJ[Edge,ElementVertices]+LocMJ[i]
        #    i = i+1
        #for j in range(len(Element)):
            #MJ[ElementVertices[j],Element] = MJ[ElementVertices[j],Element]+LocMJ[j]
    #return ME,WE,MV,WV
    
    ME = ME.tocsr()
    MV = MV.tocsr()
    MJ = MJ.tocsr()
    
    return ME,MV,MJ
def GIAssembly(J,Nodes,EdgeNodes,ElementEdges,Orientations):
    #This routine takes a mesh and assembles the global mass matrices and their inverses
    
    Basis  = [Poly1,Poly2,Poly]
    NumberElements = len(ElementEdges)
    NumberEdges = len(EdgeNodes)
    NumberNodes = len(Nodes)
    
    
    ME = lil_matrix((NumberEdges,NumberEdges))
    #ME = np.zeros((NumberEdges,NumberEdges))
    #WE=np.zeros((NumberEdges,NumberEdges))
    MV = lil_matrix((NumberNodes,NumberNodes))
    #MV = np.zeros((NumberNodes,NumberNodes))
    #WV=np.zeros((NumberNodes,NumberNodes))
    MJ = lil_matrix((NumberNodes,NumberEdges))
    #MJ = np.zeros((NumberNodes,NumberEdges))
    
    #loop over the elements
    k = 0
    for Element in ElementEdges: 
        #Compute the local mass and stiffness matrices
        #LocME,LocWE,LocMV,LocWV,Edges=LocalMEWEMVWV(Element,EdgeNodes,Nodes) 
        Ori = Orientations[k]
        k = k+1
        LocME,LocMV,LocMJ,Edges = PieceWiseLocalMEWEMVWV(J,Basis,Element,EdgeNodes,Nodes,Ori) 
        
        #The assembly of the edge-based functions is easier since Element
        #Contains the order in which the edges ought to be assembled
       
        for j in range(len(Element)):
            ME[Element[j],Element] = ME[Element[j],Element] + LocME[j]
            #WE[Element[j],Element] = WE[Element[j],Element] + LocWE[j]
        n=len(Edges)-1
        ElementVertices = [0]*n
        
        for i in range(n):
            ElementVertices[i] = Edges[i][0]
         
        for j in range(len(Element)):
            MV[ElementVertices[j],ElementVertices] = MV[ElementVertices[j],ElementVertices]+LocMV[j]
            MJ[ElementVertices[j],Element] = MJ[ElementVertices[j],Element]+LocMJ[j]
            #WV[ElementVertices[j],ElementVertices]=WV[ElementVertices[j],ElementVertices]+LocWV[j]
        
        #print(MJ)
        #i = 0
        #print(np.linalg.norm(LocMJ))
        #for Edge in Element:
        #    MJ[Edge,ElementVertices] = MJ[Edge,ElementVertices]+LocMJ[i]
        #    i = i+1
        #for j in range(len(Element)):
            #MJ[ElementVertices[j],Element] = MJ[ElementVertices[j],Element]+LocMJ[j]
    #return ME,WE,MV,WV
    
    ME = ME.tocsr()
    MV = MV.tocsr()
    MJ = MJ.tocsr()
    
    return ME,MV,MJ
#Interpolators 

def projV(func,Nodes):
    #This function will, provided a function and a set of nodes, compute the projection onto 
    #the space of node-based functions.
    n    = len(Nodes)
    proj = np.zeros((n,1))
    for i in range(n):
        proj[i] = func(Nodes[i][0],Nodes[i][1])
    return proj

def projE(Func,EdgeNodes,Nodes):
    #This function will, provided a function a set of nodes and edges, compute the 
    #projection onto the space of edge-based functions. The direction of the unit normal
    #will be assumed to be the clockwise rotation of the tangential vector.
    
    N=len(EdgeNodes)
    proj=np.zeros((N,1))
    for i in range(N):
        x1      = Nodes[EdgeNodes[i][0]][0]
        y1      = Nodes[EdgeNodes[i][0]][1]
        x2      = Nodes[EdgeNodes[i][1]][0]
        y2      = Nodes[EdgeNodes[i][1]][1]
        lengthe = math.sqrt((x2-x1)**2+(y2-y1)**2)
        xmid    = 0.5*(x1+x2)
        ymid    = 0.5*(y1+y2)
        etimesnormal = [y2-y1,x1-x2]
        Fx,Fy        = Func(xmid,ymid)
        proj[i]      = (etimesnormal[0]*Fx+etimesnormal[1]*Fy)*lengthe**-1  #midpoint rule
    return proj

def HighOrder3projE(Func,EdgeNodes,Nodes):
    #This function will, provided a function a set of nodes and edges, compute the 
    #projection onto the space of edge-based functions. The direction of the unit normal
    #will be assumed to be the clockwise rotation of the tangential vector.
    
    N=len(EdgeNodes)
    proj=np.zeros((N,1))
    for i in range(N):
        x1      = Nodes[EdgeNodes[i][0]][0]
        y1      = Nodes[EdgeNodes[i][0]][1]
        x2      = Nodes[EdgeNodes[i][1]][0]
        y2      = Nodes[EdgeNodes[i][1]][1]
        lengthe = math.sqrt((x2-x1)**2+(y2-y1)**2)
        xmid    = 0.5*(x1+x2)
        ymid    = 0.5*(y1+y2)
        etimesnormal       = [y2-y1,x1-x2]
        Fxmid,Fymid        = Func(xmid,ymid)
        Fx1,Fy1            = Func(x1,y1)
        Fx2,Fy2            = Func(x2,y2)

        proj[i] = (Fx1+4*Fxmid+Fx2)*etimesnormal[0]
        proj[i] = proj[i] +(Fy1+4*Fymid+Fy2)*etimesnormal[1]
        proj[i] = proj[i]/(6*lengthe)   
    return proj

def HighOrder5projE(Func,EdgeNodes,Nodes):
    #This function will, provided a function a set of nodes and edges, compute the 
    #projection onto the space of edge-based functions. The direction of the unit normal
    #will be assumed to be the clockwise rotation of the tangential vector.
    
    N=len(EdgeNodes)
    proj=np.zeros((N,1))
    for i in range(N):
        x1      = Nodes[EdgeNodes[i][0]][0]
        y1      = Nodes[EdgeNodes[i][0]][1]
        x2      = Nodes[EdgeNodes[i][1]][0]
        y2      = Nodes[EdgeNodes[i][1]][1]
        lengthe = math.sqrt((x2-x1)**2+(y2-y1)**2)
        etimesnormal = [y2-y1,x1-x2]
        pt1     = -math.sqrt(3/7)
        pt3     = math.sqrt(3/7)
        Fx0,Fy0 = Func(x1,y1)
        Fx1,Fy1 = Func((x1*(1-pt1)+x2*(1+pt1))/2,(y1*(1-pt1)+y2*(1+pt1))/2)
        Fx2,Fy2 = Func(0.5*(x1+x2),0.5*(y1+y2))
        Fx3,Fy3 = Func((x1*(1-pt3)+x2*(1+pt3))/2,(y1*(1-pt3)+y2*(1+pt3))/2)
        Fx4,Fy4 = Func(x2,y2)

        proj[i] = ( (1/10)*Fx0+(49/90)*Fx1+(32/45)*Fx2+(49/90)*Fx3+(1/10)*Fx4 )*etimesnormal[0]
        proj[i] = proj[i] +( (1/10)*Fy0+(49/90)*Fy1+(32/45)*Fy2+(49/90)*Fy3+(1/10)*Fy4 )*etimesnormal[1]
        proj[i] = proj[i]/(2*lengthe)   
    return proj
#PrimaryOperators

def HighOrder7projE(Func,EdgeNodes,Nodes):
    #This function will, provided a function a set of nodes and edges, compute the 
    #projection onto the space of edge-based functions. The direction of the unit normal
    #will be assumed to be the clockwise rotation of the tangential vector.
    
    N=len(EdgeNodes)
    proj=np.zeros((N,1))
    for i in range(N):
        x1      = Nodes[EdgeNodes[i][0]][0]
        y1      = Nodes[EdgeNodes[i][0]][1]
        x2      = Nodes[EdgeNodes[i][1]][0]
        y2      = Nodes[EdgeNodes[i][1]][1]
        lengthe = math.sqrt((x2-x1)**2+(y2-y1)**2)
        etimesnormal = [y2-y1,x1-x2]

        pt0 = -1
        w0  = 1/21

        pt1 = -math.sqrt((5/11)+(2/11)*math.sqrt(5/3))
        w1  = (124-7*math.sqrt(15))/350

        pt2 = -math.sqrt((5/11)-(2/11)*math.sqrt(5/3))
        w2  = (124+7*math.sqrt(15))/350

        pt3 = 0
        w3  = 256/525

        pt4 = math.sqrt((5/11)-(2/11)*math.sqrt(5/3))
        w4  = (124+7*math.sqrt(15))/350

        pt5 = math.sqrt((5/11)+(2/11)*math.sqrt(5/3))
        w5  = (124-7*math.sqrt(15))/350

        pt6 = 1
        w6  = 1/21
        
        Fx0,Fy0 = Func(x1,y1)
        Fx1,Fy1 = Func((x1*(1-pt1)+x2*(1+pt1))/2,(y1*(1-pt1)+y2*(1+pt1))/2)
        Fx2,Fy2 = Func((x1*(1-pt2)+x2*(1+pt2))/2,(y1*(1-pt2)+y2*(1+pt2))/2)
        Fx3,Fy3 = Func(0.5*(x1+x2),0.5*(y1+y2))
        Fx4,Fy4 = Func((x1*(1-pt4)+x2*(1+pt4))/2,(y1*(1-pt4)+y2*(1+pt4))/2)
        Fx5,Fy5 = Func((x1*(1-pt5)+x2*(1+pt5))/2,(y1*(1-pt5)+y2*(1+pt5))/2)
        Fx6,Fy6 = Func(x2,y2)

        proj[i] = (w0*Fx0+w1*Fx1+w2*Fx2+w3*Fx3+w4*Fx4+w5*Fx5+w6*Fx6)*etimesnormal[0]
        proj[i] = proj[i] +(w0*Fy0+w1*Fy1+w2*Fy2+w3*Fy3+w4*Fy4+w5*Fy5+w6*Fy6)*etimesnormal[1]
        proj[i] = proj[i]/(2*lengthe)   
    return proj



def primcurl(EdgeNodes,Nodes):
    #This routine computes the primary curl as a matrix
    nN = len(Nodes)
    nE = len(EdgeNodes)
    #curl = np.zeros((nE,nN))
    curl = lil_matrix((nE,nN))
    for i in range(nE):
        Node1 = EdgeNodes[i][0]
        Node2 = EdgeNodes[i][1]
        x1 = Nodes[Node1][0]
        y1 = Nodes[Node1][1]
        x2 = Nodes[Node2][0]
        y2 = Nodes[Node2][1]
        lengthe = math.sqrt((x2-x1)**2+(y2-y1)**2)
        
        curl[i,Node2] = 1/lengthe #These formulas are derived in the pdf document
        curl[i,Node1] = -1/lengthe
    curl = curl.tocsr()
    return curl

def primdiv(ElementEdges,EdgeNodes,Nodes,Orientations):
    #this routine computes the primary divergence matrix
    NEl = len(ElementEdges)
    NE  = len(EdgeNodes)
    div = np.zeros((NEl,NE))
    for i in range(NEl):
        Element = ElementEdges[i]
        N       = len(Element)
        Ori     = Orientations[i]
        for j in range(N):
            Node1 = EdgeNodes[Element[j]][0]
            Node2 = EdgeNodes[Element[j]][1]

            x1    = Nodes[Node1][0]
            y1    = Nodes[Node1][1] #these formulas are derived in the pdf document
            x2    = Nodes[Node2][0]
            y2    = Nodes[Node2][1]

            lengthe                = math.sqrt((x2-x1)**2+(y2-y1)**2)
            xP,yP,A,Vertices,Edges = Centroid(Element,EdgeNodes,Nodes,Ori)
            div[i,Element[j]]      = Ori[j]*lengthe/A
    return div


#Solver
def MFDSolver(J,Nodes,EdgeNodes,ElementEdges,BoundaryNodes,Orientations,EssentialBoundaryCond,InitialCond,ExactE,ExactB,T,dt,theta):
    #This routine will, provided a mesh, final time and time step, return the values of the electric and magnetic field at
    #the given time.
    #The boundary conditions are given above 
    time = np.arange(0,T,dt)
    InternalNodes,NumberInternalNodes = InternalObjects(BoundaryNodes,Nodes)
    ME,MV,MJ = MFDAssembly(J,Nodes,EdgeNodes,ElementEdges,Orientations) #compute the mass matrices
    
    
    #Let us construct the required matrices
    
    curl = primcurl(EdgeNodes,Nodes) #the primary curl
    #D = np.zeros((len(Nodes),len(Nodes))) #this matrix will is explained in the pdf
    D = lil_matrix((len(Nodes),len(Nodes)))
    for i in InternalNodes:
        D[i,i]=1
    D = D.tocsr()
   
    Aprime = MV+theta*dt*( ( np.transpose(curl) ).dot(ME)+MJ ).dot(curl)#MV.dot(MJ) ).dot(curl)
    Aprime = D.dot(Aprime)
    #A = np.zeros((NumberInternalNodes,NumberInternalNodes))
    A = lil_matrix((NumberInternalNodes,NumberInternalNodes))
    

    for i in range(NumberInternalNodes):
        A[i,:] = Aprime[InternalNodes[i],InternalNodes]
    A = A.tocsr()
    
    b = np.transpose(curl).dot(ME)+MJ#+MV.dot(MJ)
    b = D.dot(b)
    
    #Bh = projE(InitialCond,EdgeNodes,Nodes)
    #Bh = HighOrder3projE(InitialCond,EdgeNodes,Nodes)
    #Bh = HighOrder5projE(InitialCond,EdgeNodes,Nodes)
    Bh = HighOrder7projE(InitialCond,EdgeNodes,Nodes)
    Bh = np.transpose(Bh)[0]
  
    Eh = np.zeros(len(Nodes))
    
    
    EhInterior = np.zeros(len(Nodes)) #This is Eh in the interior
    
    EhBoundary = np.zeros(len(Nodes))   #This is Eh on the boundary
   
    for t in time:
        
        #We update the time dependant boundary conditions
        #i.e. The boundary values of the electric field
        for NodeNumber in BoundaryNodes:
            Node = Nodes[NodeNumber]
            EhBoundary[NodeNumber] = EssentialBoundaryCond(Node[0],Node[1],t+theta*dt)
        
        #Solve  for the internal values of the electric field
        
        W1 = b.dot(Bh)
        W2 = Aprime.dot(EhBoundary)
        
        #EhInterior[InternalNodes] = np.linalg.solve(A,W1[InternalNodes]-W2[InternalNodes]) 
        
        #EhInterior[InternalNodes] = spsolve(A,W1[InternalNodes]-W2[InternalNodes]) 
        #f = spsolve(A,W1[InternalNodes]-W2[InternalNodes]) 
        #EhInterior[InternalNodes] = f
        
        EhInterior[InternalNodes] = spsolve(A,W1[InternalNodes]-W2[InternalNodes]) 
        
        Eh = EhInterior+EhBoundary
        
        
        #Update the magnetic field
        Bh = Bh-dt*curl.dot(Eh) 
          
     
    #Now we compute the error
    def ContB(x,y):
        return ExactB(x,y,T)
    def ContE(x,y):
        return ExactE(x,y,T-(1-theta)*dt)
    
    Bex = HighOrder7projE(ContB,EdgeNodes,Nodes)
    Eex = projV(ContE,Nodes)
    
    Bex = np.transpose(Bex)[0]
    Eex = np.transpose(Eex)[0]
    
    B = Bh-Bex
    E = Eh-Eex
    #Eintex = projV(InteriorCond,Nodes)
    #Eintex = np.transpose(Eintex)[0]
    #EI = EhInterior - Eintex
    #InteriorError = (MV.dot(EI).dot(EI))/(MV.dot(Eintex).dot(Eintex))
    #InteriorError = math.sqrt(InteriorError)  
    #MagneticError = np.transpose(B).dot(ME).dot(B)
    #ElectricError = np.transpose(E).dot(MV).dot(E)
    
    MagneticError = (ME.dot(B).dot(B))/(ME.dot(Bex).dot(Bex))
    ElectricError = (MV.dot(E).dot(E))/(MV.dot(Eex).dot(Eex))
    
    #MagneticError = math.sqrt(MagneticError[0,0])
    #ElectricError = math.sqrt(ElectricError[0,0])
    
    MagneticError = math.sqrt(MagneticError)
    ElectricError = math.sqrt(ElectricError)
    
    return Bh,Eh,MagneticError,ElectricError

def ESolver(J,Nodes,EdgeNodes,ElementEdges,BoundaryNodes,Orientations,EssentialBoundaryCond,InitialCond,ExactE,ExactB,T,dt,theta):
    #This routine will, provided a mesh, final time and time step, return the values of the electric and magnetic field at
    #the given time.
    #The boundary conditions are given above 
    time = np.arange(0,T,dt)
    InternalNodes,NumberInternalNodes = InternalObjects(BoundaryNodes,Nodes)
    ME,MV,MJ = EAssembly(J,Nodes,EdgeNodes,ElementEdges,Orientations) #compute the mass matrices
    
    
    #Let us construct the required matrices
    
    curl = primcurl(EdgeNodes,Nodes) #the primary curl
    #D = np.zeros((len(Nodes),len(Nodes))) #this matrix will is explained in the pdf
    D = lil_matrix((len(Nodes),len(Nodes)))
    for i in InternalNodes:
        D[i,i]=1
    D = D.tocsr()
   
    Aprime = MV+theta*dt*( ( np.transpose(curl) ).dot(ME)+MJ ).dot(curl)#MV.dot(MJ) ).dot(curl)
    Aprime = D.dot(Aprime)
    #A = np.zeros((NumberInternalNodes,NumberInternalNodes))
    A = lil_matrix((NumberInternalNodes,NumberInternalNodes))
    

    for i in range(NumberInternalNodes):
        A[i,:] = Aprime[InternalNodes[i],InternalNodes]
    A = A.tocsr()
    
    b = np.transpose(curl).dot(ME)+MJ#+MV.dot(MJ)
    b = D.dot(b)
    
    #Bh = projE(InitialCond,EdgeNodes,Nodes)
    #Bh = HighOrder3projE(InitialCond,EdgeNodes,Nodes)
    #Bh = HighOrder5projE(InitialCond,EdgeNodes,Nodes)
    Bh = HighOrder7projE(InitialCond,EdgeNodes,Nodes)
    Bh = np.transpose(Bh)[0]
  
    Eh = np.zeros(len(Nodes))
    
    
    EhInterior = np.zeros(len(Nodes)) #This is Eh in the interior
    
    EhBoundary = np.zeros(len(Nodes))   #This is Eh on the boundary
   
    
    a = 0
    for t in time:
        
        #We update the time dependant boundary conditions
        #i.e. The boundary values of the electric field
        for NodeNumber in BoundaryNodes:
            Node = Nodes[NodeNumber]
            EhBoundary[NodeNumber] = EssentialBoundaryCond(Node[0],Node[1],t+theta*dt)
        
        #Solve  for the internal values of the electric field
        
        W1 = b.dot(Bh)
        W2 = Aprime.dot(EhBoundary)
        
        #EhInterior[InternalNodes] = np.linalg.solve(A,W1[InternalNodes]-W2[InternalNodes]) 
        
        #EhInterior[InternalNodes] = spsolve(A,W1[InternalNodes]-W2[InternalNodes]) 
        #f = spsolve(A,W1[InternalNodes]-W2[InternalNodes]) 
        #EhInterior[InternalNodes] = f
        
        EhInterior[InternalNodes] = spsolve(A,W1[InternalNodes]-W2[InternalNodes]) 
        
        Eh = EhInterior+EhBoundary
        
        
        #Update the magnetic field
        Bh = Bh-dt*curl.dot(Eh) 
        a = a+1   
    print(a)    
    #Now we compute the error
    def ContB(x,y):
        return ExactB(x,y,T)
    def ContE(x,y):
        return ExactE(x,y,T-(1-theta)*dt)
    
    Bex = HighOrder7projE(ContB,EdgeNodes,Nodes)
    Eex = projV(ContE,Nodes)
    
    Bex = np.transpose(Bex)[0]
    Eex = np.transpose(Eex)[0]
    
    B = Bh-Bex
    E = Eh-Eex
    #Eintex = projV(InteriorCond,Nodes)
    #Eintex = np.transpose(Eintex)[0]
    #EI = EhInterior - Eintex
    #InteriorError = (MV.dot(EI).dot(EI))/(MV.dot(Eintex).dot(Eintex))
    #InteriorError = math.sqrt(InteriorError)  
    #MagneticError = np.transpose(B).dot(ME).dot(B)
    #ElectricError = np.transpose(E).dot(MV).dot(E)
    
    MagneticError = (ME.dot(B).dot(B))/(ME.dot(Bex).dot(Bex))
    ElectricError = (MV.dot(E).dot(E))/(MV.dot(Eex).dot(Eex))
    
    #MagneticError = math.sqrt(MagneticError[0,0])
    #ElectricError = math.sqrt(ElectricError[0,0])
    
    MagneticError = math.sqrt(MagneticError)
    ElectricError = math.sqrt(ElectricError)
    
    return Bh,Eh,MagneticError,ElectricError

def LSSolver(J,Nodes,EdgeNodes,ElementEdges,BoundaryNodes,Orientations,EssentialBoundaryCond,InitialCond,ExactE,ExactB,T,dt,theta):
    #This routine will, provided a mesh, final time and time step, return the values of the electric and magnetic field at
    #the given time.
    #The boundary conditions are given above 
    time = np.arange(0,T,dt)
    InternalNodes,NumberInternalNodes = InternalObjects(BoundaryNodes,Nodes)
    ME,MV,MJ = LSAssembly(J,Nodes,EdgeNodes,ElementEdges,Orientations) #compute the mass matrices
    
    #Let us construct the required matrices
    
    curl = primcurl(EdgeNodes,Nodes) #the primary curl
    #D = np.zeros((len(Nodes),len(Nodes))) #this matrix will is explained in the pdf
    D = lil_matrix((len(Nodes),len(Nodes)))
    for i in InternalNodes:
        D[i,i]=1
    D = D.tocsr()
   
    Aprime = MV+theta*dt*( ( np.transpose(curl) ).dot(ME)+MJ ).dot(curl)#MV.dot(MJ) ).dot(curl)
    Aprime = D.dot(Aprime)
    #A = np.zeros((NumberInternalNodes,NumberInternalNodes))
    A = lil_matrix((NumberInternalNodes,NumberInternalNodes))
    

    for i in range(NumberInternalNodes):
        A[i,:] = Aprime[InternalNodes[i],InternalNodes]
    A = A.tocsr()
    
    b = np.transpose(curl).dot(ME)+MJ#+MV.dot(MJ)
    b = D.dot(b)
    
    #Bh = projE(InitialCond,EdgeNodes,Nodes)
    #Bh = HighOrder3projE(InitialCond,EdgeNodes,Nodes)
    #Bh = HighOrder5projE(InitialCond,EdgeNodes,Nodes)
    Bh = HighOrder7projE(InitialCond,EdgeNodes,Nodes)
    Bh = np.transpose(Bh)[0]
  
    Eh = np.zeros(len(Nodes))
    
    
    EhInterior = np.zeros(len(Nodes)) #This is Eh in the interior
    
    EhBoundary = np.zeros(len(Nodes))   #This is Eh on the boundary
   
    
    a = 0 
    for t in time:
        
        #We update the time dependant boundary conditions
        #i.e. The boundary values of the electric field
        for NodeNumber in BoundaryNodes:
            Node = Nodes[NodeNumber]
            EhBoundary[NodeNumber] = EssentialBoundaryCond(Node[0],Node[1],t+theta*dt)
        
        #Solve  for the internal values of the electric field
        
        W1 = b.dot(Bh)
        W2 = Aprime.dot(EhBoundary)
        
        #EhInterior[InternalNodes] = np.linalg.solve(A,W1[InternalNodes]-W2[InternalNodes]) 
        
        #EhInterior[InternalNodes] = spsolve(A,W1[InternalNodes]-W2[InternalNodes]) 
        #f = spsolve(A,W1[InternalNodes]-W2[InternalNodes]) 
        #EhInterior[InternalNodes] = f
        
        EhInterior[InternalNodes] = spsolve(A,W1[InternalNodes]-W2[InternalNodes]) 
        
        Eh = EhInterior+EhBoundary
        
        
        #Update the magnetic field
        Bh = Bh-dt*curl.dot(Eh) 
        a = a+1   
    print(a)    
    #Now we compute the error
    def ContB(x,y):
        return ExactB(x,y,T)
    def ContE(x,y):
        return ExactE(x,y,T-(1-theta)*dt)
    
    Bex = HighOrder7projE(ContB,EdgeNodes,Nodes)
    Eex = projV(ContE,Nodes)
    
    Bex = np.transpose(Bex)[0]
    Eex = np.transpose(Eex)[0]
    
    B = Bh-Bex
    E = Eh-Eex
    #MagneticError = np.transpose(B).dot(ME).dot(B)
    #ElectricError = np.transpose(E).dot(MV).dot(E)
    
    MagneticError = (ME.dot(B).dot(B))/(ME.dot(Bex).dot(Bex))
    ElectricError = (MV.dot(E).dot(E))/(MV.dot(Eex).dot(Eex))
    
    #MagneticError = math.sqrt(MagneticError[0,0])
    #ElectricError = math.sqrt(ElectricError[0,0])
    
    MagneticError = math.sqrt(MagneticError)
    ElectricError = math.sqrt(ElectricError)
    
    Eintex = projV(InteriorCond,Nodes)
    Eintex = np.transpose(Eintex)[0]
    EI = EhInterior - Eintex
    InteriorError = (MV.dot(EI).dot(EI))/(MV.dot(Eintex).dot(Eintex))
    InteriorError = math.sqrt(InteriorError)  
    return Bh,Eh,EhInterior,MagneticError,ElectricError,InteriorError,t

def GISolver(J,Nodes,EdgeNodes,ElementEdges,BoundaryNodes,Orientations,EssentialBoundaryCond,InitialCond,ExactE,ExactB,T,dt,theta):
    #This routine will, provided a mesh, final time and time step, return the values of the electric and magnetic field at
    #the given time.
    #The boundary conditions are given above 
    time = np.arange(0,T,dt)
    InternalNodes,NumberInternalNodes = InternalObjects(BoundaryNodes,Nodes)
    ME,MV,MJ = GIAssembly(J,Nodes,EdgeNodes,ElementEdges,Orientations) #compute the mass matrices
    
    #Let us construct the required matrices
    
    curl = primcurl(EdgeNodes,Nodes) #the primary curl
    #D = np.zeros((len(Nodes),len(Nodes))) #this matrix will is explained in the pdf
    D = lil_matrix((len(Nodes),len(Nodes)))
    for i in InternalNodes:
        D[i,i]=1
    D = D.tocsr()
   
    Aprime = MV+theta*dt*( ( np.transpose(curl) ).dot(ME)+MJ ).dot(curl)#MV.dot(MJ) ).dot(curl)
    Aprime = D.dot(Aprime)
    #A = np.zeros((NumberInternalNodes,NumberInternalNodes))
    A = lil_matrix((NumberInternalNodes,NumberInternalNodes))
    

    for i in range(NumberInternalNodes):
        A[i,:] = Aprime[InternalNodes[i],InternalNodes]
    A = A.tocsr()
    
    b = np.transpose(curl).dot(ME)+MJ#+MV.dot(MJ)
    b = D.dot(b)
    
    #Bh = projE(InitialCond,EdgeNodes,Nodes)
    #Bh = HighOrder3projE(InitialCond,EdgeNodes,Nodes)
    #Bh = HighOrder5projE(InitialCond,EdgeNodes,Nodes)
    Bh = HighOrder7projE(InitialCond,EdgeNodes,Nodes)
    Bh = np.transpose(Bh)[0]
  
    Eh = np.zeros(len(Nodes))
    
    
    EhInterior = np.zeros(len(Nodes)) #This is Eh in the interior
    
    EhBoundary = np.zeros(len(Nodes))   #This is Eh on the boundary
   
    
    a = 0
    for t in time:
        
        #We update the time dependant boundary conditions
        #i.e. The boundary values of the electric field
        for NodeNumber in BoundaryNodes:
            Node = Nodes[NodeNumber]
            EhBoundary[NodeNumber] = EssentialBoundaryCond(Node[0],Node[1],t+theta*dt)
        
        #Solve  for the internal values of the electric field
        
        W1 = b.dot(Bh)
        W2 = Aprime.dot(EhBoundary)
        
        #EhInterior[InternalNodes] = np.linalg.solve(A,W1[InternalNodes]-W2[InternalNodes]) 
        
        #EhInterior[InternalNodes] = spsolve(A,W1[InternalNodes]-W2[InternalNodes]) 
        #f = spsolve(A,W1[InternalNodes]-W2[InternalNodes]) 
        #EhInterior[InternalNodes] = f
        
        EhInterior[InternalNodes] = spsolve(A,W1[InternalNodes]-W2[InternalNodes]) 
        
        Eh = EhInterior+EhBoundary
        
        
        #Update the magnetic field
        Bh = Bh-dt*curl.dot(Eh) 
        a = a+1   
    print(a)    
    #Now we compute the error
    def ContB(x,y):
        return ExactB(x,y,T)
    def ContE(x,y):
        return ExactE(x,y,T-(1-theta)*dt)
    
    Bex = HighOrder7projE(ContB,EdgeNodes,Nodes)
    Eex = projV(ContE,Nodes)
    
    Bex = np.transpose(Bex)[0]
    Eex = np.transpose(Eex)[0]
    
    B = Bh-Bex
    E = Eh-Eex
    #MagneticError = np.transpose(B).dot(ME).dot(B)
    #ElectricError = np.transpose(E).dot(MV).dot(E)
    
    MagneticError = (ME.dot(B).dot(B))/(ME.dot(Bex).dot(Bex))
    ElectricError = (MV.dot(E).dot(E))/(MV.dot(Eex).dot(Eex))
    
    #MagneticError = math.sqrt(MagneticError[0,0])
    #ElectricError = math.sqrt(ElectricError[0,0])
    
    MagneticError = math.sqrt(MagneticError)
    ElectricError = math.sqrt(ElectricError)
    Eintex = projV(InteriorCond,Nodes)
    Eintex = np.transpose(Eintex)[0]
    EI = EhInterior - Eintex
    InteriorError = (MV.dot(EI).dot(EI))/(MV.dot(Eintex).dot(Eintex)) 
    InteriorError = math.sqrt(InteriorError)   
    return Bh,Eh,EhInterior,MagneticError,ElectricError,InteriorError,t



def SetEnergy(T,dt,theta,Pfile,task):
    #In this function we will compute and save the norms of the electric, magnetic fields
    Nodes,EdgeNodes,ElementEdges,BoundaryNodes,Orientations = ProcessedMesh(Pfile)  
    InternalNodes,NumberInternalNodes = InternalObjects(BoundaryNodes,Nodes)
    if task == 'E':
        ME,MV,MJ = EAssembly(J,Nodes,EdgeNodes,ElementEdges,Orientations) #compute the mass matrices
    elif task == 'LS':
        ME,MV,MJ = LSAssembly(J,Nodes,EdgeNodes,ElementEdges,Orientations)
    elif task == 'GI':
        ME,MV,MJ = GIAssembly(J,Nodes,EdgeNodes,ElementEdges,Orientations)
    
    time = np.arange(0,T,dt)
    N    = len(time)

    curl = primcurl(EdgeNodes,Nodes)
    D    = lil_matrix((len(Nodes),len(Nodes)))
    for i in InternalNodes:
        D[i,i] = 1
    D = D.tocsr()

    Aprime = MV+theta*dt*( ( np.transpose(curl) ).dot(ME)+MJ ).dot(curl)#MV.dot(MJ) ).dot(curl)
    Aprime = D.dot(Aprime)
    A      = lil_matrix((NumberInternalNodes,NumberInternalNodes))

    for i in range(NumberInternalNodes):
        A[i,:] = Aprime[InternalNodes[i],InternalNodes]
    A = A.tocsr()

    b = np.transpose(curl).dot(ME)+MJ#+MV.dot(MJ)
    b = D.dot(b)

    Bh   = HighOrder7projE(InitialCond,EdgeNodes,Nodes)
    Bh   = np.transpose(Bh)[0]
    R1   = Bh.dot(ME.dot(Bh))

    Eh         = np.zeros(len(Nodes))
    EhInterior = np.zeros(len(Nodes))
    EhBoundary = np.zeros(len(Nodes))
    
    L2         = np.arange(0,N-1)*0
    R2         = np.arange(0,N-1)*0
    step       = 0
    for t in time:

        for NodeNumber in BoundaryNodes:
            Node = Nodes[NodeNumber]
            EhBoundary[NodeNumber] = EssentialBoundaryCond(Node[0],Node[1],t+theta*dt)
        
        
        # RHS2 = RHS2+\
        #        (beta**step)*gamma*dt*\
        #        ( EhBoundary.dot(MV.dot(EhBoundary))+\
        #        (curl.dot(EhBoundary)).dot(ME.dot(curl.dot(EhBoundary))) )

        W1 = b.dot(Bh)
        W2 = Aprime.dot(EhBoundary)

        EhInterior[InternalNodes] = spsolve(A,W1[InternalNodes]-W2[InternalNodes])

        Eh = EhInterior+EhBoundary
        #LHS2 = LHS2+0.5*gamma*dt(beta**(step))*Eh.dot(MV.dot(Eh))

        Bh = Bh-dt*curl.dot(Eh)

        if step == N-1:
            L1 =  Bh.dot(ME.dot(Bh))
        else:

            R2[step] = ( EhBoundary.dot(MV.dot(EhBoundary))+\
                       (curl.dot(EhBoundary)).dot(ME.dot(curl.dot(EhBoundary))) )

            L2[step] = EhInterior.dot(MV.dot(EhInterior))

        step = step+1

    En = Energy(theta,dt,N,L1,L2,R1,R2)
    return En

def EnergyPlot(Q,N,dt,theta,Pfile,task):
    #In this function we will compute and save the norms of the electric, magnetic fields
    Nodes,EdgeNodes,ElementEdges,BoundaryNodes,Orientations = ProcessedMesh(Pfile)  
    InternalNodes,NumberInternalNodes = InternalObjects(BoundaryNodes,Nodes)
    if task == 'E':
        ME,MV,MJ = EAssembly(J,Nodes,EdgeNodes,ElementEdges,Orientations) #compute the mass matrices
    elif task == 'LS':
        ME,MV,MJ = LSAssembly(J,Nodes,EdgeNodes,ElementEdges,Orientations)
    elif task == 'GI':
        ME,MV,MJ = GIAssembly(J,Nodes,EdgeNodes,ElementEdges,Orientations)
    T    =  dt*N
    time = np.arange(0,T,dt)
    #N    = len(time)
    MV   = MV/C
    MJ   = MJ/C
    curl = primcurl(EdgeNodes,Nodes)
    D    = lil_matrix((len(Nodes),len(Nodes)))
    for i in InternalNodes:
        D[i,i] = 1
    D = D.tocsr()

    Aprime = MV+theta*dt*( ( np.transpose(curl) ).dot(ME)+MJ ).dot(curl)#MV.dot(MJ) ).dot(curl)
    Aprime = D.dot(Aprime)
    A      = lil_matrix((NumberInternalNodes,NumberInternalNodes))

    for i in range(NumberInternalNodes):
        A[i,:] = Aprime[InternalNodes[i],InternalNodes]
    A = A.tocsr()

    b = np.transpose(curl).dot(ME)+MJ
    b = D.dot(b)

    Bh   = HighOrder7projE(InitialCond,EdgeNodes,Nodes)
    Bh   = np.transpose(Bh)[0]

    Eh         = np.zeros(len(Nodes))
    EhInterior = np.zeros(len(Nodes))
    EhBoundary = np.zeros(len(Nodes))
    
    beta  = (1-Q*theta)/(1+Q*(1-theta))
    gamma = 1/(1-Q*theta)
    step  = 1
    R1    = Bh.dot(ME.dot(Bh))
    L1    = 0
    R2    = 0
    L2    = 0  
    F     = time*0
    R     = time*0 
    L     = time*0
    print(beta)
    print(gamma) 
    for t in time[0:N-1]:

        for NodeNumber in BoundaryNodes:
            Node = Nodes[NodeNumber]
            EhBoundary[NodeNumber] = EssentialBoundaryCond(Node[0],Node[1],t+theta*dt)
        
        
        R2 = R2+\
                (beta**step)*gamma*dt*\
                ( EhBoundary.dot(MV.dot(EhBoundary))+\
                (curl.dot(EhBoundary)).dot(ME.dot(curl.dot(EhBoundary))) )

        W1 = b.dot(Bh)
        W2 = Aprime.dot(EhBoundary)

        EhInterior[InternalNodes] = spsolve(A,W1[InternalNodes]-W2[InternalNodes])
        
        L2 = L2+0.5*gamma*dt*(beta**(step))*EhInterior.dot(MV.dot(EhInterior))
        Eh = EhInterior+EhBoundary
        

        Bh      = Bh-dt*curl.dot(Eh)
        L1      = (beta**step)*Bh.dot(ME.dot(Bh))
        #print(R2)
        #print(L2)
        #print(L1)
        F[step] = (R2+R1-L2-L1)/R1
        R[step] = (R2+R1)/R1
        L[step] = (L2+L1)/R1
        step    = step+1
    #print(F)
    # for i in range(10):
    #     print(L[i])
    # fig, (ax1, ax2,ax3) = plt.subplots(3, sharex=True)
    # ax1.plot(time, R)
    # ax2.plot(time, L)
    # ax3.plot(time, F)
    # plt.show()
    return L,R,F

# def VisualizeB(Bh,ElementEdges):
#     Bx = np.zeros(())
#     By = np.zeros(())
#     for Element in ElementEdges:


def SaveInmFile(name,aname,array):
    #This function will save the array with name aname in a matlab m-file with the name.
    with open("../Matlab/"+name+".m",'w') as file:
        file.writelines(aname+" = [")
        for e in range(len(array)):
            if e == 0:
                file.writelines( str( array[e] ) )
            else:
                file.writelines( ","+str( array[e] ) )
        
        file.writelines('];')

def tderivB(t,Bh,curl,A1,A,B,Nodes,BoundaryNodes,InternalNodes,EhBoundary,EhInterior):
    #this function will compute the time derivative of B at t, Bh
    for NodeNumber in BoundaryNodes:
        Node                   = Nodes[NodeNumber]
        EhBoundary[NodeNumber] = EssentialBoundaryCond(Node[0],Node[1],t)
    
    W1 = B.dot(Bh)
    W2 = A1.dot(EhBoundary)
    EhInterior[InternalNodes] = spsolve(A,W1[InternalNodes]+W2[InternalNodes]) 
    
    Eh = EhInterior+EhBoundary
    return -curl.dot(Eh)

def HartSolver(J,Nodes,EdgeNodes,ElementEdges,BoundaryNodes,Orientations,EssentialBoundaryCond,InitialCond,ExactE,ExactB,T,dt,theta):
    #This routine will, provided a mesh, final time and time step, return the values of the electric and magnetic field at
    #the given time.
    #The boundary conditions are given above 
    time = np.arange(0,T,dt)
    InternalNodes,NumberInternalNodes = InternalObjects(BoundaryNodes,Nodes)
    ME,MV,MJ = EAssembly(J,Nodes,EdgeNodes,ElementEdges,Orientations) #compute the mass matrices
    
    #Let us construct the required matrices
    
    curl = primcurl(EdgeNodes,Nodes) #the primary curl
    #D = np.zeros((len(Nodes),len(Nodes))) #this matrix will is explained in the pdf
    D = lil_matrix((len(Nodes),len(Nodes)))
    for i in InternalNodes:
        D[i,i]=1
    D = D.tocsr()
   
    A1 = D.dot(MV)
    A = lil_matrix((NumberInternalNodes,NumberInternalNodes))
    for i in range(NumberInternalNodes):
        A[i,:] = A1[InternalNodes[i],InternalNodes]
    A = A.tocsr()
    
    B = np.transpose(curl).dot(ME)-MJ
    B = D.dot(B)

    Bh = HighOrder7projE(InitialCond,EdgeNodes,Nodes)
    Bh = np.transpose(Bh)[0]
  
    Eh = np.zeros(len(Nodes))
    
    
    EhInterior = np.zeros(len(Nodes)) #This is Eh in the interior
    
    EhBoundary = np.zeros(len(Nodes))   #This is Eh on the boundary
   
    
    
    for t in time:
        
        k1 = tderivB(t,Bh,curl,A1,A,B,Nodes,BoundaryNodes,InternalNodes,EhBoundary,EhInterior)
        k2 = tderivB(t+0.5*dt,Bh+0.5*k1,curl,A1,A,B,Nodes,BoundaryNodes,InternalNodes,EhBoundary,EhInterior)
        k3 = tderivB(t+0.5*dt,Bh+0.5*k2,curl,A1,A,B,Nodes,BoundaryNodes,InternalNodes,EhBoundary,EhInterior)
        k4 = tderivB(t+dt,Bh+k3,curl,A1,A,B,Nodes,BoundaryNodes,InternalNodes,EhBoundary,EhInterior)
        #Update the magnetic field
        Bh = Bh+(dt/6)*(k1+k2+k3+k4) 
           
        
    #Now we compute the error
    def ContB(x,y):
        return ExactB(x,y,T)
    #def ContE(x,y):
    #    return ExactE(x,y,T-(1-theta)*dt)
    
    Bex = HighOrder7projE(ContB,EdgeNodes,Nodes)
    #Eex = projV(ContE,Nodes)
    
    Bex = np.transpose(Bex)[0]
    #Eex = np.transpose(Eex)[0]
    
    B = Bh-Bex
    #E = Eh-Eex
    #MagneticError = np.transpose(B).dot(ME).dot(B)
    #ElectricError = np.transpose(E).dot(MV).dot(E)
    
    MagneticError = (ME.dot(B).dot(B))/(ME.dot(Bex).dot(Bex))
    #ElectricError = (MV.dot(E).dot(E))/(MV.dot(Eex).dot(Eex))
    
    #MagneticError = math.sqrt(MagneticError[0,0])
    #ElectricError = math.sqrt(ElectricError[0,0])
    
    MagneticError = math.sqrt(MagneticError)
    #ElectricError = math.sqrt(ElectricError)
    
    return Bh,MagneticError