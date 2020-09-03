from PDEClass import PDEFullMHD
from MeshHelios import HeliosMesh
import numpy as np
import math
#1
def test_init():
    Nodes         = [[-1,-1],[0,-1],[1,-1],[-1,0],[0,0],[1,0],[-1,1],[0,1],[1,1]]
    EdgeNodes     = [[0,1],[4,1],[8,5],[4,7],[7,8],[6,7],[3,6],[0,3],[5,2],[1,2],[3,4],[4,5]]
    ElementEdges  = [[9,8,11,1],[0,1,10,7],[10,3,5,6],[11,2,4,3]]                                     
    NumBoundaryNodes = [0,1,2,3,5,6,7,8] 
    Orientations  = [[1,-1,-1,1],[1,-1,-1,-1],[1,1,-1,-1],[1,-1,-1,-1]]
    TestMesh      = HeliosMesh(Nodes,EdgeNodes,ElementEdges,Orientations)

    def Inu(xv):
        return np.array([1,1])
    def InB(xv):
        return np.array([1,1])
    Re, Rm, dt, theta = 1, 1, 0.5, 0.5
    TestPDE = PDEFullMHD(TestMesh,Re,Rm,Inu,InB,dt,theta)
    assert ( np.all(TestPDE.unx == np.array([1,1,1,1,1,1,1,1,1])) )
#2
def test_concatenate():
    Nodes         = [[-1,-1],[0,-1],[1,-1],[-1,0],[0,0],[1,0],[-1,1],[0,1],[1,1]]
    EdgeNodes     = [[0,1],[4,1],[8,5],[4,7],[7,8],[6,7],[3,6],[0,3],[5,2],[1,2],[3,4],[4,5]]
    ElementEdges  = [[9,8,11,1],[0,1,10,7],[10,3,5,6],[11,2,4,3]]                                     
    NumBoundaryNodes = [0,1,2,3,5,6,7,8] 
    Orientations  = [[1,-1,-1,1],[1,-1,-1,-1],[1,1,-1,-1],[1,-1,-1,-1]]
    TestMesh      = HeliosMesh(Nodes,EdgeNodes,ElementEdges,Orientations)

    def Inu(xv):
        return np.array([1,1])
    def InB(xv):
        return np.array([1,1])
    Re, Rm, dt, theta = 1, 1, 0.5, 0.5
    TestPDE = PDEFullMHD(TestMesh,Re,Rm,Inu,InB,dt,theta)

    x = np.array([1,1,1,1,1,1,1,1,1,1,-1,-1,-1,1,-1,-1,1,1,-1,-1,-1,-1,0,0,0,0,0])
    y = TestPDE.DirichletConcatenate()
    assert(np.all(x==y))
#3
def test_UpdateInterior():
    Nodes         = [[-1,-1],[0,-1],[1,-1],[-1,0],[0,0],[1,0],[-1,1],[0,1],[1,1]]
    EdgeNodes     = [[0,1],[4,1],[8,5],[4,7],[7,8],[6,7],[3,6],[0,3],[5,2],[1,2],[3,4],[4,5]]
    ElementEdges  = [[9,8,11,1],[0,1,10,7],[10,3,5,6],[11,2,4,3]]                                     
    NumBoundaryNodes = [0,1,2,3,5,6,7,8] 
    Orientations  = [[1,-1,-1,1],[1,-1,-1,-1],[1,1,-1,-1],[1,-1,-1,-1]]
    TestMesh      = HeliosMesh(Nodes,EdgeNodes,ElementEdges,Orientations)

    def Inu(xv):
        return np.array([1,1])
    def InB(xv):
        return np.array([1,1])
    Re, Rm, dt, theta = 1, 1, 0.5, 0.5
    TestPDE = PDEFullMHD(TestMesh,Re,Rm,Inu,InB,dt,theta)

    x = np.array([5,6,1,3,4,5,-1,-2,-3,4,-5,-6,7,8,-9,-10,-11,-1,3,1,2,3,4])
    TestPDE.DirichletUpdateInterior(x)
    y = TestPDE.DirichletConcatenate()
    
    assert(np.all(x==y))
#4
def test_Divu1():
    Nodes         = [[-1,-1],[0,-1],[1,-1],[-1,0],[0,0],[1,0],[-1,1],[0,1],[1,1]]
    EdgeNodes     = [[0,1],[4,1],[8,5],[4,7],[7,8],[6,7],[3,6],[0,3],[5,2],[1,2],[3,4],[4,5]]
    ElementEdges  = [[9,8,11,1],[0,1,10,7],[10,3,5,6],[11,2,4,3]]                                     
    NumBoundaryNodes = [0,1,2,3,5,6,7,8] 
    Orientations  = [[1,-1,-1,1],[1,-1,-1,-1],[1,1,-1,-1],[1,-1,-1,-1]]
    TestMesh      = HeliosMesh(Nodes,EdgeNodes,ElementEdges,Orientations)

    def Inu(xv):
        return np.array([1,1])
    def InB(xv):
        return np.array([1,1])
    Re, Rm, dt, theta = 1, 1, 0.5, 0.5
    TestPDE = PDEFullMHD(TestMesh,Re,Rm,Inu,InB,dt,theta)

    testans = TestPDE.DIVu(0,TestPDE.unx,TestPDE.uny,TestPDE.umx,TestPDE.umy)
    assert (abs(testans)<1E-5)

def test_Divu2():
    Nodes         = [[-1,-1],[0,-1],[1,-1],[-1,0],[0,0],[1,0],[-1,1],[0,1],[1,1]]
    EdgeNodes     = [[0,1],[4,1],[8,5],[4,7],[7,8],[6,7],[3,6],[0,3],[5,2],[1,2],[3,4],[4,5]]
    ElementEdges  = [[9,8,11,1],[0,1,10,7],[10,3,5,6],[11,2,4,3]]                                     
    NumBoundaryNodes = [0,1,2,3,5,6,7,8] 
    Orientations  = [[1,-1,-1,1],[1,-1,-1,-1],[1,1,-1,-1],[1,-1,-1,-1]]
    TestMesh      = HeliosMesh(Nodes,EdgeNodes,ElementEdges,Orientations)

    def Inu(xv):
        return np.array([xv[1]**2,xv[0]**2])
    def InB(xv):
        return np.array([1,1])
    Re, Rm, dt, theta = 1, 1, 0.5, 0.5
    TestPDE = PDEFullMHD(TestMesh,Re,Rm,Inu,InB,dt,theta)
    lunx,luny,lumx,lumy = TestPDE.GetLocalTVhDOF(0,TestPDE.unx,TestPDE.uny,TestPDE.umx,TestPDE.umy)
    testans = TestPDE.DIVu(0,lunx,luny,lumx,lumy)
    assert (abs(testans)<1E-5)
#5
def test_EhVhInProd():
    Nodes         = [[-1,-1],[0,-1],[1,-1],[-1,0],[0,0],[1,0],[-1,1],[0,1],[1,1]]
    EdgeNodes     = [[0,1],[4,1],[8,5],[4,7],[7,8],[6,7],[3,6],[0,3],[5,2],[1,2],[3,4],[4,5]]
    ElementEdges  = [[9,8,11,1],[0,1,10,7],[10,3,5,6],[11,2,4,3]]                                     
    NumBoundaryNodes = [0,1,2,3,5,6,7,8] 
    Orientations  = [[1,-1,-1,1],[1,-1,-1,-1],[1,1,-1,-1],[1,-1,-1,-1]]
    TestMesh      = HeliosMesh(Nodes,EdgeNodes,ElementEdges,Orientations)

    def Inu(xv):
        return np.array([1,1])
    def InB(xv):
        return np.array([1,1])
    Re, Rm, dt, theta = 1, 1, 0.5, 0.5
    TestPDE = PDEFullMHD(TestMesh,Re,Rm,Inu,InB,dt,theta)

    LocB = TestPDE.GetLocalEhDOF(0,TestPDE.B)
    A1 = LocB.dot(TestPDE.MEList[0].dot(LocB))

    E    = [1 for i in range(len(TestMesh.Nodes))]
    LocE = TestPDE.GetLocalVhDOF(0,E)
    A2   = LocE.dot( TestPDE.MVList[0].dot(LocE) )

    def Bv(xv):
        return np.array([1,2])
    
    Barr    = TestPDE.MagDOFs(Bv)
    LocBarr = TestPDE.GetLocalEhDOF(0,Barr)
    A3      = LocBarr.dot( TestPDE.MEList[0].dot(LocBarr))

    def Es(xv):
        return xv[0]

    Earr    = TestPDE.NodalDOFs(Es,TestMesh.Nodes)
    LocEarr = TestPDE.GetLocalVhDOF(0,Earr)
    A4      = LocEarr.dot( TestPDE.MVList[0].dot(LocE))

    assert( abs(A1-2)<1E-5 and abs(A2-1)<1E-5 and abs(A3-5)<1E-5 and abs(A4-1/2)<1E-5)

# def test_TVSGandH():
#     Nodes         = [[-1,-1],[0,-1],[1,-1],[-1,0],[0,0],[1,0],[-1,1],[0,1],[1,1]]
#     EdgeNodes     = [[0,1],[4,1],[8,5],[4,7],[7,8],[6,7],[3,6],[0,3],[5,2],[1,2],[3,4],[4,5]]
#     ElementEdges  = [[9,8,11,1],[0,1,10,7],[10,3,5,6],[11,2,4,3]]                                     
#     NumBoundaryNodes = [0,1,2,3,5,6,7,8] 
#     Orientations  = [[1,-1,-1,1],[1,-1,-1,-1],[1,1,-1,-1],[1,-1,-1,-1]]
#     TestMesh      = HeliosMesh(Nodes,EdgeNodes,ElementEdges,Orientations)

#     def Inu(xv):
#         return np.array([xv[0],0])
#     def InB(xv):
#         return np.array([1,1])
#     Re, Rm, dt, theta = 1, 1, 0.5, 0.5

#     TestPDE = PDEFullMHD(TestMesh,Re,Rm,Inu,InB,dt,theta)
#     G      = TestPDE.GISTVList[0]
#     for i in range(12):
#         print('Row'+str(i))
#         print(G[i])
#     assert (1==1)
def testProjector():
    Nodes         = [[-1,-1],[0,-1],[1,-1],[-1,0],[0,0],[1,0],[-1,1],[0,1],[1,1]]
    EdgeNodes     = [[0,1],[4,1],[8,5],[4,7],[7,8],[6,7],[3,6],[0,3],[5,2],[1,2],[3,4],[4,5]]
    ElementEdges  = [[9,8,11,1],[0,1,10,7],[10,3,5,6],[11,2,4,3]]                                     
    NumBoundaryNodes = [0,1,2,3,5,6,7,8] 
    Orientations  = [[1,-1,-1,1],[1,-1,-1,-1],[1,1,-1,-1],[1,-1,-1,-1]]
    TestMesh      = HeliosMesh(Nodes,EdgeNodes,ElementEdges,Orientations)

    def Inu(xv):
        return np.array([xv[1]+xv[0]+1,xv[1]+xv[0]+1])
    def InB(xv):
        return np.array([1,1])
    Re, Rm, dt, theta = 1, 1, 0.5, 0.5
    TestPDE = PDEFullMHD(TestMesh,Re,Rm,Inu,InB,dt,theta)
    GI      = TestPDE.GISTVList[0]
    Element = TestMesh.ElementEdges[0]
    lunx,luny,lumx,lumy = TestPDE.GetLocalTVhDOF(0,TestPDE.unx,TestPDE.uny,TestPDE.umx,TestPDE.umy)
    xP,yP,A,V,E         = TestMesh.Centroid(Element,TestMesh.Orientations[0])
    Bu                  = TestPDE.TVhSemiInProdColumn(Element,0,lunx,luny,lumx,lumy,xP,yP,A,E)
    Ans     = GI.dot(Bu)
    RealAns = np.array([1,1,1,1,1,1,0,0,0,0,0,0])
    diff    = Ans-RealAns
    norm    = 0
    for v in diff:
        norm = norm+v**2
    assert (math.sqrt(norm)<1E-5)

def test_TVhSemiInnerProd1():
    Nodes         = [[-1,-1],[0,-1],[1,-1],[-1,0],[0,0],[1,0],[-1,1],[0,1],[1,1]]
    EdgeNodes     = [[0,1],[4,1],[8,5],[4,7],[7,8],[6,7],[3,6],[0,3],[5,2],[1,2],[3,4],[4,5]]
    ElementEdges  = [[9,8,11,1],[0,1,10,7],[10,3,5,6],[11,2,4,3]]                                     
    NumBoundaryNodes = [0,1,2,3,5,6,7,8] 
    Orientations  = [[1,-1,-1,1],[1,-1,-1,-1],[1,1,-1,-1],[1,-1,-1,-1]]
    TestMesh      = HeliosMesh(Nodes,EdgeNodes,ElementEdges,Orientations)

    def Inu(xv):
        return np.array([xv[0],xv[1]])
    def InB(xv):
        return np.array([1,1])
    Re, Rm, dt, theta = 1, 1, 0.5, 0.5
    TestPDE = PDEFullMHD(TestMesh,Re,Rm,Inu,InB,dt,theta)
    
    lunx,luny,lumx,lumy = TestPDE.GetLocalTVhDOF(0,TestPDE.unx,TestPDE.uny,TestPDE.umx,TestPDE.umy)
    Ans                 = TestPDE.TVhSemiInProd( 0,lunx,luny,lumx,lumy,lunx,luny,lumx,lumy)

    assert (abs(Ans-2)<1E-5)

def test_TVhInProd1():
    Nodes         = [[-1,-1],[0,-1],[1,-1],[-1,0],[0,0],[1,0],[-1,1],[0,1],[1,1]]
    EdgeNodes     = [[0,1],[4,1],[8,5],[4,7],[7,8],[6,7],[3,6],[0,3],[5,2],[1,2],[3,4],[4,5]]
    ElementEdges  = [[9,8,11,1],[0,1,10,7],[10,3,5,6],[11,2,4,3]]                                     
    NumBoundaryNodes = [0,1,2,3,5,6,7,8] 
    Orientations  = [[1,-1,-1,1],[1,-1,-1,-1],[1,1,-1,-1],[1,-1,-1,-1]]
    TestMesh      = HeliosMesh(Nodes,EdgeNodes,ElementEdges,Orientations)

    def Inu(xv):
        return np.array([1,1])
    def InB(xv):
        return np.array([1,1])
    Re, Rm, dt, theta = 1, 1, 0.5, 0.5
    TestPDE = PDEFullMHD(TestMesh,Re,Rm,Inu,InB,dt,theta)
    
    lunx,luny,lumx,lumy = TestPDE.GetLocalTVhDOF(0,TestPDE.unx,TestPDE.uny,TestPDE.umx,TestPDE.umy)
    Ans                 = TestPDE.TVhInProd( 0,lunx,luny,lumx,lumy,lunx,luny,lumx,lumy)

    assert (abs(Ans-2)<1E-5)

def test_Quadrature():
    Nodes         = [[-1,-1],[0,-1],[1,-1],[-1,0],[0,0],[1,0],[-1,1],[0,1],[1,1]]
    EdgeNodes     = [[0,1],[4,1],[8,5],[4,7],[7,8],[6,7],[3,6],[0,3],[5,2],[1,2],[3,4],[4,5]]
    ElementEdges  = [[9,8,11,1],[0,1,10,7],[10,3,5,6],[11,2,4,3]]                                     
    NumBoundaryNodes = [0,1,2,3,5,6,7,8] 
    Orientations  = [[1,-1,-1,1],[1,-1,-1,-1],[1,1,-1,-1],[1,-1,-1,-1]]
    TestMesh      = HeliosMesh(Nodes,EdgeNodes,ElementEdges,Orientations)

    def Inu(xv):
        return np.array([1,1])
    def InB(xv):
        return np.array([1,1])
    Re, Rm, dt, theta = 1, 1, 0.5, 0.5
    TestPDE = PDEFullMHD(TestMesh,Re,Rm,Inu,InB,dt,theta)
    xxxxP,yyyyP,xxyyP,xyyyP,xxxyP = 0,0,0,0,0
    xxxP,yyyP,xyyP,xxyP           = 0,0,0,0
    for u in range(6):
        x,y   = TestPDE.xs[u],TestPDE.ys[u]

        xxxxP = xxxxP+TestPDE.ws[u]*(x**4)
        yyyyP = yyyyP+TestPDE.ws[u]*(y**4)
        xxyyP = xxyyP+TestPDE.ws[u]*(x**2)*(y**2)
        xyyyP = xyyyP+TestPDE.ws[u]*(x)*(y**3)
        xxxyP = xxxyP+TestPDE.ws[u]*(x**3)*(y)
        xxxP  = xxxP+TestPDE.ws[u]*(x**3)
        yyyP  = yyyP+TestPDE.ws[u]*(y**3)
        xyyP  = xyyP+TestPDE.ws[u]*(x)*(y**2)
        xxyP  = xxyP+TestPDE.ws[u]*(x**2)*(y)
    assert (abs(xxxxP-1/30) <1E-5)
    assert (abs(yyyyP-1/30) <1E-5)
    assert (abs(xxyyP-1/180)<1E-5)
    assert (abs(xxxyP-1/120)<1E-5)
    assert (abs(xyyyP-1/120)<1E-5)
    assert (abs(xxxP-1/20)<1E-5)
    assert (abs(yyyP-1/20)<1E-5)
    assert (abs(xxyP-1/60)<1E-5)
    assert (abs(xyyP-1/60)<1E-5)

def test_TVhInProd2():
    Nodes         = [[-1,-1],[0,-1],[1,-1],[-1,0],[0,0],[1,0],[-1,1],[0,1],[1,1]]
    EdgeNodes     = [[0,1],[4,1],[8,5],[4,7],[7,8],[6,7],[3,6],[0,3],[5,2],[1,2],[3,4],[4,5]]
    ElementEdges  = [[9,8,11,1],[0,1,10,7],[10,3,5,6],[11,2,4,3]]                                     
    NumBoundaryNodes = [0,1,2,3,5,6,7,8] 
    Orientations  = [[1,-1,-1,1],[1,-1,-1,-1],[1,1,-1,-1],[1,-1,-1,-1]]
    TestMesh      = HeliosMesh(Nodes,EdgeNodes,ElementEdges,Orientations)

    def Inu(xv):
        return np.array([1,1])
    def InB(xv):
        return np.array([1,1])
    Re, Rm, dt, theta = 1, 1, 0.5, 0.5
    TestPDE = PDEFullMHD(TestMesh,Re,Rm,Inu,InB,dt,theta)
    
    K = TestPDE.KTVList[0]
    assert(abs(K[0,0]-1)    <1E-5) #Area
    assert(abs(K[2,2]-1/3)  <1E-5) #xxP
    assert(abs(K[4,4]-1/3)  <1E-5) #yyP
    assert(abs(K[6,6]-0.2)  <1E-5) #xxxxP
    assert(abs(K[8,8]-0.2)  <1E-5) #yyyyP
    assert(abs(K[10,10]-1/9)<1E-5) #xxyyP

    assert(abs(K[8,10]+1/8) <1E-5) #xyyyP
    assert(abs(K[6,10]+1/8) <1E-5) #xxxyP
    assert(abs(K[5,7]+1/6)  <1E-5) #xxyP
    assert(abs(K[4,10]-1/6) <1E-5) #xyyP
    assert(abs(K[3,7]-1/4)  <1E-5) #xxxP
    assert(abs(K[4,8]+1/4)  <1E-5) #yyyP