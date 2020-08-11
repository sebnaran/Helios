from PDEClass import PDEFullMHD
from MeshHelios import HeliosMesh
import numpy as np
#1
def test_init():
    Nodes         = [[-1,-1],[0,-1],[1,-1],[-1,0],[0,0],[1,0],[-1,1],[0,1],[1,1]]
    EdgeNodes     = [[0,1],[4,1],[8,5],[4,7],[7,8],[6,7],[3,6],[0,3],[5,2],[1,2],[3,4],[4,5]]
    ElementEdges  = [[9,8,11,1],[0,1,10,7],[10,3,5,6],[11,2,4,3]]                                     
    NumBoundaryNodes = [0,1,2,3,5,6,7,8] 
    Orientations  = [[1,-1,-1,1],[1,-1,-1,-1],[1,1,-1,-1],[1,-1,-1,-1]]
    TestMesh      = HeliosMesh(Nodes,EdgeNodes,ElementEdges,NumBoundaryNodes,Orientations)

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
    TestMesh      = HeliosMesh(Nodes,EdgeNodes,ElementEdges,NumBoundaryNodes,Orientations)

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
    TestMesh      = HeliosMesh(Nodes,EdgeNodes,ElementEdges,NumBoundaryNodes,Orientations)

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
def test_Divu():
    Nodes         = [[-1,-1],[0,-1],[1,-1],[-1,0],[0,0],[1,0],[-1,1],[0,1],[1,1]]
    EdgeNodes     = [[0,1],[4,1],[8,5],[4,7],[7,8],[6,7],[3,6],[0,3],[5,2],[1,2],[3,4],[4,5]]
    ElementEdges  = [[9,8,11,1],[0,1,10,7],[10,3,5,6],[11,2,4,3]]                                     
    NumBoundaryNodes = [0,1,2,3,5,6,7,8] 
    Orientations  = [[1,-1,-1,1],[1,-1,-1,-1],[1,1,-1,-1],[1,-1,-1,-1]]
    TestMesh      = HeliosMesh(Nodes,EdgeNodes,ElementEdges,NumBoundaryNodes,Orientations)

    def Inu(xv):
        return np.array([1,1])
    def InB(xv):
        return np.array([1,1])
    Re, Rm, dt, theta = 1, 1, 0.5, 0.5
    TestPDE = PDEFullMHD(TestMesh,Re,Rm,Inu,InB,dt,theta)

    testans = TestPDE.DIVu(0,TestPDE.unx,TestPDE.uny,TestPDE.umx,TestPDE.umy)
    assert (abs(testans)<1E-5)
#5
def test_EhVhInProd():
    Nodes         = [[-1,-1],[0,-1],[1,-1],[-1,0],[0,0],[1,0],[-1,1],[0,1],[1,1]]
    EdgeNodes     = [[0,1],[4,1],[8,5],[4,7],[7,8],[6,7],[3,6],[0,3],[5,2],[1,2],[3,4],[4,5]]
    ElementEdges  = [[9,8,11,1],[0,1,10,7],[10,3,5,6],[11,2,4,3]]                                     
    NumBoundaryNodes = [0,1,2,3,5,6,7,8] 
    Orientations  = [[1,-1,-1,1],[1,-1,-1,-1],[1,1,-1,-1],[1,-1,-1,-1]]
    TestMesh      = HeliosMesh(Nodes,EdgeNodes,ElementEdges,NumBoundaryNodes,Orientations)

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
# def testProjector():
#     Nodes         = [[-1,-1],[0,-1],[1,-1],[-1,0],[0,0],[1,0],[-1,1],[0,1],[1,1]]
#     EdgeNodes     = [[0,1],[4,1],[8,5],[4,7],[7,8],[6,7],[3,6],[0,3],[5,2],[1,2],[3,4],[4,5]]
#     ElementEdges  = [[9,8,11,1],[0,1,10,7],[10,3,5,6],[11,2,4,3]]                                     
#     NumBoundaryNodes = [0,1,2,3,5,6,7,8] 
#     Orientations  = [[1,-1,-1,1],[1,-1,-1,-1],[1,1,-1,-1],[1,-1,-1,-1]]
#     TestMesh      = HeliosMesh(Nodes,EdgeNodes,ElementEdges,NumBoundaryNodes,Orientations)

#     def Inu(xv):
#         return np.array([xv[0],0])
#     def InB(xv):
#         return np.array([1,1])
#     Re, Rm, dt, theta = 1, 1, 0.5, 0.5
#     TestPDE = PDEFullMHD(TestMesh,Re,Rm,Inu,InB,dt,theta)
#     GI      = TestPDE.GISTVList[0]
#     Element = TestMesh.ElementEdges[0]
#     lunx,luny,lumx,lumy = TestPDE.GetLocalTVhDOF(0,TestPDE.unx,TestPDE.uny,TestPDE.umx,TestPDE.umy)
#     xP,yP,A,V,E         = TestMesh.Centroid(Element,TestMesh.Orientations[0])
#     Bu                  = TestPDE.TVhSemiInProdColumn(Element,0,lunx,luny,lumx,lumy,xP,yP,A,E)
#     print(Bu)
#     print(GI.dot(Bu))
#     print(lunx)
#     print(luny)
#     print(lumx)
#     print(lumy)
#     assert (1==1)
#6
# def test_TVhSemiInnerProd():
#     Nodes         = [[-1,-1],[0,-1],[1,-1],[-1,0],[0,0],[1,0],[-1,1],[0,1],[1,1]]
#     EdgeNodes     = [[0,1],[4,1],[8,5],[4,7],[7,8],[6,7],[3,6],[0,3],[5,2],[1,2],[3,4],[4,5]]
#     ElementEdges  = [[9,8,11,1],[0,1,10,7],[10,3,5,6],[11,2,4,3]]                                     
#     NumBoundaryNodes = [0,1,2,3,5,6,7,8] 
#     Orientations  = [[1,-1,-1,1],[1,-1,-1,-1],[1,1,-1,-1],[1,-1,-1,-1]]
#     TestMesh      = HeliosMesh(Nodes,EdgeNodes,ElementEdges,NumBoundaryNodes,Orientations)

#     def Inu(xv):
#         return np.array([xv[0],xv[1]])
#     def InB(xv):
#         return np.array([1,1])
#     Re, Rm, dt, theta = 1, 1, 0.5, 0.5
#     TestPDE = PDEFullMHD(TestMesh,Re,Rm,Inu,InB,dt,theta)
    
#     lunx,luny,lumx,lumy = TestPDE.GetLocalTVhDOF(0,TestPDE.unx,TestPDE.uny,TestPDE.umx,TestPDE.umy)
#     ap,st               = TestPDE.TVhSemiInProd( 0,lunx,luny,lumx,lumy,lunx,luny,lumx,lumy)
#     print(ap)
#     print(st)
#     assert(1 == 1)
#     #assert (abs(Ans-2)<1E-5)