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

    testans = TestPDE.DIVu(0)
    assert (abs(testans)<1E-5)
#5
def test_TVhSemiInnerProd():
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

    B = TestPDE.TVhSemiInProd(1)
    #print("mat\n")
    #print(B)
    assert (1==1)