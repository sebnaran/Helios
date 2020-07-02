from PDEClass import PDEFullMHD
from MeshHelios import HeliosMesh
import numpy as np
def test_init():
    Nodes         = [[-1,-1],[0,-1],[1,-1],[-1,0],[0,0],[1,0],[-1,1],[0,1],[1,1]]
    EdgeNodes     = [[0,1],[4,1],[8,5],[4,7],[7,8],[6,7],[3,6],[0,3],[5,2],[1,2],[3,4],[4,5]]
    ElementEdges  = [[9,8,11,1],[0,1,10,7],[10,3,5,6],[11,2,4,3]]                                     
    BoundaryNodes = [0,1,2,3,5,6,7,8] 
    Orientations  = [[1,-1,-1,1],[1,-1,-1,-1],[1,1,-1,-1],[1,-1,-1,-1]]
    TestMesh      = HeliosMesh(Nodes,EdgeNodes,ElementEdges,BoundaryNodes,Orientations)
    def Inu(x,y):
        return 1,1
    def InB(x,y):
        return 1,1
    Re = 1
    Rm = 1
    dt = 0.5
    testPDE = PDEFullMHD(TestMesh,Re,Rm,Inu,InB,dt)
    assert (np.all(testPDE.uxb == np.array([1,1,1,1,1,1,1,1])))