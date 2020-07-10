from PDEClass import PDEFullMHD
from MeshHelios import HeliosMesh
import numpy as np
from Solver import InexactNewtonTimeInt
Nodes            = [[-1,-1],[0,-1],[1,-1],[-1,0],[0,0],[1,0],[-1,1],[0,1],[1,1]]
EdgeNodes        = [[0,1],[4,1],[8,5],[4,7],[7,8],[6,7],[3,6],[0,3],[5,2],[1,2],[3,4],[4,5]]
ElementEdges     = [[9,8,11,1],[0,1,10,7],[10,3,5,6],[11,2,4,3]]                                     
NumBoundaryNodes = [0,1,2,3,5,6,7,8] 
Orientations     = [[1,-1,-1,1],[1,-1,-1,-1],[1,1,-1,-1],[1,-1,-1,-1]]
TestMesh         = HeliosMesh(Nodes,EdgeNodes,ElementEdges,NumBoundaryNodes,Orientations)
Re,Rm,dt,theta   = 1, 1, 0.5, 0.5
T                = 10

def Inu(xv):
    return np.array([1,1])
def InB(xv):
    return np.array([1,1])
def f(xvt):
    return np.array([1,1])
def g(xvt):
    return np.array([1,1])
def h(xvt):
    return 1
def ub(xvt):
    return np.array([1,1])
def Eb(xvt):
    return 1
PDE    = PDEFullMHD(TestMesh,Re,Rm,Inu,InB,dt,theta)
PDE.SetConvTestBCAndSource(f,g,h,ub,Eb)
Solver = InexactNewtonTimeInt()
time   = np.arange(0,T,dt)

for t in time:
    tempx = Solver.Newtoniter(PDE.GDirichlet,PDE.DirichletConcatenate(),PDE.NumDirichletDOF(),1E-5,50)
    PDE.DirichletUpdateInterior(tempx)
    PDE.DirichletupdateBC()
