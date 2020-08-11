from PDEClass import PDEFullMHD
from MeshHelios import HeliosMesh
import numpy as np
import math
from Solver import InexactNewtonTimeInt
import pickle

def ProcessedMesh(Pfile):
    with open(Pfile, "rb") as fp:   # Unpickling
        N,E,EE,B,O = pickle.load(fp)
    return N,E,EE,B,O
    
Re,Rm,theta   = 1, 1, 0.5
T                = 0.25
MTypes = ['Quad','Vor']
def InB(xv):
    return np.array([0,math.cos(xv[0])])
def Inu(xv):
    return np.array([1,1])
def h(xvt):
    return math.cos(xvt[2]+xvt[0])+math.sin(xvt[2]+xvt[0])
def Eb(xvt):
    return math.cos(xvt[0]+xvt[2])

for MType in MTypes:
    if MType == 'Trig':
        ProcessedFiles = ['PTh=0.101015.txt','PTh=0.051886.txt','PTh=0.0251418.txt','PTh=0.0125255.txt',\
                      'PTh=0.0062613.txt']

        dx = [0.10101525445522107, 0.05018856132284956, 0.025141822757713456, 0.012525468249897755,\
         0.006261260829309998]
    
    if MType == 'Quad':
        ProcessedFiles = ['PertPQh=0.166666.txt','PertPQh=0.0833333.txt','PertPQh=0.043478.txt',\
                      'PertPQh=0.021739.txt','PertPQh=0.010989.txt']

        dx = [0.16666666666666666, 0.08333333333333333, 0.043478260869565216, 0.021739130434782608,\
         0.010989010989010988]
    if MType == 'Vor':
        ProcessedFiles = ['PVh=0.128037.txt','PVh=0.0677285.txt','PVh=0.0345033.txt','PVh=0.0174767.txt',\
                      'PVh=0.0087872.txt']

        dx = [0.12803687993289598, 0.06772854614785964, 0.03450327796711771, 0.017476749542968805,\
        0.008787156237382746]

    i = 0
    for Pfile in ProcessedFiles:
        print(Pfile)
        Nodes,EdgeNodes,ElementEdges,BoundaryNodes,Orientations = ProcessedMesh(Pfile)
        Mesh = HeliosMesh(Nodes,EdgeNodes,ElementEdges,Orientations)
        dt = dx[i]**2
        PDE    = PDEFullMHD(Mesh,Re,Rm,Inu,InB,dt,theta)
        PDE.SetElectroBCAndSource(h,Eb)
        Solver = InexactNewtonTimeInt()
        time   = np.arange(0,T,dt)
        for t in time:
            PDE.ElectroComputeBC(t)
            PDE.Electroupdateh(t)
            tempx = Solver.Newtoniter(PDE.ElectroG,PDE.ElectroConcatenate(),PDE.NumElectroDOF(),1E-5,50)
            PDE.ElectroUpdateUnknownDOFs(tempx)
            PDE.E = PDE.ElectroupdateBC(PDE.E)

        def exactB(xv):
            Bx = 0
            By = math.cos(xv[0]+T)
            return np.array([Bx,By])

        def exactE(xv):
            return math.cos(xv[0]+T-dt)

        Earr = PDE.NodalDOFs(exactE,PDE.Mesh.Nodes)
        Eerr = Earr-PDE.E
        L2E  = PDE.VhL2Norm(Eerr)

        Barr = PDE.MagDOFs(exactB)
        Berr = Barr-PDE.B
        L2B  = PDE.EhL2Norm(Berr) 

        print('ElectricErr = '+str(L2E))
        print('MagneticErr = '+str(L2B))
        i = i+1