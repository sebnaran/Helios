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
MTypes = ['Trig','Quad','Vor']
def InB(xv):
    return np.array([0,math.cos(xv[0])])
def Inu(xv):
    return np.array([0.95*xv[1]**2,0.95*xv[0]**2])
def ub(xv):
    return np.array([xv[1]**2,xv[0]**2])

def Eb(xvt):
    return math.cos(xvt[0]+xvt[2])

for MType in MTypes:
    if MType == 'Trig':
        ProcessedFiles = ['PTh=0.101015.txt','PTh=0.051886.txt','PTh=0.0251418.txt']#,'PTh=0.0125255.txt',\
                      #'PTh=0.0062613.txt']

        dx = [0.10101525445522107, 0.05018856132284956, 0.025141822757713456, 0.012525468249897755,\
         0.006261260829309998]
    
    if MType == 'Quad':
        ProcessedFiles = ['PertPQh=0.166666.txt','PertPQh=0.0833333.txt']#,'PertPQh=0.043478.txt',\
                      #'PertPQh=0.021739.txt','PertPQh=0.010989.txt']

        dx = [0.16666666666666666, 0.08333333333333333, 0.043478260869565216, 0.021739130434782608,\
         0.010989010989010988]
    if MType == 'Vor':
        ProcessedFiles = ['PVh=0.128037.txt','PVh=0.0677285.txt','PVh=0.0345033.txt']#,'PVh=0.0174767.txt',\
                      #'PVh=0.0087872.txt']

        dx = [0.12803687993289598, 0.06772854614785964, 0.03450327796711771, 0.017476749542968805,\
        0.008787156237382746]

    i = 0
    for Pfile in ProcessedFiles:
        print(Pfile)
        Nodes,EdgeNodes,ElementEdges,BoundaryNodes,Orientations = ProcessedMesh(Pfile)
        Mesh = HeliosMesh(Nodes,EdgeNodes,ElementEdges,Orientations)
        dt = dx[i]**2
        PDE    = PDEFullMHD(Mesh,Re,Rm,Inu,InB,dt,theta)
        PDE.SetFlowBC(ub)
        Solver = InexactNewtonTimeInt()
        t = 0
        PDE.FlowComputeBC(t)
        print('Computing Solution')
        tempx = Solver.Newtoniter(PDE.FlowG,PDE.FlowConcatenate(),PDE.NumFlowDOF(),1E-1,50)
        PDE.FlowUpdateUnknownDOFs(tempx)
        PDE.unx,PDE.uny,PDE.umx,PDE.umy = PDE.FlowupdateBC(PDE.unx,PDE.uny,PDE.umx,PDE.umy)

        def exactu(xv):
            return np.array([xv[0]**2,xv[1]**2])

        tempun     = PDE.NodalDOFs(exactu,PDE.Mesh.Nodes)
        eunx,euny  = PDE.DecompIntoCoord(tempun)
        tempum     = PDE.NodalDOFs(exactu,PDE.Mesh.MidNodes)
        eumx, eumy = PDE.DecompIntoCoord(tempum)
        errunx = np.array(eunx)-np.array(PDE.unx)
        erruny = np.array(euny)-np.array(PDE.uny)
        errumx = np.array(eumx)-np.array(PDE.umx)
        errumy = np.array(eumy)-np.array(PDE.umy)

        uL2err = PDE.TVhL2Norm(errunx,erruny,errumx,errumy)
        def exactp(xv):
            return 2*xv[0]+2*xv[1]

        exph = PDE.PhDOF(exactp)

        parr = exph-PDE.p
        perr = PDE.PhL2Norm(parr)
        print('pErr = '+str(perr))
        print('uErr = '+str(uL2err))