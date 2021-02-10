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
    #return np.array([0.95*xv[1]**2,0.95*xv[0]**2])
    return np.array([0,0])
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
        PDE.nSetFlowBC(ub)
        Solver = InexactNewtonTimeInt()
        t = 0
        PDE.nFlowComputeBC(t)
        #tempx = PDE.nFlowConcatenate()
        #for i in range(len(tempx)):
        #    tempx[i] = i+1
        #tempx[0] = 0.024031434275685715
        tempx = Solver.FlowSolve(PDE.nFlowG,PDE.nFlowConcatenate(),PDE.nNumFlowDOF(),50,1E-5)
        PDE.unx,PDE.uny,PDE.umx,PDE.umy,PDE.p = PDE.nFlowUpdateUnknownDOFs(tempx,PDE.unx,PDE.uny,PDE.umx,PDE.umy,PDE.p)
        #print(f'G(tempx)={PDE.nFlowG(tempx)}')
        tempx2 = PDE.nFlowConcatenate()
        #print(f'G(tempx2)={PDE.nFlowG(tempx2)}')
        #print(f'lentempx ={len(tempx)}')
        #print(f'lentempx2={len(tempx2)}')
        for k in range(len(tempx)):
            if abs(tempx[k]-tempx2[k])>1e-5:
                print(k)
        print(f'numintnodes={len(Mesh.NumInternalNodes)}')
        print(f'Els={len(Mesh.ElementEdges)}')
        print(f'nummidintnodes={len(Mesh.NumInternalMidNodes)}')

        #print(f'tempx[0]={tempx[0]}')
        #print(f'tempx1={tempx}')
        #print(tempx)
        #PDE.unx,PDE.uny,PDE.umx,PDE.umy       = PDE.nFlowupdateBC(PDE.unx,PDE.PDE.uny,PDE.umx,PDE.umy)
        #
        #print(f'AfterFuncunx[j=0]={PDE.unx[Mesh.NumInternalNodes[0]]}')
        #print(f'tempx2={tempx}')
        #print(f'tempx[0]={tempx[0]}')
        #print(f'unx[int]={PDE.unx[Mesh.NumInternalNodes[0]]}')

       #if np.all(tempx,PDE.nFlowConcatenate()):
        t = 0
        PDE.unx,PDE.uny,PDE.umx,PDE.umy = PDE.nFlowupdateBC(PDE.unx,PDE.uny,PDE.umx,PDE.umy)

        def exactu(xv):
            return np.array([xv[1]**2,xv[0]**2])

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
        i = i+1 