from PDEClass import PDEFullMHD
from MeshHelios import HeliosMesh
import numpy as np
import math
from Solver import InexactNewtonTimeInt
import pickle
from numpy.linalg import norm as n2

def ProcessedMesh(Pfile):
    with open(Pfile, "rb") as fp:   # Unpickling
        N,E,EE,B,O = pickle.load(fp)
    return N,E,EE,B,O
    
Re,Rm,theta   = 1, 1, 0.5
T                = 0.001
MTypes = ['Trig','Quad','Vor']
#MTypes = ['Small']
def exactu(xv,t):
    return np.array([math.exp(t)*math.cos(xv[1]),0])
def exactB(xv,t):
    return np.array([0,math.cos(xv[0]+t)])
def exactE(xv,t):
    return math.cos(xv[0]+t)
def exactp(xv,t):
    return -xv[0]*math.cos(xv[1])
#WithElectromagnetics
def f(xv,t):
    y1 = math.cos(t+xv[0])**2-math.cos(xv[1])+0.5*math.cos(t+xv[0])\
         +0.5*math.exp(1/1000)*math.cos(y)*math.cos(xv[0]+t)**2
    y2 = xv[0]*math.sin(xv[1])
    return np.array([y1,y2])
def Inu(xv):
    return exactu(xv,0)
def InB(xv):
    return exactu(xv,0)
for MType in MTypes:
    if MType == 'Trig':
        ProcessedFiles = ['PTh=0.2.txt','PTh=0.101015.txt','PTh=0.051886.txt']#,'PTh=0.0251418.txt']#,'PTh=0.0125255.txt',\
                      #'PTh=0.0062613.txt']

        dx = [0.2,0.10101525445522107, 0.05018856132284956, 0.025141822757713456, 0.012525468249897755,\
         0.006261260829309998]
    
    if MType == 'Quad':
        ProcessedFiles = ['PertPQh=0.166666.txt','PertPQh=0.0833333.txt','PertPQh=0.043478.txt']#,\
                      #'PertPQh=0.021739.txt','PertPQh=0.010989.txt']

        dx = [0.16666666666666666, 0.08333333333333333, 0.043478260869565216, 0.021739130434782608,\
         0.010989010989010988]
    if MType == 'Vor':
        ProcessedFiles = ['PVh=0.333333.txt','PVh=0.128037.txt','PVh=0.0677285.txt']#,'PVh=0.0345033.txt']#,'PVh=0.0174767.txt',\
                      #'PVh=0.0087872.txt']

        dx = [0.3333333333,0.12803687993289598, 0.06772854614785964, 0.03450327796711771, 0.017476749542968805,\
        0.008787156237382746]
    if MType == 'Small':
        ProcessedFiles = ['PTh=0.408248.txt','PTh=0.2.txt','PTh=0.101015.txt']
        #ProcessedFiles = ['PTh=0.2.txt','PTh=0.101015.txt']
        dx = [0.408024,0.2,0.101015]
        #dx = [0.408248,0.2,0.101015]
    #i = 0
    for Pfile in ProcessedFiles:
        print(Pfile)
        Nodes,EdgeNodes,ElementEdges,BoundaryNodes,Orientations = ProcessedMesh(Pfile)
        Mesh = HeliosMesh(Nodes,EdgeNodes,ElementEdges,Orientations)
        #dt   = 0.05*dx[i]**2
        #i    = i+1
        #print(dt)
        dt = 0.001
        PDE    = PDEFullMHD(Mesh,Re,Rm,Inu,InB,dt,theta)
        PDE.SetFlowBCandSource(exactu,f)
        Solver = InexactNewtonTimeInt()
        time = [0]
        T    = dt
        #time   = np.arange(0,T,dt)
        print(2*len(Mesh.NumInternalNodes)+2*len(Mesh.NumInternalMidNodes)+len(Mesh.ElementEdges)-1)
        #print(f'IntNodes={len(Mesh.NumInternalNodes)},IntMid={len(Mesh.NumInternalMidNodes)},Cells={len(Mesh.ElementEdges)}')
        #print(len(time))
        tempx = PDE.FlowConcatenate()
        for t in time:
            PDE.FlowComputeBC(t)
            PDE.Flowupdatef(t)
            tempx = Solver.FlowSolve(PDE.FlowG,tempx,PDE.NumFlowDOF(),50,1E-5)
            #print(f'tempx[0]={tempx[0]}')
            #print(f'G(tempx)={PDE.FlowG(tempx)}')
            #print(f'n2G(tempx)={n2(PDE.FlowG(tempx))}')
            PDE.unx,PDE.uny,PDE.umx,PDE.umy,PDE.p = PDE.FlowUpdateInt(tempx,PDE.unx,PDE.uny,PDE.umx,PDE.umy,PDE.p)
            #print(f'asasan2G(tempx)={n2(PDE.FlowG(tempx))}')
            PDE.unx,PDE.uny,PDE.umx,PDE.umy       = PDE.FlowUpdateBC(PDE.unx,PDE.uny,PDE.umx,PDE.umy)
            #print(f'PDE.unx[j=0]={PDE.unx[Mesh.NumInternalNodes[0]]}')
            j = 0 
            tempx2 = PDE.FlowConcatenate()
            intN   = len(Mesh.NumInternalNodes)
            intMN  = len(Mesh.NumInternalMidNodes)
            #for i in range(len(tempx)):
            #    print(tempx[i])
            #    print(tempx2[i])
            #    if abs(tempx[i]-tempx2[i])>1e-5:
            #        print(i)
            #for i in range(len(Mesh.ElementEdges)):
            #    if abs(PDE.p-tempx[2])
            #for i in range(len(tempx)):
            #    if abs(tempx[i]-tempx2[i])>1e-5:
            #        print(f'i={i}')
            #        print(f'tempx[i]={tempx[i]}')
            #        print(f'tempx2[i]={tempx2[i]}')
            #print(f'G(tempx2)={PDE.FlowG(tempx2)}')
            #print(f'n2G(tempx2)={n2(PDE.FlowG(tempx2))}')
            #print(f'n2G(tempx)={n2(PDE.FlowG(tempx))}')
        
        def exu(xv):
            return exactu(xv,time[len(time)-1])
        def exp(xv):
            return exactp(xv,0.5*time[len(time)-1]+0.5*time[len(time)-2])
        tempun     = PDE.NodalDOFs(exu,PDE.Mesh.Nodes)
        eunx,euny  = PDE.DecompIntoCoord(tempun)
        tempum     = PDE.NodalDOFs(exu,PDE.Mesh.MidNodes)
        eumx, eumy = PDE.DecompIntoCoord(tempum)
        errunx = np.array(eunx)-np.array(PDE.unx)
        erruny = np.array(euny)-np.array(PDE.uny)
        errumx = np.array(eumx)-np.array(PDE.umx)
        errumy = np.array(eumy)-np.array(PDE.umy)

        uL2err = PDE.TVhL2Norm(errunx,erruny,errumx,errumy)
        exph = PDE.PhDOF(exp)

        parr = exph-PDE.p
        perr = PDE.PhL2Norm(parr)
        print('pErr = '+str(perr))
        print('uErr = '+str(uL2err))