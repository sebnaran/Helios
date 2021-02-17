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
    
def RetrieveAMRMesh(Pfile):
    with open(Pfile, "rb") as fp:   # Unpickling
        N,E,EE,B,O,BT,LR,C = pickle.load(fp)
    return N,E,EE,B,O,BT,LR,C

def SaveInmFile(name,aname,array):
    #This function will save the array with name aname in a matlab m-file with the name.
    with open("../MATLAB/CavityProblem/"+name+".m",'w') as file:
        file.writelines(aname+" = [")
        for e in range(len(array)):
            if e == 0:
                file.writelines( str( array[e] ) )
            else:
                file.writelines( ","+str( array[e] ) )
        
        file.writelines('];')
Re,Rm,theta   = 1, 1, 0.5
#T             = 0.1
#MTypes = ['Trig','Quad','Vor']
#MTypes = ['OnlyOne']
MTypes = ['Small']

def f(xv,t):
    return np.array([0,0])

def Inu(xv):
    if abs(xv[1]-1)<1E-8:
        v = 1
    else:
        v = 0
    return np.array([0,v])

def InB(xv):
    return np.array([1,0])

def B(xv,t):
    return np.array([1,0])

def E(xv,t):
    return 0

def ub(xv,t):
    if abs(xv[1]-1)<1E-8:
        v = 1
    else:
        v = 0
    return np.array([v,0])

def Eb(xv,t):
    return 0

for MType in MTypes:
    if MType =='OnlyOne':
        ProcessedFiles = ['PVh=0.0345033.txt']
        dx = [0.0345033]
    if MType == 'Trig':
        ProcessedFiles = ['PTh=0.101015.txt','PTh=0.051886.txt','PTh=0.0251418.txt']#,'PTh=0.0125255.txt',\
                      #'PTh=0.0062613.txt']

        dx = [0.10101525445522107, 0.05018856132284956, 0.025141822757713456, 0.012525468249897755,\
         0.006261260829309998]
    
    if MType == 'Quad':
        ProcessedFiles = ['PertPQh=0.166666.txt','PertPQh=0.0833333.txt','PertPQh=0.043478.txt']#,\
                      #'PertPQh=0.021739.txt','PertPQh=0.010989.txt']

        dx = [0.16666666666666666, 0.08333333333333333, 0.043478260869565216, 0.021739130434782608,\
         0.010989010989010988]
    if MType == 'Vor':
        ProcessedFiles = ['PVh=0.128037.txt','PVh=0.0677285.txt','PVh=0.0345033.txt']#,'PVh=0.0174767.txt',\
                      #'PVh=0.0087872.txt']

        dx = [0.12803687993289598, 0.06772854614785964, 0.03450327796711771, 0.017476749542968805,\
        0.008787156237382746]
    if MType == 'Small':
        ProcessedFiles = ['PTh=0.2.txt']#,'PTh=0.101015.txt']
        dx = [0.2]#,0.101015]
    i = 0
    for Pfile in ProcessedFiles:
        print(Pfile)

        Nodes,EdgeNodes,ElementEdges,BoundaryNodes,Orientations = ProcessedMesh(Pfile)
        Mesh = HeliosMesh(Nodes,EdgeNodes,ElementEdges,Orientations)
        dt   = 0.001
        PDE    = PDEFullMHD(Mesh,Re,Rm,Inu,InB,dt,theta)
        PDE.MHDSetFlowBCandSource(ub,f)
        PDE.MHDsetElecMagField(B,E)
        PDE.ComputeElecMagDOF(0)
        PDE.MHDFlowupdatef(0)
        Solver = InexactNewtonTimeInt()
        T      = 10*dt
        time   = np.arange(0,T,dt)
        tempx  = PDE.FlowConcatenate()
        for t in time:
            PDE.MHDFlowComputeBC(t)
            PDE.unx,PDE.uny,PDE.umx,PDE.umy       = PDE.FlowUpdateBC(PDE.unx,PDE.uny,PDE.umx,PDE.umy)
            print('unx='+str(PDE.unx))
            print('uny='+str(PDE.uny))

            tempx = Solver.FlowSolve(PDE.FlowG,tempx,PDE.NumFlowDOF(),50,1E-5)
            PDE.unx,PDE.uny,PDE.umx,PDE.umy,PDE.p = PDE.FlowUpdateInt(tempx,PDE.unx,PDE.uny,PDE.umx,PDE.umy,PDE.p)

        i = i+1
        SaveInmFile('funx','unx',PDE.unx)
        SaveInmFile('funy','uny',PDE.uny)
        x = [pt[0] for pt in Mesh.Nodes]
        y = [pt[1] for pt in Mesh.Nodes]
        SaveInmFile('fx','x',x)
        SaveInmFile('fy','y',y)