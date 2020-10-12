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
T                = 0.5
#MTypes = ['Trig','Quad','Vor']
#MTypes = ['OnlyOne']
MTypes = ['Small']

def f(xv,t):
    return np.array([0,0])

def h(xv,t):
    return 0

def Inu(xv):
    if abs(xv[1]-1)<1E-8:
        v = 1
    else:
        v = 0
    return np.array([0,v])

def InB(xv):
    return np.array([1,0])

def ub(xv,t):
    if abs(xv[1]-1)<1E-8:
        v = 1
    else:
        v = 0
    return np.array([0,v])

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
        dt = 0.05*dx[i]**2
        PDE    = PDEFullMHD(Mesh,Re,Rm,Inu,InB,dt,theta)
        PDE.SetMHDBCandSource(ub,Eb,f,h)
        Solver = InexactNewtonTimeInt()
        time   = np.arange(0,T,dt)
        for t in time:
            PDE.MHDComputeBC(t)
            PDE.MHDComputeSources(t)
            tempx = Solver.Newtoniter(PDE.MHDG,PDE.MHDConcatenate(PDE.unx,PDE.uny,PDE.umx,PDE.umy,PDE.B,PDE.E,PDE.p),PDE.SetNumMHDDof(),1E-5,50,PDE)
            PDE.unx,PDE.uny,PDE.umx,PDE.umy,PDE.B,PDE.E,PDE.p = PDE.MHDUpdateInt(tempx,PDE.unx,PDE.uny,PDE.umx,PDE.umy,PDE.B,PDE.E,PDE.p)
            PDE.unx,PDE.uny,PDE.umx,PDE.umy,PDE.E             = PDE.MHDUpdateBC(PDE.unx,PDE.uny,PDE.umx,PDE.umy,PDE.E)
        i = i+1
        SaveInmFile('funx','unx',PDE.unx)
        SaveInmFile('funy','uny',PDE.uny)
        x = [u[0] for u in Mesh.Nodes]
        y = [u[1] for u in Mesh.Nodes]
        SaveInmFile('fx','x',x)
        SaveInmFile('fy','y',y)