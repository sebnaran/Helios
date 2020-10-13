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
Re,Rm,theta   = 1, 1, 0.5
T                = 0.005
#MTypes = ['Trig','Quad','Vor']
#MTypes = ['OnlyOne']
MTypes = ['Small']
def exactu(xv,t):
    return np.array([math.exp(t)*math.cos(xv[1]),0])
def exactB(xv,t):
    return np.array([0,math.cos(xv[0]+t)])
def exactE(xv,t):
    return math.cos(xv[0]+t)
def exactp(xv,t):
    return -xv[0]*math.cos(xv[1])

uL = []
BL = []
EL = []
pL = []
def f(xv,t):
    y1 = math.cos(t+xv[0])**2-math.cos(xv[1])+2*math.exp(t)*math.cos(xv[1])\
         +math.exp(t)*math.cos(xv[1])*math.cos(t+xv[0])**2
    y2 = xv[0]*math.sin(xv[1])
    return np.array([y1,y2])
def h(xv,t):
    return math.cos(t+xv[0])*(1+math.exp(t)*math.cos(xv[1]))+math.sin(t+xv[0])

def Inu(xv):
    return exactu(xv,0)
def InB(xv):
    return exactB(xv,0)

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
        ProcessedFiles = ['PTh=0.408248.txt','PTh=0.2.txt']#,'PTh=0.101015.txt']
        dx = [0.408248,0.2]#,0.101015]
    i = 0
    for Pfile in ProcessedFiles:
        print(Pfile)
        uL.append(Pfile)
        BL.append(Pfile)
        EL.append(Pfile)
        pL.append(Pfile)
        Nodes,EdgeNodes,ElementEdges,BoundaryNodes,Orientations = ProcessedMesh(Pfile)
        Mesh = HeliosMesh(Nodes,EdgeNodes,ElementEdges,Orientations)
        dt = 0.05*dx[i]**2
        PDE    = PDEFullMHD(Mesh,Re,Rm,Inu,InB,dt,theta)
        PDE.SetMHDBCandSource(exactu,exactE,f,h)
        Solver = InexactNewtonTimeInt()
        time   = np.arange(0,T,dt)
        for t in time:
            PDE.MHDComputeBC(t)
            PDE.MHDComputeSources(t)
            tempx = Solver.Newtoniter(PDE.MHDG,PDE.MHDConcatenate(PDE.unx,PDE.uny,PDE.umx,PDE.umy,PDE.B,PDE.E,PDE.p),PDE.SetNumMHDDof(),1E-5,50,PDE)
            PDE.unx,PDE.uny,PDE.umx,PDE.umy,PDE.B,PDE.E,PDE.p = PDE.MHDUpdateInt(tempx,PDE.unx,PDE.uny,PDE.umx,PDE.umy,PDE.B,PDE.E,PDE.p)
            PDE.unx,PDE.uny,PDE.umx,PDE.umy,PDE.E             = PDE.MHDUpdateBC(PDE.unx,PDE.uny,PDE.umx,PDE.umy,PDE.E)
            

        def exu(xv):
            return exactu(xv,T)

        def exB(xv):
            return exactB(xv,T)

        def exE(xv):
            return exactE(xv,T-dt*theta)

        def exp(xv):
            return exactp(xv,T-dt*theta)
        
        tempxun   = PDE.NodalDOFs(exu,Mesh.Nodes)
        exunx,exuny = PDE.DecompIntoCoord(tempxun)
        tempxum   = PDE.NodalDOFs(exu,Mesh.MidNodes)
        exumx,exumy = PDE.DecompIntoCoord(tempxum)
        unxerr,unyerr = PDE.unx-exunx,PDE.uny-exuny
        umxerr,umyerr = PDE.umx-exumx,PDE.umy-exumy
        L2u           = PDE.TVhL2Norm(unxerr,unyerr,umxerr,umyerr)

        Barr = PDE.MagDOFs(exB)
        Berr = Barr-PDE.B
        L2B  = PDE.EhL2Norm(Berr)
        
        Earr = PDE.NodalDOFs(exE,PDE.Mesh.Nodes)
        Eerr = Earr-PDE.E
        L2E  = PDE.VhL2Norm(Eerr)

        parr = PDE.PhDOF(exp)
        perr = parr-PDE.p
        L2p  = PDE.PhL2Norm(perr)

        print('VelErr = '+str(L2u))
        print('ElectricErr = '+str(L2E))
        print('MagneticErr = '+str(L2B))
        print('PresErr = '+str(L2p))
        uL.append(L2u)
        BL.append(L2B)
        EL.append(L2E)
        pL.append(L2p)
        
        i = i+1
    print('uL'+str(uL))
    print('BL'+str(BL))
    print('EL'+str(EL))
    print('pL'+str(pL))