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
theta,dt,T = 1/2,0.001,0.5
Re,Rm,theta   = 1, 1, 0.5
T                = 0.25
MTypes = ['Trig','Quad','Vor']
#MTypes = ['Quad','Vor']
#MTypes = ['Quad']
#MTypes  = ['Onlyone']
def InB(xv):
    return np.array([0,math.cos(xv[0])])

def exactu(xv):
    return np.array([xv[1]**2,xv[0]**2])

def exactp(xv):
    return 2*xv[0]+2*xv[1]
for MType in MTypes:
    if MType == 'Onlyone':
        ProcessedFiles = ['PTh=0.101015.txt']
        dx = [0.10101525445522107]
        #ProcessedFiles = ['PertPQh=0.021739.txt']
        #ProcessedFiles = ['AMRMesh']
        #dx = [0.021739]
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
    print(MType)
    for Pfile in ProcessedFiles:
        print('\n dx='+str(dx[i]))
        Nodes,EdgeNodes,ElementEdges,BoundaryNodes,Orientations = ProcessedMesh(Pfile)
        #Nodes,EdgeNodes,ElementEdges,BoundaryNodes,Orientations,BottomToTop,LeftToRight,Corners = RetrieveAMRMesh("AMRmesh.txt")
        Mesh = HeliosMesh(Nodes,EdgeNodes,ElementEdges,Orientations)
        dt = dx[i]**2
        PDE    = PDEFullMHD(Mesh,Re,Rm,exactu,InB,dt,theta)
        PDE.nSetFlowBC(exactu)
        PDE.nFlowComputeBC(0)
        #PDE.unx,PDE.uny,PDE.umx,PDE.umy = PDE.FlowupdateBC(PDE.unx,PDE.uny,PDE.umx,PDE.umy)
        PDE.p = PDE.PhDOF(exactp)
        y     = PDE.nFlowG(PDE.FlowConcatenate())
        ynx,yny,ymx,ymy,yp = PDE.nSplity(y)
        print(yp)
        #err = PDE.TVhL2Norm(ynx,yny,ymx,ymy) v
        H1err = PDE.TVhH1Norm(ynx,yny,ymx,ymy)

        L2err = PDE.TVhL2Norm(ynx,yny,ymx,ymy)

        perr = PDE.PhL2Norm(yp)
        print('H1Err='+str(H1err**2))
        print('L2Err='+str(L2err**2))
        print('perr='+str(perr**2))
        # y    = np.concatenate((ynx,yny,ymx,ymy), axis=None)
        # ave = np.sum(y)/len(y)
        # l2norm = np.linalg.norm(y)
        # maxy   = max(abs(y))
        # #print('y =' +str(y))
        # print('max|y|='+str(maxy))
        # j = 0
        # k = 0
        # for c in abs(y):
        #     if abs(c-maxy)<1E-5 and k==0:
        #         print('MaxHappensat'+str(j))
        #         k = 1
        #     j +=1

        # print('ave='+str(ave))
        # print('EuclideanNorm='+str(l2norm))
        # print('NumEntries='+str(len(y)))
        # #print(y)
        # print(len(PDE.Mesh.NumInternalNodes))
        # print(len(PDE.Mesh.NumInternalMidNodes))
        i = i+1