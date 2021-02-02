from PDEClass import PDEFullMHD
from MeshHelios import HeliosMesh
import numpy as np
import math
from Solver import InexactNewtonTimeInt
import pickle
import time

def ProcessedMesh(Pfile):
    with open(Pfile, "rb") as fp:   # Unpickling
        N,E,EE,B,O = pickle.load(fp)
    return N,E,EE,B,O
def RetrieveAMRMesh(Pfile):
    with open(Pfile, "rb") as fp:   # Unpickling
        N,E,EE,B,O,BT,LR,C = pickle.load(fp)
    return N,E,EE,B,O,BT,LR,C
Re,Rm,theta   = 1, 1, 0.5
T                = 0.25
MTypes = ['Trig','Quad','Vor']
#MTypes = ['Vor']
#MTypes = ['Quad','Vor']
#MTypes = ['Quad']
#MTypes  = ['Onlyone']
def exactu(xv,t):
    return np.array([math.exp(t)*math.cos(xv[1]),0])
#def exactB(xv,t):
#    return np.array([math.cos(xv[0]+t),0])
def exactB(xv,t):
    return np.array([0,math.cos(xv[0]+t)])
def exactE(xv,t):
    return math.cos(xv[0]+t)
def exactp(xv,t):
    return -xv[0]*math.cos(xv[1])
def f(xv,t):
    y1 = math.cos(t+xv[0])**2-math.cos(xv[1])+2*math.exp(t)*math.cos(xv[1])\
         +math.exp(t)*math.cos(xv[1])*math.cos(t+xv[0])**2
    #y1 = -math.cos(xv[1])+2*math.exp(t)*math.cos(xv[1])
    y2 = xv[0]*math.sin(xv[1])
    return np.array([y1,y2])
def h(xv,t):
    return math.cos(t+xv[0])*(1+math.exp(t)*math.cos(xv[1]))+math.sin(t+xv[0])

def initu(xv):
    return exactu(xv,0)
def initB(xv):
    return exactB(xv,0)
def initE(xv):
    return exactE(xv,0)
def initp(xv):
    return exactp(xv,0)


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
        PDE    = PDEFullMHD(Mesh,Re,Rm,initu,initB,dt,theta)
        PDE.SetMHDBCandSource(exactu,exactE,f,h)
        PDE.MHDComputeBC(0)
        PDE.MHDComputeSources(0)
        def ut(xv):
            return exactu(xv,dt)
        def Bt(xv):
            return exactB(xv,dt)
        def Et(xv):
            return exactE(xv,theta*dt)
        def pt(xv):
            return exactp(xv,theta*dt)
        tempxun   = PDE.NodalDOFs(ut,Mesh.Nodes)
        xunx,xuny = PDE.DecompIntoCoord(tempxun)
        tempxum   = PDE.NodalDOFs(ut,Mesh.MidNodes)
        xumx,xumy = PDE.DecompIntoCoord(tempxum)
        xB        = PDE.MagDOFs(Bt)
        xE        = PDE.NodalDOFs(Et,Mesh.Nodes)
        xp        = PDE.PhDOF(pt)
        start = time.time()
        y         = PDE.MHDG(PDE.MHDConcatenate(xunx,xuny,xumx,xumy,xB,xE,xp))
        #y         = PDE.pMHDG(PDE.MHDConcatenate(xunx,xuny,xumx,xumy,xB,xE,xp))
        end = time.time()
        print('time ='+str(end-start))
        fnx,fny,fmx,fmy,farf,elecf,divf = PDE.MHDSplity(y)
        #ynx,yny,ymx,ymy,yp = PDE.Splity(y)
        momerr  = PDE.TVhL2Norm(fnx,fny,fmx,fmy)
        Farerr  = PDE.EhL2Norm(farf)
        Elecerr = PDE.VhL2Norm(elecf)
        Masserr = 0
        for j in range(len(Mesh.ElementEdges)):
            lunx,luny,lumx,lumy = PDE.GetLocalTVhDOF(j,PDE.unx,PDE.uny,PDE.umx,PDE.umy)
            divu, A = PDE.DIVu(j,lunx,luny,lumx,lumy)
            Masserr = Masserr+PDE.PhInProd(i,divu*A,divu*A)
        Masserr = math.sqrt(Masserr)

        print('momerr='+str(momerr))
        print('Farerr='+str(Farerr))
        print('Elecerr='+str(Elecerr))
        print('Masserr='+str(Masserr))
    
        # #H1err = PDE.TVhH1Norm(ynx,yny,ymx,ymy)

        # #L2err = PDE.TVhL2Norm(ynx,yny,ymx,ymy)

        # #perr = PDE.PhL2Norm(yp)
        # #print('H1Err='+str(H1err**2))
        # #print('L2Err='+str(L2err**2))
        # #print('perr='+str(perr**2))
        # # y    = np.concatenate((ynx,yny,ymx,ymy), axis=None)
        # print('NumIntNodes='+str(len(Mesh.NumInternalNodes)))
        # print('NumIntMidNodes='+str(len(Mesh.NumInternalMidNodes)))
        # print('NumCells='+str(len(Mesh.ElementEdges)))

        # ave = np.sum(y)/len(y)
        # l2norm = np.linalg.norm(y)
        #maxy   = max(abs(y))
        # miny   = min(abs(y))
        # #print('y =' +str(y))
        #print('max|y|='+str(maxy))
        # print('min|y|='+str(miny))
        #j = 0
        #k = 0
        #for c in abs(y):
            #if abs(c-maxy)<1E-5 and k==0:
            #    print('MaxHappensat'+str(j))
                #k = 1
            #j +=1

        # print('ave='+str(ave))
        # print('EuclideanNorm='+str(l2norm))
        # print('NumEntries='+str(len(y)))
        # #print(y)
        #print(len(PDE.Mesh.NumInternalNodes))
        #print(len(PDE.Mesh.NumInternalMidNodes))
        i = i+1