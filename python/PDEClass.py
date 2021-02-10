from MeshHelios import HeliosMesh
import multiprocessing as mp
from scipy.sparse import csr_matrix
from scipy.sparse import lil_matrix
import math
import numpy as np
from numpy.linalg import norm as n2

class PDEFullMHD(object):
    def __init__(self,Mesh,Re,Rm,Inu,InB,dt,theta):
        #The Following values are useful for the implementation of some quadrature rules 
        self.pt0, self.w0  = -1, 1/21
        self.pt1, self.w1 = -math.sqrt((5/11)+(2/11)*math.sqrt(5/3)), (124-7*math.sqrt(15))/350
        self.pt2, self.w2 = -math.sqrt((5/11)-(2/11)*math.sqrt(5/3)), (124+7*math.sqrt(15))/350
        self.pt3, self.w3 = 0, 256/525
        self.pt4, self.w4 = math.sqrt((5/11)-(2/11)*math.sqrt(5/3)), (124+7*math.sqrt(15))/350
        self.pt5, self.w5  = math.sqrt((5/11)+(2/11)*math.sqrt(5/3)), (124-7*math.sqrt(15))/350
        self.pt6, self.w6 = 1, 1/21

        self.xs = [0.0915762135098,0.8168475729805,0.0915762135098,0.1081030181681,0.4459484909160,0.4459484909160]
        self.ys = [0.0915762135098,0.0915762135098,0.8168475729805,0.4459484909160,0.1081030181681,0.4459484909160]
        self.ws = [0.2199034873106/4,0.2199034873106/4,0.2199034873106/4,0.4467631793560/4,0.4467631793560/4,0.4467631793560/4]

        #self.xs = [-0.1081,-0.1081,-0.7838,-0.8168,-0.8168,0.6337]
        #self.ys = [-0.1081,-0.7838,-0.1081,-0.8168,0.6337,-0.8168]
        #self.ws = [0.4468,0.4468,0.4468,0.2199,0.2199,0.2199]
        self.Mesh, self.Re, self.Rm, self.dt, self.theta  = Mesh, Re, Rm, dt, theta
        #We initialize the dofs
        #The boundary values on the vel and elec fields are decoupled from the
        #internal values. Thus, we will keep track of 4 arrays, they are:
        #dofs for elec field,
        #dof of the pressure and magnetic field.
        #the dofs for the vel fields are stored in two arrays, x and y comp
        tempun             = self.NodalDOFs(Inu,self.Mesh.Nodes)
        self.unx,self.uny  = self.DecompIntoCoord(tempun)
        tempum             = self.NodalDOFs(Inu,self.Mesh.MidNodes)
        self.umx, self.umy = self.DecompIntoCoord(tempum)
        self.B             = self.MagDOFs(InB)
        self.p             = np.zeros(len(Mesh.ElementEdges),dtype = float)
        self.E             = np.zeros(len(self.Mesh.Nodes),  dtype = float)
        
        self.evalcount = 0
        self.MRot    = self.Rot(self.Mesh.EdgeNodes,self.Mesh.Nodes)
        self.MEList    = []
        self.MVList    = []
        self.HSTVList  = []
        self.GISTVList = []
        self.DTVList   = []
        self.KTVList   = []
        self.RTKIList  = []
        for i in range(len(self.Mesh.ElementEdges)):
            tempME,tempMV             = self.ElecMagStandMassMat(self.Mesh.ElementEdges[i],self.Mesh.Orientations[i])
            TempSH,TempGI,TempD,TempK,TempRTKI = self.TVhInnerPreCompute(i)
            self.MEList.append(tempME)
            self.MVList.append(tempMV)
            self.HSTVList.append(TempSH)
            self.GISTVList.append(TempGI)
            self.DTVList.append(TempD)
            self.KTVList.append(TempK)
            self.RTKIList.append(TempRTKI)
        
    ##################################################################################
    ##################################################################################    
    #Compute DOFs from func
    def NodalDOFs(self,Func,Nodes):
        #This function computes the dof of the init cond on the vel field.
        return np.array([Func(Node) for Node in Nodes])

    def DecompIntoCoord(self,Array):
        #This function will, given a list of pairs return two arrays.
        #The first and second components.
        arrayx = [x[0] for x in Array]
        Arrayy = [x[1] for x in Array]
        return np.array(arrayx),np.array(Arrayy)

    def PhDOF(self,p):
        ph = np.zeros((len(self.Mesh.ElementEdges)), dtype=float)
        j  = 0
        for Element in self.Mesh.ElementEdges:
            xP,yP,A,V,E = self.Mesh.Centroid(Element,self.Mesh.Orientations[j])
            ph[j]       = A*p([xP,yP])
            j = j+1
            # for i in range(len(V)-1):
            #     Node1 = V[i]
            #     Node2 = V[i+1]
            #     x1,y1 = Node1[0],Node1[1]
            #     x2,y2 = Node2[0],Node2[1]
            
            #     xh1,yh1 = (x1+x2)/2,(y1+y2)/2
            #     xh2,yh2 = (x1+xP)/2,(y1+yP)/2
            #     xh3,yh3 = (x2+xP)/2,(y2+yP)/2

            #     AT    = 0.5*abs( (x2-xP)*(y1-yP)-(y2-yP)*(x1-xP) )
            #     ph[j] = ph[j]+(AT/3)*(p([xh1,yh1])+p([xh2,yh2])+p([xh3,yh3]))
        return ph

    def MagDOFs(self,Func):
    #This computes the dofs of the initial magnetic field
    #This routine could be sped up by vectorizing. For our purposes this
    #this is not necessary
        N    = len(self.Mesh.EdgeNodes)
        proj = np.zeros(N)

        for i in range(N):
            x1      = self.Mesh.Nodes[self.Mesh.EdgeNodes[i][0]][0]
            y1      = self.Mesh.Nodes[self.Mesh.EdgeNodes[i][0]][1]
            x2      = self.Mesh.Nodes[self.Mesh.EdgeNodes[i][1]][0]
            y2      = self.Mesh.Nodes[self.Mesh.EdgeNodes[i][1]][1]
            lengthe = math.sqrt((x2-x1)**2+(y2-y1)**2)
            etimesnormal = [y2-y1,x1-x2]

            [Fx0,Fy0] = Func([x1,y1])
            [Fx1,Fy1] = Func([(x1*(1-self.pt1)+x2*(1+self.pt1))/2,(y1*(1-self.pt1)+y2*(1+self.pt1))/2])
            [Fx2,Fy2] = Func([(x1*(1-self.pt2)+x2*(1+self.pt2))/2,(y1*(1-self.pt2)+y2*(1+self.pt2))/2])
            [Fx3,Fy3] = Func([0.5*(x1+x2),0.5*(y1+y2)])
            [Fx4,Fy4] = Func([(x1*(1-self.pt4)+x2*(1+self.pt4))/2,(y1*(1-self.pt4)+y2*(1+self.pt4))/2])
            [Fx5,Fy5] = Func([(x1*(1-self.pt5)+x2*(1+self.pt5))/2,(y1*(1-self.pt5)+y2*(1+self.pt5))/2])
            [Fx6,Fy6] = Func([x2,y2])

            proj[i] = (self.w0*Fx0+self.w1*Fx1+self.w2*Fx2+self.w3*Fx3+self.w4*Fx4+self.w5*Fx5+self.w6*Fx6)*etimesnormal[0]
            proj[i] = proj[i] +(self.w0*Fy0+self.w1*Fy1+self.w2*Fy2+self.w3*Fy3+self.w4*Fy4+self.w5*Fy5+self.w6*Fy6)*etimesnormal[1]
            proj[i] = proj[i]/(2*lengthe)   
        return proj

    ##################################################################################
    ##################################################################################    
    #These functions work as an interface with the solver class. These are for Stokes flow
    # def StokesConcatenate(self,unx,uny,umx,umy,p):
    #     #This function returns an array that concatenates all the unknowns
    #     intunx = [unx[i] for i in self.Mesh.NumInternalNodes]
    #     intuny = [uny[i] for i in self.Mesh.NumInternalNodes]
    #     intumx = [umx[i] for i in self.Mesh.NumInternalMidNodes]
    #     intumy = [umy[i] for i in self.Mesh.NumInternalMidNodes]

    #     return np.concatenate((intunx,intuny,intumx,intumy,p[0:len(p)-1]), axis=None)
    
    # def StokesSetNumMHDDof(self):
    #     a = len(self.Mesh.NumInternalNodes)
    #     b = len(self.Mesh.NumInternalMidNodes)
    #     d = len(self.Mesh.ElementEdges)
    #     return 2*a + 2*b + d - 1
    
    # # def MHDSplitdelx(self,delx):
    # #     intn  = len(self.Mesh.NumInternalNodes)
    # #     intmn = len(self.Mesh.NumInternalMidNodes)
    # #     intE  = len(self.Mesh.EdgeNodes)
    # #     intC  = len(self.Mesh.ElementEdges)

    # #     cut1 = intn
    # #     cut2 = cut1+intn
    # #     cut3 = cut2+intmn
    # #     cut4 = cut3+intmn
    # #     cut5 = cut4+intE
    # #     cut6 = cut5+intC

    # #     Cutdelx = np.split(delx,[cut1,cut2,cut3,cut4,cut5,cut6])
    # #     delunx  = np.zeros(intn, dtype = float)
    # #     deluny  = np.zeros(intn, dtype = float)
    # #     delE    = np.zeros(intn, dtype = float)
    # #     delumx  = np.zeros(intmn,dtype = float)
    # #     delumy  = np.zeros(intmn,dtype = float)
    # #     delB    = np.zeros(intE, dtype = float)
    # #     delp    = np.zeros(intC, dtype = float)

    # def StokesUpdateInt(self,x,unx,uny,umx,umy,B,E,p):
    #     cut1 = len(self.Mesh.NumInternalNodes) #Number of internal dofs for ux
    #     cut2 = 2*cut1                          #Number of internal dofs for uy
    #     cut3 = cut2+len(self.Mesh.NumInternalMidNodes)                #Number of internal dofs for umx
    #     cut4 = cut3+len(self.Mesh.NumInternalMidNodes)                #Number of internal dofs for umy
    #     Cutx = np.split(x,[cut1,cut2,cut3,cut4])
    #     runx = unx*0
    #     runy = uny*0
    #     rumx = umx*0
    #     rumy = umy*0
    #     rp   = p*0
        
    #     for i in self.Mesh.NumBoundaryNodes:
    #         runx[i] = unx[i]
    #         runy[i] = uny[i]
    #     for i in self.Mesh.NumBMidNodes:
    #         rumx[i] = umx[i]
    #         rumy[i] = umy[i]
    #     j = 0
    #     for i in self.Mesh.NumInternalNodes:
    #         runx[i] = Cutx[0][j]
    #         runy[i] = Cutx[1][j]
    #         j = j+1
    #     j = 0
    #     for i in self.Mesh.NumInternalMidNodes:
    #         rumx[i] = Cutx[2][j]
    #         rumy[i] = Cutx[3][j]
    #         j = j+1
    #     rB = Cutx[4]
    #     rp[0:len(rp)-1]   = Cutx[4]
    #     rp[len(rp)-1]     = -np.sum(p[0:len(p)-1])
    #     return runx,runy,rumx,rumy,rp

    # def StokesUpdateBC(self,unx,uny,umx,umy,E):
    #     runx = unx*0
    #     runy = uny*0
    #     rumx = umx*0
    #     rumy = umy*0
        
    #     for i in self.Mesh.NumInternalNodes:
    #         runx[i] = unx[i]
    #         runy[i] = uny[i]

    #     for i in self.Mesh.NumInternalMidNodes:
    #         rumx[i] = umx[i]
    #         rumy[i] = umy[i]

    #     j = 0
    #     for i in self.Mesh.NumBoundaryNodes:
    #         runx[i] = self.ubnx[j]
    #         runy[i] = self.ubny[j]
    #         j      = j+1
    #     j = 0
    #     for i in self.Mesh.NumBMidNodes:
    #         rumx[i] = self.ubmx[j]
    #         rumy[i] = self.ubmy[j]
    #         j      = j+1
    #     return runx,runy,rumx,rumy

    # def SetStokesBCandSource(self,ub,f):
    #     self.ub,self.f = ub,f

    # def StokesComputeBC(self,t):
    #     def dummyub(xv):
    #         return self.ub(xv,t+self.dt)
    #     tempubn              = self.NodalDOFs(dummyub,self.Mesh.BNodes)
    #     self.ubnx, self.ubny = self.DecompIntoCoord(tempubn)
    #     tempum               = self.NodalDOFs(dummyub,self.Mesh.BMidNodes)
    #     self.ubmx, self.ubmy = self.DecompIntoCoord(tempum)

    
    # def StokesComputeSources(self,t):
    #     def dummyf(xv):
    #         return self.f(xv,t+self.theta*self.dt)
        
    #     tempfn             = self.NodalDOFs(dummyf,self.Mesh.Nodes)
    #     self.fnx, self.fny = self.DecompIntoCoord(tempfn)
    #     tempfm             = self.NodalDOFs(dummyf,self.Mesh.MidNodes)
    #     self.fmx, self.fmy = self.DecompIntoCoord(tempfm)

    # def StokesG(self,x):
    #     #The x is passed because the Scipy Linear Function class requires it.
    #     #It will use the current values of the internal variables
    #     #self.evalcount = self.evalcount+1
    #     unp1x,unp1y = np.zeros((len(self.Mesh.Nodes)),dtype =float),np.zeros((len(self.Mesh.Nodes)),dtype =float)
    #     ump1x,ump1y = np.zeros((len(self.Mesh.MidNodes)),dtype =float),np.zeros((len(self.Mesh.MidNodes)),dtype =float)

    #     p           = np.zeros((len(self.Mesh.ElementEdges)),dtype =float)
    #     unp1x,unp1y,ump1x,ump1y,p = self.StokesUpdateInt(x,unp1x,unp1y,ump1x,ump1y,p)
    #     unp1x,unp1y,ump1x,ump1y = self.StokesUpdateBC(unp1x,unp1y,ump1x,ump1y)

    #     y       = np.zeros(len(x))
    #     nx      = (unp1x-self.unx)/self.dt - self.fnx
    #     ny      = (unp1y-self.uny)/self.dt - self.fny
    #     mx      = (ump1x-self.umx)/self.dt - self.fmx
    #     my      = (ump1y-self.umy)/self.dt - self.fmy
    #     unthetax = (1-self.theta)*self.unx+self.theta*unp1x 
    #     unthetay = (1-self.theta)*self.uny+self.theta*unp1y
    #     umthetax = (1-self.theta)*self.umx+self.theta*ump1x 
    #     umthetay = (1-self.theta)*self.umy+self.theta*ump1y
        
    #     intN,intNM,k = len(self.Mesh.NumInternalNodes),len(self.Mesh.NumInternalMidNodes),0 
    #     for i in self.Mesh.NumInternalNodes:
    #         Cells = self.Mesh.NodestoCells[i]
    #         v1nx,v1ny = np.zeros((len(self.Mesh.Nodes))),np.zeros((len(self.Mesh.Nodes)))
    #         v1mx,v1my = np.zeros((len(self.Mesh.MidNodes))),np.zeros((len(self.Mesh.MidNodes)))
    #         v1nx[i]   = 1
            
    #         v2nx,v2ny = np.zeros((len(self.Mesh.Nodes))),np.zeros((len(self.Mesh.Nodes)))
    #         v2mx,v2my = np.zeros((len(self.Mesh.MidNodes))),np.zeros((len(self.Mesh.MidNodes)))
    #         v2ny[i]   = 1
    #         #Momentum, nodal DOFs
    #         for Cell in Cells:
    #             lnx,lny,lmx,lmy         = self.GetLocalTVhDOF(Cell,nx,ny,mx,my)
    #             lv1nx,lv1ny,lv1mx,lv1my = self.GetLocalTVhDOF(Cell,v1nx,v1ny,v1mx,v1my)
    #             lv2nx,lv2ny,lv2mx,lv2my = self.GetLocalTVhDOF(Cell,v2nx,v2ny,v2mx,v2my)
    #             divv1,A                 = self.DIVu(Cell,lv1nx,lv1ny,lv1mx,lv1my)
    #             divv2,A                 = self.DIVu(Cell,lv2nx,lv2ny,lv2mx,lv2my)
    #             #UseWithOldDisc
    #             #v1xB                     = self.Cross2Dto1D(lv1nx,lv1ny,RTBthetanx,RTBthetany)
    #             #v2xB                     = self.Cross2Dto1D(lv2nx,lv2ny,RTBthetanx,RTBthetany)

    #             y[k] = y[k]\
    #                 +self.TVhInProd(Cell,lnx,lny,lmx,lmy,lv1nx,lv1ny,lv1mx,lv1my)\
    #                 +(1/self.Re)*self.TVhSemiInProd(Cell,locunthetax,locunthetay,locumthetax,locumthetay,lv1nx,lv1ny,lv1mx,lv1my)\
    #                 -self.PhInProd(Cell,p[Cell],divv1*A)
    #                 #+Jn.dot(self.MVList[Cell].dot(v1xB))
    #                 #UseWithNewDisc
                    
    #             y[k+intN]  = y[k+intN]\
    #                 +self.TVhInProd(Cell,lnx,lny,lmx,lmy,lv2nx,lv2ny,lv2mx,lv2my)\
    #                 +(1/self.Re)*self.TVhSemiInProd(Cell,locunthetax,locunthetay, locumthetax,locumthetay,lv2nx,lv2ny,lv2mx,lv2my)\
    #                 -self.PhInProd(Cell,p[Cell],divv2*A)
    #                 #UseWitholdDisc
    #                 #+Jn.dot( self.MVList[Cell].dot(v2xB) )
    #                 #UseWithNewDisc
                    
                    
               
    #         k = k+1
    #     #Momentum Midnodes
    #     k = 0
    #     for i in self.Mesh.NumInternalMidNodes:
    #         Cells = self.Mesh.EdgestoCells[i]
    #         v1nx,v1ny = np.zeros((len(self.Mesh.Nodes))),np.zeros((len(self.Mesh.Nodes)))
    #         v1mx,v1my = np.zeros((len(self.Mesh.MidNodes))),np.zeros((len(self.Mesh.MidNodes)))
    #         v1mx[i]   = 1
            
    #         v2nx,v2ny = np.zeros((len(self.Mesh.Nodes))),np.zeros((len(self.Mesh.Nodes)))
    #         v2mx,v2my = np.zeros((len(self.Mesh.MidNodes))),np.zeros((len(self.Mesh.MidNodes)))
    #         v2my[i]   = 1
    #         for Cell in Cells:
    #             lv1nx,lv1ny,lv1mx,lv1my = self.GetLocalTVhDOF(Cell,v1nx,v1ny,v1mx,v1my)
    #             lv2nx,lv2ny,lv2mx,lv2my = self.GetLocalTVhDOF(Cell,v2nx,v2ny,v2mx,v2my)
    #             lnx,lny,lmx,lmy         = self.GetLocalTVhDOF(Cell,nx,ny,mx,my)
    #             locunthetax,locunthetay, locumthetax,locumthetay = self.GetLocalTVhDOF(Cell,unthetax,unthetay,umthetax,umthetay)
    #             divv1,A       = self.DIVu(Cell,lv1nx,lv1ny,lv1mx,lv1my)
    #             locEn = self.GetLocalVhDOF(Cell,E)

    #             y[k+2*intN] = y[k+2*intN]\
    #                 +self.TVhInProd(Cell,lnx,lny,lmx,lmy,lv1nx,lv1ny,lv1mx,lv1my)\
    #                +(1/self.Re)*self.TVhSemiInProd(Cell,locunthetax,locunthetay, locumthetax,locumthetay,lv1nx,lv1ny,lv1mx,lv1my)\
    #                 -self.PhInProd(Cell,p[Cell],divv1*A)
                    

    #             divv2,A             = self.DIVu(Cell,lv2nx,lv2ny,lv2mx,lv2my)
    #             y[k+2*intN+intNM] = y[k+2*intN+intNM]\
    #                             +self.TVhInProd(Cell,lnx,lny,lmx,lmy,lv2nx,lv2ny,lv2mx,lv2my)\
    #                             +(1/self.Re)*self.TVhSemiInProd(Cell,locunthetax,locunthetay, locumthetax,locumthetay,lv2nx,lv2ny,lv2mx,lv2my)\
    #                             -self.PhInProd(Cell,p[Cell],divv2*A)
                 
    #         k = k+1

    #     k    = 0
    #     Nump = 2*intN+2*intNM
    #     for i in range(NumE-1):
    #         lunxi ,lunyi ,lumxi ,lumyi = self.GetLocalTVhDOF(i,unthetax,unthetay,umthetax,umthetay)
    #         lunxl ,lunyl ,lumxl ,lumyl = self.GetLocalTVhDOF(NumE-1,unthetax,unthetay,umthetax,umthetay)
    #         Divui,Ai                   = self.DIVu(i,lunxi,lunyi,lumxi,lumyi)
            
    #         Divul,Al                   = self.DIVu(NumE-1,lunxl,lunyl,lumxl,lumyl)
    #         y[k+Nump]     = y[k+Nump]+self.PhInProd(i,Divui*Ai,1)-self.PhInProd(NumE-1,Divul*Al,1)
    #         k = k+1
    #     return y
    
    # def MHDSplity(self,y):
    #     fnx,fny = np.zeros((len(self.Mesh.Nodes)),dtype=float),np.zeros((len(self.Mesh.Nodes)),dtype=float)
    #     fmx,fmy = np.zeros((len(self.Mesh.MidNodes)),dtype=float),np.zeros((len(self.Mesh.MidNodes)),dtype=float)
    #     divf    = np.zeros((len(self.Mesh.ElementEdges)),dtype=float)
    #     farf    = np.zeros((len(self.Mesh.EdgeNodes)),dtype=float)
    #     h       = np.zeros((len(self.Mesh.Nodes)),dtype=float)
    #     return self.MHDUpdateInt(y,fnx,fny,fmx,fmy,farf,h,divf)

    ##################################################################################
    ##################################################################################    
    #These functions work as an interface with the solver class. These are for Full MHD
    def MHDConcatenate(self,unx,uny,umx,umy,B,E,p):
        #This function returns an array that concatenates all the unknowns
        intunx = [unx[i] for i in self.Mesh.NumInternalNodes]
        intuny = [uny[i] for i in self.Mesh.NumInternalNodes]
        intumx = [umx[i] for i in self.Mesh.NumInternalMidNodes]
        intumy = [umy[i] for i in self.Mesh.NumInternalMidNodes]
        intE   = [E[i] for i in self.Mesh.NumInternalNodes]

        return np.concatenate((intunx,intuny,intumx,intumy,B,intE,p[0:len(p)-1]), axis=None)
    
    def SetNumMHDDof(self):
        a = len(self.Mesh.NumInternalNodes)
        b = len(self.Mesh.NumInternalMidNodes)
        c = len(self.Mesh.EdgeNodes)
        d = len(self.Mesh.ElementEdges)
        return 3*a + 2*b + c + d - 1
    
    # def MHDSplitdelx(self,delx):
    #     intn  = len(self.Mesh.NumInternalNodes)
    #     intmn = len(self.Mesh.NumInternalMidNodes)
    #     intE  = len(self.Mesh.EdgeNodes)
    #     intC  = len(self.Mesh.ElementEdges)

    #     cut1 = intn
    #     cut2 = cut1+intn
    #     cut3 = cut2+intmn
    #     cut4 = cut3+intmn
    #     cut5 = cut4+intE
    #     cut6 = cut5+intC

    #     Cutdelx = np.split(delx,[cut1,cut2,cut3,cut4,cut5,cut6])
    #     delunx  = np.zeros(intn, dtype = float)
    #     deluny  = np.zeros(intn, dtype = float)
    #     delE    = np.zeros(intn, dtype = float)
    #     delumx  = np.zeros(intmn,dtype = float)
    #     delumy  = np.zeros(intmn,dtype = float)
    #     delB    = np.zeros(intE, dtype = float)
    #     delp    = np.zeros(intC, dtype = float)

    def MHDUpdateInt(self,x,unx,uny,umx,umy,B,E,p):
        cut1 = len(self.Mesh.NumInternalNodes) #Number of internal dofs for ux
        cut2 = 2*cut1                          #Number of internal dofs for uy
        cut3 = cut2+len(self.Mesh.NumInternalMidNodes)                #Number of internal dofs for umx
        cut4 = cut3+len(self.Mesh.NumInternalMidNodes)                #Number of internal dofs for umy
        cut5 = cut4+len(self.B)                #Number of dofs for B
        cut6 = cut5+cut1                       #The number of internal dofs for E is the same as for the vel field.
        Cutx = np.split(x,[cut1,cut2,cut3,cut4,cut5,cut6])
        runx = unx*0
        runy = uny*0
        rumx = umx*0
        rumy = umy*0
        rB   = B*0
        rE   = E*0
        rp   = p*0
        
        for i in self.Mesh.NumBoundaryNodes:
            runx[i] = unx[i]
            runy[i] = uny[i]
            rE[i]   = E[i]
        for i in self.Mesh.NumBMidNodes:
            rumx[i] = umx[i]
            rumy[i] = umy[i]
        j = 0
        for i in self.Mesh.NumInternalNodes:
            runx[i] = Cutx[0][j]
            runy[i] = Cutx[1][j]
            rE[i]   = Cutx[5][j]
            j = j+1
        j = 0
        for i in self.Mesh.NumInternalMidNodes:
            rumx[i] = Cutx[2][j]
            rumy[i] = Cutx[3][j]
            j = j+1
        rB = Cutx[4]
        rp[0:len(rp)-1]   = Cutx[6]
        rp[len(rp)-1]     = -np.sum(p[0:len(p)-1])
        return runx,runy,rumx,rumy,rB,rE,rp

    def MHDUpdateBC(self,unx,uny,umx,umy,E):
        runx = unx*0
        runy = uny*0
        rumx = umx*0
        rumy = umy*0
        rE   = E*0
        
        for i in self.Mesh.NumInternalNodes:
            runx[i] = unx[i]
            runy[i] = uny[i]
            rE[i]   = E[i]

        for i in self.Mesh.NumInternalMidNodes:
            rumx[i] = umx[i]
            rumy[i] = umy[i]

        j = 0
        for i in self.Mesh.NumBoundaryNodes:
            runx[i] = self.ubnx[j]
            runy[i] = self.ubny[j]
            rE[i]   = self.Ebarr[j]
            j      = j+1
        j = 0
        for i in self.Mesh.NumBMidNodes:
            rumx[i] = self.ubmx[j]
            rumy[i] = self.ubmy[j]
            j      = j+1
        return runx,runy,rumx,rumy,rE

    def SetMHDBCandSource(self,ub,Eb,f,h):
        self.ub,self.Eb,self.f,self.h = ub,Eb,f,h

    def MHDComputeBC(self,t):
        def dummyub(xv):
            return self.ub(xv,t+self.dt)
        def dummyEb(xv):
            return self.Eb(xv,t+self.theta*self.dt)
        tempubn              = self.NodalDOFs(dummyub,self.Mesh.BNodes)
        self.ubnx, self.ubny = self.DecompIntoCoord(tempubn)
        tempum               = self.NodalDOFs(dummyub,self.Mesh.BMidNodes)
        self.ubmx, self.ubmy = self.DecompIntoCoord(tempum)
        self.Ebarr           = self.NodalDOFs(dummyEb,self.Mesh.BNodes)
    
    def MHDComputeSources(self,t):
        def dummyf(xv):
            return self.f(xv,t+self.theta*self.dt)
        def dummyh(xv):
            return self.h(xv,t+self.theta*self.dt)
        
        tempfn             = self.NodalDOFs(dummyf,self.Mesh.Nodes)
        self.fnx, self.fny = self.DecompIntoCoord(tempfn)
        tempfm             = self.NodalDOFs(dummyf,self.Mesh.MidNodes)
        self.fmx, self.fmy = self.DecompIntoCoord(tempfm)
        self.hdof          = self.NodalDOFs(dummyh,self.Mesh.Nodes)

    def MomentumyNode(self,node,unx,uny,umx,umy,B,E,p):
        Cells = self.Mesh.NodestoCells[i]
        v1nx,v1ny = np.zeros((len(self.Mesh.Nodes))),np.zeros((len(self.Mesh.Nodes)))
        v1mx,v1my = np.zeros((len(self.Mesh.MidNodes))),np.zeros((len(self.Mesh.MidNodes)))
        v1nx[i]   = 1
            
        v2nx,v2ny = np.zeros((len(self.Mesh.Nodes))),np.zeros((len(self.Mesh.Nodes)))
        v2mx,v2my = np.zeros((len(self.Mesh.MidNodes))),np.zeros((len(self.Mesh.MidNodes)))
        v2ny[i]   = 1
        #Momentum, nodal DOFs
        for Cell in Cells:
            locBtheta           = self.GetLocalEhDOF(Cell,Bntheta)
            #UseForOldDisc
            #RTBthetanx,RTBthetany = self.PiRTBn(locBtheta,Cell)

            locunthetax,locunthetay, locumthetax,locumthetay = self.GetLocalTVhDOF(Cell,unthetax,unthetay,umthetax,umthetay)
            #UseForNewDisc
            RTBthetanx,RTBthetany,RTBthetamx,RTBthetamy,locEm = self.PiRTBnm(locBtheta,E,Cell)
                
            locEn = self.GetLocalVhDOF(Cell,E) 
            Jn    = locEn+self.Cross2Dto1D(locunthetax,locunthetay,RTBthetanx,RTBthetany)
            #UseWithNewDisc
            Jm    = locEm+self.Cross2Dto1D(locumthetax,locumthetay,RTBthetamx,RTBthetamy)

            #UseWithNewDisc
            JxBnx,JxBny = self.Cross1Dto2D(Jn,RTBthetanx,RTBthetany)
            JxBmx,JxBmy = self.Cross1Dto2D(Jm,RTBthetamx,RTBthetamy)
                
            lnx,lny,lmx,lmy         = self.GetLocalTVhDOF(Cell,nx,ny,mx,my)
            lv1nx,lv1ny,lv1mx,lv1my = self.GetLocalTVhDOF(Cell,v1nx,v1ny,v1mx,v1my)
            lv2nx,lv2ny,lv2mx,lv2my = self.GetLocalTVhDOF(Cell,v2nx,v2ny,v2mx,v2my)
            divv1,A                 = self.DIVu(Cell,lv1nx,lv1ny,lv1mx,lv1my)
            divv2,A                 = self.DIVu(Cell,lv2nx,lv2ny,lv2mx,lv2my)
            #UseWithOldDisc
            #v1xB                     = self.Cross2Dto1D(lv1nx,lv1ny,RTBthetanx,RTBthetany)
            #v2xB                     = self.Cross2Dto1D(lv2nx,lv2ny,RTBthetanx,RTBthetany)

            y[k] = y[k]\
                +self.TVhInProd(Cell,lnx,lny,lmx,lmy,lv1nx,lv1ny,lv1mx,lv1my)\
                +(1/self.Re)*self.TVhSemiInProd(Cell,locunthetax,locunthetay,locumthetax,locumthetay,lv1nx,lv1ny,lv1mx,lv1my)\
                -self.PhInProd(Cell,p[Cell],divv1*A)\
                -self.TVhInProd(Cell,JxBnx,JxBny,JxBmx,JxBmy,lv1nx,lv1ny,lv1mx,lv1my)
                #UseWitholdDisc
                #+Jn.dot(self.MVList[Cell].dot(v1xB))
                #UseWithNewDisc
                    
            y[k+intN]  = y[k+intN]\
                +self.TVhInProd(Cell,lnx,lny,lmx,lmy,lv2nx,lv2ny,lv2mx,lv2my)\
                +(1/self.Re)*self.TVhSemiInProd(Cell,locunthetax,locunthetay, locumthetax,locumthetay,lv2nx,lv2ny,lv2mx,lv2my)\
                -self.PhInProd(Cell,p[Cell],divv2*A)\
                -self.TVhInProd(Cell,JxBnx,JxBny,JxBmx,JxBmy,lv2nx,lv2ny,lv2mx,lv2my)
                #UseWitholdDisc
                #+Jn.dot( self.MVList[Cell].dot(v2xB) )
                #UseWithNewDisc
    
    def pMHDG(self,x,Gunx):
        self.evalcount = self.evalcount+1
        unp1x,unp1y = np.zeros((len(self.Mesh.Nodes)),dtype =float),np.zeros((len(self.Mesh.Nodes)),dtype =float)
        ump1x,ump1y = np.zeros((len(self.Mesh.MidNodes)),dtype =float),np.zeros((len(self.Mesh.MidNodes)),dtype =float)
        Bp1         = np.zeros((len(self.Mesh.EdgeNodes)),dtype =float)
        E           = np.zeros((len(self.Mesh.Nodes)),dtype =float)
        p           = np.zeros((len(self.Mesh.ElementEdges)),dtype =float)
        unp1x,unp1y,ump1x,ump1y,Bp1,E,p = self.MHDUpdateInt(x,unp1x,unp1y,ump1x,ump1y,Bp1,E,p)
        unp1x,unp1y,ump1x,ump1y,E       = self.MHDUpdateBC(unp1x,unp1y,ump1x,ump1y,E)

        y       = np.zeros(len(x))
        nx      = (unp1x-self.unx)/self.dt - self.fnx
        ny      = (unp1y-self.uny)/self.dt - self.fny
        mx      = (ump1x-self.umx)/self.dt - self.fmx
        my      = (ump1y-self.umy)/self.dt - self.fmy
        unthetax = (1-self.theta)*self.unx+self.theta*unp1x 
        unthetay = (1-self.theta)*self.uny+self.theta*unp1y
        umthetax = (1-self.theta)*self.umx+self.theta*ump1x 
        umthetay = (1-self.theta)*self.umy+self.theta*ump1y
        Bntheta  = (1-self.theta)*self.B+self.theta*Bp1
        
        input_lst    = self.Mesh.NumInternalNodes

        def fGunx(Node):
            return Gunx(Node,nx,ny,mx,my,unthetax,unthetay,umthetax,umthetay,Bntheta,E,p) 
        y1 = self.pool.map(fGunx,input_lst)
        return y1

    def MHDG(self,x):
        #The x is passed because the Scipy Linear Function class requires it.
        #It will use the current values of the internal variables
        self.evalcount = self.evalcount+1
        unp1x,unp1y = np.zeros((len(self.Mesh.Nodes)),dtype =float),np.zeros((len(self.Mesh.Nodes)),dtype =float)
        ump1x,ump1y = np.zeros((len(self.Mesh.MidNodes)),dtype =float),np.zeros((len(self.Mesh.MidNodes)),dtype =float)
        Bp1         = np.zeros((len(self.Mesh.EdgeNodes)),dtype =float)
        E           = np.zeros((len(self.Mesh.Nodes)),dtype =float)
        p           = np.zeros((len(self.Mesh.ElementEdges)),dtype =float)
        unp1x,unp1y,ump1x,ump1y,Bp1,E,p = self.MHDUpdateInt(x,unp1x,unp1y,ump1x,ump1y,Bp1,E,p)
        unp1x,unp1y,ump1x,ump1y,E       = self.MHDUpdateBC(unp1x,unp1y,ump1x,ump1y,E)
       
        y       = np.zeros(len(x))
        nx      = (unp1x-self.unx)/self.dt - self.fnx
        ny      = (unp1y-self.uny)/self.dt - self.fny
        mx      = (ump1x-self.umx)/self.dt - self.fmx
        my      = (ump1y-self.umy)/self.dt - self.fmy
        unthetax = (1-self.theta)*self.unx+self.theta*unp1x 
        unthetay = (1-self.theta)*self.uny+self.theta*unp1y
        umthetax = (1-self.theta)*self.umx+self.theta*ump1x 
        umthetay = (1-self.theta)*self.umy+self.theta*ump1y
        Bntheta  = (1-self.theta)*self.B+self.theta*Bp1
        
        intN,intNM,k = len(self.Mesh.NumInternalNodes),len(self.Mesh.NumInternalMidNodes),0 
        print('intN='+str(intN))
        for i in self.Mesh.NumInternalNodes:
            print('here')
            Cells = self.Mesh.NodestoCells[i]
            v1nx,v1ny = np.zeros((len(self.Mesh.Nodes))),np.zeros((len(self.Mesh.Nodes)))
            v1mx,v1my = np.zeros((len(self.Mesh.MidNodes))),np.zeros((len(self.Mesh.MidNodes)))
            v1nx[i]   = 1
            
            v2nx,v2ny = np.zeros((len(self.Mesh.Nodes))),np.zeros((len(self.Mesh.Nodes)))
            v2mx,v2my = np.zeros((len(self.Mesh.MidNodes))),np.zeros((len(self.Mesh.MidNodes)))
            v2ny[i]   = 1
            #Momentum, nodal DOFs
            for Cell in Cells:
                locBtheta           = self.GetLocalEhDOF(Cell,Bntheta)
                #UseForOldDisc
                #RTBthetanx,RTBthetany = self.PiRTBn(locBtheta,Cell)

                locunthetax,locunthetay, locumthetax,locumthetay = self.GetLocalTVhDOF(Cell,unthetax,unthetay,umthetax,umthetay)
                #UseForNewDisc
                RTBthetanx,RTBthetany,RTBthetamx,RTBthetamy,locEm = self.PiRTBnm(locBtheta,E,Cell)
                
                locEn = self.GetLocalVhDOF(Cell,E) 
                Jn    = locEn+self.Cross2Dto1D(locunthetax,locunthetay,RTBthetanx,RTBthetany)
                #UseWithNewDisc
                Jm    = locEm+self.Cross2Dto1D(locumthetax,locumthetay,RTBthetamx,RTBthetamy)

                #UseWithNewDisc
                JxBnx,JxBny = self.Cross1Dto2D(Jn,RTBthetanx,RTBthetany)
                JxBmx,JxBmy = self.Cross1Dto2D(Jm,RTBthetamx,RTBthetamy)
                
                lnx,lny,lmx,lmy         = self.GetLocalTVhDOF(Cell,nx,ny,mx,my)
                lv1nx,lv1ny,lv1mx,lv1my = self.GetLocalTVhDOF(Cell,v1nx,v1ny,v1mx,v1my)
                lv2nx,lv2ny,lv2mx,lv2my = self.GetLocalTVhDOF(Cell,v2nx,v2ny,v2mx,v2my)
                divv1,A                 = self.DIVu(Cell,lv1nx,lv1ny,lv1mx,lv1my)
                divv2,A                 = self.DIVu(Cell,lv2nx,lv2ny,lv2mx,lv2my)
                #UseWithOldDisc
                #v1xB                     = self.Cross2Dto1D(lv1nx,lv1ny,RTBthetanx,RTBthetany)
                #v2xB                     = self.Cross2Dto1D(lv2nx,lv2ny,RTBthetanx,RTBthetany)

                y[k] = y[k]\
                    +(1/self.Re)*self.TVhSemiInProd(Cell,locunthetax,locunthetay,locumthetax,locumthetay,lv1nx,lv1ny,lv1mx,lv1my)\
                    -self.PhInProd(Cell,p[Cell],divv1*A)\
                    +Jn.dot(self.MVList[Cell].dot(v1xB))
                    #-self.TVhInProd(Cell,JxBnx,JxBny,JxBmx,JxBmy,lv1nx,lv1ny,lv1mx,lv1my)
                    #UseWitholdDisc
                   
                    #UseWithNewDisc
                    #+self.TVhInProd(Cell,lnx,lny,lmx,lmy,lv1nx,lv1ny,lv1mx,lv1my)
                    
                y[k+intN]  = y[k+intN]\
                    +(1/self.Re)*self.TVhSemiInProd(Cell,locunthetax,locunthetay, locumthetax,locumthetay,lv2nx,lv2ny,lv2mx,lv2my)\
                    -self.PhInProd(Cell,p[Cell],divv2*A)\
                    +Jn.dot( self.MVList[Cell].dot(v2xB) )
                    #-self.TVhInProd(Cell,JxBnx,JxBny,JxBmx,JxBmy,lv2nx,lv2ny,lv2mx,lv2my)
                    #UseWitholdDisc
                    
                    #UseWithNewDisc
                    #+self.TVhInProd(Cell,lnx,lny,lmx,lmy,lv2nx,lv2ny,lv2mx,lv2my)
                    
                    
               
            k = k+1
        #Momentum Midnodes
        k = 0
        for i in self.Mesh.NumInternalMidNodes:
            Cells = self.Mesh.EdgestoCells[i]
            v1nx,v1ny = np.zeros((len(self.Mesh.Nodes))),np.zeros((len(self.Mesh.Nodes)))
            v1mx,v1my = np.zeros((len(self.Mesh.MidNodes))),np.zeros((len(self.Mesh.MidNodes)))
            v1mx[i]   = 1
            
            v2nx,v2ny = np.zeros((len(self.Mesh.Nodes))),np.zeros((len(self.Mesh.Nodes)))
            v2mx,v2my = np.zeros((len(self.Mesh.MidNodes))),np.zeros((len(self.Mesh.MidNodes)))
            v2my[i]   = 1
            for Cell in Cells:
                #UseWithNewDisc
                locBtheta           = self.GetLocalEhDOF(Cell,Bntheta)
                lv1nx,lv1ny,lv1mx,lv1my = self.GetLocalTVhDOF(Cell,v1nx,v1ny,v1mx,v1my)
                lv2nx,lv2ny,lv2mx,lv2my = self.GetLocalTVhDOF(Cell,v2nx,v2ny,v2mx,v2my)
                lnx,lny,lmx,lmy         = self.GetLocalTVhDOF(Cell,nx,ny,mx,my)
                locunthetax,locunthetay, locumthetax,locumthetay = self.GetLocalTVhDOF(Cell,unthetax,unthetay,umthetax,umthetay)
                divv1,A       = self.DIVu(Cell,lv1nx,lv1ny,lv1mx,lv1my)
                locEn = self.GetLocalVhDOF(Cell,E)
                #UseWithNewDisc
                RTBthetanx,RTBthetany,RTBthetamx,RTBthetamy,locEm = self.PiRTBnm(locBtheta,E,Cell)
                 
                Jn    = locEn+self.Cross2Dto1D(locunthetax,locunthetay,RTBthetanx,RTBthetany)
                Jm    = locEm+self.Cross2Dto1D(locumthetax,locumthetay,RTBthetamx,RTBthetamy)
                JxBnx,JxBny = self.Cross1Dto2D(Jn,RTBthetanx,RTBthetany)
                JxBmx,JxBmy = self.Cross1Dto2D(Jm,RTBthetamx,RTBthetamy)

                y[k+2*intN] = y[k+2*intN]\
                   +(1/self.Re)*self.TVhSemiInProd(Cell,locunthetax,locunthetay, locumthetax,locumthetay,lv1nx,lv1ny,lv1mx,lv1my)\
                    -self.PhInProd(Cell,p[Cell],divv1*A)\
                    #-self.TVhInProd(Cell,JxBnx,JxBny,JxBmx,JxBmy,lv1nx,lv1ny,lv1mx,lv1my)
                    #+self.TVhInProd(Cell,lnx,lny,lmx,lmy,lv1nx,lv1ny,lv1mx,lv1my)
                    

                divv2,A             = self.DIVu(Cell,lv2nx,lv2ny,lv2mx,lv2my)
                y[k+2*intN+intNM] = y[k+2*intN+intNM]\
                                +(1/self.Re)*self.TVhSemiInProd(Cell,locunthetax,locunthetay, locumthetax,locumthetay,lv2nx,lv2ny,lv2mx,lv2my)\
                                -self.PhInProd(Cell,p[Cell],divv2*A)\
                                #-self.TVhInProd(Cell,JxBnx,JxBny,JxBmx,JxBmy,lv2nx,lv2ny,lv2mx,lv2my)
                                #+self.TVhInProd(Cell,lnx,lny,lmx,lmy,lv2nx,lv2ny,lv2mx,lv2my)
                 
            k = k+1

        Faraday = (self.MRot).dot(E)+(Bp1-self.B)/self.dt

        MagnN  = 2*intN+2*intNM
        # Faraday
        for k in range(len(self.Mesh.EdgeNodes)):
            Cells = self.Mesh.EdgestoCells[i]
            for Cell in Cells:
                Element      = self.Mesh.ElementEdges[Cell]
                ind          = Element.index(i)
                TestFar      = np.zeros(len(Element))
                TestFar[ind] = 1
                locFar       = self.GetLocalEhDOF(Cell,Faraday)
                y[k+MagnN]         = y[k+MagnN]+locFar.dot(self.MEList[Cell].dot(TestFar))
        #Ampere-Ohm
        k,ElecN = 0,MagnN+len(self.Mesh.EdgeNodes)
        for i in self.Mesh.NumInternalNodes:
            D      = np.zeros((len(self.Mesh.Nodes)))
            D[i]   = 1
            RotD   = (self.MRot).dot(D)
            Cells  = self.Mesh.NodestoCells[i]
            for Cell in Cells:
                LocRotD   = self.GetLocalEhDOF(Cell,RotD)
                LocthetaB = self.GetLocalEhDOF(Cell,Bntheta)
                LocE      = self.GetLocalVhDOF(Cell,E)
                LocD      = self.GetLocalVhDOF(Cell,D)
                Loch      = self.GetLocalVhDOF(Cell,self.hdof)
                RTBthetax,RTBthetay = self.PiRTBn(LocthetaB,Cell)
                
                locunthetax,locunthetay, nouse1,nouse2 = self.GetLocalTVhDOF(Cell,unthetax,unthetay,mx,my)

                J    = LocE+self.Cross2Dto1D(locunthetax,locunthetay,RTBthetax,RTBthetay)
                y[k+ElecN]    = y[k+ElecN]\
                               +(J-Loch).dot(self.MVList[Cell].dot(LocD))\
                               -(1/self.Rm)*LocthetaB.dot(self.MEList[Cell].dot(LocRotD))
            k = k+1

        k = 0
        Nump = ElecN+intN
        NumE = len(self.Mesh.ElementEdges)
        # Div = np.zeros((NumE),dtype=float)
        # As  = []
        # for i in range(NumE-1):
        #     lunxi ,lunyi ,lumxi ,lumyi  = self.GetLocalTVhDOF(i,unthetax,unthetay,umthetax,umthetay)
        #     Divui,Ai                 = self.DIVu(i,lunxi,lunyi,lumxi,lumyi)
        #     Div[i] = Divui
        #     As.append(Ai)
        
        # Div[NumE-1] = -np.sum(Divui)
        # print('divu'+str(Div[NumE-1]))
        # Ai,V,En = self.Mesh.Area(self.Mesh.ElementEdges[NumE-1],self.Mesh.Orientations[NumE-1])
        # As.append(Ai)
        for i in range(NumE-1):
            lunxi ,lunyi ,lumxi ,lumyi = self.GetLocalTVhDOF(i,unthetax,unthetay,umthetax,umthetay)
            lunxl ,lunyl ,lumxl ,lumyl = self.GetLocalTVhDOF(NumE-1,unthetax,unthetay,umthetax,umthetay)
            Divui,Ai                   = self.DIVu(i,lunxi,lunyi,lumxi,lumyi)
            
            Divul,Al                   = self.DIVu(NumE-1,lunxl,lunyl,lumxl,lumyl)
            y[k+Nump]     = y[k+Nump]+self.PhInProd(i,Divui*Ai,1)-self.PhInProd(NumE-1,Divul*Al,1)
            k = k+1
        return y
    
    def MHDSplity(self,y):
        fnx,fny = np.zeros((len(self.Mesh.Nodes)),dtype=float),np.zeros((len(self.Mesh.Nodes)),dtype=float)
        fmx,fmy = np.zeros((len(self.Mesh.MidNodes)),dtype=float),np.zeros((len(self.Mesh.MidNodes)),dtype=float)
        divf    = np.zeros((len(self.Mesh.ElementEdges)),dtype=float)
        farf    = np.zeros((len(self.Mesh.EdgeNodes)),dtype=float)
        h       = np.zeros((len(self.Mesh.Nodes)),dtype=float)
        return self.MHDUpdateInt(y,fnx,fny,fmx,fmy,farf,h,divf)

    ##########################################################################################
    ##########################################################################################
    def nFlowConcatenate(self):
        #This function returns an array that concatenates all the unknowns
        intunx = [self.unx[i] for i in self.Mesh.NumInternalNodes]
        intuny = [self.uny[i] for i in self.Mesh.NumInternalNodes]
        intumx = [self.umx[i] for i in self.Mesh.NumInternalMidNodes]
        intumy = [self.umy[i] for i in self.Mesh.NumInternalMidNodes]
        tempp  = self.p[0:len(self.p)-1]
        return np.concatenate((intunx,intuny,intumx,intumy,tempp), axis=None)

    def nNumFlowDOF(self):
        return 2*len(self.Mesh.NumInternalNodes)+2*len(self.Mesh.NumInternalMidNodes)+len(self.Mesh.ElementEdges)-1

    def nSetFlowBC(self,ub):
        self.ub = ub  #source terms and BC

    def nFlowComputeBC(self,t):
        def dummyub(xv):
            return self.ub([xv[0],xv[1],t+self.dt])
        tempubn              = self.NodalDOFs(dummyub,self.Mesh.BNodes)
        self.ubnx,self.ubny  = self.DecompIntoCoord(tempubn)
        tempum               = self.NodalDOFs(dummyub,self.Mesh.BMidNodes)
        self.ubmx, self.ubmy = self.DecompIntoCoord(tempum)

    def nFlowupdateBC(self,unx,uny,umx,umy):
        j = 0
        for i in self.Mesh.NumBoundaryNodes:
            unx[i]  = self.ubnx[j]
            uny[i]  = self.ubny[j]
            j = j+1
        j = 0
        for i in self.Mesh.NumBMidNodes:
            umx[i]  = self.ubmx[j]
            umy[i]  = self.ubmy[j]
            j = j+1
        return unx,uny,umx,umy

    def nFlowUpdateUnknownDOFs(self,x,unx,uny,umx,umy,p):
        cut1 = len(self.Mesh.NumInternalNodes) #Number of internal dofs for ux
        cut2 = 2*cut1                          #Number of internal dofs for uy
        cut3 = cut2+len(self.Mesh.NumInternalMidNodes)                #Number of internal dofs for umx
        cut4 = cut3+len(self.Mesh.NumInternalMidNodes)                #Number of internal dofs for umy 
        Cutx = np.split(x,[cut1,cut2,cut3,cut4])

        runx,runy = np.zeros((len(self.Mesh.Nodes))),np.zeros((len(self.Mesh.Nodes)))
        rumx,rumy = np.zeros((len(self.Mesh.MidNodes))),np.zeros((len(self.Mesh.MidNodes)))
        rp        = np.zeros((len(self.Mesh.ElementEdges)))
        
        for i in self.Mesh.NumBoundaryNodes:
            runx[i] = unx[i]
            runy[i] = uny[i]
        for i in self.Mesh.NumBMidNodes:
            rumx[i] = umx[i]
            rumy[i] = umy[i]
        j = 0
        for i in self.Mesh.NumInternalNodes:
            runx[i] = Cutx[0][j]
            runy[i] = Cutx[1][j]
            j = j+1
        j = 0
        for i in self.Mesh.NumInternalMidNodes:
            rumx[i] = Cutx[2][j]
            rumy[i] = Cutx[3][j]
            j = j+1
        rp[0:len(rp)-1]   = Cutx[4]
        rp[len(rp)-1]     = -np.sum(rp[0:len(rp)-1])
        return runx,runy,rumx,rumy,rp

    def nFlowG(self,x):
        unx,uny = np.zeros((len(self.Mesh.Nodes))),np.zeros((len(self.Mesh.Nodes)))
        umx,umy = np.zeros((len(self.Mesh.MidNodes))),np.zeros((len(self.Mesh.MidNodes)))
        p       = np.zeros((len(self.Mesh.ElementEdges)))

        intN  = len(self.Mesh.NumInternalNodes) 
        intMN = len(self.Mesh.NumInternalMidNodes)
        Cutx  = np.split(x,[intN,2*intN,2*intN+intMN,2*intN+2*intMN])
        j = 0
        for i in self.Mesh.NumInternalNodes:
            unx[i] = Cutx[0][j]
            uny[i] = Cutx[1][j]
            j = j+1
        j = 0
        for i in self.Mesh.NumInternalMidNodes:
            umx[i] = Cutx[2][j]
            umy[i] = Cutx[3][j]
            j = j+1

        p[0:len(p)-1]   = Cutx[4]

        p[len(p)-1]     = -np.sum(p[0:len(p)-1])
        unx,uny,umx,umy = self.nFlowupdateBC(unx,uny,umx,umy)
        y = np.zeros((len(x)),dtype=float) 
        k = 0 
        for i in self.Mesh.NumInternalNodes:
            Cells = self.Mesh.NodestoCells[i]
            v1nx,v1ny = np.zeros((len(self.Mesh.Nodes))),np.zeros((len(self.Mesh.Nodes)))
            v1mx,v1my = np.zeros((len(self.Mesh.MidNodes))),np.zeros((len(self.Mesh.MidNodes)))
            v1nx[i]   = 1
            
            v2nx,v2ny = np.zeros((len(self.Mesh.Nodes))),np.zeros((len(self.Mesh.Nodes)))
            v2mx,v2my = np.zeros((len(self.Mesh.MidNodes))),np.zeros((len(self.Mesh.MidNodes)))
            v2ny[i]   = 1
            for Cell in Cells:
                lunx ,luny ,lumx ,lumy  = self.GetLocalTVhDOF(Cell,unx,uny,umx,umy)
                lv1nx,lv1ny,lv1mx,lv1my = self.GetLocalTVhDOF(Cell,v1nx,v1ny,v1mx,v1my)
                lv2nx,lv2ny,lv2mx,lv2my = self.GetLocalTVhDOF(Cell,v2nx,v2ny,v2mx,v2my)

                divv1,A = self.DIVu(Cell,lv1nx,lv1ny,lv1mx,lv1my)
                y[k]  = y[k]+(1/self.Re)*self.TVhSemiInProd(Cell,lunx,luny,lumx,lumy,lv1nx,lv1ny,lv1mx,lv1my)-self.PhInProd(Cell,p[Cell],divv1*A)

                divv2,A      = self.DIVu(Cell,lv2nx,lv2ny,lv2mx,lv2my)
                y[k+intN]  = y[k+intN]+(1/self.Re)*self.TVhSemiInProd(Cell,lunx,luny,lumx,lumy,lv2nx,lv2ny,lv2mx,lv2my)-self.PhInProd(Cell,p[Cell],divv2*A)
            k = k+1
        k = 0
        for i in self.Mesh.NumInternalMidNodes:
            Cells = self.Mesh.EdgestoCells[i]
            v1nx,v1ny = np.zeros((len(self.Mesh.Nodes))),np.zeros((len(self.Mesh.Nodes)))
            v1mx,v1my = np.zeros((len(self.Mesh.MidNodes))),np.zeros((len(self.Mesh.MidNodes)))
            v1mx[i]   = 1
            
            v2nx,v2ny = np.zeros((len(self.Mesh.Nodes))),np.zeros((len(self.Mesh.Nodes)))
            v2mx,v2my = np.zeros((len(self.Mesh.MidNodes))),np.zeros((len(self.Mesh.MidNodes)))
            v2my[i]   = 1
            for Cell in Cells:
                
                lunx ,luny ,lumx ,lumy  = self.GetLocalTVhDOF(Cell,unx,uny,umx,umy)
                lv1nx,lv1ny,lv1mx,lv1my = self.GetLocalTVhDOF(Cell,v1nx,v1ny,v1mx,v1my)
                lv2nx,lv2ny,lv2mx,lv2my = self.GetLocalTVhDOF(Cell,v2nx,v2ny,v2mx,v2my)
                
                divv1,A       = self.DIVu(Cell,lv1nx,lv1ny,lv1mx,lv1my)
                y[k+2*intN] = y[k+2*intN]+(1/self.Re)*self.TVhSemiInProd(Cell,lunx,luny,lumx,lumy,lv1nx,lv1ny,lv1mx,lv1my)-self.PhInProd(Cell,p[Cell],divv1*A)

                divv2,A             = self.DIVu(Cell,lv2nx,lv2ny,lv2mx,lv2my)
                y[k+2*intN+intMN] = y[k+2*intN+intMN]+(1/self.Re)*self.TVhSemiInProd(Cell,lunx,luny,lumx,lumy,lv2nx,lv2ny,lv2mx,lv2my)-self.PhInProd(Cell,p[Cell],divv2*A)
            k = k+1
        k = 0
        last = len(self.Mesh.ElementEdges)-1
        lunx2 ,luny2 ,lumx2 ,lumy2  = self.GetLocalTVhDOF(last,unx,uny,umx,umy)
        Divu2,A2                    = self.DIVu(last,lunx2,luny2,lumx2,lumy2)
        for i in range(last):
            lunx ,luny ,lumx ,lumy  = self.GetLocalTVhDOF(i,unx,uny,umx,umy)
            Divu,A                  = self.DIVu(i,lunx,luny,lumx,lumy)
            
            y[k+2*intN+2*intMN]     = y[k+2*intN+2*intMN]+self.PhInProd(i,Divu*A,1)-self.PhInProd(last,Divu2*A2,1)
            k = k+1
        return y

    def nSplity(self,y):
        intN  = len(self.Mesh.NumInternalNodes) #Number of internal dofs for ux
        intBN = len(self.Mesh.NumInternalMidNodes)
        Cuty  = np.split(y,[intN,2*intN,2*intN+intBN,2*intN+2*intBN])
        ynx   = np.zeros(len(self.Mesh.Nodes))
        yny   = np.zeros(len(self.Mesh.Nodes))
        ymx   = np.zeros(len(self.Mesh.MidNodes))
        ymy   = np.zeros(len(self.Mesh.MidNodes))
        yp    = np.zeros(len(self.Mesh.ElementEdges))
        j = 0
        for i in self.Mesh.NumInternalNodes:
            ynx[i] = Cuty[0][j]
            yny[i] = Cuty[1][j]
            j = j+1
        j = 0
        for i in self.Mesh.NumInternalMidNodes:
            ymx[i] = Cuty[2][j]
            ymy[i] = Cuty[3][j]
            j = j+1

        yp[0:len(yp)-1] = Cuty[4]
        yp[len(yp)-1]   = -np.sum(yp[0:len(yp)-1])
        return ynx,yny,ymx,ymy,yp
    ##########################################################################################
    ##########################################################################################
    def FlowConcatenate(self):
        #This function returns an array that concatenates all the unknowns
        intunx = [self.unx[i] for i in self.Mesh.NumInternalNodes]
        intuny = [self.uny[i] for i in self.Mesh.NumInternalNodes]
        intumx = [self.umx[i] for i in self.Mesh.NumInternalMidNodes]
        intumy = [self.umy[i] for i in self.Mesh.NumInternalMidNodes]
        tempp  = self.p[0:len(self.p)-1]
        return np.concatenate((intunx,intuny,intumx,intumy,tempp), axis=None)

    def NumFlowDOF(self):
        return 2*len(self.Mesh.NumInternalNodes)+2*len(self.Mesh.NumInternalMidNodes)+len(self.Mesh.ElementEdges)-1

    def Flowupdatef(self,t):
        def dummyf(xv):
            return self.f(xv,t+self.theta*self.dt)
        tempfn             = self.NodalDOFs(dummyf,self.Mesh.Nodes)
        self.fnx,self.fny  = self.DecompIntoCoord(tempfn)
        tempfm             = self.NodalDOFs(dummyf,self.Mesh.MidNodes)
        self.fmx, self.fmy = self.DecompIntoCoord(tempfm)

    def FlowUpdateInt(self,x,unx,uny,umx,umy,p):
        cut1 = len(self.Mesh.NumInternalNodes) #Number of internal dofs for ux
        cut2 = 2*cut1                          #Number of internal dofs for uy
        cut3 = cut2+len(self.Mesh.NumInternalMidNodes)                #Number of internal dofs for umx
        cut4 = cut3+len(self.Mesh.NumInternalMidNodes)                #Number of internal dofs for umy 
        Cutx = np.split(x,[cut1,cut2,cut3,cut4])
        runx = unx*0
        runy = uny*0
        rumx = umx*0
        rumy = umy*0
        rp   = p*0
        
        for i in self.Mesh.NumBoundaryNodes:
            runx[i] = unx[i]
            runy[i] = uny[i]
        for i in self.Mesh.NumBMidNodes:
            rumx[i] = umx[i]
            rumy[i] = umy[i]
        j = 0
        for i in self.Mesh.NumInternalNodes:
            runx[i] = Cutx[0][j]
            runy[i] = Cutx[1][j]
            j = j+1
        j = 0
        for i in self.Mesh.NumInternalMidNodes:
            rumx[i] = Cutx[2][j]
            rumy[i] = Cutx[3][j]
            j = j+1
        rp[0:len(rp)-1]   = Cutx[4]
        rp[len(rp)-1]     = -np.sum(p[0:len(p)-1])
        return runx,runy,rumx,rumy,rp

    def FlowUpdateBC(self,unx,uny,umx,umy):
        runx = unx*0
        runy = uny*0
        rumx = umx*0
        rumy = umy*0
        
        for i in self.Mesh.NumInternalNodes:
            runx[i] = unx[i]
            runy[i] = uny[i]

        for i in self.Mesh.NumInternalMidNodes:
            rumx[i] = umx[i]
            rumy[i] = umy[i]

        j = 0
        for i in self.Mesh.NumBoundaryNodes:
            runx[i] = self.ubnx[j]
            runy[i] = self.ubny[j]
            j      = j+1

        j = 0
        for i in self.Mesh.NumBMidNodes:
            rumx[i] = self.ubmx[j]
            rumy[i] = self.ubmy[j]
            j      = j+1
        return runx,runy,rumx,rumy

    def SetFlowBCandSource(self,ub,f):
        self.ub,self.f = ub,f

    def FlowComputeBC(self,t):
        def dummyub(xv):
            return self.ub(xv,t+self.dt)
        tempubn              = self.NodalDOFs(dummyub,self.Mesh.BNodes)
        self.ubnx, self.ubny = self.DecompIntoCoord(tempubn)
        tempum               = self.NodalDOFs(dummyub,self.Mesh.BMidNodes)
        self.ubmx, self.ubmy = self.DecompIntoCoord(tempum)

    def FlowG(self,x):
        unp1x,unp1y = np.zeros((len(self.Mesh.Nodes)),dtype =float),np.zeros((len(self.Mesh.Nodes)),dtype =float)
        ump1x,ump1y = np.zeros((len(self.Mesh.MidNodes)),dtype =float),np.zeros((len(self.Mesh.MidNodes)),dtype =float)
        p           = np.zeros((len(self.Mesh.ElementEdges)),dtype =float)
        unp1x,unp1y,ump1x,ump1y,p = self.FlowUpdateInt(x,unp1x,unp1y,ump1x,ump1y,p)
        unp1x,unp1y,ump1x,ump1y   = self.FlowUpdateBC(unp1x,unp1y,ump1x,ump1y)
        #p = self.dt*p
        intN  = len(self.Mesh.NumInternalNodes) 
        intMN = len(self.Mesh.NumInternalMidNodes)

        nx      = (unp1x-self.unx) - self.fnx*self.dt
        ny      = (unp1y-self.uny) - self.fny*self.dt
        mx      = (ump1x-self.umx) - self.fmx*self.dt
        my      = (ump1y-self.umy) - self.fmy*self.dt
        unthetax = (1-self.theta)*self.unx+self.theta*unp1x 
        unthetay = (1-self.theta)*self.uny+self.theta*unp1y
        umthetax = (1-self.theta)*self.umx+self.theta*ump1x 
        umthetay = (1-self.theta)*self.umy+self.theta*ump1y
        
        y = np.zeros((len(x)),dtype=float) 
        k = 0 
        for i in self.Mesh.NumInternalNodes:
            Cells = self.Mesh.NodestoCells[i]
            v1nx,v1ny = np.zeros((len(self.Mesh.Nodes))),np.zeros((len(self.Mesh.Nodes)))
            v1mx,v1my = np.zeros((len(self.Mesh.MidNodes))),np.zeros((len(self.Mesh.MidNodes)))
            v1nx[i]   = 1
            
            v2nx,v2ny = np.zeros((len(self.Mesh.Nodes))),np.zeros((len(self.Mesh.Nodes)))
            v2mx,v2my = np.zeros((len(self.Mesh.MidNodes))),np.zeros((len(self.Mesh.MidNodes)))
            v2ny[i]   = 1
            for Cell in Cells:
                luntx ,lunty ,lumtx ,lumty  = self.GetLocalTVhDOF(Cell,unthetax,unthetay,umthetax,umthetay)
                lnx ,lny ,lmx ,lmy          = self.GetLocalTVhDOF(Cell,nx,ny,mx,my)
                lv1nx,lv1ny,lv1mx,lv1my     = self.GetLocalTVhDOF(Cell,v1nx,v1ny,v1mx,v1my)
                lv2nx,lv2ny,lv2mx,lv2my     = self.GetLocalTVhDOF(Cell,v2nx,v2ny,v2mx,v2my)

                divv1,A = self.DIVu(Cell,lv1nx,lv1ny,lv1mx,lv1my)
                y[k]  = y[k]+\
                    self.TVhInProd(Cell,lnx,lny,lmx,lmy,lv1nx,lv1ny,lv1mx,lv1my)+\
                    self.dt*(1/self.Re)*self.TVhSemiInProd(Cell,luntx,lunty,lumtx,lumty,lv1nx,lv1ny,lv1mx,lv1my)-self.PhInProd(Cell,p[Cell],divv1*A)

                divv2,A      = self.DIVu(Cell,lv2nx,lv2ny,lv2mx,lv2my)
                y[k+intN]  = y[k+intN]+\
                    self.TVhInProd(Cell,lnx,lny,lmx,lmy,lv1nx,lv1ny,lv1mx,lv1my)+\
                        self.dt*(1/self.Re)*self.TVhSemiInProd(Cell,luntx,lunty,lumtx,lumty,lv2nx,lv2ny,lv2mx,lv2my)-self.PhInProd(Cell,p[Cell],divv2*A)
            k = k+1
        k = 0
        for i in self.Mesh.NumInternalMidNodes:
            Cells = self.Mesh.EdgestoCells[i]
            v1nx,v1ny = np.zeros((len(self.Mesh.Nodes))),np.zeros((len(self.Mesh.Nodes)))
            v1mx,v1my = np.zeros((len(self.Mesh.MidNodes))),np.zeros((len(self.Mesh.MidNodes)))
            v1mx[i]   = 1
            
            v2nx,v2ny = np.zeros((len(self.Mesh.Nodes))),np.zeros((len(self.Mesh.Nodes)))
            v2mx,v2my = np.zeros((len(self.Mesh.MidNodes))),np.zeros((len(self.Mesh.MidNodes)))
            v2my[i]   = 1
            for Cell in Cells:
                
                luntx ,lunty ,lumtx ,lumty  = self.GetLocalTVhDOF(Cell,unthetax,unthetay,umthetax,umthetay)
                lnx ,lny ,lmx ,lmy          = self.GetLocalTVhDOF(Cell,nx,ny,mx,my)
                lv1nx,lv1ny,lv1mx,lv1my = self.GetLocalTVhDOF(Cell,v1nx,v1ny,v1mx,v1my)
                lv2nx,lv2ny,lv2mx,lv2my = self.GetLocalTVhDOF(Cell,v2nx,v2ny,v2mx,v2my)
                
                divv1,A       = self.DIVu(Cell,lv1nx,lv1ny,lv1mx,lv1my)
                y[k+2*intN] = y[k+2*intN]+\
                    self.TVhInProd(Cell,lnx,lny,lmx,lmy,lv1nx,lv1ny,lv1mx,lv1my)+\
                    self.dt*(1/self.Re)*self.TVhSemiInProd(Cell,luntx,lunty,lumtx,lumty,lv1nx,lv1ny,lv1mx,lv1my)-self.PhInProd(Cell,p[Cell],divv1*A)

                divv2,A             = self.DIVu(Cell,lv2nx,lv2ny,lv2mx,lv2my)
                y[k+2*intN+intMN] = y[k+2*intN+intMN]+\
                    self.TVhInProd(Cell,lnx,lny,lmx,lmy,lv1nx,lv1ny,lv1mx,lv1my)+\
                    self.dt*(1/self.Re)*self.TVhSemiInProd(Cell,luntx,lunty,lumtx,lumty,lv2nx,lv2ny,lv2mx,lv2my)-self.PhInProd(Cell,p[Cell],divv2*A)
            k = k+1
        k = 0
        for i in range(len(self.Mesh.ElementEdges)-1):
            last = len(self.Mesh.ElementEdges)-1
            luntx ,lunty ,lumtx ,lumty  = self.GetLocalTVhDOF(i,unthetax,unthetay,umthetax,umthetay)
            Divu,A                      = self.DIVu(i,luntx,lunty,lumtx,lumty)
            luntx2 ,lunty2 ,lumtx2 ,lumty2  = self.GetLocalTVhDOF(last,unthetax,unthetay,umthetax,umthetay)
            Divu2,A2                        = self.DIVu(last,luntx2,lunty2,lumtx2,lumty2)
            y[k+2*intN+2*intMN]             = y[k+2*intN+2*intMN]+self.PhInProd(i,Divu*A,1)-self.PhInProd(last,Divu2*A2,1)
            k = k+1
        return y

    def Splity(self,y):
        intN  = len(self.Mesh.NumInternalNodes) #Number of internal dofs for ux
        intBN = len(self.Mesh.NumInternalMidNodes)
        Cuty  = np.split(y,[intN,2*intN,2*intN+intBN,2*intN+2*intBN])
        ynx   = np.zeros(len(self.Mesh.Nodes))
        yny   = np.zeros(len(self.Mesh.Nodes))
        ymx   = np.zeros(len(self.Mesh.MidNodes))
        ymy   = np.zeros(len(self.Mesh.MidNodes))
        yp    = np.zeros(len(self.Mesh.ElementEdges))
        j = 0
        for i in self.Mesh.NumInternalNodes:
            ynx[i] = Cuty[0][j]
            yny[i] = Cuty[1][j]
            j = j+1
        j = 0
        for i in self.Mesh.NumInternalMidNodes:
            ymx[i] = Cuty[2][j]
            ymy[i] = Cuty[3][j]
            j = j+1

        yp[0:len(yp)-1] = Cuty[4]
        yp[len(yp)-1]   = -np.sum(yp[0:len(yp)-1])
        return ynx,yny,ymx,ymy,yp
    ##########################################################################################
    ##########################################################################################
    #Here are the interface to Solver functions regarding a tests involving only the electromagnetics.
    #         
    def ElectroConcatenate(self):
        #This function returns an array that concatenates all the unknowns
        intE = [self.E[i] for i in self.Mesh.NumInternalNodes]
        return np.concatenate((self.B,intE), axis=None)
    
    def NumElectroDOF(self):
        return len(self.Mesh.NumInternalNodes)+len(self.Mesh.EdgeNodes)

    def Electroupdateh(self,t):
        def dummyh(xv):
            return self.h([xv[0],xv[1],t+self.theta*self.dt])
        self.hdof = self.NodalDOFs(dummyh,self.Mesh.Nodes)

    def ElectroComputeBC(self,t):
        def dummyEb(xv):
            return self.Eb([xv[0],xv[1],t+self.theta*self.dt])
        self.ElectroBC = self.NodalDOFs(dummyEb,self.Mesh.BNodes)

    def ElectroupdateBC(self,E):
        j = 0
        for i in self.Mesh.NumBoundaryNodes:
            E[i]  = self.ElectroBC[j]
            j = j+1
        return E

    def SetElectroBCAndSource(self,h,Eb):
        self.h, self.Eb = h, Eb  #source terms and BC

    def ElectroUpdateUnknownDOFs(self,x):
        cut1 = len(self.Mesh.EdgeNodes) #Number of internal dofs for ux
        Cutx = np.split(x,[cut1])
        j = 0
        for i in self.Mesh.NumInternalNodes:
            self.E[i]   = Cutx[1][j]
            j = j+1
        self.B = Cutx[0]

    def ElectroG(self,x):
        cut1 = len(self.Mesh.EdgeNodes) #Number of internal dofs for ux
        Cutx = np.split(x,[cut1])
        Bnp1 = Cutx[0]
        E    = np.zeros((len(self.Mesh.Nodes)))
        j = 0
        for i in self.Mesh.NumInternalNodes:
            E[i] = Cutx[1][j]
            j    = j+1
        j = 0
        for i in self.Mesh.NumBoundaryNodes:
            E[i] = self.ElectroBC[j]
            j    = j+1
        y       = np.zeros(len(x))
        MRot    = self.Rot(self.Mesh.EdgeNodes,self.Mesh.Nodes)
        Faraday = MRot.dot(E)+(Bnp1-self.B)/self.dt
        N       = len(self.Mesh.EdgeNodes)
        for i in range(N):
            Cells = self.Mesh.EdgestoCells[i]
            for Cell in Cells:
                Element      = self.Mesh.ElementEdges[Cell]
                ind          = Element.index(i)
                TestFar      = np.zeros(len(Element))
                TestFar[ind] = 1
                locFar       = self.GetLocalEhDOF(Cell,Faraday)
                y[i]         = y[i]+locFar.dot(self.MEList[Cell].dot(TestFar))
        i = 0
        for v in self.Mesh.NumInternalNodes:
            D      = np.zeros((len(self.Mesh.Nodes)))
            D[v]   = 1
            RotD   = MRot.dot(D)
            Cells  = self.Mesh.NodestoCells[v]
            thetaB = (1-self.theta)*self.B+self.theta*Bnp1
            for Cell in Cells:
                LocRotD   = self.GetLocalEhDOF(Cell,RotD)
                LocthetaB = self.GetLocalEhDOF(Cell,thetaB)
                LocE      = self.GetLocalVhDOF(Cell,E)
                LocD      = self.GetLocalVhDOF(Cell,D)
                Loch      = self.GetLocalVhDOF(Cell,self.hdof)
                y[i+N]    = y[i+N]+(LocE-Loch).dot(self.MVList[Cell].dot(LocD))-(1/self.Rm)*LocthetaB.dot(self.MEList[Cell].dot(LocRotD))
            i = i+1
        return y

    #########################################################################################
    #########################################################################################
    #The following routines will, given a cell compute each of the bilinear forms in the var form.
    #Construct MFD-type matrices
    def LocalMassMatrix(self,N,R,n,A):
        #Given the matrices N,R as defined in Ch.4 of MFD book and the dimension
        #of the reconstruction space this function assembles the local mass matrix
        #The formula is M=M0+M1 where M0=R(N^T R)^-1R^T and M1=lamb*DD^T where the 
        #columns of D span the null-space of N^T and lamb=2*trace(M0)/n 
        #n is the dimension of the reconstruction space
        #nu is the average, over the element, of the diffusion coefficient
        #A is the area of the element
    
        #These commands compute M0
        M0    = np.matmul(np.transpose(N),R) 
        M0    = np.linalg.inv(M0)
        M0    = np.matmul(R,M0)
        M0    = np.matmul(M0,np.transpose(R))
    
        M1    = np.linalg.inv(np.transpose(N).dot(N))
        M1    = np.identity(n)-N.dot(M1).dot(np.transpose(N))
    
        gamma = np.trace(R.dot(np.transpose(R)))/(n*A)
        #And finally we put the two matrices together
        return M0+M1*gamma

    ############Electromagnetics
    def ElecMagStandMassMat(self,Element,Ori):
        n                = len(Element)
        NE               = np.zeros((n,2))
        RE               = np.zeros((n,2))
        xP,yP,A,Vertices,Edges = self.Mesh.Centroid(Element,Ori)
        for i in range(n):
            x1 = Vertices[i][0]
            y1 = Vertices[i][1]
            x2 = Vertices[i+1][0]
            y2 = Vertices[i+1][1]
            lengthEdge = math.sqrt((x2-x1)**2+(y2-y1)**2)
            NE[i][0] = (y2-y1)*Ori[i]*lengthEdge**-1
            NE[i][1] = (x1-x2)*Ori[i]*lengthEdge**-1
            RE[i][0] = (0.5*(x1+x2)-xP)*Ori[i]*lengthEdge #These formulas are derived in the tex-document
            RE[i][1] = (0.5*(y1+y2)-yP)*Ori[i]*lengthEdge
        
        NV = np.ones( (n,1))
        RV = np.zeros((n,1))
        
        x1n = Vertices[n-1][0] #first vertex of n-1th-edge
        y1n = Vertices[n-1][1]
    
        x2n = Vertices[0][0]
        y2n = Vertices[0][1] #second vertex of n-1th edge
    
        x11 = x2n #first vertex of first edge
        y11 = y2n
    
        x21 = Vertices[1][0]
        y21 = Vertices[1][1]  #second vertex of first edge
    
        omegan2 = (x2n-x1n)*((yP-y2n)+(2*yP-y1n-y2n))/6
        omega11 = (x21-x11)*((yP-y11)+(2*yP-y11-y21))/6
        RV[0] = omegan2+omega11
   
        for i in range(1,n):
        
            x1iminusone = Vertices[i-1][0] #first vertex of i-1th-edge
            y1iminusone = Vertices[i-1][1]
               
            x2iminusone = Vertices[i][0]    
            y2iminusone = Vertices[i][1] #second vertex of i-1th edge

            x1i = x2iminusone #first vertex of i+1 edge
            y1i = y2iminusone
    
            x2i = Vertices[i+1][0] #second vertex of i+1 edge
            y2i = Vertices[i+1][1]  
        
            omega2iminusone = (x2iminusone-x1iminusone)*((yP-y2iminusone)+\
                                                   (2*yP-y1iminusone-y2iminusone))/6
            omega1i = (x2i-x1i)*((yP-y1i)+(2*yP-y2i-y1i))/6   
        
            RV[i] = omega2iminusone+omega1i
        ME = self.LocalMassMatrix(NE,RE,n,A)
        MV = self.LocalMassMatrix(NV,RV,n,A)
        return ME,MV

    def BDivSquared(self,B):
        #This function computes the divergence of B.
        D      = 0
        DivMat = self.BDiv()
        divB   = DivMat.dot(B)
        for k in range(len(divB)):
            A,V,E = self.Mesh.Area(self.Mesh.ElementEdges[k],self.Mesh.Orientations[k])
            D     = D + A*divB[k]**2 
        return D
    
    def BDiv(self):
        NEl = len(self.Mesh.ElementEdges)
        NE  = len(self.Mesh.EdgeNodes)
        div = np.zeros((NEl,NE))
        for i in range(NEl):
            Element = self.Mesh.ElementEdges[i]
            N       = len(Element)
            Ori     = self.Mesh.Orientations[i]
            for j in range(N):
                Node1 = self.Mesh.EdgeNodes[Element[j]][0]
                Node2 = self.Mesh.EdgeNodes[Element[j]][1]

                x1    = self.Mesh.Nodes[Node1][0]
                y1    = self.Mesh.Nodes[Node1][1] #these formulas are derived in the pdf document
                x2    = self.Mesh.Nodes[Node2][0]
                y2    = self.Mesh.Nodes[Node2][1]

                lengthe           = math.sqrt((x2-x1)**2+(y2-y1)**2)
                div[i,Element[j]] = Ori[j]*lengthe
        return div
    
    def Rot(self,EdgeNodes,Nodes):
    #This routine computes the primary curl as a matrix
        nN = len(Nodes)
        nE = len(EdgeNodes)
    #curl = np.zeros((nE,nN))
        curl = lil_matrix((nE,nN))
        for i in range(nE):
            Node1 = EdgeNodes[i][0]
            Node2 = EdgeNodes[i][1]
            x1 = Nodes[Node1][0]
            y1 = Nodes[Node1][1]
            x2 = Nodes[Node2][0]
            y2 = Nodes[Node2][1]
            lengthe = math.sqrt((x2-x1)**2+(y2-y1)**2)
        
            curl[i,Node2] = 1/lengthe #These formulas are derived in the pdf document
            curl[i,Node1] = -1/lengthe
        curl = curl.tocsr()
        return curl

    def GetLocalEhDOF(self,ElementNum,arr):
        Element = self.Mesh.ElementEdges[ElementNum]
        Locarr  = np.zeros( (len(Element)) )
        i = 0
        for e in Element:
            Locarr[i] = arr[e]
            i         = i+1
        return Locarr

    def GetLocalVhDOF(self,ElementNum,arr):
        Element = self.Mesh.ElementEdges[ElementNum]
        Locarr  = np.zeros( (len(Element)) )
        V,E     = self.Mesh.StandardElement(Element,self.Mesh.Orientations[ElementNum])
        i = 0
        for e in range(len(E)-1):
            Edge      = E[e]
            v         = Edge[0]
            Locarr[i] = arr[v]
            i         = i+1
        return Locarr

    def VhL2Norm(self,Arr):
        Norm = 0
        for i in range(len(self.Mesh.ElementEdges)):
            LocArr = self.GetLocalVhDOF(i,Arr)
            Norm   = Norm+LocArr.dot( self.MVList[i].dot(LocArr))
        return math.sqrt(Norm)
    
    def EhL2Norm(self,Arr):
        Norm = 0
        for i in range(len(self.Mesh.ElementEdges)):
            LocArr = self.GetLocalEhDOF(i,Arr)
            Norm   = Norm+LocArr.dot( self.MEList[i].dot(LocArr))
        return math.sqrt(Norm)

    def PiRTBB(self,B,Element,ElementNum,E):
        #The DOF must be local.
        N  = len(Element)
        BB = np.zeros((3),dtype=float)
        ori = self.Mesh.Orientations[ElementNum]
        for i in range(N):
            Edge = E[i]
            Node1,Node2  = self.Mesh.Nodes[Edge[0]],self.Mesh.Nodes[Edge[1]]
            x1,y1,x2,y2  = Node1[0],Node1[1],Node2[0],Node2[1]
            xh,yh        = (x1+x2)/2,(y1+y2)/2
            ell          = math.sqrt( (x2-x1)**2+(y2-y1)**2 )
            BB[0] = BB[0]+ori[i]*B[i]*(ell/6)*(x1+4*xh+x2)
            BB[1] = BB[1]+ori[i]*B[i]*(ell/6)*(y1+4*yh+y2) 
            BB[2] = BB[2]+ori[i]*B[i]*(ell/12)*((x1**2+y1**2)+4*(xh**2+yh**2)+(x2**2+y2**2))
        
        return BB
            
    def PiRTBn(self,LocB,ElementNum):
        #The DOF
        Element = self.Mesh.ElementEdges[ElementNum]
        N       = len(Element)
        V,E     = self.Mesh.StandardElement(Element,self.Mesh.Orientations[ElementNum])
        BB      = self.PiRTBB(LocB,Element,ElementNum,E)
        coeffs  = self.RTKIList[ElementNum].dot(BB)
        PiRTBx  = np.zeros((N),dtype=float)
        PiRTBy  = np.zeros((N),dtype=float)

        for i in range(N):
            Edge = E[i]
            Node = self.Mesh.Nodes[Edge[0]]
            x,y  = Node[0],Node[1]
            PiRTBx[i] = coeffs[0]+coeffs[2]*x 
            PiRTBy[i] = coeffs[1]+coeffs[2]*y
        
        return PiRTBx,PiRTBy

    def PiRTBnm(self,LocB,El,ElementNum):
        #The DOF
        Element = self.Mesh.ElementEdges[ElementNum]
        N       = len(Element)
        V,E     = self.Mesh.StandardElement(Element,self.Mesh.Orientations[ElementNum])
        BB      = self.PiRTBB(LocB,Element,ElementNum,E)
        coeffs  = self.RTKIList[ElementNum].dot(BB)
        PiRTBnx  = np.zeros((N),dtype=float)
        PiRTBny  = np.zeros((N),dtype=float)
        PiRTBmx  = np.zeros((N),dtype=float)
        PiRTBmy  = np.zeros((N),dtype=float)
        Em       = np.zeros((N),dtype=float)
        for i in range(N):
            Edge = E[i]
            n1,n2 = Edge[0],Edge[1]
            Node1 = self.Mesh.Nodes[n1]
            Node2 = self.Mesh.Nodes[n2]
            x1,y1  = Node1[0],Node1[1]
            x2,y2  = Node2[0],Node2[1]
            PiRTBnx[i] = coeffs[0]+coeffs[2]*x1 
            PiRTBny[i] = coeffs[1]+coeffs[2]*y1
            PiRTBmx[i] = coeffs[0]+coeffs[2]*0.5*(x1+x2)
            PiRTBmy[i] = coeffs[1]+coeffs[2]*0.5*(y1+y2)
            Em[i]      = 0.5*(El[n1]+El[n2])
        return PiRTBnx,PiRTBny,PiRTBmx,PiRTBmy,Em

    ######################################################################################
    #Fluid Flow
    
    def q0(self,x,y):
        return np.array([1,0])

    def q1(self,x,y):
        return np.array([0,1])
    
    def q2(self,x,y):
        return np.array([x,0])
    
    def q3(self,x,y):
        return np.array([0,x])

    def q4(self,x,y):
        return np.array([y,0])

    def q5(self,x,y):
        return np.array([0,y])
    
    def q6(self,x,y):
        return np.array([x**2,0])

    def q7(self,x,y):
        return np.array([0,x**2])

    def q8(self,x,y):
        return np.array([y**2,0])

    def q9(self,x,y):
        return np.array([0,y**2])

    def q10(self,x,y):
        return np.array([x*y,0])

    def q11(self,x,y):
        return np.array([0,x*y])

    def PhInProd(self,ElementNumber,ph,qh):
        #This function integrates two functions in Ph over the provided element.
        Element = self.Mesh.ElementEdges[ElementNumber]
        A,V,E = self.Mesh.Area(Element,self.Mesh.Orientations[ElementNumber])
        return ph*qh/A
    
    def PhL2Norm(self,ph):
        Norm, j = 0,0
        for Element in self.Mesh.ElementEdges:
            A,V,E = self.Mesh.Area(Element,self.Mesh.Orientations[j])
            Norm  = Norm+ph[j]*ph[j]/A
            j = j+1
        return math.sqrt(Norm)
    
    def DIVu(self,ElementNumber,unx,uny,umx,umy):
        #This routine computes the divergence of u over the element provided.
        Element = self.Mesh.ElementEdges[ElementNumber]
        A,V,E   = self.Mesh.Area(Element,self.Mesh.Orientations[ElementNumber])
        unx     = np.append(unx,unx[0])
        uny     = np.append(uny,uny[0])
        S, k = 0, 0

        for i in range(len(Element)):
            n1, n2 = E[k][0], E[k][1]
            v1, v2 = self.Mesh.Nodes[n1], self.Mesh.Nodes[n2]
            x1, x2 = v1[0], v2[0]
            y1, y2 = v1[1], v2[1]
            etimesnormal = [y2-y1,x1-x2]
            S = S+(unx[i]+unx[i+1]+4*umx[i])*etimesnormal[0]+(uny[i]+uny[i+1]+4*umy[i])*etimesnormal[1]
            k = k+1
        return S/(6*A),A     
    
    def L(self,x0,y0,x1,y1,x2,y2,x,y):
        # D = (x1-x0)*(y2-y0)-(x2-x0)*(y1-y0)
        # xtemp = x-x0
        # ytemp = y-y0

        # l1 = (y2-y0)*xtemp+(x0-x2)*ytemp
        # l2 = (y0-y1)*xtemp+(x1-x0)*ytemp
        # return l1/D,l2/D
        l1 = (x1-x0)*x+(x2-x0)*y+x0
        l2 = (y1-y0)*x+(y2-y0)*y+y0
        return l1,l2
        

    def TVhInnerPreCompute(self, ElementNumber):
        #This function will compute one of the matrices that make the Semi-inner product
        Element     = self.Mesh.ElementEdges[ElementNumber]
        xP,yP,A,V,E = self.Mesh.Centroid(Element,self.Mesh.Orientations[ElementNumber])
        xyP,xxP,yyP                   = 0,0,0
        xxxP,yyyP,xxyP,xyyP           = 0,0,0,0
        xxxxP,yyyyP,xyyyP,xxyyP,xxxyP = 0,0,0,0,0
        for i in range(len(V)-1):
            Node1 = V[i]
            Node2 = V[i+1]
            x1,y1 = Node1[0],Node1[1]
            x2,y2 = Node2[0],Node2[1]
        
            xh1,yh1 = (x1+x2)/2,(y1+y2)/2
            xh2,yh2 = (x1+xP)/2,(y1+yP)/2
            xh3,yh3 = (x2+xP)/2,(y2+yP)/2

            AT  = 0.5*abs( (x2-xP)*(y1-yP)-(y2-yP)*(x1-xP) )
            xyP = xyP+(AT/3)*(xh1*yh1+xh2*yh2+xh3*yh3)
            xxP = xxP+(AT/3)*(xh1**2+xh2**2+xh3**2)
            yyP = yyP+(AT/3)*(yh1**2+yh2**2+yh3**2)

            Jac = abs( (x1-xP)*(y2-yP)-(x2-xP)*(y1-yP) )
            for u in range(6):
                x,y   = self.xs[u],self.ys[u]
                lx,ly = self.L(xP,yP,x1,y1,x2,y2,x,y)
                   
                xxxP  = xxxP +Jac*self.ws[u]*lx**3
                yyyP  = yyyP +Jac*self.ws[u]*ly**3
                xxyP  = xxyP +Jac*self.ws[u]*(lx**2)*ly
                xyyP  = xyyP +Jac*self.ws[u]*(lx)*(ly**2)
                xxxxP = xxxxP+Jac*self.ws[u]*(lx**4)
                yyyyP = yyyyP+Jac*self.ws[u]*(ly**4)
                xyyyP = xyyyP+Jac*self.ws[u]*(lx)*(ly**3)
                xxyyP = xxyyP+Jac*self.ws[u]*(lx**2)*(ly**2)
                xxxyP = xxxyP+Jac*self.ws[u]*(lx**3)*(ly)

        H   = np.zeros((12,12),dtype=float)
        G   = np.zeros((12,12),dtype=float)
        K   = np.zeros((12,12),dtype=float)
        RTK = np.zeros((3,3),  dtype=float)


        RTK[0,0]          = A
        RTK[0,2],RTK[2,0] = xP*A,xP*A

        RTK[1,1]          = A
        RTK[1,2],RTK[2,1] = yP*A,yP*A

        RTK[2,2]          = xxP+yyP  

        H[2,6], H[6,2]  = 2*xP*A,2*xP*A
        H[2,10],H[10,2] = yP*A,yP*A
        H[2,2]          = A

        H[3,7], H[7,3]  = 2*xP*A,2*xP*A
        H[3,11],H[11,3] = yP*A,yP*A
        H[3,3]          = A

        H[4,8], H[8,4]  = 2*yP*A,2*yP*A
        H[4,10],H[10,4] = xP*A,xP*A
        H[4,4]          = A 

        H[5,9],H[9,5]   = 2*yP*A,2*yP*A
        H[5,11],H[11,5] = xP*A,xP*A
        H[5,5]          = A
    
        H[6,10],H[10,6] = 2*xyP,2*xyP
        H[6,6]          = 4*xxP

        H[7,11],H[11,7] = 2*xyP,2*xyP
        H[7,7]          = 4*xxP

        H[8,10],H[10,8] = 2*xyP,2*xyP
        H[8,8]          = 4*yyP

        H[9,11],H[11,9] = 2*xyP,2*xyP
        H[9,9]          = 4*yyP

        H[10,10]        = yyP+xxP

        H[11,11]        = yyP+xxP

        G[2,6], G[6,2]  = 2*xP*A,2*xP*A
        G[2,10],G[10,2] = yP*A,yP*A
        G[2,2]          = A

        G[3,7], G[7,3]  = 2*xP*A,2*xP*A
        G[3,11],G[11,3] = yP*A,yP*A
        G[3,3]          = A

        G[4,8], G[8,4]  = 2*yP*A,2*yP*A
        G[4,10],G[10,4] = xP*A,xP*A 
        G[4,4]          = A

        G[5,9],G[9,5]   = 2*yP*A,2*yP*A
        G[5,11],G[11,5] = xP*A,xP*A
        G[5,5]          = A

        G[6,10],G[10,6] = 2*xyP,2*xyP
        G[6,6]          = 4*xxP

        G[7,11],G[11,7] = 2*xyP,2*xyP
        G[7,7]          = 4*xxP

        G[8,10],G[10,8] = 2*xyP,2*xyP
        G[8,8]          = 4*yyP

        G[9,11],G[11,9] = 2*xyP,2*xyP
        G[9,9]          = 4*yyP

        G[10,10]        = yyP+xxP

        G[11,11]        = yyP+xxP

        K[0,0]          = A
        K[0,2],K[2,0]   = A*xP, A*xP
        K[0,4],K[4,0]   = A*yP, A*yP
        K[0,6],K[6,0]   = xxP,xxP
        K[0,8],K[8,0]   = yyP,yyP
        K[0,10],K[10,0] = xyP,xyP

        K[1,1]          = A
        K[1,3],K[3,1]   = A*xP,A*xP
        K[1,5],K[5,1]   = A*yP,A*yP
        K[1,7],K[7,1]   = xxP,xxP
        K[1,9],K[9,1]   = yyP,yyP
        K[1,11],K[11,1] = xyP,xyP

        K[2,2]          = xxP
        K[2,4],K[4,2]   = xyP,xyP
        K[2,6],K[6,2]   = xxxP,xxxP
        K[2,8],K[8,2]   = xyyP,xyyP
        K[2,10],K[10,2] = xxyP,xxyP

        K[3,3]          = xxP
        K[3,5],K[5,3]   = xyP,xyP
        K[3,7],K[7,3]   = xxxP,xxxP 
        K[3,9],K[9,3]   = xyyP,xyyP
        K[3,11],K[11,3] = xxyP,xxyP

        K[4,4]          = yyP
        K[4,6],K[6,4]   = xxyP,xxyP
        K[4,8],K[8,4]   = yyyP,yyyP
        K[4,10],K[10,4] = xyyP,xyyP

        K[5,5]          = yyP
        K[5,7],K[7,5]   = xxyP,xxyP
        K[5,9],K[9,5]   = yyyP,yyyP
        K[5,11],K[11,5] = xyyP,xyyP

        K[6,6]          = xxxxP
        K[6,8],K[8,6]   = xxyyP,xxyyP
        K[6,10],K[10,6] = xxxyP,xxxyP

        K[7,7]          = xxxxP
        K[7,9],K[9,7]   = xxyyP,xxyyP
        K[7,11],K[11,7] = xxxyP,xxxyP

        K[8,8]          = yyyyP
        K[8,10],K[10,8] = xyyyP,xyyyP

        K[9,9]          = yyyyP
        K[9,11],K[11,9] = xyyyP,xyyyP

        K[10,10]        = xxyyP

        K[11,11]        = xxyyP   
 
        Basis = [self.q0,self.q1,self.q2,self.q3,self.q4,self.q5,self.q6,self.q7,self.q8,self.q9,self.q10,self.q11]
        Num   = len(Element)
        D = np.zeros((4*Num,12),dtype=float)
        j = 0
        for q in Basis:
            oqx,oqy = 0,0
            for i in range(len(V)-1):
                Node       = V[i]
                x,y        = Node[0],Node[1]
                qx,qy      = q(x,y)
                D[i,j]     = qx
                D[i+Num,j] = qy
                oqx        = oqx+qx
                oqy        = oqy+qy

            for i in range(len(V)-1):
                Edge = E[i]
                Node1,Node2 = self.Mesh.Nodes[Edge[0]],self.Mesh.Nodes[Edge[1]]
                x1,y1,x2,y2 = Node1[0],Node1[1],Node2[0],Node2[1]
                
                qx,qy        = q(0.5*(x1+x2),0.5*(y1+y2))
                D[2*Num+i,j] = qx
                D[3*Num+i,j] = qy
                oqx          = oqx+qx
                oqy          = oqy+qy
            G[0,j] = oqx
            G[1,j] = oqy
            j= j+1 
        return H,np.linalg.inv(G),D,K,np.linalg.inv(RTK) 
        
    def TVhSemiInProdColumn(self,Element,ElementNumber,unx,uny,umx,umy,xP,yP,A,E):
        #This function will create a column vector, the first two entries 
        #in this vector will be the sum of the x coordinates and y coordinates
        #of the evaluation of a function at the nodes and midedges.
        #The rest will be the semi-inner product of the function against each
        #basis function.
        N   = len(Element)
        B   = np.zeros((12),dtype=float)

        Div,A = self.DIVu(ElementNumber,unx,uny,umx,umy)
        B[0] = np.sum(unx)+np.sum(umx)
        B[1] = np.sum(uny)+np.sum(umy)
        unx = np.append(unx,unx[0])
        uny = np.append(uny,uny[0])
        
        for i in range(N):
            Edge = E[i]
            Node1,Node2  = self.Mesh.Nodes[Edge[0]],self.Mesh.Nodes[Edge[1]]
            x1,y1,x2,y2  = Node1[0],Node1[1],Node2[0],Node2[1]
            xh,yh = 0.5*(x1+x2),0.5*(y1+y2)
            en    = [y2-y1,x1-x2]
            #Compute T1
            B[2]  = B[2]+(en[0]/6)*(unx[i]+4*umx[i]+unx[i+1]) #q3
            B[3]  = B[3]+(en[0]/6)*(uny[i]+4*umy[i]+uny[i+1]) #q4
            B[4]  = B[4]+(en[1]/6)*(unx[i]+4*umx[i]+unx[i+1]) #q5
            B[5]  = B[5]+(en[1]/6)*(uny[i]+4*umy[i]+uny[i+1]) #q6
            B[6]  = B[6]+(en[0]/3)*(x1*unx[i]+4*xh*umx[i]+x2*unx[i+1]) #q7
            B[7]  = B[7]+(en[0]/3)*(x1*uny[i]+4*xh*umy[i]+x2*uny[i+1]) #q8
            B[8]  = B[8]+(en[1]/3)*(y1*unx[i]+4*yh*umx[i]+y2*unx[i+1]) #q9
            B[9]  = B[9]+(en[1]/3)*(y1*uny[i]+4*yh*umy[i]+y2*uny[i+1]) #q10
            B[10] = B[10]+(1/6)*((y1*en[0]+x1*en[1])*unx[i]+4*(yh*en[0]+xh*en[1])*umx[i]+(y2*en[0]+x2*en[1])*unx[i+1]) #q11
            B[11] = B[11]+(1/6)*((y1*en[0]+x1*en[1])*uny[i]+4*(yh*en[0]+xh*en[1])*umy[i]+(y2*en[0]+x2*en[1])*uny[i+1]) #q12

            #Compute T3
            udn1 = unx[i]*en[0]+uny[i]*en[1]
            udnh = umx[i]*en[0]+umy[i]*en[1]
            udn2 = unx[i+1]*en[0]+uny[i+1]*en[1]

            B[6] = B[6]-(1/3)*(x1*udn1+4*xh*udnh+x2*udn2)#q7
            B[7] = B[7]-(1/3)*(y1*udn1+4*yh*udnh+y2*udn2)#q8
            B[8] = B[8]-(1/3)*(x1*udn1+4*xh*udnh+x2*udn2)#q9
            B[9] = B[9]-(1/3)*(y1*udn1+4*yh*udnh+y2*udn2)#q10
            
        #Compute T2
        B[6] = B[6]+2*A*xP*Div#q7
        B[7] = B[7]+2*A*yP*Div#q8
        B[8] = B[8]+2*A*xP*Div#q9
        B[9] = B[9]+2*A*yP*Div#q10
        
        return B

    def TVhSemiInProd(self,ElementNumber,unx,uny,umx,umy,vnx,vny,vmx,vmy):
        #This function computes the semi-inner product in TVh of u against v over the selected element.
        #The inputed DOF must be local.
        H, GI, D    = self.HSTVList[ElementNumber], self.GISTVList[ElementNumber], self.DTVList[ElementNumber]
        Element     = self.Mesh.ElementEdges[ElementNumber]
        xP,yP,A,V,E = self.Mesh.Centroid(Element,self.Mesh.Orientations[ElementNumber])
        Bu          = self.TVhSemiInProdColumn(Element,ElementNumber,unx,uny,umx,umy,xP,yP,A,E)
        Bv          = self.TVhSemiInProdColumn(Element,ElementNumber,vnx,vny,vmx,vmy,xP,yP,A,E)

        Pistaru ,Pistarv = GI.dot(Bu), GI.dot(Bv)
        Piu,     Piv     = D.dot(Pistaru), D.dot(Pistarv)
        u                = np.concatenate((unx,uny,umx,umy), axis=None)
        v                = np.concatenate((vnx,vny,vmx,vmy), axis=None)
        return np.transpose(Pistaru).dot(H.dot(Pistarv))+A*(u-Piu).dot(v-Piv)

    def TVhInProd(self,ElementNumber,unx,uny,umx,umy,vnx,vny,vmx,vmy):
        #This function computes the inner product in TVh of u against v over the selected element.
        #The inputed DOF must be local.
        K, GI, D    = self.KTVList[ElementNumber], self.GISTVList[ElementNumber], self.DTVList[ElementNumber]
        Element     = self.Mesh.ElementEdges[ElementNumber]
        xP,yP,A,V,E = self.Mesh.Centroid(Element,self.Mesh.Orientations[ElementNumber])
        Bu          = self.TVhSemiInProdColumn(Element,ElementNumber,unx,uny,umx,umy,xP,yP,A,E)
        Bv          = self.TVhSemiInProdColumn(Element,ElementNumber,vnx,vny,vmx,vmy,xP,yP,A,E)

        Pistaru ,Pistarv = GI.dot(Bu), GI.dot(Bv)
        Piu,     Piv     = D.dot(Pistaru), D.dot(Pistarv)
        u                = np.concatenate((unx,uny,umx,umy), axis=None)
        v                = np.concatenate((vnx,vny,vmx,vmy), axis=None)
        ans1             = np.transpose(Pistaru).dot( K.dot(Pistarv) )
        ans2             = A*(u-Piu).dot(v-Piv)
        return np.transpose(Pistaru).dot(K.dot(Pistarv))+A*(u-Piu).dot(v-Piv)
    
    def GetLocalTVhDOF(self,ElementNumber,Gunx,Guny,Gumx,Gumy):
        #This function will, provided 
        Element = self.Mesh.ElementEdges[ElementNumber]
        V,E     = self.Mesh.StandardElement(Element,self.Mesh.Orientations[ElementNumber])
        N       = len(Element)
        lunx    = np.zeros((N),dtype=float)
        luny    = np.zeros((N),dtype=float)
        lumx    = np.zeros((N),dtype=float)
        lumy    = np.zeros((N),dtype=float)
        for i in range(len(E)-1):
            Edge    = E[i]
            v1      = Edge[0]

            lunx[i] = Gunx[v1]
            luny[i] = Guny[v1]

        for i in range(len(E)-1):
            Edge    = Element[i]
            lumx[i] = Gumx[Edge]
            lumy[i] = Gumy[Edge]

        return lunx,luny,lumx,lumy

    def TVhL2Norm(self,unx,uny,umx,umy):
        Norm = 0
        for i in range(len(self.Mesh.ElementEdges)):
            lunx,luny,lumx,lumy = self.GetLocalTVhDOF(i,unx,uny,umx,umy)
            Norm = Norm+self.TVhInProd(i,lunx,luny,lumx,lumy,lunx,luny,lumx,lumy)
        return math.sqrt(Norm)

    def TVhH1Norm(self,unx,uny,umx,umy):
        Norm = 0
        for i in range(len(self.Mesh.ElementEdges)):
            lunx,luny,lumx,lumy = self.GetLocalTVhDOF(i,unx,uny,umx,umy)
            Norm = Norm+self.TVhSemiInProd(i,lunx,luny,lumx,lumy,lunx,luny,lumx,lumy)
        return math.sqrt(Norm) 

    def Cross2Dto1D(self,Ax,Ay,Bx,By):
        #This routine takes the evalution of two vector valued functions over the nodes of a cell
        #and returns A cross B:
        return Ax*By-Ay*Bx
    
    def Cross1Dto2D(self,J,Bx,By):
        return -J*By,J*Bx