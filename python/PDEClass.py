from MeshHelios import HeliosMesh
import math
import numpy as np

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

        self.Mesh, self.Re, self.Rm, self.dt, self.theta  = Mesh, Re, Rm, dt, theta
        #We initialize the dofs
        #The boundary values on the vel and elec fields are decoupled from the
        #internal values. Thus, we will keep track of 4 arrays, they are:
        #dofs for elec field,
        #dof of the pressure and magnetic field.
        #the dofs for the vel fields are stored in two arrays, x and y comp
        tempu  = self.NodalDOFs(Inu,self.Mesh.Nodes)
        self.ux,self.uy = self.DecompIntoCoord(tempu)
        self.B          = self.MagDOFs(InB)
        self.p          = np.zeros(len(Mesh.ElementEdges))
        self.E          = np.zeros(len(self.Mesh.Nodes))
    ##################################################################################
    ##################################################################################    
    #Initiation of Boundary Conditions/DifferentTypes of Simulations and their updates
    def SetConvTestBCAndSource(self,f,g,h,ub,Eb):
        self.f, self.g, self.h, self.ub, self.Eb = f, g, h, ub, Eb  #source terms and BC
        self.DirichletUpdateSources( self.theta*self.dt) #Initiatiates, assumes that this was done at the init time
        #The expectation is that all these functions return np.arrays

    #The following routines update the dofs of the bc and sources
    def DirichletUpdateSources(self,t):
        self.updatef(t+self.theta*self.dt)
        self.updateg(t+self.theta*self.dt) 
        self.updateh(t+self.theta*self.dt)

    def updatef(self,t):
        def dummyf(xv):
            return self.f([xv[0],xv[1],t])
        self.fdof = self.NodalDOFs(dummyf,self.Mesh.Nodes)
    
    def updateg(self,t):
        def dummyg(xv):
            return self.g([xv[0],xv[1],t])
        self.gdof = self.MagDOFs(dummyg)
    
    def updateh(self,t):
        def dummyh(xv):
            return self.h([xv[0],xv[1],t])
        self.hdof = self.NodalDOFs(dummyh,self.Mesh.Nodes)

    def DirichletupdateBC(self,t):
        def dummyubn(xv):
            return self.ub([xv[0],xv[1],t])
        def dummyubnp1(xv):
            return self.ub([xv[0],xv[1],t+self.dt*self.theta])
        def dummyEb(xv):
            return self.Eb([xv[0],xv[1],t])
        tempubn   = self.NodalDOFs(dummyubn,self.Mesh.BNodes)
        tempEb    = self.NodalDOFs(dummyEb,self.Mesh.BNodes)
        tempubnp1 = self.NodalDOFs(dummyubn,self.Mesh.BNodes)
        j = 0
        for i in self.Mesh.NumBoundaryNodes:
            self.ux[i] = tempub[j][0]
            self.uy[i] = tempub[j][1]
            self.E[i]  = tempEb[j]
            j = j+1
    ##################################################################################
    ##################################################################################    
    #Initiation of Boundary Conditions/DifferentTypes of Simulations
    def NodalDOFs(self,Func,Nodes):
        #This function computes the dof of the init cond on the vel field.
        return np.array([Func(Node) for Node in Nodes])

    def DecompIntoCoord(self,Array):
        #This function will, given a list of pairs return two arrays.
        #The first and second components.
        arrayx = [x[0] for x in Array]
        Arrayy = [x[1] for x in Array]
        return arrayx,Arrayy

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
    #These functions work as an interface with the solver class.
    def DirichletConcatenate(self):
        #This function returns an array that concatenates all the unknowns
        intux = [self.ux[i] for i in self.Mesh.NumInternalNodes]
        intuy = [self.uy[i] for i in self.Mesh.NumInternalNodes]
        intE  = [self.E[i] for i in self.Mesh.NumInternalNodes]
        return np.concatenate((intux,intuy,self.B,intE,self.p), axis=None)
    
    def NumDirichletDOF(self):
        return 3*len(self.Mesh.NumInternalNodes)+len(self.Mesh.EdgeNodes)+len(self.Mesh.ElementEdges)

    def DirichletUpdateInterior(self,x):
        cut1 = len(self.Mesh.NumInternalNodes) #Number of internal dof for ux
        cut2 = 2*cut1                          #Number of internal dof for uy
        cut3 = cut2+len(self.B)                #Number of dofs for B
        cut4 = cut3+cut1                       #The number of internal dofs for E is the same as for the vel field.
        Cutx = np.split(x,[cut1,cut2,cut3,cut4])
        j = 0
        for i in self.Mesh.NumInternalNodes:
            self.ux[i] = Cutx[0][j]
            self.uy[i] = Cutx[1][j]
            self.E[i]  = Cutx[3][j]
            j = j+1
        
        self.B = Cutx[2]
        self.p = Cutx[4]
        
    def GDirichlet(self,x):
        return x
    