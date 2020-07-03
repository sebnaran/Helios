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
        #Finally the array of dofs for the vel field is has length 2,
        #self.u[0] is the set of x values and self.u[1] is the set of y values
        [self.ux,self.uy] = self.NodalDOFs(Inu,self.Mesh.Nodes)
        self.B            = self.MagDOFs(InB)
        self.p            = np.zeros(len(Mesh.ElementEdges))
        self.E            = np.zeros(len(self.Mesh.Nodes))
    ##################################################################################
    ##################################################################################    
    #Initiation of Boundary Conditions/DifferentTypes of Simulations and their updates
    def SetConvTestBCAndSource(self,f,g,h,ub,Eb):
        self.f, self.g, self.h, self.ub, self.Eb = f, g, h, ub, Eb  #momentum eq source 
        self.updatef(self.theta*self.dt)
        self.updateg(self.theta*self.dt)
        self.updateh(self.theta*self.dt)
        #The expectation is that all these functions return np.arrays

    #The following routines update the dofs of the bc and sources
    def updatef(self,t):
        def dummyf(xv):
            return self.f([xv[0],xv[1],t])
        [self.fdofx,self.fdofy] = self.NodalDOFs(dummyf,self.Mesh.Nodes)
    
    def updateg(self,t):
        def dummyg(xv):
            return self.g([xv[0],xv[1],t])
        self.gdof = self.MagDOFs(dummyg)
    
    def updateh(self,t):
        def dummyh(xv):
            return self.h([xv[0],xv[1],t])
        self.hdof = self.NodalDOFs(dummyh,self.Mesh.Nodes)

    def updateub(self,t):
        def dummyub(xv):
            return self.ub([xv[0],xv[1],t])
    
    def updateEb(self,t):
        def dummyEb(xv):
            return self.Eb([xv[0],xv[1],t])

    ##################################################################################
    ##################################################################################    
    #Initiation of Boundary Conditions/DifferentTypes of Simulations
    def NodalDOFs(self,Func,Nodes):
        #This function computes the dof of the init cond on the vel field.
        return Func(Nodes)

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
    def Concatenate(self):
        x = np.concatenate(self.ux,self.uy,self.B,self.E,self.p)
        return x
    
    def NumUnknownDOFDirichlet(self):
        return (3*len(self.ux)+3*len(self.Eb)+len(self.p)+len(self.B))

    def GDirichlet(self,x):
        return x
