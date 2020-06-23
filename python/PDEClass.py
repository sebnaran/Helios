from MeshHelios import HeliosMesh
import math
import numpy as np

class PDEFullMHD(object):
    def __init__(self,Mesh,Re,Rm,Inu,InB,dt,T):
        self.Mesh = Mesh
        self.Mesh.MakeDictionaries()
        self.Re   = Re
        self.Rm   = Rm
        self.u    = Inu #initial velocity
        self.T    = T   #Final time
        self.dt   = dt  #time step size

        #We initialize the dofs
        N1              = len(Mesh.Nodes)
        N2              = len(Mesh.ElementEdges)
        self.p          = np.zeros(N1) 
        self.E          = np.zeros(N2)
        self.B          = self.CIMagDOFs(InB) #initial mag field
        self.ux,self.uy = self.CIVelDOFs(Inu)         #initial vel field
        

    def ConvTestBoundarySourceCond(self,f,g,h,Bb,Eb):
        self.f  = f  #momentum eq source 
        self.g  = g  #Faraday source
        self.h  = h  #Ampere-Ohm source
        self.Bb = Bb #Boundary Cond
        self.Eb = Eb

    def CIMagDOFs(self,Func):
    #This computes the dofs of the initial magnetic field
        N    = len(self.Mesh.EdgeNodes)
        proj = np.zeros(N)

        for i in range(N):
            x1      = self.Mesh.Nodes[self.Mesh.EdgeNodes[i][0]][0]
            y1      = self.Mesh.Nodes[self.Mesh.EdgeNodes[i][0]][1]
            x2      = self.Mesh.Nodes[self.Mesh.EdgeNodes[i][1]][0]
            y2      = self.Mesh.Nodes[self.Mesh.EdgeNodes[i][1]][1]
            lengthe = math.sqrt((x2-x1)**2+(y2-y1)**2)
            etimesnormal = [y2-y1,x1-x2]

            pt0 = -1
            w0  = 1/21

            pt1 = -math.sqrt((5/11)+(2/11)*math.sqrt(5/3))
            w1  = (124-7*math.sqrt(15))/350

            pt2 = -math.sqrt((5/11)-(2/11)*math.sqrt(5/3))
            w2  = (124+7*math.sqrt(15))/350

            pt3 = 0
            w3  = 256/525

            pt4 = math.sqrt((5/11)-(2/11)*math.sqrt(5/3))
            w4  = (124+7*math.sqrt(15))/350

            pt5 = math.sqrt((5/11)+(2/11)*math.sqrt(5/3))
            w5  = (124-7*math.sqrt(15))/350

            pt6 = 1
            w6  = 1/21
        
            Fx0,Fy0 = Func(x1,y1)
            Fx1,Fy1 = Func((x1*(1-pt1)+x2*(1+pt1))/2,(y1*(1-pt1)+y2*(1+pt1))/2)
            Fx2,Fy2 = Func((x1*(1-pt2)+x2*(1+pt2))/2,(y1*(1-pt2)+y2*(1+pt2))/2)
            Fx3,Fy3 = Func(0.5*(x1+x2),0.5*(y1+y2))
            Fx4,Fy4 = Func((x1*(1-pt4)+x2*(1+pt4))/2,(y1*(1-pt4)+y2*(1+pt4))/2)
            Fx5,Fy5 = Func((x1*(1-pt5)+x2*(1+pt5))/2,(y1*(1-pt5)+y2*(1+pt5))/2)
            Fx6,Fy6 = Func(x2,y2)

            proj[i] = (w0*Fx0+w1*Fx1+w2*Fx2+w3*Fx3+w4*Fx4+w5*Fx5+w6*Fx6)*etimesnormal[0]
            proj[i] = proj[i] +(w0*Fy0+w1*Fy1+w2*Fy2+w3*Fy3+w4*Fy4+w5*Fy5+w6*Fy6)*etimesnormal[1]
            proj[i] = proj[i]/(2*lengthe)   
        return proj

    def CIVelDOFs(self,Func):
        N  = len(self.Mesh.Nodes)
        ux = np.zeros(N)
        uy = np.zeros(N)
        for i in range(N):
            Node  = self.Mesh.Nodes[i]
            x     = Node[0]
            y     = Node[1]

            insux,insuy = Func(x,y)
            ux[i] = insux
            uy[i] = insuy
        return ux,uy
    
    def Concatenate(self):
        x = np.concatenate(self.ux,self.uy,self.B,self.E,self.p)
        return x