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
        tempun  = self.NodalDOFs(Inu,self.Mesh.Nodes)
        self.unx,self.uny  = self.DecompIntoCoord(tempun)
        tempum  = self.NodalDOFs(Inu,self.Mesh.MidNodes)
        self.umx, self.umy = self.DecompIntoCoord(tempum)
        self.B             = self.MagDOFs(InB)
        self.p             = np.zeros(len(Mesh.ElementEdges))
        self.E             = np.zeros(len(self.Mesh.Nodes))

        self.MEList     = []
        self.MVList     = []
        self.HSTVList   = []
        self.GISTVList  = []
        self.DTVList    = []
        for i in range(len(self.Mesh.ElementEdges)):
            tempME,tempMV       = self.ElecMagStandMassMat(self.Mesh.ElementEdges[i],self.Mesh.Orientations[i])
            TempSH,TempGI,TempD = self.TVhSemiInnerPreCompute(i)
            self.MEList.append(tempME)
            self.MVList.append(tempMV)
            self.HSTVList.append(TempSH)
            self.GISTVList.append(TempGI)
            self.DTVList.append(TempD)
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
            return self.Eb([xv[0],xv[1],t+self.theta*self.dt])
   
        tempubn               = self.NodalDOFs(dummyubn,self.Mesh.BNodes)
        tempubnx,tempubny     = self.DecompIntoCoord(tempubn)
        #tempubnp1             = self.NodalDOFs(dummyubnp1,self.Mesh.BNodes)
        #tempubnp1x,tempubnp1y = self.DecompIntoCoord(temp1ubnp1)
        tempEb                = self.NodalDOFs(dummyEb,self.Mesh.BNodes)
        
        j = 0
        for i in self.Mesh.NumBoundaryNodes:
            self.unx[i] = tempubnx[j]
            self.uny[i] = tempubny[j]
            self.E[i]  = tempEb[j]
            j = j+1
        
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
        intunx = [self.unx[i] for i in self.Mesh.NumInternalNodes]
        intuny = [self.uny[i] for i in self.Mesh.NumInternalNodes]
        intumx = [self.umx[i] for i in self.Mesh.NumInternalMidNodes]
        intumy = [self.umy[i] for i in self.Mesh.NumInternalMidNodes]
        intE   = [self.E[i] for i in self.Mesh.NumInternalNodes]
        return np.concatenate((intunx,intuny,intumx,intumy,self.B,intE,self.p), axis=None)
    
    def NumDirichletDOF(self):
        return 3*len(self.Mesh.NumInternalNodes)+3*len(self.Mesh.EdgeNodes)+len(self.Mesh.ElementEdges)

    def DirichletUpdateInterior(self,x):
        cut1 = len(self.Mesh.NumInternalNodes) #Number of internal dofs for ux
        cut2 = 2*cut1                          #Number of internal dofs for uy
        cut3 = cut2+len(self.Mesh.NumInternalMidNodes)                #Number of internal dofs for umx
        cut4 = cut3+len(self.Mesh.NumInternalMidNodes)                #Number of internal dofs for umy
        cut5 = cut4+len(self.B)                #Number of dofs for B
        cut6 = cut5+cut1                       #The number of internal dofs for E is the same as for the vel field.
        Cutx = np.split(x,[cut1,cut2,cut3,cut4,cut5,cut6])
        j = 0
        for i in self.Mesh.NumInternalNodes:
            self.unx[i] = Cutx[0][j]
            self.uny[i] = Cutx[1][j]
            self.E[i]   = Cutx[5][j]
            j = j+1
        j = 0
        for i in self.Mesh.NumInternalMidNodes:
            self.umx[i] = Cutx[2][j]
            self.umy[i] = Cutx[3][j]
            j = j+1
        self.B = Cutx[4]
        self.p = Cutx[6]
        
    def GDirichlet(self,x):
        #The x is passed because the Scipy Linear Function class requires it.
        #It will use the current values of the internal variables
        for i in self.Mesh.NumInternalNodes:
            CellNums = self.Mesh.NodestoCells[i]
            
            

    
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

    def BDivSquared(self):
        #This function computes the divergence of B.
        D      = 0
        DivMat = self.BDiv()
        divB   = DivMat.dot(self.B)
        for k in range(len(divB)):
            D = D + k 
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

    def PhInProd(self,ElementNumber,dof):
        #This function integrates two function in Ph over the provided element.
        #The first function is p and the second is has dof value over this cell as provided.
        Element = self.Mesh.ElementEdges[ElementNumber]
        A,V,E = self.Mesh.Area(Element,self.Mesh.Orientations[ElementNumber])
        return dof*A*self.p[ElementNumber]
    
    def DIVu(self,ElementNumber,unx,uny,umx,umy):
        #This routine computes the divergence of u over the element provided.
        Element = self.Mesh.ElementEdges[ElementNumber]
        A,V,E   = self.Mesh.Area(Element,self.Mesh.Orientations[ElementNumber])
        unx     = np.append(unx,unx[0])
        uny     = np.append(uny,uny[0])
        S, k = 0, 0
        for i in range(len(Element)):
            n1, n2 = E[k][0], E[k+1][1]
            v1, v2 = self.Mesh.Nodes[n1], self.Mesh.Nodes[n2]
            x1, x2 = v1[0], v2[0]
            y1, y2 = v1[1], v2[1]
            
            etimesnormal = [y2-y1,x1-x2]
            
            S = S+(unx[i]+unx[i+1]+4*umx[i])*etimesnormal[0]
            S = S+(uny[i]+uny[i+1]+4*umy[i])*etimesnormal[1]
            k = k+1
        return S/(6*A)     
    
    def TVhSemiInnerPreCompute(self, ElementNumber):
        #This function will compute one of the matrices that make the Semi-inner product
        Element     = self.Mesh.ElementEdges[ElementNumber]
        xP,yP,A,V,E = self.Mesh.Centroid(Element,self.Mesh.Orientations[ElementNumber])
        xyP         = 0
        for i in range(len(V)-1):
            Node1 = V[i]
            Node2 = V[i+1]
            x1,y1 = Node1[0],Node1[1]
            x2,y2 = Node2[0],Node2[1]
            
            xh1,yh1 = (x1+x2)/2,(y1+y2)/2
            xh2,yh2 = (x1+xP)/2,(y1+yP)/2
            xh3,yh3 = (x2+xP)/2,(y2+yP)/2

            AT  = 0.5*abs( (x2-xP)*(y1-yP)-(y2-yP)*(x2-xP) )
            xyP = xyP+(AT/3)*(xh1*yh1+xh2*yh2+xh3*yh3)

        H = np.zeros((12,12),dtype=float)
        G = np.zeros((12,12),dtype=float)
        H[2,6], H[6,2]  = 2*xP*A,2*xP*A
        H[2,10],H[10,2] = yP*A,yP*A

        H[3,7], H[7,3]  = 2*xP*A,2*xP*A
        H[3,11],H[11,3] = yP*A,yP*A

        H[4,8], H[8,4]  = 2*yP*A,2*yP*A
        H[4,10],H[10,4] = xP*A,xP*A 

        H[5,9],H[9,5]   = 2*yP*A,2*yP*A
        H[5,11],H[11,5] = xP*A,xP*A

        H[6,10],H[10,6] = 2*xyP,2*xyP

        H[7,11],H[11,7] = 2*xyP,2*xyP

        H[8,10],H[10,8] = 2*xyP,2*xyP

        H[9,11],H[11,9] = 2*xyP,2*xyP

        G[2,6], G[6,2]  = 2*xP*A,2*xP*A
        G[2,10],G[10,2] = yP*A,yP*A

        G[3,7], G[7,3]  = 2*xP*A,2*xP*A
        G[3,11],G[11,3] = yP*A,yP*A

        G[4,8], G[8,4]  = 2*yP*A,2*yP*A
        G[4,10],G[10,4] = xP*A,xP*A 

        G[5,9],G[9,5]   = 2*yP*A,2*yP*A
        G[5,11],G[11,5] = xP*A,xP*A

        G[6,10],G[10,6] = 2*xyP,2*xyP

        G[7,11],G[11,7] = 2*xyP,2*xyP

        G[8,10],G[10,8] = 2*xyP,2*xyP

        G[9,11],G[11,9] = 2*xyP,2*xyP

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
        return H,np.linalg.inv(G),D 

    def TVhSemiInProdColumn(self,Element,ElementNumber,unx,uny,umx,umy,xP,yP,A,E):
        #This function will create a column vector, the first two entries 
        #in this vector will be the sum of the x coordinates and y coordinates
        #of the evaluation of a function at the nodes and midedges.
        #The rest will be the semi-inner product of the function against each
        #basis function.
        N   = len(Element)
        B   = np.zeros((12),dtype=float)

        Div = self.DIVu(ElementNumber,unx,uny,umx,umy)
        unx = np.append(unx,unx[0])
        uny = np.append(uny,uny[0])

        B[0] = np.sum(unx)+np.sum(umx)
        B[1] = np.sum(uny)+np.sum(umy)
        
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

            B[6] = B[6]-(1/3)*(x1*udn1+4*xh*udnh+x2*udn2)
            B[7] = B[7]-(1/3)*(y1*udn1+4*yh*udnh+y2*udn2)
            B[8] = B[8]-(1/3)*(x1*udn1+4*xh*udnh+x2*udn2)
            B[9] = B[9]-(1/3)*(y1*udn1+4*yh*udnh+y2*udn2)
        
        #Compute T2
        B[6] = B[6]+2*A*xP*Div
        B[7] = B[7]+2*A*yP*Div
        B[8] = B[8]+2*A*xP*Div
        B[9] = B[9]+2*A*yP*Div
        return B

    def TVhSemiInProd(self,ElementNumber,unx,uny,umx,umy,vnx,vny,vmx,vmy):
        #This function computes the semi-inner product in TVh of u against v over the selected element.
        H, GI, D    = self.HSTVList[ElementNumber], self.GISTVList[ElementNumber], self.DTVList[ElementNumber]
        Element     = self.Mesh.ElementEdges[ElementNumber]
        xP,yP,A,V,E = self.Mesh.Centroid(Element,self.Mesh.Orientations[ElementNumber])
        Bu          = self.TVhSemiInProdColumn(Element,ElementNumber,unx,uny,umx,umy,xP,yP,A,E)
        Bv          = self.TVhSemiInProdColumn(Element,ElementNumber,vnx,vny,vmx,vmy,xP,yP,A,E)
        Element     = self.Mesh.ElementEdges[ElementNumber]

        Pistaru ,Pistarv = GI.dot(Bu), GI.dot(Bv)
        Piu,     Piv     = D.dot(Pistaru), D.dot(Pistarv)
        u                = np.concatenate((unx,uny,umx,umy), axis=None)
        v                = np.concatenate((vnx,vny,vmx,vmy), axis=None)
        aprox = np.transpose(Pistaru).dot(H).dot(Pistarv)
        stab  = A*(u-Piu).dot(v-Piv)
        print(Pistaru)
        print('printed')
        #return np.transpose(Pistaru).dot(H).dot(Pistarv)+A*(u-Piu).dot(v-Piv)
        return aprox,stab
    
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