from scipy.sparse.linalg import LinearOperator
from scipy.sparse.linalg import gmres
from scipy.sparse.linalg import lgmres
from scipy.linalg import det
from numpy.linalg import matrix_rank as rank
import multiprocessing as mp
from numpy.linalg import norm as n2
import numpy as np
import math
import time
from scipy.sparse.linalg import spsolve
from scipy import linalg

def outside_func(par):
    obj,G,Gxm,xm,ndof,i = par
    return obj.ithCol(G,Gxm,xm,ndof,i)

class InexactNewtonTimeInt(object):
    def __init__(self):
        #Series of constants for the error in the GMRES
        self.eps    = 1E-7
        self.etamax = 0.8
        self.alpha  = 1.5
        self.gamma  = 0.9
        self.epsr   = 1E-4
        #self.pool   = mp.Pool(12)
    def J(self,cols):
        ndof = len(cols[0])
        J = np.zeros((ndof,ndof),dtype = float)
        for i in range(ndof):
            J[:,i] = cols[i]
        return J

    def ithCol(self,G,Gxm,xm,ndof,i):
        unitnormal = np.zeros(ndof,dtype=float)
        unitnormal[i] = 1
        return (G(xm+self.eps*unitnormal)-Gxm)/self.eps

    def FlowSolve(self,G,x0,ndof,maxiter,tol):
        xm     = x0
        delxm = 0.0 * xm
        Gxm    = G(x0)
        nGxm   = n2(Gxm)
        
        j = 0

        while nGxm>tol and j<maxiter:
            #print('ngxm='+str(nGxm))
            Cols   = []
            for i in range(ndof):
               col = self.ithCol(G,Gxm,xm,ndof,i)
               #print(f'norm of{i}-th Node={n2(col)}')
               Cols.append(col)
            J    = self.J(Cols)
            #for i in range(ndof):
            #    col = J[i,:]
            #    print(f'norm of{i}-th row={n2(col)}')

            cond = np.linalg.cond(J)
            print('Cond #='+str(cond))
            #print('shape='+str(J.shape))
            #print('det ='+str(det(J))) 
            #print('rank='+str(rank(J)))

            #delxm = spsolve(J,-Gxm)

            delxm = linalg.solve(J, -Gxm)
            #def fDGxm(delx):
            #    return (G(xm+self.eps*delx)-Gxm)/(self.eps)
            #DGxm  = LinearOperator((ndof,ndof), matvec = fDGxm)
            #delxm, exitcode = gmres(DGxm,-Gxm,tol=tol/10.0,atol=tol/10.0,x0=delxm)
            #print(delxm)
            #delxm, exitcode = gmres(J,-Gxm,tol=tol/10.0,atol=tol/10.0,x0=delxm)
            #print(exitcode)
            xm              = xm + delxm
            Gxm             = G(xm)
            nGxm            = n2(Gxm)
            j               = j+1

            #print('ngxm='+str(nGxm))

            #if exitcode>1E-5:
            #    print('error ocurred, exitcode='+str(exitcode))
        return xm

    #Third Attepmt
    def Newtoniter(self,G,x0,ndof,tol,maxiter,PDE,unx, uny,umx,umy,B,E,p):
        #This function will perform a Newton Iteration
        #To find the zeroes of the function G provided the initial guess x0
        #Within the tolerance tol in the 2-norm. ndof is the number of unknowns.
        xm     = x0
        delxm = 0.0 * xm
        print('here2')
        Gxm    = G(x0)
        print('here3')
        etamm1 = self.etamax
        nGxmm1 = 1
        epsa   = math.sqrt(ndof)*(10**(-15))
        epst   = epsa+self.epsr*n2(G(x0))
        nGxm = n2(G(xm))
        print('Entered')
        j = 0
        maxiter = 1
        while nGxm>tol and j<maxiter:
            print('err before gmres='+str(nGxm))  
            

            etamA = self.gamma*(nGxm/nGxmm1)**(self.alpha)
            etamB = min([self.etamax,max([etamA,self.gamma*etamm1**self.alpha])])
            etam  = min([self.etamax,max([etamB,self.gamma*(epst/nGxm)])])
            Cols  = []
            for i in range(ndof):
               col = self.ithCol(G,Gxm,xm,ndof,i)
               Cols.append(col)
            
            J               = self.J(Cols)
            #print(J)
            #print('Number of Newton Iterations='+str(j))
            print('shape='+str(J.shape))
            print('det ='+str(det(J))) 
            print('rank='+str(rank(J)))
            j = j+1
            
            # delxm, exitcode = gmres(J,-Gxm,tol=tol/10.0,atol=tol/10.0,x0=delxm)
            # #delxm, exitcode = gmres(DGxm,-Gxm,tol=tol/10.0,atol=tol/10.0,x0=delxm)
            # xm              = xm + delxm
            # Gxm             = G(xm)
                
            # etamm1 = etam
            # nGxmm1 = nGxm
            # nGxm = n2(Gxm)
            # print('err after gmres='+str(nGxm))
            # if exitcode>1E-5:
            #     print('GMRES finished without reaching tolerance')
            #     print('Num of GMRES iterations='+str(exitcode))
                
            # if nGxm<tol:
            #     y = Gxm
            #     fnx,fny,fmx,fmy,farf,elecf,divf = PDE.MHDSplity(y)
        
            #     momerr  = PDE.TVhL2Norm(fnx,fny,fmx,fmy)
            #     Farerr  = PDE.EhL2Norm(farf)
            #     Elecerr = PDE.VhL2Norm(elecf)
            #     #print('momerr='+str(momerr))
            #     #print('Farerr='+str(Farerr))
            #     #print('Elecerr='+str(Elecerr))



            #     #print('err before div cleaning'+str(nGxm))
            #     uny,uny,umx,umy,B,E,p = PDE.MHDUpdateInt(xm,unx,uny,umx,umy,B,E,p)
            #     divB = PDE.BDivSquared(B)
            #     #print('Div before cleaning'+str(divB))
            #     B1   = B
            #     B    = PDE.B-PDE.dt*PDE.MRot.dot(E)
            #     u    = (PDE.MRot).dot(E)+(B-PDE.B)/PDE.dt

            #     divB = PDE.BDivSquared(B)

            #     xm                    = PDE.MHDConcatenate(unx,uny,umx,umy,B,E,p)
            #     uny,uny,umx,umy,B,E,p = PDE.MHDUpdateInt(xm,unx,uny,umx,umy,B,E,p)

            #     errB = n2(B1-B)

            #     y    = G(xm)
            #     nGxm = n2(y)
            #     Gxm  = y

            #     #print('err after div cleaning'+str(nGxm))
                
            #     fnx,fny,fmx,fmy,farf,elecf,divf = PDE.MHDSplity(y)
        
            #     momerr  = PDE.TVhL2Norm(fnx,fny,fmx,fmy)
            #     Farerr  = PDE.EhL2Norm(farf)
            #     Elecerr = PDE.VhL2Norm(elecf)
            #     #print('momerr='+str(momerr))
            #     #print('Farerr='+str(Farerr))
            #     #print('Elecerr='+str(Elecerr))
            # i += 1
        #if i == maxiter:
        #    print('Surpassed max number of iter without arriving at sol')
        #else:
        #    print('Successfully completed Newton iterations')
        #return xm  
    #Third attempt

    #Second Attempt
    # def Newtoniter(self,G,x0,ndof,tol,maxiter,PDE,unx, uny,umx,umy,B,E,p):
    #     #This function will perform a Newton Iteration
    #     #To find the zeroes of the function G provided the initial guess x0
    #     #Within the tolerance tol in the 2-norm. ndof is the number of unknowns.
    #     xm     = x0
    #     delxm = 0.0 * xm
    #     Gxm    = G(x0)
    #     etamm1 = self.etamax
    #     nGxmm1 = 1
    #     epsa   = math.sqrt(ndof)*(10**(-15))
    #     epst   = epsa+self.epsr*n2(G(x0))
    #     nGxm = n2(G(xm))
    #     i = 0
    #     while nGxm>tol and i<maxiter:
    #         print('err before gmres='+str(nGxm))  
            
    #         # def fDGxm(delx):
    #         #     sum = 0.0*Gxm
    #         #     for i in range(ndof):
    #         #         unitnormal = np.zeros(ndof)
    #         #         unitnormal[i] = 1.0
    #         #         sum = sum + delx[i]*( Gxm-G(xm-self.eps*unitnormal) ) / (self.eps)
    #         #     return sum
    #         #def fDGxm(delx):
    #         #    return (G(xm+self.eps*delx)-Gxm)/(self.eps)
    #         # def fDGxm(delx):
    #         #     ndelx = n2(delx)
    #         #     if ndelx<1E-5:
    #         #         return delx*0
    #         #     else:
    #         #         eps   = 1E-7/ndelx
    #         #         return (G(xm+eps*delx)-Gxm)/(eps)
    #         #DGxm  = LinearOperator((ndof,ndof), matvec = fDGxm)

    #         etamA = self.gamma*(nGxm/nGxmm1)**(self.alpha)
    #         etamB = min([self.etamax,max([etamA,self.gamma*etamm1**self.alpha])])
    #         etam  = min([self.etamax,max([etamB,self.gamma*(epst/nGxm)])])
    #         Cols  = []
    #         for i in range(ndof):
    #            col = self.ithCol(G,Gxm,xm,ndof,i)
    #            Cols.append(col)
    #         # input_lst = []
    #         # for i in range(ndof):
    #         #     par = (self,G,Gxm,xm,ndof,i)
    #         #     input_lst.append(par) 
    #         # Cols            = self.pool.map(outside_func,input_lst)
    #         J               = self.J(Cols)
    #         delxm, exitcode = gmres(J,-Gxm,tol=tol/10.0,atol=tol/10.0,x0=delxm)
    #         #delxm, exitcode = gmres(DGxm,-Gxm,tol=tol/10.0,atol=tol/10.0,x0=delxm)
    #         xm              = xm + delxm
    #         Gxm             = G(xm)
                
    #         etamm1 = etam
    #         nGxmm1 = nGxm
    #         nGxm = n2(Gxm)
    #         print('err after gmres='+str(nGxm))
    #         if exitcode>1E-5:
    #             print('GMRES finished without reaching tolerance')
    #             print('Num of GMRES iterations='+str(exitcode))
                
    #         if nGxm<tol:
    #             y = Gxm
    #             fnx,fny,fmx,fmy,farf,elecf,divf = PDE.MHDSplity(y)
        
    #             momerr  = PDE.TVhL2Norm(fnx,fny,fmx,fmy)
    #             Farerr  = PDE.EhL2Norm(farf)
    #             Elecerr = PDE.VhL2Norm(elecf)
    #             #print('momerr='+str(momerr))
    #             #print('Farerr='+str(Farerr))
    #             #print('Elecerr='+str(Elecerr))



    #             #print('err before div cleaning'+str(nGxm))
    #             uny,uny,umx,umy,B,E,p = PDE.MHDUpdateInt(xm,unx,uny,umx,umy,B,E,p)
    #             divB = PDE.BDivSquared(B)
    #             #print('Div before cleaning'+str(divB))
    #             B1   = B
    #             B    = PDE.B-PDE.dt*PDE.MRot.dot(E)
    #             u    = (PDE.MRot).dot(E)+(B-PDE.B)/PDE.dt

    #             divB = PDE.BDivSquared(B)

    #             xm                    = PDE.MHDConcatenate(unx,uny,umx,umy,B,E,p)
    #             uny,uny,umx,umy,B,E,p = PDE.MHDUpdateInt(xm,unx,uny,umx,umy,B,E,p)

    #             errB = n2(B1-B)

    #             y    = G(xm)
    #             nGxm = n2(y)
    #             Gxm  = y

    #             #print('err after div cleaning'+str(nGxm))
                
    #             fnx,fny,fmx,fmy,farf,elecf,divf = PDE.MHDSplity(y)
        
    #             momerr  = PDE.TVhL2Norm(fnx,fny,fmx,fmy)
    #             Farerr  = PDE.EhL2Norm(farf)
    #             Elecerr = PDE.VhL2Norm(elecf)
    #             #print('momerr='+str(momerr))
    #             #print('Farerr='+str(Farerr))
    #             #print('Elecerr='+str(Elecerr))
    #         i += 1
    #     if i == maxiter:
    #         print('Surpassed max number of iter without arriving at sol')
    #     else:
    #         print('Successfully completed Newton iterations')
    #     return xm  
        
    # def Newtoniter(self,G,x0,ndof,tol,maxiter,PDE,unx, uny,umx,umy,B,E,p):
    #     #This function will perform a Newton Iteration
    #     #To find the zeroes of the function G provided the initial guess x0
    #     #Within the tolerance tol in the 2-norm. ndof is the number of unknowns.
    #     xm     = x0
    #     #print(x0)
    #     #print('eval G')
    #     Gxm    = G(x0)
    #     etamm1 = self.etamax
    #     nGxmm1 = 1
    #     epsa   = math.sqrt(ndof)*(10**(-15))
    #     epst   = epsa+self.epsr*n2(G(x0))
    #     for i in range(maxiter):
    #         nGxm = n2(G(xm))
    #         if abs(nGxm)<tol:
    #             return xm
    #         else:
    #             def fDGxm(delx):
    #                 return (G(xm+self.eps*delx)-Gxm)/(self.eps)
    #             DGxm  = LinearOperator((ndof,ndof), matvec = fDGxm)
    #             etamA = self.gamma*(nGxm/nGxmm1)**(self.alpha)
    #             etamB = min([self.etamax,max([etamA,self.gamma*etamm1**self.alpha])])
    #             etam  = min([self.etamax,max([etamB,self.gamma*(epst/nGxm)])])
                
    #             PDE.evalcount   = 0
    #             uny,uny,umx,umy,B,E,p = PDE.MHDUpdateInt(xm,unx,uny,umx,umy,B,E,p)
    #             divB = PDE.BDivSquared(B)
    #             print('Number of Newton iter '+str(i))
    #             print('Current err on Newton Iter before GMRES='+str(nGxm))
    #             print('Current square of the L2 Norm on DivB before GMRES ='+str(divB))
    #             #print('Tolerance on GMRES = '+str(etam*nGxm))
    #             start = time.time()
    #             delxm, exitcode = gmres(DGxm,-Gxm,tol=etam*nGxm)#tol = 1e-3,maxiter=5000000)
    #             end   = time.time() 
    #             #delxm, exitcode = lgmres(DGxm,-Gxm, tol = 1E-8)#tol=etam*nGxm)
                
    #             #print('TimeFor1gmres'+str(end-start))
    #             #print('evalcount='+str(PDE.evalcount))
    #             uny,uny,umx,umy,B,E,p                      = PDE.MHDUpdateInt(xm,unx,uny,umx,umy,B,E,p)
    #             delunx,deluny,delumx,delumy,delB,delE,delp = PDE.MHDUpdateInt(delxm,unx,uny,umx,umy,B,E,p)
                
    #             delB  = PDE.B-B-PDE.dt*((PDE.MRot).dot(E+delE))
    #             delxm = PDE.MHDConcatenate(delunx,deluny,delumx,delumy,delB,delE,delp)
    #             xm    = xm+delxm
    #             nGxm  = n2(G(xm))
    #             uny,uny,umx,umy,B,E,p = PDE.MHDUpdateInt(xm,unx,uny,umx,umy,B,E,p)

    #             divB = PDE.BDivSquared(B)
    #             print('Current err on Newton Iter after GMRES='+str(nGxm))
    #             print('Current square of theL2 Norm on DivB after GMRES ='+str(divB))
    #             Gxm             = G(xm)
    #             etamm1 = etam
    #             nGxmm1 = nGxm

                
                
    #             if exitcode>1E-5:
    #                 print('GMRES finished without reaching tolerance')
    #                 print('Num of GMRES iterations='+str(exitcode))
    #     print('Surpassed max number of iter without arriving at sol')