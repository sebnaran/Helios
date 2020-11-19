from scipy.sparse.linalg import LinearOperator
from scipy.sparse.linalg import gmres
from scipy.sparse.linalg import lgmres
from numpy.linalg import norm as n2
import numpy as np
import math
import time

class InexactNewtonTimeInt(object):
    def __init__(self):
        #Series of constants for the error in the GMRES
        self.eps    = 1E-7
        self.etamax = 0.8
        self.alpha  = 1.5
        self.gamma  = 0.9
        self.epsr   = 1E-4

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

    def Newtoniter(self,G,x0,ndof,tol,maxiter,PDE,unx, uny,umx,umy,B,E,p):
        #This function will perform a Newton Iteration
        #To find the zeroes of the function G provided the initial guess x0
        #Within the tolerance tol in the 2-norm. ndof is the number of unknowns.
        xm     = x0

        Gxm    = G(x0)
        etamm1 = self.etamax
        nGxmm1 = 1
        epsa   = math.sqrt(ndof)*(10**(-15))
        epst   = epsa+self.epsr*n2(G(x0))
        nGxm = n2(G(xm))
        i = 0
        while nGxm>tol and i<maxiter:
            print('err before gmres='+str(nGxm))  
            def fDGxm(delx):
                return (G(xm+self.eps*delx)-Gxm)/(self.eps)
            DGxm  = LinearOperator((ndof,ndof), matvec = fDGxm)
            etamA = self.gamma*(nGxm/nGxmm1)**(self.alpha)
            etamB = min([self.etamax,max([etamA,self.gamma*etamm1**self.alpha])])
            etam  = min([self.etamax,max([etamB,self.gamma*(epst/nGxm)])])
                
                
            delxm, exitcode = gmres(DGxm,-Gxm,tol=self.eps)#tol=etam*nGxm)
            xm              = xm + delxm
            Gxm             = G(xm)
                
            etamm1 = etam
            nGxmm1 = nGxm
            nGxm = n2(Gxm)
            print('err after gmres='+str(nGxm))
            if exitcode>1E-5:
                print('GMRES finished without reaching tolerance')
                print('Num of GMRES iterations='+str(exitcode))
                
            if nGxm<tol:
                y = Gxm
                fnx,fny,fmx,fmy,farf,elecf,divf = PDE.MHDSplity(y)
        
                momerr  = PDE.TVhL2Norm(fnx,fny,fmx,fmy)
                Farerr  = PDE.EhL2Norm(farf)
                Elecerr = PDE.VhL2Norm(elecf)
                #print('momerr='+str(momerr))
                #print('Farerr='+str(Farerr))
                #print('Elecerr='+str(Elecerr))



                #print('err before div cleaning'+str(nGxm))
                uny,uny,umx,umy,B,E,p = PDE.MHDUpdateInt(xm,unx,uny,umx,umy,B,E,p)
                divB = PDE.BDivSquared(B)
                #print('Div before cleaning'+str(divB))
                B1   = B
                B    = PDE.B-PDE.dt*PDE.MRot.dot(E)
                u    = (PDE.MRot).dot(E)+(B-PDE.B)/PDE.dt

                divB = PDE.BDivSquared(B)

                xm                    = PDE.MHDConcatenate(unx,uny,umx,umy,B,E,p)
                uny,uny,umx,umy,B,E,p = PDE.MHDUpdateInt(xm,unx,uny,umx,umy,B,E,p)

                errB = n2(B1-B)

                y    = G(xm)
                nGxm = n2(y)
                Gxm  = y

                #print('err after div cleaning'+str(nGxm))
                
                fnx,fny,fmx,fmy,farf,elecf,divf = PDE.MHDSplity(y)
        
                momerr  = PDE.TVhL2Norm(fnx,fny,fmx,fmy)
                Farerr  = PDE.EhL2Norm(farf)
                Elecerr = PDE.VhL2Norm(elecf)
                #print('momerr='+str(momerr))
                #print('Farerr='+str(Farerr))
                #print('Elecerr='+str(Elecerr))
            i += 1
        if i == maxiter:
            print('Surpassed max number of iter without arriving at sol')
        else:
            print('Successfully completed Newton iterations')
        return xm  
            