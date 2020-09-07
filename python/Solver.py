from scipy.sparse.linalg import LinearOperator
from scipy.sparse.linalg import gmres
from numpy.linalg import norm as n2
import numpy as np
import math

class InexactNewtonTimeInt(object):
    def __init__(self):
        #Series of constants for the error in the GMRES
        self.eps    = 1E-7
        self.etamax = 0.8
        self.alpha  = 1.5
        self.gamma  = 0.9
        self.epsr   = 1E-4

    def Newtoniter(self,G,x0,ndof,tol,maxiter):
        #This function will perform a Newton Iteration
        #To find the zeroes of the function G provided the initial guess x0
        #Within the tolerance tol in the 2-norm. ndof is the number of unknowns.
        xm     = x0
        #print(x0)
        #print('eval G')
        Gxm    = G(x0)
        etamm1 = self.etamax
        nGxmm1 = 1
        epsa   = math.sqrt(ndof)*(10**(-15))
        epst   = epsa+self.epsr*n2(G(x0))
        for i in range(maxiter):
            nGxm = n2(G(xm))
            #print('err='+str(nGxm))
            if abs(nGxm)<tol:
                return xm
            else:
                def fDGxm(delx):
                    return (G(xm+self.eps*delx)-Gxm)/(self.eps)
                DGxm = LinearOperator((ndof,ndof), matvec = fDGxm)
                etamA = self.gamma*(nGxm/nGxmm1)**(self.alpha)
                etamB = min([self.etamax,max([etamA,self.gamma*etamm1**self.alpha])])
                etam  = min([self.etamax,max([etamB,self.gamma*(epst/nGxm)])])
                delxm, exitcode = gmres(DGxm,-Gxm,atol=etam*nGxm)
                
                xm              = xm+delxm
                Gxm             = G(xm)
                etamm1 = etam
                nGxmm1 = nGxm
                if exitcode>1E-5:
                    print('Error Ocurred in the GMRES Iteration')
        print('Surpassed max number of iter without arriving at sol')