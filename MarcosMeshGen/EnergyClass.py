import numpy as np
import matplotlib.pyplot as plt

class Energy(object):

    def __init__(self,theta,dt,N,L1,L2,R1,R2):
        self.Np1    = N  #This is should be equal to N+1 where N is the last step taken in the time integration
        self.itheta = theta
        self.idt    = dt
        self.iL1    = L1
        self.iL2    = L2 
        self.iR1    = R1 
        self.iR2    = R2
    
    def SearchQ(self):
        
        Qs = np.arange(0,1/self.itheta,0.01)
        F  = Qs*0
        step = 0
        for Q in Qs:
            beta  = (1-Q*self.itheta)/(1+Q*(1-self.itheta))
            gamma = 1/(1-Q*self.itheta)
            LHS   = 0
            RHS   = 0
            for n in range(self.Np1-1):
                Ei = self.iL2[n] 
                Eb = self.iR2[n]

                LHS = beta**(n)*Ei + LHS
                RHS = beta**(n)*Eb + RHS

            LHS = 0.5*LHS*gamma*self.idt
            RHS = RHS*gamma*self.idt
        
            F[step] = self.iR1+RHS-LHS-(self.iL1)*beta**(self.Np1)
            F[step] = F[step]/self.iR1
            step    = step+1

        return F