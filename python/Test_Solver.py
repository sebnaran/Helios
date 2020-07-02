import numpy as np
from Solver import InexactNewtonTimeInt
from numpy.linalg import norm as n2
from scipy.sparse.linalg import LinearOperator
from scipy.sparse.linalg import gmres
from scipy.sparse import csc_matrix
class TestPDE(object):
    #This class is made for the sole purpose of testing the solver.
    #It represents the system
    # y^2-3x-2y = 0
    # x^2-3y-2x = 0
    #The exact solutions are x=y=5 or x=y=0
    def __init__(self):
        self.x = np.array([0,0])
    def G(self,x):
        y1 = x[1]**2-3*x[0]-2*x[1]
        y2 = x[0]**2-3*x[1]-2*x[0]
        return np.array([y1,y2])
    
def test_NewtonSolver():
    PDE    = TestPDE()
    Solver = InexactNewtonTimeInt(10,PDE)
    guess  = np.array([4.5,5.5])
    sol    = Solver.Newtoniter(PDE.G,guess,2,1E-5,50)
    print(sol)
    assert np.allclose(sol,np.array([5,5]))

def test_UsingLinOpWithGMRES():
    def G(x):
        return np.array([2*x[0],6*x[1]])
    
    DG = LinearOperator( (2,2), matvec = G)
    delxm, exitcode = gmres(DG,np.array([1,1]),atol=1E-5)

    assert ( np.allclose(delxm,np.array([1/2,1/6])) and exitcode <1E-5)

