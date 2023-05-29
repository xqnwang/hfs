import numpy as np
import gurobipy as gp
from gurobipy import GRB, quicksum
# Gurobi Optimizer version 10.0.1 build v10.0.1rc0

def miqp(y, S, W, l0 = 0, M = 10):
    """
    Solve the OP problem: min_{A, G} (y - SGy)' W^{-1} (y - SGy) + l0sum(A)  
                          s.t. G = (S'W^{-1}AS)^{-1}S'W^{-1}A

    Parameters
    ----------
    y : np.array
        1-d numpy array of base forecasts with sisze n.
    S : np.array
        n x nb numpy array describing the hierarchy structure.
    W : np.array
        n x n numpy array. The covariance matrix of the base forecast errors.
    l0 : float, optional
        lagrange multiplier.
    M : float, optional
        bound of G matrix elements.

    Returns
    -------
    1-d numpy array of diagonal elements of A.

    """
    n = S.shape[0] # number of samples. n = S.shape[0]
    nb = S.shape[1]
    p = nb * n
    I = np.identity(nb)
    inv_W = np.linalg.inv(W)
    
    """ MODEL """
    model = gp.Model('MIP') # the optimization model
    
    """ PARAMETERS """
    ## G matrix
    G = model.addMVar(shape=(nb, n), vtype=GRB.CONTINUOUS,
                      ub=np.repeat(M, p).reshape((nb, n)), lb=np.repeat(-M, p).reshape((nb, n)))
    ## A matrix
    A = model.addMVar(shape=(n, n), vtype=GRB.BINARY,
                      ub=np.identity(n), lb=np.diag(np.repeat(0, n)))
    ## Error vetor used to enable Quadratic objective function
    E = model.addMVar(shape=(n, ), vtype=GRB.CONTINUOUS,
                      ub=GRB.INFINITY, lb=-GRB.INFINITY)
    ## Matrix introduced to enforce the inverse equality constraint (inverse cannot be directly inluded in objective function and constraints) 
    C = model.addMVar(shape=(nb, nb), vtype=GRB.CONTINUOUS,
                      ub=GRB.INFINITY, lb=-GRB.INFINITY)
    model.update()

    """ OBJECTIVE """
    model.setObjective(0.5 * E.T @ inv_W @ E + l0 * quicksum(A@np.repeat(1, n)), GRB.MINIMIZE)

    """ CONSTRAINTS """
    SW = S.T @ inv_W
    model.addConstr(E == y - S@G@y)
    model.addConstr(C @ SW @ A @ S == I) # invertible
    model.addConstr(G == C @ SW @ A)
    
    """ OPTIMIZE """
    model.setParam('OutputFlag', False)
    # model.setParam('BarConvTol', 1e-16)
    # model.setParam('BarIterLimit', 100000)
    # model.setParam('BarQCPConvTol', 1e-16)
    model.optimize()
    
    out = np.diag(A.x)
    return out
    
    