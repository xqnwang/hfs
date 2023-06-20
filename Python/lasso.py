import numpy as np
from numpy import linalg as LA
import gurobipy as gp
from gurobipy import GRB, quicksum
# Gurobi Optimizer version 10.0.1 build v10.0.1rc0
# import math

def socp(y, S, W, l1 = 0, m = None, M = None, weight = True, unbiased = True, TimeLimit = 0, Quick = 1, LogToConsole = 0, OutputFlag = 0):
    """
    Solve the OP problem: min_{G} 0.5 * (y - SGy)' W^{-1} (y - SGy) + l1 * sum_{j}(||G_{.j}||_2)
                          s.t. GS = I
                          

    Parameters
    ----------
    y : np.array
        1-d numpy array of base forecasts with size n.
    S : np.array
        n x nb numpy array describing the hierarchy structure.
    W : np.array
        n x n numpy array. The covariance matrix of the base forecast errors.
    l1 : float, optional
        lagrange multiplier.
    weight: int
        Enables weighted group-lasso when it is True.
    unbiased: int
        Includes unbiasedness constraint
    TimeLimit: float, optional
        set a timeout for gurobi.
    LogToConsole: int, optional
        Enables or disables console logging. Use OutputFlag to shut off all logging.
    OutputFlag: int, optional
        Enables or disables solver output. Use LogFile and LogToConsole for finer-grain control. Setting OutputFlag to 0 is equivalent to setting LogFile to "" and LogToConsole to 0.
        
        
    Returns
    -------
    1-d numpy array of diagonal elements of A.

    """
            
    n = S.shape[0]
    nb = S.shape[1]
    p = nb * n
    
    y = y.reshape((n,)) # reshape imported R object from (n, 1) to (n,)
    I = np.identity(nb)
    inv_W = np.linalg.inv(W)
    
    """ MinT solution """
    R = S.T @ inv_W
    G_mint = np.linalg.inv(R @ S) @ R
    
    """ POSSIBLE MAX LAMBDA
    obj_init = 0.5 * (y - S@G_mint@y).T @ inv_W @ (y - S@G_mint@y)
    ndigits = math.floor(math.log10(abs(obj_init))) + 2
    lambda_max = 10**ndigits
    """
    
    if weight:
        w = 1/LA.norm(G_mint, axis=0)
    else:
        w = np.repeat(1, n)
    
    if m is None:
        m = np.amax(abs(G_mint)) + 1
    if M is None:
        M = nb
    emax = np.amax(abs(y))
    
    """ SUPPRESS ALL OUTPUT """
    env = gp.Env(empty=True)
    env.setParam("OutputFlag",OutputFlag)
    env.start()
    
    """ SOCP MODEL """
    model = gp.Model('SOCP', env=env) # the optimization model
    
    """ PARAMETERS """
    # G matrix
    G = model.addMVar(shape=(p, ), vtype=GRB.CONTINUOUS,
                      ub=np.repeat(m, p), lb=np.repeat(-m, p))
    # Error
    E = model.addMVar(shape=(n, ), vtype=GRB.CONTINUOUS,
                      ub=np.repeat(emax, n), lb=np.repeat(-emax, n))
    # Auxiliary variables for l2 norm
    AUX = model.addMVar(shape=(n, ), vtype=GRB.CONTINUOUS,
                        ub=GRB.INFINITY, lb=np.repeat(0, n))
    model.update()

    """ OBJECTIVE """
    model.setObjective(0.5 * E.T @ inv_W @ E + l1 * quicksum(w * AUX), GRB.MINIMIZE)

    """ CONSTRAINTS """
    model.addConstr(y == E + np.kron(y.T, S) @ G)
    # SOC constraints
    for j in range(n):
        model.addGenConstrNorm(AUX[j], G[(j*nb):((j+1)*nb)], 2)
    # Unbiasedness constraint
    if unbiased:
        model.addConstr(I.reshape(-1) == np.kron(S.T, I) @ G)
    model.update()
    
    """ OUTPUT THE MODEL TO A FILE """ 
    # model.write("myfile.lp")
    
    """ OPTIMIZE """
    model.Params.OutputFlag = OutputFlag
    model.Params.LogToConsole = LogToConsole
    if Quick:
        model.Params.NumericFocus = 1
        model.Params.OptimalityTol = 1e-4
        model.Params.FeasibilityTol = 1e-4
        model.Params.BarConvTol = 1e-4
        model.Params.BarQCPConvTol = 1e-4
    if TimeLimit > 0:
        model.params.TimeLimit = TimeLimit
    model.optimize()
    # model.Params.Threads = 1
    
    g = G.X
    G = g.reshape(n, nb).T
    Z = 1 - (~(abs(G) > 1e-8).any(axis=0))*1
    obj = model.objval
    
    return G, Z, obj