import numpy as np
from numpy import linalg as LA
import gurobipy as gp
from gurobipy import GRB, quicksum
# Gurobi Optimizer version 10.0.1 build v10.0.1rc0
# import math

def eglasso(Y, Y_hat, S, l1 = 0, m = None, M = None, weight = True, TimeLimit = 0, LogToConsole = 0, OutputFlag = 0):
    """
    Solve the SOCP problem: min_{G} 0.5 * 1/N * ||vec(Y) - np.kron(S, Y) vec(G')||_2^2 + l1 * sum_{j}(w_{j} ||G_{.j}||_2)
                          

    Parameters
    ----------
    Y : np.array
        N x n numpy array of true values.
    Y_hat: np.array
        N x n numpy array of base forecasts.
    S : np.array
        n x nb numpy array describing the hierarchy structure.
    l1 : float, optional
        lagrange multiplier.
    weight: int
        Enables weighted group-lasso when it is True.
    TimeLimit: float, optional
        set a timeout for gurobi.
    LogToConsole: int, optional
        Enables or disables console logging. Use OutputFlag to shut off all logging.
    OutputFlag: int, optional
        Enables or disables solver output. Use LogFile and LogToConsole for finer-grain control. Setting OutputFlag to 0 is equivalent to setting LogFile to "" and LogToConsole to 0.
        
        
    Returns
    -------
    G, Z, obj

    """
    
    N = Y.shape[0]
    n = S.shape[0]
    nb = S.shape[1]
    p = nb * n
    Nn = N * n
    
    y = Y.T.reshape(Nn) # vec(Y)
    
    """ OLS solution """
    W = np.identity(n)
    inv_W = np.linalg.inv(W)
    R = S.T @ inv_W
    G_ols = np.linalg.inv(R @ S) @ R
    
    """ Penalty factor """
    if weight:
        w = 1/LA.norm(G_ols, axis=0)
    else:
        w = np.repeat(1, n)
    
    # w = w/w.sum() # normalize weight vector to sum to 1 
    
    """ Bound """
    if m is None:
        m = np.amax(abs(G_ols)) + 1
    if M is None:
        M = nb
    emax = np.amax(abs(y))
    
    """ SUPPRESS ALL OUTPUT """
    env = gp.Env(empty=True)
    env.setParam("OutputFlag",OutputFlag)
    env.start()
    
    """ EGLASSO MODEL """
    model = gp.Model('EGLASSO', env=env) # the optimization model
    
    """ PARAMETERS """
    # G matrix: vec(G')
    G = model.addMVar(shape=(p, ), vtype=GRB.CONTINUOUS,
                      ub=np.repeat(m, p), lb=np.repeat(-m, p))
    # Error
    E = model.addMVar(shape=(Nn, ), vtype=GRB.CONTINUOUS,
                      ub=np.repeat(emax, Nn), lb=np.repeat(-emax, Nn))
    # Auxiliary variables for l2 norm
    AUX = model.addMVar(shape=(n, ), vtype=GRB.CONTINUOUS,
                        ub=GRB.INFINITY, lb=np.repeat(0, n))
    model.update()

    """ OBJECTIVE """
    model.setObjective(0.5 * 1/N * E.T @ E + l1 * quicksum(w * AUX), GRB.MINIMIZE)

    """ CONSTRAINTS """
    model.addConstr(y == E + np.kron(S, Y_hat) @ G)
    # SOC constraints
    for j in range(n):
        model.addGenConstrNorm(AUX[j], G[n * np.arange(nb) + j], 2)
    model.update()
    
    """ OUTPUT THE MODEL TO A FILE """ 
    # model.write("myfile.lp")
    
    """ OPTIMIZE """
    model.Params.OutputFlag = OutputFlag
    model.Params.LogToConsole = LogToConsole
    if n > 50:
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
    G = g.reshape(nb, n)
    Z = 1 - (~(abs(G) > 1e-5).any(axis=0))*1
    obj = model.objval
    
    return G, Z, obj
