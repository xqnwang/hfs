import numpy as np
import pandas as pd
import math
import gurobipy as gp
from gurobipy import *
# Gurobi Optimizer version 10.0.1 build v10.0.1rc0

def miqp_AS(y, S, W, l0 = 0, m = None, MIPGap = None, TimeLimit = 600, LogToConsole = 0, OutputFlag = 0, WarmStart = 1, MIPFocus = 0, Cuts = -1):
    """
    Solve the OP problem: min_{A, G_bar, C} 0.5 * (y - S G_bar y)' W^{-1} (y - S G_bar y) + l0 * sum(A)
                          s.t. C(S'A'W^{-1}AS) = I ...invertible
                               G_bar = (S'A'W^{-1}AS)^{-1}S'A'W^{-1} ...G_bar
                               G_bar S = I ...unbiasedness

    Parameters
    ----------
    y : np.array
        1-d numpy array of base forecasts with size n.
    S : np.array
        n x nb numpy array describing the hierarchy structure.
    W : np.array
        n x n numpy array. The covariance matrix of the base forecast errors.
    l0 : float, optional
        lagrange multiplier. 
    m : float, optional
        bound of G matrix elements.
    MIPGap: float, optional
        the MIP solver will terminate (with an optimal result) when the gap between the lower and upper objective bound is less than MIPGap times the absolute value of the incumbent objective value
    TimeLimit: float, optional
        set a timeout for gurobi.
    LogToConsole: int, optional
        Enables or disables console logging. Use OutputFlag to shut off all logging.
    OutputFlag: int, optional
        Enables or disables solver output. Use LogFile and LogToConsole for finer-grain control. Setting OutputFlag to 0 is equivalent to setting LogFile to "" and LogToConsole to 0.
    WarmStart: int, optional
        bottom-up warm start. Default True.
    MIPFocus: int, optional
        If you are more interested in finding feasible solutions quickly, you can select MIPFocus=1. If you believe the solver is having no trouble finding good quality solutions, and wish to focus more attention on proving optimality, select MIPFocus=2. If the best objective bound is moving very slowly (or not at all), you may want to try MIPFocus=3 to focus on the bound.
    Cuts: int, optional
        Global cut aggressiveness setting. Use value 0 to shut off cuts, 1 for moderate cut generation, 2 for aggressive cut generation, and 3 for very aggressive cut generation. This parameter is overridden by the parameters that control individual cut types (e.g., CliqueCuts).
    
    
    Returns
    -------
    G, Z (1-d numpy array of diagonal elements of A), obj, gap, opt

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
    obj_guess = 0.5 * (y - S@G_mint@y).T @ inv_W @ (y - S@G_mint@y) + l0 * nb
    
    if m is None:
        m = np.amax(abs(G_mint)) + 1
        
    if MIPGap is None:
        if n <= 50:
            MIPGap = 0.0001
        else:
            MIPGap = 0.001
    
    emax = np.amax(abs(y))
    
    """ SUPPRESS ALL OUTPUT """
    env = gp.Env(empty=True)
    env.setParam("OutputFlag",OutputFlag)
    env.start()
    
    """ MODEL """
    model = gp.Model('MIP_AS', env=env)
    
    """ PARAMETERS """
    ## A matrix
    A = model.addMVar(shape=(n, n), vtype=GRB.BINARY,
                      ub=np.identity(n), lb=np.diag(np.repeat(0, n)))
    ## G matrix
    G = model.addMVar(shape=(nb, n), vtype=GRB.CONTINUOUS,
                      ub=np.repeat(m, p).reshape((nb, n)), lb=np.repeat(-m, p).reshape((nb, n)))
    ## Error vetor used to enable Quadratic objective function
    E = model.addMVar(shape=(n, ), vtype=GRB.CONTINUOUS,
                      ub=np.repeat(emax, n), lb=np.repeat(-emax, n))
    ## Matrix introduced to enforce the inverse equality constraint (inverse cannot be directly inluded in objective function and constraints) 
    C = model.addMVar(shape=(nb, nb), vtype=GRB.CONTINUOUS,
                      ub=GRB.INFINITY, lb=-GRB.INFINITY)
    model.update()

    """ OBJECTIVE """
    model.setObjective(0.5 * E.T @ inv_W @ E + l0 * quicksum(A@np.repeat(1, n)), GRB.MINIMIZE)

    """ CONSTRAINTS """
    AUX = S.T @ A.T @ inv_W
    model.addConstr(E == y - S@G@y)
    model.addConstr(G @ A @ S == I) # Avoid C(S'A'W^{-1}AS) = I as Gurobi doesn't support multiplication of three MVars
    model.addConstr(C @ AUX == G)
    model.addConstr(G @ S == I)
    model.update()
    
    """ WARM START """
    if WarmStart:
        model.NumStart = 2
        model.update()
        
        # 1. Bottom-up
        model.params.StartNumber = 0 # set StartNumber
        rS = np.sum(S, axis=1)
        z_bu = np.diag(np.where(rS == 1, 1, 0))
        A.Start = z_bu
        
        # 2. All
        model.params.StartNumber = 1
        z_all = np.identity(n)
        A.start = z_all
    
    """ OPTIMIZE """
    model.Params.MIPGap = MIPGap
    model.Params.OutputFlag = OutputFlag
    model.Params.LogToConsole = LogToConsole
    model.Params.MIPFocus = MIPFocus
    model.Params.Cuts = Cuts
    if TimeLimit > 0:
        model.params.TimeLimit = TimeLimit
    model.optimize()
    
    G = G.X
    Z = np.diag(A.X)
    obj = model.objval
    gap = model.MIPGap
    if gap <= MIPGap or obj < obj_guess or abs(obj - obj_guess)/abs(obj_guess) < 0.01:
        opt = 1
    else:
        opt = 0
    
    return G, Z, obj, gap, opt
    
