import numpy as np
from numpy import linalg as LA
import gurobipy as gp
from gurobipy import GRB, quicksum
# Gurobi Optimizer version 10.0.1 build v10.0.1rc0

def miqp(y, S, W, l0 = 0, l2 = 0, m = None, M = None, MIPGap = None, TimeLimit = 600, LogToConsole = 0, OutputFlag = 0, WarmStart = 1, MIPFocus = 0, Cuts = -1):
    """
    Solve the OP problem: min_{G} 0.5 * (y - SGy)' W^{-1} (y - SGy) + l0 * Z(G) + l2 * ||vec(G)||^2
                          s.t. GS = I
    where Z(G) counts the number of nonzero columns of G


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
    l2 : float, optional
        lagrange multiplier.   
    m : float, optional
        bound of G matrix elements.
    M : float, optional
        bound of sum of absolute values of a column of G.
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
    1-d numpy array of diagonal elements of A.

    """
    
    """ Change parameters in a callback
    def my_callback(model, where):
        if where == GRB.Callback.MIP:
            run_time = model.cbGet(GRB.Callback.RUNTIME)
            # objbst = model.cbGet(GRB.Callback.MIP_OBJBST)
            # objbnd = model.cbGet(GRB.Callback.MIP_OBJBND)
            # gap = abs((objbst - objbnd) / objbst)
            
            if run_time > TimeChange:
                model._changeParam = True
                model.terminate()
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
    obj_guess = 0.5 * (y - S@G_mint@y).T @ inv_W @ (y - S@G_mint@y) + l0 * nb + l2 * LA.norm(G_mint, 2)**2
    
    if m is None:
        m = np.amax(G_mint) + 1
    
    if M is None:
        M = nb
    
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
    
    """ MIP MODEL """
    model = gp.Model('MIP', env=env) # the optimization model
    
    """ PARAMETERS """
    # G matrix
    G = model.addMVar(shape=(p, ), vtype=GRB.CONTINUOUS,
                  ub=np.repeat(m, p), lb=np.repeat(-m, p))
    # Positive G matrix
    PG = model.addMVar(shape=(p, ), vtype=GRB.CONTINUOUS,
                  ub=np.repeat(m, p), lb=np.repeat(0, p))
    # Error
    E = model.addMVar(shape=(n, ), vtype=GRB.CONTINUOUS,
                  ub=np.repeat(emax, n), lb=np.repeat(-emax, n))
                  # ub=[i * 3 for i in np.absolute(y)], lb=[-i * 3 for i in np.absolute(y)]
    # Binary
    Z = model.addMVar(shape=(n, ), vtype=GRB.BINARY)
    model.update()

    """ OBJECTIVE """
    model.setObjective(0.5 * E.T @ inv_W @ E + l0 * quicksum(Z) + l2 * PG.T @ PG, GRB.MINIMIZE)

    """ CONSTRAINTS """
    model.addConstr(y == E + np.kron(y.T, S) @ G)
    model.addConstr(I.reshape(-1) == np.kron(S.T, I) @ G)
    model.addConstrs(quicksum(PG[(j*nb):((j+1)*nb)]) <= M*Z[j] for j in range(n))
    model.addConstr(PG >= G)
    model.addConstr(PG >= -G)
    model.update()
    
    """ WARM START """
    if WarmStart:
        model.NumStart = 3
        model.update()
        
        # 1. Bottom-up
        model.params.StartNumber = 0 # set StartNumber
        rS = np.sum(S, axis=1)
        z_bu = np.where(rS == 1, 1, 0)
        for j in range(n):
            Z[j].Start = z_bu[j]
        
        # 2. All
        model.params.StartNumber = 1
        for j in range(n):
            Z[j].Start = 1
        
        # 3. Optimal
        """ QP MODEL """
        qp = gp.Model('MIP', env=env) # the optimization model
        
        """ PARAMETERS """
        # G matrix
        G_qp = qp.addMVar(shape=(p, ), vtype=GRB.CONTINUOUS,
                         ub=np.repeat(m, p), lb=np.repeat(-m, p))
        # Positive G matrix
        PG_qp = qp.addMVar(shape=(p, ), vtype=GRB.CONTINUOUS,
                          ub=np.repeat(m, p), lb=np.repeat(0, p))
        # Error
        E_qp = qp.addMVar(shape=(n, ), vtype=GRB.CONTINUOUS,
                         ub=np.repeat(emax, n), lb=np.repeat(-emax, n))
        # Relaxed binary
        Z_qp = qp.addMVar(shape=(n, ), vtype=GRB.CONTINUOUS,
                         ub=np.repeat(1, n), lb=np.repeat(0, n))
        qp.update()
        
        """ OBJECTIVE """
        qp.setObjective(0.5 * E_qp.T @ inv_W @ E_qp + l0 * quicksum(Z_qp) + l2 * PG_qp.T @ PG_qp, GRB.MINIMIZE)
        
        """ CONSTRAINTS """
        qp.addConstr(y == E_qp + np.kron(y.T, S) @ G_qp)
        qp.addConstr(I.reshape(-1) == np.kron(S.T, I) @ G_qp)
        qp.addConstrs(quicksum(PG_qp[(j*nb):((j+1)*nb)]) <= M*Z_qp[j] for j in range(n))
        qp.addConstr(PG_qp >= G_qp)
        qp.addConstr(PG_qp >= -G_qp)
        qp.update()
        qp.optimize()
        z_qp = np.where(Z_qp.X >= 0.001, 1, 0)
        
        model.params.StartNumber = 2
        for j in range(n):
            Z[j].Start = z_qp[j]
    
    """ OPTIMIZE """
    model.Params.MIPGap = MIPGap
    model.Params.OutputFlag = OutputFlag
    model.Params.LogToConsole = LogToConsole
    model.Params.MIPFocus = MIPFocus
    model.Params.Cuts = Cuts
    model.params.TimeLimit = TimeLimit
    model.optimize()
    # model.Params.Threads = 1
    """ Change parameters in a callback
    NewGap = 0.01
    model._changeParam = False
    model.optimize(my_callback) 
    if model._changeParam:
        model.params.MIPGap = NewGap
        model.params.TimeLimit = TimeChange
        model.optimize()
    """
    
    g = G.X
    Z = Z.X
    G = g.reshape(n, nb).T
    obj = model.objval
    gap = model.MIPGap
    if gap <= MIPGap or obj < obj_guess or abs(obj - obj_guess)/abs(obj_guess) < 0.01:
        opt = 1
    else:
        opt = 0
    
    return G, Z, obj, gap, opt
