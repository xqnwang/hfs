import numpy as np
import math
from scipy import optimize as sci_opt


def upper_bound_solve(x, y, W, S, l0, l2, m, support, z_support, group_indices):
    if len(support) != 0:
        # x_support = x[:, support]
        # if l2 > 0:
        #     x_ridge = np.sqrt(2 * l2) * np.identity(len(support))
        #     x_upper = np.concatenate((x_support, x_ridge), axis=0)
        #     y_upper = np.concatenate((y, np.zeros(len(support))), axis=0)
        # else:
        #     x_upper = x_support
        #     y_upper = y
        # res = sci_opt.lsq_linear(x_upper, y_upper, (-m, m))
        # upper_bound = res.cost + l0 * len(z_support)
        # upper_beta = res.x
        nb = S.shape[1]
        activeset = z_support
        group_indices_restricted = [group_indices[index] for index in activeset]
        group_indices_restricted_reset_indices = []
        start_index = 0
        for i in range(len(group_indices_restricted)):
            group_indices_restricted_reset_indices.append(list(range(start_index, start_index+len(group_indices_restricted[i]))))
            start_index += len(group_indices_restricted[i])
        active_coordinate_indices = []
        for group_index in activeset:
            active_coordinate_indices += group_indices[group_index]
        tS = S.transpose()
        I = np.identity(nb)
        kron_tSI = np.kron(tS, I)
        upper_bound, upper_beta = gurobi_constrained_ridge_regression(x[:, active_coordinate_indices], y, group_indices_restricted_reset_indices, W, kron_tSI[:, active_coordinate_indices], l0, l2, m)
    else:
        upper_bound = 0.5 * np.linalg.norm(y) ** 2
        upper_beta = []
    return upper_bound, upper_beta







def gurobi_constrained_ridge_regression(x, y, group_indices, W, kron_tSI, l0, l2, m):
    try:
        from gurobipy import Model, GRB, QuadExpr, MQuadExpr, LinExpr, quicksum
    except ModuleNotFoundError:
        raise Exception('Gurobi is not installed')
    model = Model()  # the optimization model
    n = x.shape[0]  # number of samples. n = S.shape[0]
    p = x.shape[1]  # number of features
    nb = math.sqrt(kron_tSI.shape[0])  # number of bottom-level series
    group_num = len(group_indices)
    
    beta = model.addMVar(shape=(p, ), vtype=GRB.CONTINUOUS,
                         name=['B' + str(feature_index) for feature_index in range(p)],
                         ub=np.repeat(m, p), lb=np.repeat(-m, p))
    
    z = model.addMVar(shape=(group_num, ), vtype=GRB.CONTINUOUS,
                      name=['z' + str(group_index) for group_index in range(group_num)],
                      ub=np.repeat(1, group_num), lb=np.repeat(1, group_num))
    
    r = model.addVar(shape=(n, ), vtype=GRB.CONTINUOUS,
                     name=['r' + str(sample_index) for sample_index in range(n)],
                     ub=GRB.INFINITY, lb=-GRB.INFINITY)
    
    model.update()
    

    """ OBJECTIVE """
    
    model.setObjective(0.5*r.T@W@r + l0*quicksum(z) + l2*beta.T@beta, GRB.MINIMIZE)
    

    """ CONSTRAINTS """
    
    model.addConstr(r == y - x@beta)
    
    for group_index in range(group_num):
        l2_sq = [beta[feature_index]*beta[feature_index] for feature_index in group_indices[group_index]]
        model.addConstr(quicksum(l2_sq) <= m * m * z[group_index]*z[group_index])
    
    I = np.identity(nb)
    vI = I.reshape((nb*nb,))
    model.addConstr(vI == kron_tSI@beta)  # kronecker(t(S), I) beta = vec(I)
    
    model.update()
    model.setParam('OutputFlag', False)
    # model.setParam('BarConvTol', 1e-16)
    # model.setParam('BarIterLimit', 100000)
    # model.setParam('BarQCPConvTol', 1e-16)
    model.optimize()

    output_beta = np.zeros(len(beta))
    output_z = np.zeros(len(z))

    for i in range(len(beta)):
        output_beta[i] = beta[i].x
    for group_index in range(group_num):
        output_z[group_index] = z[group_index].x

    return model.ObjVal, output_beta
