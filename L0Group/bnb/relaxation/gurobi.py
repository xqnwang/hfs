import numpy as np
import math

CONIC = True

def l0gurobi_activeset(x, y, initial_activeset, group_indices, inv_W, kron_tSI, l0, l2, m, lb, ub, relaxed=True):
    converged = False
    fixed_to_zero = np.where(ub == 0)[0]
    # indices of the active groups.
    activeset = initial_activeset
    while (not converged):
        group_indices_restricted = [group_indices[index] for index in activeset]
        group_indices_restricted_reset_indices = []
        start_index = 0
        for i in range(len(group_indices_restricted)):
            group_indices_restricted_reset_indices.append(list(range(start_index, start_index+len(group_indices_restricted[i]))))
            start_index += len(group_indices_restricted[i])
        active_coordinate_indices = []
        for group_index in activeset:
            active_coordinate_indices += group_indices[group_index]
        beta_restricted, z_restricted, obj, _ = \
        l0gurobi(x[:, active_coordinate_indices], y, group_indices_restricted_reset_indices, inv_W, kron_tSI[:, active_coordinate_indices], l0, l2, m, lb[activeset], ub[activeset])
        # Check the KKT conditions.
        r = y - np.dot(x[:, active_coordinate_indices], beta_restricted)
        rt_invW_x = r @ inv_W @ x
        if l2 != 0 and np.sqrt(l0/l2) <= m:
            group_norms = np.array([np.linalg.norm(rt_invW_x[group_indices[index]]) for index in range(len(group_indices))])
            violations = set(np.where(group_norms > (2*np.sqrt(l0*l2) ))[0])
        else:
            group_norms = np.array([np.linalg.norm(rt_invW_x[group_indices[index]]) for index in range(len(group_indices))])
            violations = set(np.where(group_norms > (l0/m + l2*m ))[0])
        no_check_indices = set(activeset).union(set(fixed_to_zero))
        violations = violations.difference(no_check_indices)
        # print("Number of violations: ", len(violations))
        if len(violations) == 0:
            converged = True
            # print("Objective: ", obj)
        else:
            violations_list = np.array(sorted(violations))
            if len(violations_list) > 10:
                top_violations = violations_list[np.argpartition(group_norms[violations_list], -10)[-10:]]
            else:
                top_violations = violations_list
            # activeset += list(violations)
            activeset += list(top_violations)

    beta = np.zeros(x.shape[1])
    beta[active_coordinate_indices] = beta_restricted
    z = np.zeros(len(group_indices))
    z[activeset] = z_restricted
    return beta, z, obj


def l0gurobi(x, y, group_indices, inv_W, kron_tSI, l0, l2, m, lb, ub, relaxed=True):
    try:
        from gurobipy import Model, GRB, QuadExpr, LinExpr, quicksum
    except ModuleNotFoundError:
        raise Exception('Gurobi is not installed')
    model = Model()  # the optimization model
    n = x.shape[0]  # number of samples
    p = x.shape[1]  # number of features
    # n = S.shape[0]
    nb = int(math.sqrt(kron_tSI.shape[0]))  # number of bottom-level series
    group_num = len(group_indices)
    
    beta = model.addMVar(shape=(p, ), vtype=GRB.CONTINUOUS,
                         name=['B' + str(feature_index) for feature_index in range(p)],
                         ub=np.repeat(m, p), lb=np.repeat(-m, p))
    
    if relaxed:
        z = model.addMVar(shape=(group_num, ), vtype=GRB.CONTINUOUS,
                          name=['z' + str(group_index) for group_index in range(group_num)],
                          ub=ub, lb=lb)
    else:
        z = model.addMVar(shape=(group_num, ), vtype=GRB.BINARY,
                          name=['z' + str(group_index) for group_index in range(group_num)])
    
    r = model.addMVar(shape=(n, ), vtype=GRB.CONTINUOUS,
                     name=['r' + str(sample_index) for sample_index in range(n)],
                     ub=GRB.INFINITY, lb=-GRB.INFINITY)
    
    if l2 > 0:
        if CONIC:
            s = model.addVar(shape=(group_num, ), vtype=GRB.CONTINUOUS,
                             name=['s' + str(group_index) for group_index in range(group_num)],
                             lb=0)
    
    model.update()
    

    """ OBJECTIVE """
    
    if l2 == 0:
        model.setObjective(0.5*r.T@inv_W@r + l0*quicksum(z), GRB.MINIMIZE)
    else:
        if not CONIC:
            model.setObjective(0.5*r.T@inv_W@r + l0*quicksum(z) + l2*beta.T@beta, GRB.MINIMIZE)
        else:
            model.setObjective(0.5*r.T@inv_W@r + l0*quicksum(z) + l2*s.T@s, GRB.MINIMIZE)
    
    
    """ CONSTRAINTS """
    
    model.addConstr(r == y - x@beta)
    
    for group_index in range(group_num):
        l2_sq = [beta[feature_index]*beta[feature_index] for feature_index in group_indices[group_index]]
        model.addConstr(quicksum(l2_sq) <= m * m * z[group_index]*z[group_index])
    
    if l2 > 0:
        if CONIC:
            for group_index in range(group_num):
                l2_sq = [beta[feature_index]*beta[feature_index] for feature_index in group_indices[group_index]]
                model.addConstr(s[group_index] * z[group_index] >= quicksum(l2_sq))
    
    I = np.identity(nb)
    vI = I.reshape((nb*nb,))
    model.addConstr(vI == kron_tSI@beta)  # kronecker(t(S), I) beta = vec(I)
    
    model.update()
    model.setParam('OutputFlag', False)
    # model.setParam('BarConvTol', 1e-16)
    # model.setParam('BarIterLimit', 100000)
    # model.setParam('BarQCPConvTol', 1e-16)
    model.optimize()
    
    output_beta = np.zeros(beta.shape[0])
    output_z = np.zeros(z.shape[0])

    for i in range(beta.shape[0]):
        output_beta[i] = beta[i].x
    for group_index in range(group_num):
        output_z[group_index] = z[group_index].x

    return output_beta, output_z, model.ObjVal, None
