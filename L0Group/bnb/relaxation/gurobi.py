import numpy as np

CONIC = True

def l0gurobi_activeset(x, y, initial_activeset, group_indices, l0, l2, m, lb, ub, relaxed=True):
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
        l0gurobi(x[:, active_coordinate_indices], y, group_indices_restricted_reset_indices, l0, l2, m, lb[activeset], ub[activeset])
        # Check the KKT conditions.
        r = y - np.dot(x[:, active_coordinate_indices], beta_restricted)
        r_t_x = np.dot(r.T, x)
        if l2 != 0 and np.sqrt(l0/l2) <= m:
            group_norms = np.array([np.linalg.norm(r_t_x[group_indices[index]]) for index in range(len(group_indices))])
            violations = set(np.where(group_norms > (2*np.sqrt(l0*l2) ))[0])
        else:
            group_norms = np.array([np.linalg.norm(r_t_x[group_indices[index]]) for index in range(len(group_indices))])
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


def l0gurobi(x, y, group_indices, l0, l2, m, lb, ub, relaxed=True):
    try:
        from gurobipy import Model, GRB, QuadExpr, LinExpr, quicksum
    except ModuleNotFoundError:
        raise Exception('Gurobi is not installed')
    model = Model()  # the optimization model
    n = x.shape[0]  # number of samples
    p = x.shape[1]  # number of features
    group_num = len(group_indices)

    beta = {}  # features coefficients
    z = {}
    for feature_index in range(p):
        beta[feature_index] = model.addVar(vtype=GRB.CONTINUOUS,
                                           name='B' + str(feature_index),
                                           ub=m, lb=-m)
    for group_index in range(group_num):
        if relaxed:
            z[group_index] = model.addVar(vtype=GRB.CONTINUOUS,
                                            name='z' + str(feature_index),
                                            ub=ub[group_index],
                                            lb=lb[group_index])
        else:
            z[group_index] = model.addVar(vtype=GRB.BINARY,
                                            name='z' + str(feature_index))

    r = {}
    for sample_index in range(n):
        r[sample_index] = model.addVar(vtype=GRB.CONTINUOUS,
                                       name='r' + str(sample_index),

                                       ub=GRB.INFINITY, lb=-GRB.INFINITY)
    if l2 > 0:
        s = {}
        if CONIC:
            for group_index in range(group_num):
                s[group_index] = model.addVar(vtype=GRB.CONTINUOUS,
                                                   name='s' + str(group_index), lb=0)
    model.update()

    """ OBJECTIVE """

    obj = QuadExpr()

    for sample_index in range(n):
        obj.addTerms(0.5, r[sample_index], r[sample_index])

    if l2 > 0:
        if not CONIC:
            for feature_index in range(p):
                obj.addTerms(l2, beta[feature_index], beta[feature_index])
        else:
            for group_index in range(group_num):
                obj.addTerms(l2, s[group_index])

    for group_index in range(group_num):
        obj.addTerms(l0, z[group_index])

    model.setObjective(obj, GRB.MINIMIZE)

    """ CONSTRAINTS """

    for sample_index in range(n):
        expr = LinExpr()
        expr.addTerms(x[sample_index, :], [beta[key] for key in range(p)])
        model.addConstr(r[sample_index] == y[sample_index] - expr)

    # for group_index in range(group_num):
    #     for feature_index in group_indices[group_index]:
    #         model.addConstr(beta[feature_index] <= z[group_index] * m)
    #         model.addConstr(beta[feature_index] >= -z[group_index] * m)
    for group_index in range(group_num):
        l2_sq = []
        for feature_index in group_indices[group_index]:
            l2_sq.append(beta[feature_index]*beta[feature_index])
        model.addConstr(quicksum(l2_sq) <= m * m * z[group_index]*z[group_index])

    if l2 > 0:
        if CONIC:
            for group_index in range(group_num):
                l2_sq = []
                for feature_index in group_indices[group_index]:
                    l2_sq.append(beta[feature_index]*beta[feature_index])
                model.addConstr(s[group_index] * z[group_index] >= quicksum(l2_sq))

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

    return output_beta, output_z, model.ObjVal, None
