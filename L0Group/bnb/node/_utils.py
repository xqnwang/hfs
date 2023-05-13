import numpy as np
from scipy import optimize as sci_opt


def upper_bound_solve(x, y, l0, l2, m, support, z_support, group_indices):
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
        upper_bound, upper_beta = gurobi_constrained_ridge_regression(x[:, active_coordinate_indices], y, group_indices_restricted_reset_indices, l0, l2, m)
    else:
        upper_bound = 0.5 * np.linalg.norm(y) ** 2
        upper_beta = []
    return upper_bound, upper_beta







def gurobi_constrained_ridge_regression(x, y, group_indices, l0, l2, m):
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
        z[group_index] = model.addVar(vtype=GRB.CONTINUOUS,
                                        name='z' + str(feature_index),
                                        ub=1,
                                        lb=1)

    r = {}
    for sample_index in range(n):
        r[sample_index] = model.addVar(vtype=GRB.CONTINUOUS,
                                       name='r' + str(sample_index),

                                       ub=GRB.INFINITY, lb=-GRB.INFINITY)

    model.update()

    """ OBJECTIVE """

    obj = QuadExpr()

    for sample_index in range(n):
        obj.addTerms(0.5, r[sample_index], r[sample_index])


    for feature_index in range(p):
        obj.addTerms(l2, beta[feature_index], beta[feature_index])

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
