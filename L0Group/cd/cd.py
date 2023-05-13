import numpy as np
import scipy as sc
from math import fabs, sqrt, ceil, floor, exp, inf
import itertools
import numpy as np
from scipy.linalg import block_diag
from scipy import sparse
import copy
from collections import defaultdict


def ObjectiveCDL0Group(eps,B,Lambda0,Lambda2,Indices):
    L0Groupnorm = 0
    for i in Indices:
        if np.count_nonzero(B[i]) > 0:
            L0Groupnorm += 1
    return 0.5*eps + Lambda0*L0Groupnorm + Lambda2*np.linalg.norm(B,2)**2

class cd:
    NumCalls = 0 # Class variable that counts the number of calls to CD in one session
    def __init__(self, x, y, group_indices, Params):
        cd.NumCalls += 1
        self.Params = Params
        self.X = x
        self.y = y
        self.group_indices = group_indices
        self.NumIters = Params['NumIters']
        self.Tolerance = Params['Tolerance']
        self.Logging = Params['Logging']
        self.FeaturesOrder = Params['FeaturesOrder']

        n, p = self.X.shape
        if self.FeaturesOrder == "Cyclic":
            self.order = list(range(p))
        elif self.FeaturesOrder == "Corr":
            self.order = Params["CorrOrder"]
        self.iterations = 0 # Num of iterations until convergence (used for debugging)

        self.Lg = []
        for g in group_indices:
            Xg = self.X[:,g]
            l = np.max(np.linalg.eig(np.dot(Xg.T,Xg))[0])
            self.Lg.append(l)
        self.Lg = np.array(self.Lg)

    def fit(self, lambda_0, lambda_2, warm_start=None):
        n, p = self.X.shape
        Lambda0 = lambda_0
        Lambda2 = lambda_2
        if warm_start is None:
            self.B = np.zeros(p)
        else:
            self.B = np.copy(warm_start)
        Lg = self.Lg + 2*lambda_2 # list of Lgs
        Indices = self.group_indices # [[indexset1], [indexset2], ...]
        Thresholdg = 2*Lambda0/Lg
        if self.FeaturesOrder == "Cyclic":
            self.order = list(range(len(Lg)))
        elif self.FeaturesOrder == "Corr":
            self.order = self.Params["CorrOrder"]
        self.r = self.y - self.X.dot(self.B)
        self.objective = ObjectiveCDL0Group(self.r.T.dot(self.r), self.B, lambda_0, lambda_2, self.group_indices)

        def IterationEnd(t):
            objective_old = self.objective
            se = np.dot(self.r,self.r)
            self.objective = ObjectiveCDL0Group(se, self.B, lambda_0, lambda_2, self.group_indices)
            if(self.Logging):
                print('Iteration: '+str(t+1)+'. Objective = '+str(self.objective)+'. |Supp(B)| = '+str(np.count_nonzero(self.B)))
                print('Active Set: ',np.nonzero(self.B))
            if fabs(self.objective-objective_old)/objective_old < self.Tolerance:
                self.iterations = t+1
                if(self.Logging):
                    print('Converged in '+str(t+1)+' iterations.')
                return True
            return False

        for t in range(self.NumIters):
            for g in self.order:
                Bg = np.copy(self.B[Indices[g]])
                x = Bg - (1/Lg[g]) * ( - np.dot(self.X[:,Indices[g]].T,self.r) + 2*Lambda2*Bg )
                # x = np.linalg.lstsq(self.X[:,Indices[g]], self.r + np.dot(self.X[:,Indices[g]],Bg) )[0]
                if np.dot(x,x) >= Thresholdg[g]: # ||x||^2
                    self.B[Indices[g]] = x
                    self.r += np.dot(self.X[:,Indices[g]], Bg - x)
                else:
                    self.B[Indices[g]] = 0
                    self.r += np.dot(self.X[:,Indices[g]], Bg)
            if (IterationEnd(t)): break


        return (self.B,self.objective)



class cd_swaps:
    def __init__(self, x, y, group_indices, Params):
        self.Params = Params
        self.x = x
        self.y = y
        self.group_indices = group_indices
        self.Xpinv = []
        for g in self.group_indices:
            Xgp = np.linalg.pinv(self.x[:,g])
            self.Xpinv.append(Xgp)
        self.cd_object = cd(self.x, self.y, self.group_indices, self.Params)

    def fit(self, lambda_0, lambda_2, max_swaps=1, warm_start=None):
        if lambda_2 > 0:
            # In this case, update the pseudo inverse to account for lambda_2 > 0.
            for i, g in enumerate(self.group_indices):
                Xgp = np.dot( np.linalg.inv(np.dot(self.x[:,g].T, self.x[:,g]) + 2*lambda_2*np.eye(self.x[:,g].shape[1])) , self.x[:,g].T)
                self.Xpinv[i] = Xgp
        B, objective = self.cd_object.fit(lambda_0, lambda_2, warm_start)
        Indices = self.group_indices
        for t in range(max_swaps):
            r = self.y - self.x.dot(B)
            Support = []
            for i in Indices:
                if np.count_nonzero(B[i]) > 0:
                    Support.append(i)
            found_better = False
            for g in Support:
                XgBg = np.dot(self.x[:,g],B[g])
                rnog = r + XgBg
                keep_g_objective = 0.5*np.dot(r, r) + lambda_2*np.dot(B[g], B[g])
                maxgnew = 0
                Bmaxgnew = 0
                min_objective = 1e20
                for i in range(len(Indices)):
                    gnew = Indices[i]
                    if gnew not in Support:
                        Bgnew = np.dot(self.Xpinv[i],rnog)
                        XgnewBgnew = np.dot(self.x[:,gnew],Bgnew)
                        current_residuals = rnog - XgnewBgnew
                        current_ssr = np.dot(current_residuals, current_residuals)
                        current_objective = 0.5*current_ssr + lambda_2*np.dot(Bgnew, Bgnew)
                        if current_objective < min_objective:
                            maxgnew = gnew
                            Bmaxgnew = Bgnew
                            min_objective = current_objective

                if min_objective < keep_g_objective:
                    B[maxgnew] = Bmaxgnew
                    B_g_old = B[g]
                    B[g] = np.zeros(len(g))
                    B_temp, objective_temp = self.cd_object.fit(lambda_0, lambda_2, B)
                    # Additional check in case of numerical issues.
                    if objective_temp < objective:
                        B = B_temp
                        objective = objective_temp
                        found_better = True
                        break
                    else:
                        # Revert back due to a numerical precision error.
                        print("A swap led to numerical issues. Reverting...")
                        B[maxgnew] = np.zeros(len(Bmaxgnew))
                        B[g] = B_g_old

            if not found_better:
                # At this point no swap can help (since there was no break)
                return (B,objective)
        return (B,objective)
