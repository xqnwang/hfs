import numpy as np
import scipy as sc
import math
from numpy.linalg import cholesky
from numpy.random import normal

def gen_synthetic(parameter, n, p, SuppSize, SNR, Gsize):
    np.random.seed(1)
    Indices = [[j for j in range(i,i+Gsize)] for i in range(0,p,Gsize)]
    Ind = np.linspace(0, len(Indices)-1, num=SuppSize).astype(int)
    B = np.zeros(p)
    for g in Ind:
        B[Indices[g]] = 1

    Groupnum = int(p/Gsize)
    Covariance = np.array([[parameter**np.abs(i-j) for j in range(Groupnum)] for i in range(Groupnum)])
    Xg, _, _, _ = GenGaussianDataFixed(n,Groupnum,Covariance,SNR,B[0:Groupnum], "E")
    X_training = np.empty(shape=(n,p))
    for i in range(Groupnum):
        epsilonmatrix = np.random.normal(scale=np.sqrt(0.111111),size=(n,Gsize))
        Xgi = np.reshape(Xg[:,i],newshape=[n,1])
        X_training[:,i*Gsize:(i+1)*Gsize] = np.add(epsilonmatrix, Xgi)

    mu = X_training.dot(B)
    var_XB = (np.std(mu,ddof=1))**2
    sd_epsilon = np.sqrt(var_XB/SNR)
    epsilon = normal(size=n,scale=sd_epsilon)
    y_training = mu + epsilon

    group_indices = [[Gsize*i + j for j in range(Gsize)] for i in range(int(p/Gsize))]

    return X_training, y_training, group_indices, B



def GenGaussianDataFixed(n,p,Covariance,SNR,B,D=""):
    if D == "I":
        X = normal(size=(n,p))
        var_XB = np.dot(B,B)
        mu = X.dot(B)

    elif D == "CLarge":
        rho = Covariance
        X = normal(size=(n,p)) + np.sqrt(rho/(1-rho))*normal(size=(n,1))
        mu = X.dot(B)
        var_XB = (np.std(mu,ddof=1))**2

    else:
        A = cholesky(Covariance)
        Z = normal(size=(n,p))
        X = Z.dot(A.T)
        var_XB = np.dot(np.dot(B,Covariance),B)
        mu = X.dot(B)

    # Generate epsilon
    sd_epsilon = np.sqrt(var_XB/SNR)
    epsilon = normal(size=n,scale=sd_epsilon)
    epsilondev = normal(size=n,scale=sd_epsilon)

    # Generate y
    y = mu + epsilon
    ydev = mu + epsilondev

    return X, y, ydev, B
