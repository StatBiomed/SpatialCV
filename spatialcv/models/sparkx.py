# Reimplementation of SPARK-X model

import numpy as np
import pandas as pd
from chi2comb import chi2comb_cdf, ChiSquared
import math

def sparkx(loc_mat,count_mat,X_mat = None):
    """
    Main function of sparkx
    
    Parameters
    ----------
    loc: numpy.matrix
        The matrix for spatial location information  x axis: location  y axis: coordinate
    count: numpy.matrix
        The count matrix of genes    x axis: location  y axis: genes
    Xamt: numpy.matrix
        A covariate matrix 
    
    Returns
    -------
    a dataframe with each column representing a transformation of coordinate and each row representing a gene
    """
    
    colsum = count_mat.sum(axis = 0)
    
    count = count_mat[:,np.where(colsum != 0)[1]]
    location = loc_mat[np.where(colsum != 0)[1],:]           # location  matrix S

    if(X_mat != None):
        X_mat = X_mat[np.where(colsum != 0)[1],:]                # delete some cells
        
    
    rowsum = count.sum(axis = 1)
    
    count = count[np.where(rowsum != 0)[0],:]    # delete some genes
    
    count = np.transpose(count)                      # count : n*p matrix    each row is the y vector
    
    n,p = count.shape   
    d = location.shape[1]                             
    
    loc_colmean = location.mean(axis = 0)
    S = location - np.tile(loc_colmean,(n,1))              #centering location
    
    dfdict = {}
    
    dfdict['projection'] = sparkp(S,count,X_mat)
    
    for i in range(5):
        adjust_S = np.hstack((trans_coor(S[:,0:1],i,func = "gaussian"),trans_coor(S[:,1:2],i,func = "gaussian")))   
        
        loc_colmean = adjust_S.mean(axis = 0)
        adjust_S = adjust_S - np.tile(loc_colmean,(n,1))              #centering location
        
        dfdict['gauss'+str(i+1)] = sparkp(adjust_S,count,X_mat)
        
    for i in range(5):
        adjust_S = np.hstack((trans_coor(S[:,0:1],i,func = "cosine"),trans_coor(S[:,1:2],i,func = "cosine")))
        
        loc_colmean = adjust_S.mean(axis = 0)
        adjust_S = adjust_S - np.tile(loc_colmean,(n,1))              #centering location
        
        dfdict['cosine'+str(i+1)] = sparkp(adjust_S,count,X_mat)
    
    df = pd.DataFrame(dfdict)
    
    return df

def sparkp(loc,count,X_mat = None):                           # loc: n × d matrix of spatial coordinates     count:n × p count matrix
    
    """
    sparkx algorithm function
    
    Parameters
    ----------
    loc: numpy.matrix
        The matrix for spatial location information  x axis: location  y axis: coordinate
    count: numpy.matrix
        The count matrix of genes    x axis: location  y axis: genes
    Xamt: numpy.matrix
        A covariate matrix 
    
    Returns
    -------
    a list of p-value, P (teststat > observed value) under null hypothesis for each gene
    """
    
    n,p = count.shape   
    
    loc_colmean = loc.mean(axis = 0)
    S = loc - np.tile(loc_colmean,(n,1))              #centering location
    STSinv = np.linalg.inv(S.transpose().dot(S))
    
    YTYinv = 1/(count.power(2)).sum(axis = 0)             # 1*p matrix stores the (yty)-1 for each gene
    
    if(X_mat == None):
        
        SigmaE = np.linalg.eig(STSinv.dot(S.transpose()).dot(S))[0]            # eigenvalue of Sigma
        
        YTHS = count.transpose().dot(S)                    # each row is ytHS matrix (1*d)
        
        teststat = np.matrix(list(map(lambda x: YTHS[x,:].dot(STSinv).dot(YTHS[x,:].transpose())*YTYinv[0,x], range(p))))*n  
        #teststat is n^2 times the teststat in paper, following the mixture chi dist.
        
        Ybar = count.mean(axis = 0)
        EE = list(map(lambda x: 1-n*YTYinv[0,x]*Ybar[0,x]**2,range(p)))         #eigenvalue of E, each element for each gene
        
    else:                                                # X_mat  n*q matrix
        
        XTXinv = np.linalg.inv(X_mat.transpose().dot(X_mat))
        adjust = X_mat.dot(XTXinv).dot(X_mat.transpose())
        H = np.identity(n) - adjust
        HS = H.dot(S)
        
        SigmaE = np.linalg.eig(STSinv.dot(HS.transpose()).dot(HS))[0]
        
        HY = H.dot(Y)
        YTHY = HY.power(2).sum(axis = 1)
        
        YTHS = Y.transpose().dot(HS)
        
        teststat = np.matrix(list(map(lambda x: YTHS[x,:].dot(STSinv).dot(YTHS[x,:].transpose())*YTYinv[0,x], range(p))))*n
        
        EE = list(map(lambda x: YTHY[0,x]*YTYinv[0,x],range(p)))
        
    pval = list(map(lambda x:1-chi2comb_cdf(teststat[0,x],[ChiSquared(SigmaE[i] * EE[x], 0, 1) for i in range(2)],0,atol = 1e-3)[0],range(p))) 
    
    return pval


def trans_coor(coord,iter,func = "gaussian"):          # here location is centered    one list of location coord.
    """
    Transform the coordinate
    
    Parameters
    ----------
    coord: numpy.matrix
        One column of location coordinate
    iter: integer in (0,1,2,3,4)
        
    func: "gaussian" or "cosine"
    
    Returns
    -------
    A transformed numpy.matrix as the same size of coord
    """    
    
    q = np.quantile(np.absolute(coord),0.2*(iter+1))
    
    if(func == "gaussian"):
        result = np.exp(-np.power(coord,2)/(2*q**2))

  
    if(func == "cosine"):
        result = np.cos(2*math.pi*coord/q)
        
        
    return(result)
