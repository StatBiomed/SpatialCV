# Reimplementation of SPARK-X model

def sparkx(loc,count,X_mat = None):                           # loc: n × d matrix of spatial coordinates     count:n × p count matrix
    
    loc_colmean = location.mean(axis = 0)
    S = location - np.tile(loc_colmean,(n,1))              #centering location
    STSinv = np.linalg.inv(S.transpose().dot(S))
    
    YTYinv = 1/(count.power(2)).sum(axis = 0)             # 1*p matrix stores the (yty)-1 for each gene
    
    if(X_mat == None):
        
        SigmaE = np.linalg.eig(STSinv.dot(S.transpose()).dot(S))[0]            # eigenvalue of Sigma
        
        YTHS = count.transpose().dot(S)                    # each row is ytHS matrix (1*d)
        
        teststat = np.matrix(list(map(lambda x: YTHS[x,:].dot(STSinv).dot(YTHS[x,:].transpose())*YTYinv[0,x], range(p))))*n  
        #teststat is n^2 times the teststat in paper, following the mixture chi dist.
        
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
        
    pval = list(map(lambda x:1-chi2comb_cdf(teststat[0,x],[ChiSquared(SigmaE[i] * EE[x], 0, 1) for i in range(2)],0)),range(p))     
    
    return pval
