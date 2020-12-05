# -*- coding: utf-8 -*-
"""
Created on Wed Jul 25 14:35:23 2018

a module for orthogonal linear separation analysis (OLSA)

@author: setsuo, shotaro, and tadahaya
"""
import sys
import csv
import math
import os
import numpy as np
import pandas as pd
np.seterr(divide='ignore', invalid='ignore')
import time
from sklearn.decomposition import PCA
from scipy import stats
from scipy.cluster.hierarchy import ward,leaves_list

from .data import Result


class Varimax:
    """Varimax rotation"""
    def __init__(self, X):

        self.n = np.shape(X)[0]
        self.p = np.shape(X)[1]
        self.A = np.matrix(X)
        print("ok1")

        if self.p > 1:
            self.h = []
            for i in range(0, self.n):
                sum = self.A[i].dot(self.A[i].T)[0,0]
                self.h.append(math.sqrt(sum))
                if self.h[i]!=0:
                    self.A[i] /= self.h[i]

    def rotateV(self,acceptable_error=1.0e-9):
        """varimax rotation"""
        mm = lambda x: x*x
        if (self.p < 2):
            return self.A
        while True:
            ckA = np.matrix(self.A.A)
            for i in range(0, self.p):
                for j in range(i + 1, self.p):
                    x = np.matrix(self.A.T[i].A)
                    y = np.matrix(self.A.T[j].A)
                    u = np.matrix(x.A**2 - y.A**2)
                    v = np.matrix(2.0 * x.A * y.A)
                    cA = np.sum(u)
                    cB = np.sum(v)
                    cC = u.dot(u.T)[0,0] - v.dot(v.T)[0,0]
                    cD = 2.0 * u.dot(v.T)[0,0]

                    num = cD - 2 * cA * cB / self.n
                    den = cC - (mm(cA) - mm(cB)) /self.n
                    theta4 = math.atan(num / den)
                    if (num > 0.0):
                        if (theta4 < 0.0):
                            theta4 += math.pi
                    else:
                        if (theta4 > 0.0):
                            theta4 += math.pi
                    theta = theta4 / 4.0
                    tx = self.A.T[i] * math.cos(theta) + self.A.T[j] * math.sin(theta)
                    ty = -self.A.T[i] * math.sin(theta) + self.A.T[j] * math.cos(theta)

                    self.A.T[i] = tx
                    self.A.T[j] = ty

            dif = np.sum((ckA.A-self.A.A)**2)

            print("\r" + str(dif),end="")
            if (dif < acceptable_error):
                for i in range(0, self.n):
                    self.A[i] *= self.h[i]
                break

            if math.isnan(dif):
                print("error")
                sys.exit()
        print("")
        return self.A
       

class SmirnovGrubbs:
    def __init__(self):
        self.RemainedIndexes = list()
        self.RemovedIndexesHigher = list()
        self.RemovedIndexesLower = list()
        self.NonbiasedVar = float()
        self.alpha = float()


def calcSG(TS,alpha):
    """conduct SG test to exclude the data with unusual TS"""
    res = SmirnovGrubbs()
    res.alpha = alpha
    Data = list()
    RemovedDataHigher = list()
    RemovedDataLower = list()

    for i in range(0,TS.shape[0]):
        Data.append([i,TS[i]])

    while True:
        n=len(Data)
        if n<3:
            break
        t = stats.t.isf((alpha/n) / 2, n-2)
        Gtest = (n-1)/math.sqrt(n) * math.sqrt(t**2 / (n-2 + t**2))

        mean=0.0
        for d in Data:
            mean = mean + d[1]
        mean = mean / n

        var = 0.0
        for d in Data:
            var = var + (d[1]-mean)*(d[1]-mean)
        var = var / n
        sd = math.sqrt(var)

        maxindex = Data[0][0]
        maxvalue = math.fabs(Data[0][1] - mean)
        for i in range(0,len(Data)):
            if maxvalue < math.fabs(Data[i][1] - mean):
                maxindex = i
                maxvalue =  math.fabs(Data[i][1] - mean)

        #SmirnovGrubbs
        if maxvalue / sd > Gtest:
            if (Data[maxindex][1]-mean)>0:
                RemovedDataHigher.append([Data[maxindex][0],Data[maxindex][1]])
            else:
                RemovedDataLower.append([Data[maxindex][0],Data[maxindex][1]])
            del Data[maxindex]
        else:
            break

    mean=0.0
    for d in Data:
        mean = mean + d[1]
    mean = mean / n

    ubvar = 0.0
    for d in Data:
        ubvar = ubvar + (d[1]-mean)*(d[1]-mean)
    ubvar = ubvar / (n-1.0)

    res.NonbiasedVar = ubvar

    for d in Data:
        res.RemainedIndexes.append(int(d[0]))

    for d in RemovedDataHigher:
        res.RemovedIndexesHigher.append(int(d[0]))

    for d in RemovedDataLower:
        res.RemovedIndexesLower.append(int(d[0]))

    return res


def usspca(Data,accumulation=0.95,IsSphered=True,
           UseMirror=True,UseSmirnovGrubbs=True,WardClusterSort=True):
    """
    conduct Unit Spherized Symmetric PCA
    
    Parameters
    ----------
    Data: DataClass object
        subjected to analysis
        
    accumulation: float, default 0.95
        % of cumulative contribution of the calculated vectors
    
    IsSphered: boolean, default True
        whether data is unit-sphereized before calculation
    
    UseMirror: boolean, default True
        whether the input data set is combined with the origin-symmetric set before calculation
    
    UseSmirnovGrubbs: boolean, default True
        whether outliers are excluded according to SG test 
    
    WardClusterSort: boolean, default True
        whether response score matrix is sorted according to clustering with ward method 
    
    """
    if WardClusterSort==True: Data2 = Data.wardclustersort()
    else: Data2 = Data

    res = Result()
    res.accumulation = accumulation
    res.sphered = IsSphered
    res.Name = Data2.Name[:]
    res.index = Data2.index[:]
    res.X = np.array(Data2.X)
    res.filename = Data2.filename

    print ('********************************')
    print ('The conditions of this proccess')
    print(('Contribution Accumulation Max is ' + str(accumulation)))
    print(('IsSphered Flag is ' + str(IsSphered)))
    print(('UseMirror Flag is ' + str(UseMirror)))
    print(('UseSmirnovGrubbs Flag is ' + str(UseSmirnovGrubbs)))
    print ('********************************')  

    #store Total Strength
    res.TS = np.array([np.linalg.norm(s) for s in res.X.T])

    if UseSmirnovGrubbs==True:
        SGres = calcSG(res.TS,0.05)

        print("excluded samples by SG test")
        print("too large norm：")
        print(list([res.Name[s] for s in SGres.RemovedIndexesHigher]))
        print("too small norm：")
        print(list([res.Name[s] for s in SGres.RemovedIndexesLower]))

        #data reconstruction
        remainedX = np.array([res.X.T[s] for s in SGres.RemainedIndexes]).T
    else:
        #data reconstruction
        remainedX = res.X
       
    #map the data on a unit sphere
    XS = np.array(res.X)
    remainedXS = np.array(remainedX)
    if res.sphered==True:
        XS = np.array([s/np.linalg.norm(s) for s in res.X.T]).T
        remainedXS = np.array([s/np.linalg.norm(s) for s in remainedX.T]).T

    #add a mirror data set origin-symmetric to the analyzed set
    #in order to approximate the centroid to the origin
    if UseMirror==True:
        Xpca =np.c_[remainedXS, remainedXS*-1]
    else:
        Xpca =np.c_[remainedXS]
    
    #calculate the centroid of the data set
    #np.c_[] for conversion into nx1 vector
    Centroid = np.c_[np.array([sum(s)/len(s) for s in Xpca])] 
    
    print('data preparation OK')

    ##PCA
    decomposer = PCA()
    decomposer.fit(Xpca.T)

    print('PCA calculation OK')

    #Response Profile Matrix
    SampleNum = np.shape(remainedX)[1]
    res.Rpca = decomposer.components_.T [:,0:SampleNum]

    #Contributaion Vector
    res.Cpca=decomposer.explained_variance_ratio_[0:SampleNum]
    res.ps =['%s%d' % ('P',s) for s in range(1,len(res.Cpca)+1)] #label like p1~

    AccumuC = 0.0
    res.AcceptNum=0
    for i in range(0,len(res.Cpca)):
        res.AcceptNum = res.AcceptNum+1
        AccumuC = AccumuC + res.Cpca[i]
        if AccumuC > res.accumulation:
            break

    print("The selected components number is %d" % res.AcceptNum)
    print("The accumulation value is %f" % AccumuC)

    #Normalized Weight Matrix
    #subtract the centroid from the input data in PCA
    #as for OLSA, the centroid matches the origin
    res.NWpca = res.Rpca.T.dot(XS-Centroid)

    return res


def olsa(Data,accumulation=0.60,FixSet=set(),IsSphered=True,
         UseMirror=True,UseSmirnovGrubbs=True,acceptable_error=1.0e-6,WardClusterSort=True):
    """
    conduct orthogonal linear separation analysis (OLSA) and return a Result object
    
    Parameters
    ----------
    Data: DataClass object
        subjected to analysis
        
    accumulation: float, default 0.60
        % of cumulative contribution of the vectors subjected to varimax rotation
    
    FixSet: int
        indicates the column number not subjected to varimax rotation
    
    IsSphered: boolean, default True
        whether data is unit-sphereized before calculation
    
    UseMirror: boolean, default True
        whether the input data set is combined with the origin-symmetric set before calculation
    
    UseSmirnovGrubbs: boolean, default True
        whether outliers are excluded according to SG test
    
    acceptable_error: float, default 1.0e-9
        determines acceptable error in varimax rotation
        
    WardClusterSort: boolean, default True
        whether response score matrix is sorted according to clustering with ward method 
    
    """    
    start = time.time()
    res = usspca(Data,accumulation,IsSphered,UseMirror
                 ,UseSmirnovGrubbs,WardClusterSort)
    
    #store filename
    res.filename = Data.filename
    
    XS = np.array(res.X)
    if IsSphered==True:
        XS = np.array([s/np.linalg.norm(s) for s in res.X.T]).T

    #subject the vectors accounting for *% of cumulative contribution to varimax rotation 
    #some vectors can be fixed if indicated

    Vselect = set(range(0,res.AcceptNum))
    if isinstance(FixSet, int):
        Vselect = set(Vselect)-set([FixSet])
    else:
        Vselect = set(Vselect)-set(FixSet)
    Vtotal = set(range(0,np.shape(res.Rpca)[1]))
    VResidue = Vtotal - Vselect
    print("Vselect = " , Vselect)
    print("VResidue = " , VResidue)

    ForVarimaxData = res.Rpca[:,list(Vselect)]

    for i in list(Vselect):
        res.ps[i] = "%sV" % res.ps[i]
    res.ps = [res.ps[i] for i in list(Vselect)+list(VResidue)]

    print("Call Varimax class ")
    rot = Varimax(ForVarimaxData)
    rotated = np.array(rot.rotateV(acceptable_error))

    #Sort of rotated
    SeudoNW=rotated.T.dot(XS)
    SeudoVars = np.array([s.var() for s in np.c_[SeudoNW,SeudoNW* -1]])
    Slist=list()
    for i in range(0,len(SeudoVars)):
        Slist.append([SeudoVars[i],i])
    #ascending
    Slist.sort()
    #reverse the ranking
    Slist.reverse()

    Srotated=np.array(rotated)
    for i in range(0,len(Slist)):
        Srotated[:,i]=rotated[:,Slist[i][1]]

    ###loading matrix
    res.Rpca = np.c_[Srotated,res.Rpca[:,list(VResidue)]]
    print("Calc New NWpca")

    ###loading matrix
    res.NWpca = res.Rpca.T.dot(XS)
    print("Varimax Finished")

    #calculate the contribution ratio
    Vars = np.array([s.var() for s in np.c_[res.NWpca,res.NWpca * -1]])
    res.Cpca = Vars/Vars.sum()

    #display elapsed time
    passed = time.time() - start
    h,m = divmod(passed,3600)
    m,s = divmod(m,60)
    print('elapsed time: {0} hr {1} min {2} sec'.format(int(h),int(m),round(s,-2)))
    return res
