# -*- coding: utf-8 -*-
"""
Created on Wed Jul 25 14:35:23 2018

a module for orthogonal linear separation analysis (OLSA)

@author: setsuo,shotaro, and tadahaya
"""
import sys
import csv
import math
import os
import numpy as np
import pandas as pd
np.seterr(divide='ignore', invalid='ignore')
import time
from scipy import stats
from scipy.cluster.hierarchy import ward,leaves_list

class DataClass:
    """ Handle data """
    def __init__(self):
        self.filename=""
        self.Name=list()
        self.X = np.array([[],[]])
        self.index = list()

    def load(self,filename,read_tsv=False):
        """Load data into an object"""
        self.filename=filename
        print('reading file')
        if read_tsv==False:
            csv_obj = csv.reader(open(filename, "r"))
        else:
            csv_obj = csv.reader(open(filename, "r"),delimiter='\t')
        data = [v for v in csv_obj]
        self.Name = data[0]
        #delete the 1st column
        del self.Name[0]
        print("sample names are")
        print(self.Name)
        del data[0]
        data_conved = [[float(elm) for elm in v[1:]] for v in data] #Text to float
        self.index = [v[0] for v in data]
        self.X = np.matrix(data_conved)
        print('data read OK')

    def load_df(self,dataframe):
        """Load dataframe into an object"""
        self.X = np.array(dataframe)
        self.Name = list(dataframe.columns)
        self.index = list(dataframe.index)
        print('a dataframe was loaded')

    def clone(self):
        ret = DataClass()
        ret.Name = self.Name[:] #DeepCopy
        ret.X = np.array(self.X) #DeepCopy
        ret.index = self.index[:] #DeepCopy
        ret.filename = self.filename #DeepCopy
        return ret

    def deleteC(self,index):
        """delete the selected columns"""
        a = range(0,len(self.Name))
        b = list()
        if isinstance(index, int):
            b=list([index])
        else:
            b=index[:]
        c=list(set(a)-set(b))
        return self.selectC(c)

    def selectC(self,index):
        """extract the selected columns"""
        ret = self.clone()
        ret.X = self.X[:,index]
        ret.filename = self.filename
        if isinstance(index, int):
            ret.Name = list([self.Name[index]])
        else:
            ret.Name = list([self.Name[s] for s in index])
        return ret

    def selectR(self,index):
        """extract the selected rows"""
        ret = self.clone()
        ret.X = self.X[index,:]
        ret.filename = self.filename
        if isinstance(index, int):
            ret.index = list([self.index[index]])
        else:
            ret.Name = list([self.Name[s] for s in index])
        return ret

    def bindC(X,Y):
        ret = DataClass()
        ret.Name = list(X.Name + Y.Name)
        ret.X = np.c_[X.X,Y.X]
        ret.index = list(X.index)
        return ret

    def wardclustersort(self,IsSpheredSort=True):
        """
        sort columns according to Ward clustering
        
        Parameters
        ----------    
        IsSpheredSort: boolean, default True
            apply this function to the unit-sphereized data 
        
        """
        ret = DataClass()
        ret.index = list(self.index)
        ret.filename = self.filename
        XS=self.X
        if IsSpheredSort==True:
            XS= np.array([s/np.linalg.norm(s) for s in np.array(self.X.T)]).T
        sortindex = leaves_list(ward(XS.T))
        ret.Name = list([self.Name[sortindex[0]]])
        ret.X = self.X[:,sortindex[0]]
        for i in sortindex[1:]:
            ret.Name = list(ret.Name + [self.Name[i]])
            ret.X = np.c_[ret.X,self.X[:,i]]
        return ret


class Result:
    """ Handle OLSA results """
    def __init__(self):
        self.accumulation = float()
        self.sphered = bool()
        self.filename=""
        self.Name = list() #sample name list
        self.index = list() #variable name list
        self.X = np.matrix([[],[]])
        self.TS = np.matrix([[],[]]) #total strength list
        self.Rpca = np.matrix([[],[]]) #response vector matrix
        self.Cpca=np.matrix([[],[]]) #contribution
        self.AcceptNum=int()
        self.NWpca = np.matrix([[],[]]) #response score matrix
        self.ps =list() #vector name list
    

    def export(self,savefilename='',TS=True,Contribution=True,RSM=True,RVM=True
               ,Raw=False,WxTS=False,CM=False,CMex=False,Confirmation=False):
        """
        export data into a csv file
        
        Parameters
        ----------
        savefilename: str, default ""
            a path of the output. if no description, a path is generated from the input filename
            
        TS: boolean, default True
            whether TS is exported
            
        Contribution: boolean, default True
            whether contribution is exported

        RSM: boolean, default True
            whether response score matrix is exported

        RVM: boolean, default True
            whether response vector matrix is exported
        
        Raw: boolean, default True
            whether raw is exported
        
        WxTS: boolean, default False
            whether weight x TS is exported
        
        CM: boolean, default False
            whether correlation matrix of response score is exported
        
        CMex: boolean, default False
            whether correlation matrix of response score w/o ones of
             the highest and lowest vector is exported

        Confirmation: boolean, default False
            whether data for logical confirmation is exported
        
        """
        def makelisttable(row,col,inputdata):
            data = np.array(inputdata)
            a = list(["index"] + list(col))
            if len(row)>1:
                rn = np.shape(data)[0]
                b = [a] + list([[row[i]]+list(data[i]) for i in range(0,rn)])
            else:
                b = [a] + list([list(row)+list(data)])
            return b
        
        #path preparation
        if len(savefilename) > 0:
            basefilename = savefilename.replace('.csv','')
        elif len(self.filename) > 0:
            dirname = os.path.dirname(self.filename)
            filename = os.path.basename(self.filename).split('.')[0]
            basefilename = '{0}\\{1}'.format(dirname,filename)
        else:
            print('ERROR!!: No original filename or savefilename')
            print('Enter savefilename in export()')
            sys.exit(1)
        
        #Saving
        if TS==True:
            with open('{0}_TS.csv'.format(basefilename), 'w') as f:
                writer = csv.writer(f, lineterminator='\n') #indicate the newline code （\n）
                writer.writerows(makelisttable([""],self.Name,self.TS))
        else: pass

        if Contribution==True:
            with open('{0}_Cont.csv'.format(basefilename), 'w') as f:
                writer = csv.writer(f, lineterminator='\n') #indicate the newline code （\n）
                writer.writerows(makelisttable([""],self.ps,self.Cpca))
        else: pass

        if RSM==True:
            with open('{0}_RSM.csv'.format(basefilename), 'w') as f:
                writer = csv.writer(f, lineterminator='\n') #indicate the newline code （\n）
                writer.writerows(makelisttable(self.ps,self.Name,self.NWpca))
        else: pass

        if WxTS==True:
            with open('{0}_WxTS.csv'.format(basefilename), 'w') as f:
                writer = csv.writer(f, lineterminator='\n') #indicate the newline code （\n）
                wt = self.NWpca.dot(np.diag(self.TS))
                writer.writerows(makelisttable(self.ps,self.Name,wt))
        else: pass

        if CM==True:
            with open('{0}_CM.csv'.format(basefilename), 'w') as f:
                writer = csv.writer(f, lineterminator='\n') #indicate the newline code （\n）
                writer.writerows(makelisttable(self.Name,self.Name,np.corrcoef(self.NWpca.T)))
        else: pass

        if CMex==True:
            with open('{0}_CMex.csv'.format(basefilename), 'w') as f:
                writer = csv.writer(f, lineterminator='\n') #indicate the newline code （\n）
                writer.writerows(makelisttable(self.Name,self.Name,np.corrcoef(self.NWpca.T[:,1:self.AcceptNum])))
        else: pass

        if RVM==True:
            with open(basefilename + '_RVM.csv', 'w') as f:
                writer = csv.writer(f, lineterminator='\n') #indicate the newline code （\n）
                writer.writerows(makelisttable(self.index,self.ps,self.Rpca))
        else: pass

        if Raw==True:
            with open('{0}_Raw.csv'.format(basefilename), 'w') as f:
                writer = csv.writer(f, lineterminator='\n') #indicate the newline code （\n）
                writer.writerows(makelisttable(self.index,self.Name,self.X))
        else: pass

        if Confirmation==True:
            with open('{0}_Confirmation.csv'.format(basefilename), 'w') as f:
                writer = csv.writer(f, lineterminator='\n') #indicate the newline code （\n）
                mat = self.Rpca.dot(self.NWpca.dot(np.diag(self.TS)))
                writer.writerows(makelisttable(self.index,self.Name,mat))
        else: pass
        print("data save finished")

    def rsm(self,varimax_only=True):
        """
        return response score matrix as a dataframe        
        """
        df = pd.DataFrame(self.NWpca)
        df.index = self.ps
        df.columns = self.Name
        if varimax_only==True:
            df2 = df.iloc[:self.AcceptNum,:]
        else:
            df2 = df
        return df2
    
    def rvm(self,varimax_only=True):
        """
        return response vector matrix as a dataframe        
        """
        df = pd.DataFrame(self.Rpca)
        df.index = self.index
        df.columns = self.ps
        if varimax_only==True:
            df2 = df.iloc[:,:self.AcceptNum]
            print("return response vector matrix")
        else:
            df2 = df
        return df2
        
    def ts(self):
        """
        return total strength as a dataframe        
        """
        df = pd.DataFrame(self.TS)
        df.index = self.Name
        df.columns = ["total strength"]
        df2 = df.T
        return df2           
    
    def contribution(self):
        """
        return vector contribution as a dataframe        
        """
        df = pd.DataFrame(self.Cpca)
        df.index = self.ps
        df.columns = ["contribution"]
        df2 = df.T
        return df2  
    
    def weightedTS(self):
        """
        return weighted total strength as a dataframe        
        """
        wt = self.NWpca.dot(np.diag(self.TS))
        df = pd.DataFrame(wt)
        df.index = self.ps
        df.columns = self.Name
        return df
    
    def cm(self):
        """
        return correlation matrix as a dataframe        
        """
        df = pd.DataFrame(np.corrcoef(self.NWpca.T))
        df.index = self.Name
        df.columns = self.Name
        return df

    def export_at_once(self,savefilename='',TS=True,Contribution=True,RSM=True,RVM=True
                      ,Raw=True,WxTS=False,CM=False,CMex=False,Confirmation=False):
        """
        export all data into a csv file
        
        Parameters are the same with those of "export"
        
        """
        def makelisttable(row,col,inputdata):
            data = np.array(inputdata)
            a = list(["index"] + list(col))
            if len(row)>1:
                rn = np.shape(data)[0]
                b = [a] + list([[row[i]]+list(data[i]) for i in range(0,rn)])
            else:
                b = [a] + list([list(row)+list(data)])
            return b
        
        #path preparation
        if len(savefilename) > 0: pass
        else:
            dirname = os.path.dirname(self.filename)
            filename = os.path.basename(self.filename).split('.')[0]
            savefilename = '{0}\\{1}_res.csv'.format(dirname,filename)        
        
        #Saving
        with open(savefilename, 'w') as f:
            writer = csv.writer(f, lineterminator='\n') #indicate the newline code （\n）

            if self.sphered==False:
                writer.writerow(["This data is NOT sphered."])
            else:
                writer.writerow(["This data is sphered."])

            if TS==True:
                writer.writerow("")
                writer.writerow(["Total Strength"])
                writer.writerows(makelisttable([""],self.Name,self.TS))
            else: pass
            
            if Contribution==True:
                writer.writerow("")
                writer.writerow(["Contribution of PCA"])
                writer.writerow(['The number of principal components is {0}, accounting for 95% of cumulative contribution'.format(self.AcceptNum)])
                writer.writerows(makelisttable([""],self.ps,self.Cpca))
            else: pass
            
            if RSM==True:
                writer.writerow("")
                writer.writerow(["Response Score Matrix"])
                writer.writerows(makelisttable(self.ps,self.Name,self.NWpca))
            else: pass

            if WxTS==True:
                writer.writerow("")
                writer.writerow(["Weight * diag(TS) for graphical modeling"])
                wt = self.NWpca.dot(np.diag(self.TS))
                writer.writerows(makelisttable(self.ps,self.Name,wt))
            else: pass

            if CM==True:
                writer.writerow("")
                writer.writerow(["Correlation of Response Score Matrix"])
                writer.writerows(makelisttable(self.Name,self.Name,np.corrcoef(self.NWpca.T)))
            else: pass

            if CMex==True:           
                writer.writerow("")
                writer.writerow(["Correlation of Response Score Matrix Excluding First and Minor Components"])
                writer.writerow(self.Name)
                writer.writerows(makelisttable(self.Name,self.Name,np.corrcoef(self.NWpca.T[:,1:self.AcceptNum])))
            else: pass

            if RVM==True:
                writer.writerow("")
                writer.writerow(["Response Vector Matrix"])
                writer.writerows(makelisttable(self.index,self.ps,self.Rpca))
            else: pass

            if Raw==True:
                writer.writerow("")
                writer.writerow(["Raw Data"])
                writer.writerows(makelisttable(self.index,self.Name,self.X))
            else: pass

            if Confirmation==True:
                writer.writerow("")
                writer.writerow(["Reponse * Weight * diag(TS) for logical confirmation"])
                mat = self.Rpca.dot(self.NWpca.dot(np.diag(self.TS)))
                writer.writerows(makelisttable(self.index,self.Name,mat))
            else: pass

        print("data save finished")