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
from sklearn.decomposition import PCA
from scipy import stats
from scipy.cluster.hierarchy import ward,leaves_list

from .data import DataClass,Result
from .plot import RadarChart
from . import calculator as calc


class OLSA():
    def __init__(self):
        self.data = DataClass()
        self.res = Result()


    ### setter & getter ###
    def load_df(self,df):
        """
        load data
        
        Parameters
        ----------
        df: dataframe
            feature x sample dataframe
        
        """
        self.data.load_df(df)


    def get_input(self):
        """ get loaded data """
        if len(self.data.index)==0:
            raise ValueError("!! No data loaded !!")
        print("return data, sample list, and feature list")
        return self.data.X,self.data.Name,self.data.index


    def get_res(self):
        """ get main output """
        if self.res.X.shape[0]==0:
            raise ValueError("!! No result: use olsa or usspca before this !!")
        print("return Response Score Matrix, Response Vector Matrix, Total Strength, and Contribution")
        rsm = self.res.rsm()
        rvm = self.res.rvm()
        ts = self.res.ts()
        cont = self.res.contribution()
        return rsm,rvm,ts,cont


    ### calculation ###
    def olsa(self,accumulation:float=0.60,acceptable_error:float=1.0e-6,FixSet:set=set(),IsSphered:bool=True,
            UseMirror:bool=True,UseSmirnovGrubbs:bool=True,WardClusterSort:bool=True):
        """
        conduct orthogonal linear separation analysis (OLSA) and return a Result object
        
        Parameters
        ----------
        accumulation: float, default 0.60
            % of cumulative contribution of the vectors subjected to varimax rotation
        
        acceptable_error: float, default 1.0e-6
            determines acceptable error in varimax rotation

        FixSet: int
            indicates the column number not subjected to varimax rotation
        
        IsSphered: boolean, default True
            whether data is unit-sphereized before calculation
        
        UseMirror: boolean, default True
            whether the input data set is combined with the origin-symmetric set before calculation
        
        UseSmirnovGrubbs: boolean, default True
            whether outliers are excluded according to SG test
                    
        WardClusterSort: boolean, default True
            whether response score matrix is sorted according to clustering with ward method 
        
        """
        if len(self.data.index)==0:
            raise ValueError("!! No data: load data before calculation !!")
        self.res = calc.olsa(self.data,accumulation=accumulation,acceptable_error=acceptable_error,FixSet=FixSet,
                             IsSphered=IsSphered,UseMirror=UseMirror,UseSmirnovGrubbs=UseSmirnovGrubbs,
                             WardClusterSort=WardClusterSort)


    def usspca(self,accumulation:float=0.95,IsSphered:bool=True,
            UseMirror:bool=True,UseSmirnovGrubbs:bool=True,WardClusterSort:bool=True):
        """
        conduct Unit Spherized Symmetric PCA
        
        Parameters
        ----------
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
        if len(self.data.index)==0:
            raise ValueError("!! No data: load data before calculation !!")
        self.res = calc.usspca(self.data,accumulation=accumulation,IsSphered=IsSphered,UseMirror=UseMirror,
                               UseSmirnovGrubbs=UseSmirnovGrubbs,WardClusterSort=WardClusterSort)


    def export(self,savefilename:str='',TS:bool=True,Contribution:bool=True,RSM:bool=True,
               RVM:bool=True,Raw:bool=False,WxTS:bool=True,CM:bool=False,CMex:bool=False,
               Confirmation:bool=False):
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
        
        WxTS: boolean, default True
            whether weight x TS is exported
        
        CM: boolean, default False
            whether correlation matrix of response score is exported
        
        CMex: boolean, default False
            whether correlation matrix of response score w/o ones of
             the highest and lowest vector is exported

        Confirmation: boolean, default False
            whether data for logical confirmation is exported
        
        """
        self.res.export(savefilename,TS,Contribution,RSM,RVM,Raw,WxTS,CM,CMex,Confirmation)


    ### visualization ###
    def plot(self,score=None,focus:list=[],absolute:bool=True,overlay:bool=False,
             color:str="mediumblue",n_labeled:int=None,plot_all:bool=False,
             figsize:tuple=None,savefig:str="",dpi:int=100,**kwargs):
        """
        plot scores as a radarchart
        
        Parameters
        ----------
        score: dataframe
            data to be visualized (contracted x samples) such as response score matrix

        absolute: boolean
            indicates whether data is converted into absolute values before visualization

        focus: list
            a list of sample names to be visualized

        overlay: boolean
            whether plots are overlaid or not

        color: str or list
            indicates the color of plot

        n_labeled: int
            indicates the No. of labels to be plotted

        plot_all: boolean
            whether all samples are visualized or not

        limit: tuple
            determine min and max values of plot: (min,max)
        
        color: str
            determine the color of plot
        
        figsize: tupple
            determine figsize
        
        markersize: int
            determine size of markers
            
        labelsize: int
            determine size of axis-labels
            
        ticksize: int
            deterimne font size of polar axis
                        
        axes_label: boolean
            whether labels of axes are neeeded
        
        alpha: float, default 0.25
            indicates transparency
            
        savefig: str
            path for exported file
            
        """    
        if score is None:
            if self.res.X.shape==0:
                raise ValueError("!! Indicate score such as responase score !!")
            else:
                score = self.res.rsm()
        if len(focus)==0:
            focus = list(score.columns)
            if len(focus) > 10:
                if plot_all==False:
                    raise ValueError("!! Too many samples: turn on plot_all to plot all data !!")
        feature = list(score.index)
        abs_score = np.abs(score)
        if absolute:
            score = abs_score
        if overlay:
            if type(color)!=list:
                color = []
            dat = RadarChart()
            axes = []
            for i,f in enumerate(focus):
                if len(color)==0:
                    color = [""]*len(focus)
                else:
                    div,mod = divmod(len(focus),len(color))
                    color = color*div + color[:mod]
                temp = score[f].values
                temp2 = abs_score[f].sort_values(ascending=False)
                if n_labeled is None:
                    l2d = dat.plot(values=temp,features=feature,label=f,color=color[i],**kwargs)
                else:
                    shown = list(temp2.index)[:n_labeled]
                    feature2 = [v if v in shown else "" for v in feature]
                    l2d = dat.plot(values=temp,features=feature2,label=f,color=color[i],**kwargs)
                axes.append(l2d)
            dat.close(savefig=savefig,dpi=dpi,handles=axes,labels=focus)
        else:
            if type(color)!=str:
                raise TypeError("!! Color should be string when overlay is False !!")
            for f in focus:
                if figsize is None:                
                    figsize = (4,6) # hard coding to align the sizes of plots
                dat = RadarChart(figsize=figsize)
                temp = score[f].values
                temp2 = abs_score[f].sort_values(ascending=False)
                if n_labeled is None:
                    dat.plot(values=temp,features=feature,label=f,color=color,**kwargs)
                else:
                    shown = list(temp2.index)[:n_labeled]
                    feature2 = [v if v in shown else "" for v in feature]
                    dat.plot(values=temp,features=feature2,label=f,color=color,**kwargs)
                dat.close(savefig,dpi)
