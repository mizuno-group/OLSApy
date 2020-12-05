# -*- coding: utf-8 -*-
"""
Created on Mon Sep 23 13:05:34 2019

module for plot

@author: tadahaya

"""

import pandas as pd
import numpy as np
from scipy.stats import rankdata
import matplotlib.pyplot as plt
import seaborn as sns
import csv
from pathlib import Path

class RadarChart():
    def __init__(self,figsize=None):
        """ prep polar chart """
        if figsize is not None:
            self.fig = plt.figure(figsize=figsize)
        else:
            self.fig = plt.figure()
        self.ax = self.fig.add_subplot(111,polar=True)
        plt.gca().spines['polar'].set_color('lightgray')
        plt.gca().spines['start'].set_color('lightgray')
        plt.gca().spines['end'].set_color('lightgray')
        plt.gca().spines['inner'].set_color('lightgray')
        plt.grid(color='lightgray')
        self.feature = None # list
        self.rlim = None # tuple


    def plot(self,values,features,label="",limit=None,color="",
            alpha=0.25,markersize=5,axes_label=True,labelsize=10):
        """
        radar chart
        
        Parameters
        ----------
        values: array
            values to be plotted

        features: list
            list of features to be plotted
        
        label: str
            sample name
            
        limit: tuple
            determine min and max values of plot: (min,max)
        
        color: str
            determine the color of plot
                
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
                        
        """    
        # prep
        if len(values)!=len(features):
            raise ValueError("!! len(values) and len(features) should be equal !!")
        angles = np.linspace(0,2 * np.pi, len(features) + 1,endpoint=True)
        maxvalue = np.max(values)
        minvalue = np.min(values)
        values = np.concatenate((values,[values[0]])) # since it rounds once

        # limit
        if limit==None:
            limit = (minvalue + minvalue/len(features),maxvalue + maxvalue/len(features))
        else:
            limit = limit
        if self.rlim is None:
            self.rlim = limit
        else:
            if self.rlim[0] > limit[0]:
                self.rlim = (limit[0],self.rlim[1])
            if self.rlim[1] < limit[1]:
                self.rlim = (self.rlim[0],limit[1])

        # feature
        if self.feature is None:
            self.feature = features
        else:
            self.feature = [v if len(v) > 0 else self.feature[i] for i,v in enumerate(features)]

        # visualize
        self.ax.plot(angles,np.zeros((len(angles),)),"-",color="black",linewidth=.5) # 0 line
        if len(color)==0:
            l2d, = self.ax.plot(angles,values,"o-",markersize=markersize) # determines the outer frame
            self.ax.fill(angles,values,alpha=alpha) # fill the frame
        else:
            l2d, = self.ax.plot(angles,values,"o-",color=color,label=label,markersize=markersize) # determines the outer frame
            self.ax.fill(angles,values,alpha=alpha,color=color) # fill the frame
        self.ax.set_theta_zero_location('N', offset=0) # determines the location of zero 'North'
        self.ax.set_theta_direction(-1) # determines the rotation: clockwise
        if axes_label:
            self.ax.set_thetagrids(angles[:-1] * 180 / np.pi, self.feature, fontsize=labelsize) #label of axis
        else:
            self.ax.set_thetagrids(angles[:-1] * 180 / np.pi, [""]*len(self.feature)) #label of axis
        self.ax.set_rlim(self.rlim)
        return l2d


    def close(self,savefig="",dpi=100,handles=[],labels=[]):
        """
        close the canvas

        Parameters
        ----------
        savefig: str
            path for exported file

        dpi: int
            indicates the dpi of the output

        """
        if (len(handles) > 0) and (len(labels) > 0):
            plt.legend(handles=handles,labels=labels,loc="best")
        else:
            plt.legend(loc="lower center",bbox_to_anchor=(0.5,-0.25))
        plt.tight_layout()
        if len(savefig) > 1:
            plt.savefig(fname=savefig,dpi=dpi,format="tif",bbox_inches="tight")
        plt.show()
        plt.close(self.fig)