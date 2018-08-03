# -*- coding: utf-8 -*-
"""
Created on Wed Aug  1 11:44:52 2018

to test whether olsa correctly runs

@author: tadahaya
"""
import pytest
import sys
sys.path.append(".\\main") #to make a path for the modules
sys.path.append(".\\tests")
from olsapy import DataClass,olsa,usspca
import pandas as pd

### fixures ###
@pytest.fixture
def data_path():
    return "tests\\data\\test_data.txt" #CMap derived signatures, tsv format

@pytest.fixture
def df_test():
    df = pd.read_csv("tests\\data\\test_data.txt",delimiter="\t",index_col=0)
    return df

@pytest.fixture
def dat_test():
    return DataClass()

@pytest.fixture
def name_test():
    return "tests\\data\\test_data_res.txt"

@pytest.fixture
def res_test():
    dat = DataClass()
    dat.load("tests\\data\\test_data.txt",read_tsv=True)
    res = olsa(dat)
    return res

@pytest.fixture
def dat_test2():
    dat = DataClass()
    dat.load("tests\\data\\test_data.txt",read_tsv=True)
    return dat


### tests ###
def test_load1(data_path,dat_test): #normal run
    dat_test.load(data_path,read_tsv=True)
    res = olsa(dat_test)
    res.export()
"""
1. generate a DataClass object
2. load a file path of data (usually csv is supposed but tsv is employed considering the size)
3. run OLSA and obtain a Result object
4. export result files as dataframes
"""

def test_load2(df_test,dat_test,name_test): #run after loading dataframe
    dat_test.load_df(df_test) #.load_df: load dataframe into a DataClass object
    res = olsa(dat_test)
    res.export(savefilename=name_test) #savefilename is necessary
                                       #because self.filename does not exist in this case

def test_run1(dat_test2): #normal run
    res = olsa(dat_test2)
    res.export()

@pytest.mark.parametrize("accum,fix,issph,usemi,usesm,accept,ward",
        [
        (0.5,"",True,True,True,1.0e-9,True), #accumulation modification
        (0.8,0,True,True,True,1.0e-9,True), #vector fixation
        (0.8,"",False,True,True,1.0e-9,True), #sphererization option
        (0.8,"",True,False,True,1.0e-9,True), #mirror data option
        (0.8,"",True,True,False,1.0e-9,True), #SG test option
        (0.8,"",True,True,True,1.0e-6,True), #error rate change
        (0.8,"",True,True,True,1.0e-9,False) #sorting option
        ])

def test_run2(dat_test2,accum,fix,issph,usemi,usesm,accept,ward): #check olsa() options
    res = olsa(dat_test2,accumulation=accum,FixSet=fix,IsSphered=issph,UseMirror=usemi,
               UseSmirnovGrubbs=usesm,acceptable_error=accept,WardClusterSort=ward)
    res.export()  

def test_run3(dat_test2): #check usspca (without varimax rotation)
    res = usspca(dat_test2)
    res.export()

@pytest.mark.parametrize("raw,wxts,cm,cmex,conf",
        [
        (True,False,False,False,False), #Raw activated
        (False,True,False,False,False), #WxTS activated
        (False,False,True,False,False), #CM activated
        (False,False,False,True,False), #CMex activated
        (False,False,False,False,True), #Confirmation activated
        ])

def test_result1(res_test,raw,wxts,cm,cmex,conf): #Raw activated
    res_test.export(Raw=raw,WxTS=wxts,CM=cm,CMex=cmex,Confirmation=conf)
    
def test_result2(res_test): #extract each result
    res_test.rsm()
    res_test.rvm()
    res_test.ts()
    res_test.contribution()
    res_test.cm()
    res_test.cm()
    res_test.weightedTS()

def test_result3(res_test): #check export_at_once
    res_test.export_at_once()

def test_result4(res_test,name_test): #check both orifinal file name and savefilename exist
    res_test.export(savefilename=name_test)

    
if __name__ == '__main__':
    pytest.main()