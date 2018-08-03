# -*- coding: utf-8 -*-
"""
Created on Mon Jul 30 22:44:40 2018

to test whether data is correctly loaded

@author: tadahaya
"""
import pytest
import sys
sys.path.append("C:\\Scripts\\olsa_package\\main") #to make a path for the modules
from olsapy import DataClass,olsa
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

                                       
if __name__ == '__main__':
    pytest.main()
