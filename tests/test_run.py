# -*- coding: utf-8 -*-
"""
Created on Mon Jul 30 22:44:40 2018

to test whether olsa correctly runs

@author: tadahaya
"""
import pytest
import sys
sys.path.append("C:\\Scripts\\olsa_package\\main") #to make a path for the modules
from olsapy import DataClass,olsa,usspca

### fixtures ###
@pytest.fixture
def dat_test2():
    dat = DataClass()
    dat.load("tests\\data\\test_data.txt",read_tsv=True)
    return dat

### tests ###
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

    
if __name__ == '__main__':
    pytest.main()