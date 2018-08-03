# -*- coding: utf-8 -*-
"""

Created on Mon Jul 30 22:44:40 2018

to test whether result objects correctly work

@author: tadahaya
"""
import pytest
import sys
sys.path.append("C:\\Scripts\\olsa_package\\main") #to make a path for the modules
from olsapy import DataClass,olsa

@pytest.fixture
def res_test():
    dat = DataClass()
    dat.load("tests\\data\\test_data.txt",read_tsv=True)
    res = olsa(dat)
    return res

@pytest.fixture
def name_test():
    return "tests\\data\\test_data_res.txt"


### tests ###
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