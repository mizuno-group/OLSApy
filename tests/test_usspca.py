# -*- coding: utf-8 -*-
"""
Created on Mon Sep 23 13:05:34 2019

@author: tadahaya
"""
import unittest
import pandas as pd
import numpy as np
import os
import sys
import math

from olsapy import olsa

class SampleTest(unittest.TestCase):
    CLS_VAL = 'none'

    # called when test class initialization
    @classmethod
    def setUpClass(cls):
        if sys.flags.debug:
            print('> setUpClass method is called.')
        cls.CLS_VAL = '> setUpClass : initialized!'
        if sys.flags.debug:
            print(cls.CLS_VAL)

    # called when test class end
    @classmethod
    def tearDownClass(cls):
        if sys.flags.debug:
            print('> tearDownClass method is called.')
        cls.CLS_VAL = '> tearDownClass : released!'
        if sys.flags.debug:
            print(cls.CLS_VAL)

    # called when a test method runs
    def setUp(self):
        if sys.flags.debug:
            print(os.linesep + '> setUp method is called.')
        self.smpl = olsa.OLSA()
        num = list(range(10))
        col = ["sample" + str(v) for v in num]
        idx = ["feature" + str(v) for v in range(1000)]
        dat = []
        for v in num:
            dat.append(np.random.normal((-5 + v)*0.1,1,1000))
        dat = pd.DataFrame(dat,index=col,columns=idx).T
        self.smpl.load_df(dat)

    # called when a test method ends
    def tearDown(self):
        if sys.flags.debug:
            print(os.linesep + '> tearDown method is called.')

    def _df_checker(self,df):
        if type(df)!=pd.core.frame.DataFrame:
            return False
        elif df.shape[0]==0:
            return False
        else:
            head = df.head(1)
            judge = math.isnan(head.iat[0,0])
            return not judge

    def test_calc(self):
        # prepare test patterns
        test_patterns = [
            (0.6,True,True,True,True), # (arg1, arg2, ..., expected result)
            (0.8,True,True,True,True), # (arg1, arg2, ..., expected result)
            (0.6,False,True,True,True), # (arg1, arg2, ..., expected result)
            (0.6,True,False,True,True), # (arg1, arg2, ..., expected result)
            (0.6,True,True,False,True), # (arg1, arg2, ..., expected result)
            (0.6,True,True,True,False), # (arg1, arg2, ..., expected result)
            ]

        ### loop for sweeping all conditions
        for a1,a2,a3,a4,a5 in test_patterns:
            with self.subTest(accumulation=a1,IsSphered=a2,
                              UseMirror=a3,UseSmirnovGrubbs=a4,WardClusterSort=a5):
                self.smpl.usspca(accumulation=a1,IsSphered=a2,
                               UseMirror=a3,UseSmirnovGrubbs=a4,WardClusterSort=a5)
                self.assertTrue(self._df_checker(self.smpl.get_res()[0]))