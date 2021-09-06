# -*- coding: utf-8 -*-
"""
Created on Thu Jun  6 17:15:04 2019

@author: Boda Li
"""

import pandas as pd

O2data=pd.read_csv('oxygen.csv');
print(O2data.describe())