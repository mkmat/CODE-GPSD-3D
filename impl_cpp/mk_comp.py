#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Nov  3 19:11:24 2023

@author: samarth
"""

import pandas as pd
import numpy as np

df    = pd.read_csv('r_max_details.csv', index_col='id')

df    = df.iloc[:,1:]
col_mk = []

for col in df.columns:
        col_mk.append(col+'_mk')

df_mk = pd.read_csv('r_max_details-mk.csv', sep=' ', index_col=None, names=col_mk)
#df = df.join(df_mk)