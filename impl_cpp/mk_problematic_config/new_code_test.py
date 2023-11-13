#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Nov 13 17:43:11 2023

@author: samarth
"""

import pandas as pd

df_1 = pd.read_csv("../results.gpsd")
df_2 = pd.read_csv("../latest.csv")

df = df_1.join(df_2, on='id', how='inner', lsuffix='_o', rsuffix='_n')
df['diff'] = df['r_o'] - df['r_n']