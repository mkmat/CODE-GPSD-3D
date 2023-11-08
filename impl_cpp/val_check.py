#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Nov  8 17:26:13 2023

@author: samarth
"""

import pandas as pd

df_1 = pd.read_csv('./r_max_details.csv')
df_2 = pd.read_csv('./reference/r_max_details.csv')

df = df_1.join(df_2, on='id', how='inner', lsuffix = '_l', rsuffix = '_r')

df['diff'] = df['lpes_l'] - df['lpes_r']