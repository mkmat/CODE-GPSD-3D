#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Oct 30 16:15:55 2023

@author: samarth
"""

import pandas as pd

df = pd.read_csv('code_rmax.csv', names=['id', 'rmax'])
df_2 = pd.read_csv('coords_id.csv.vol', sep = ' ', names = ['id', 'r2'])

df = df.merge(df_2, on='id', how='inner')

df['ratio'] = df['rmax']/df['r2']