#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Nov  1 02:31:12 2023

@author: samarth
"""

import pandas as pd
import numpy as np

"""df = pd.read_csv('coords.csv', names=['x','y','z'])

L = 10

x = 7.25927
y = 5.26002
z = 3.81528

df['x_diff'] = df['x'] - x
df['y_diff'] = df['y'] - y
df['z_diff'] = df['z'] - z

df['x_diff'] = df['x_diff'] - L * round(df['x_diff']/L)
df['y_diff'] = df['y_diff'] - L * round(df['y_diff']/L)
df['z_diff'] = df['z_diff'] - L * round(df['z_diff']/L)

df['r'] = df['x_diff']**2 + df['y_diff']**2 + df['z_diff']**2

df = df[df['r'] < (1.77**2)]"""

df_sa = pd.read_csv('../check_2.txt', names=['r_sa'])
df_mk = pd.read_csv('results_mk.csv', sep=' ', names=['px', 'py', 'pz', 'cx', 'cy', 'cz', 'r_mk'])


df_sa = df_sa.iloc[:180]
df_sa = df_sa.join(df_mk)

df_sa['diff'] = df_sa['r_sa'] - df_sa['r_mk']