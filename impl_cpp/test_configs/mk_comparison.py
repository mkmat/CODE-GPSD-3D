#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Nov  2 15:13:53 2023

@author: samarth
"""

import pandas as pd

df_sa = pd.read_csv('lpes_vals.csv', names=['r_max_sa'])
df_mk = pd.read_csv('lpes_vals-mk.csv', names=['r_max_mk'])

df_sa['r_max_mk'] = df_mk['r_max_mk']
df_sa['diff'] = df_sa['r_max_sa'] - df_sa['r_max_mk']