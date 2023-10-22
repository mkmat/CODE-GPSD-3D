#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Oct 22 16:10:09 2023

@author: samarth
"""

import pandas as pd

df = pd.read_csv('0.csv', skiprows=19)
df_c = df[['x0', 'x1', 'x2']]

df_c.to_csv('coords.csv', index=False, header=False)