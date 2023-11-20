#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Nov 20 11:54:19 2023

@author: samarth
"""

import pandas as pd

df = pd.read_csv('config10.txt', sep=' ', names = ['id', 'x', 'y', 'z'])

df[['x','y','z']].to_csv('config10_revised.txt', sep=' ', index=False, header=False)