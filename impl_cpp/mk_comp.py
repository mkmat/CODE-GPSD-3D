#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Nov  3 19:11:24 2023

@author: samarth
"""

import pandas as pd
import numpy as np

def get_probe_centre_distance(df):
    
    L = 10
    
    df['x_diff'] = df['px'] - df['cx']
    df['y_diff'] = df['py'] - df['cy']
    df['z_diff'] = df['pz'] - df['cz']
    
    df['x_diff'] = df['x_diff'] - L*round(df['x_diff']/L)
    df['y_diff'] = df['y_diff'] - L*round(df['y_diff']/L)
    df['z_diff'] = df['z_diff'] - L*round(df['y_diff']/L)
    
    df['R2'] = np.sqrt((df['x_diff']**2) + (df['x_diff']**2) + (df['x_diff']**2))
    
    df = df.drop(columns = ['x_diff', 'y_diff', 'z_diff'])
    
    return df

def get_min_distance_from_particle(df, row):
    
    L = 10
    
    df['x_diff'] = df['x'] - row['cx']
    df['y_diff'] = df['y'] - row['cy']
    df['z_diff'] = df['z'] - row['cz']
    
    df['x_diff'] = df['x_diff'] - L*round(df['x_diff']/L)
    df['y_diff'] = df['y_diff'] - L*round(df['y_diff']/L)
    df['z_diff'] = df['z_diff'] - L*round(df['y_diff']/L)
    
    df['R2']         = np.sqrt((df['x_diff']**2) + (df['y_diff']**2) + (df['z_diff']**2)) - 1

    
    return df['R2'].min(), df['R2'].idxmin()

    
    

L = 10
coords_df = pd.read_csv('./test_configs/off_lattice_coords.csv', names=['x','y','z'])

df_sa = pd.read_csv('r_max_details.csv')
df_mk = pd.read_csv('r_max_details-mk.csv', sep=' ', index_col=False)
df    = df_sa.join(df_mk, on='id', how='inner', lsuffix='_sa', rsuffix='_mk')


df['diff'] = df['lpes_sa'] - df['lpes_mk']
test_df    = df[abs(df['diff']) > 0.01]
id_df      = test_df[['id']] 

df_sa = df_sa.merge(id_df, on='id', how='inner')
df_mk = df_mk.merge(id_df, on='id', how='inner')

df_sa = get_probe_centre_distance(df_sa)
df_mk = get_probe_centre_distance(df_mk)

size = len(df_sa)

for i in range(size):
    
    r_min, idx_min = get_min_distance_from_particle(coords_df, df_sa.iloc[i])
    df_sa.loc[i, 'R1']     = r_min
    df_sa.loc[i, 'R1_idx'] = idx_min
    
for i in range(size):
    
    r_min, idx_min = get_min_distance_from_particle(coords_df, df_mk.iloc[i])
    df_mk.loc[i, 'R1']     = r_min
    df_mk.loc[i, 'R1_idx'] = idx_min
    
#df_sa.to_csv('sa_code_anomalies_stats.csv', index=False)
#df_mk.to_csv('mk_code_anomalies_stats.csv', index=False)