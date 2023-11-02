#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Nov  2 13:15:47 2023

@author: samarth
"""

import random
import math
from itertools import product
import pandas as pd

def generate_non_overlapping_spheres(num_spheres, radius, side_length):
    spheres = []
    cell_size = 2 * radius

    def get_closest_image(coord):
        return coord % side_length

    def check_overlap(new_sphere, existing_spheres):
        for sphere in existing_spheres:
            for dx, dy, dz in product([-side_length, 0, side_length], repeat=3):
                distance = math.sqrt((new_sphere[0] - get_closest_image(sphere[0] + dx)) ** 2 +
                                     (new_sphere[1] - get_closest_image(sphere[1] + dy)) ** 2 +
                                     (new_sphere[2] - get_closest_image(sphere[2] + dz)) ** 2)
                if distance < 2 * radius:  # Check if the spheres overlap
                    return True
        return False

    while len(spheres) < num_spheres:
        new_sphere = [random.uniform(0, side_length), random.uniform(0, side_length), random.uniform(0, side_length)]

        if not check_overlap(new_sphere, spheres):
            spheres.append(new_sphere)

    return spheres

# Define the parameters
number_of_spheres = 200
sphere_radius = 0.5
cubic_side_length = 10.0

# Generate non-overlapping spheres
result = generate_non_overlapping_spheres(number_of_spheres, sphere_radius, cubic_side_length)

df = pd.DataFrame(result)
df.to_csv('off_lattice_coords.csv', index=False, header=False)