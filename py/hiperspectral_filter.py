#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Apr 29 17:27:39 2020

@author: sergiomarconi
"""

import h5py
#hdf5_file = h5py.File(refl_filename, 'r')

from pathlib import Path
import numpy as np
import h5py
import gdal, osr
import matplotlib.pyplot as plt
import sys
import ogr, os
import math
import pandas as pd
from math import pi
import random
import string

epsg = 32616
ross = "thick" #thin
li = "sparse"
wd = "/orange/ewhite/s.marconi/Chapter1/2015_Campaign/D08/TALL/L4"
wd = "//orange/idtrees-collab/hsi_brdf_corrected/"
h5_pt = "/orange/ewhite/s.marconi/Chapter1/2015_Campaign/D08/TALL/L1"
h5_pt = "/ufrc/ewhite/s.marconi/Chapter3/neonVegWrangleR/outdir/DP3.30006.001/DP3.30006.001/"
h5_pt = "/orange/ewhite/NeonData/"

result = list(Path(h5_pt).rglob("*.[h][5]"))
result = pd.Series(result).astype('str') 
pattern = pd.read_csv("//ufrc/ewhite/s.marconi/Chapter3/neonVegWrangleR/indir/all_tiles_to_get.csv")
corrected = pd.read_csv("//orange/ewhite/s.marconi/tiles_processed.csv")

tiles = result.str.contains('|'.join(pattern.final_tiles))
result = result[tiles]
tiles = result.str.contains('|'.join(corrected.tile))
result = result[~tiles]


#f = "/Users/sergiomarconi/Downloads/NEON_spectrometer-orthorectified-surface-directional-reflectance---mosaic-2/NEON_D01_HARV_DP3_732000_4713000_reflectance.h5"
epsg = 32
ross = "thick"
li = "dense"
#wd = "/Users/sergiomarconi/Documents/Data/BRDF"
#create_tiff_datasets(f,  wd = wd, epsg = epsg, ross = ross, li = li)
from joblib import Parallel, delayed
Parallel(n_jobs=11)(delayed(create_tiff_datasets)(full_path = f,  wd = wd, epsg = epsg, ross = ross, li = li) for f in result)
