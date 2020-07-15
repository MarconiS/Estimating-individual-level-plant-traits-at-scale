#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jul 13 23:45:13 2020

@author: sergiomarconi
"""



import rasterio
import geopandas as gpd
import pandas as pd
import rasterstats
import numpy as np
from rasterstats import zonal_stats
from rasterio.plot import show
import os


tr_pt = "/Users/sergiomarconi/Dropbox (UFL)/Traits_Maps/TALL/ITCs/eval/"
arr = os.listdir(tr_pt)

gpd1 = gpd.read_file("/Users/sergiomarconi/Documents/Data/TOH/OSBS_N.shp")
# we'll make a string column for the wKT geom
gpd1['WKT'] = gpd1['geometry'].apply(lambda x: str(x))
grouped_gdf = gpd1.groupby('WKT').mean().reset_index()
grouped_gdf.to_csv("/Users/sergiomarconi/Documents/Data/TOH/csv/OSBS_N.csv")

gpd1 = gpd.read_file("/Users/sergiomarconi/Documents/Data/TOH/OSBS_C.shp")
# we'll make a string column for the wKT geom
gpd1['WKT'] = gpd1['geometry'].apply(lambda x: str(x))
grouped_gdf = gpd1.groupby('WKT').mean().reset_index()
grouped_gdf.to_csv("/Users/sergiomarconi/Documents/Data/TOH/csv/OSBS_C.csv")

gpd1 = gpd.read_file("/Users/sergiomarconi/Documents/Data/TOH/OSBS_P.shp")
# we'll make a string column for the wKT geom
gpd1['WKT'] = gpd1['geometry'].apply(lambda x: str(x))
grouped_gdf = gpd1.groupby('WKT').mean().reset_index()
grouped_gdf.to_csv("/Users/sergiomarconi/Documents/Data/TOH/csv/OSBS_P.csv")

gpd1 = gpd.read_file("/Users/sergiomarconi/Documents/Data/TOH/OSBS_LMA.shp")
# we'll make a string column for the wKT geom
gpd1['WKT'] = gpd1['geometry'].apply(lambda x: str(x))
grouped_gdf = gpd1.groupby('WKT').mean().reset_index()
grouped_gdf.to_csv("/Users/sergiomarconi/Documents/Data/TOH/csv/OSBS_LMA.csv")




gpd1 = gpd.read_file("/Users/sergiomarconi/Documents/Data/TOH/TALL_N.shp")
# we'll make a string column for the wKT geom
gpd1['WKT'] = gpd1['geometry'].apply(lambda x: str(x))
grouped_gdf = gpd1.groupby('WKT').mean().reset_index()
grouped_gdf.to_csv("/Users/sergiomarconi/Documents/Data/TOH/csv/TALL_N.csv")

gpd1 = gpd.read_file("/Users/sergiomarconi/Documents/Data/TOH/TALL_C.shp")
# we'll make a string column for the wKT geom
gpd1['WKT'] = gpd1['geometry'].apply(lambda x: str(x))
grouped_gdf = gpd1.groupby('WKT').mean().reset_index()
grouped_gdf.to_csv("/Users/sergiomarconi/Documents/Data/TOH/csv/TALL_C.csv")

gpd1 = gpd.read_file("/Users/sergiomarconi/Documents/Data/TOH/TALL_P.shp")
# we'll make a string column for the wKT geom
gpd1['WKT'] = gpd1['geometry'].apply(lambda x: str(x))
grouped_gdf = gpd1.groupby('WKT').mean().reset_index()
grouped_gdf.to_csv("/Users/sergiomarconi/Documents/Data/TOH/csv/TALL_P.csv")

gpd1 = gpd.read_file("/Users/sergiomarconi/Documents/Data/TOH/TALL_LMA.shp")
# we'll make a string column for the wKT geom
gpd1['WKT'] = gpd1['geometry'].apply(lambda x: str(x))
grouped_gdf = gpd1.groupby('WKT').mean().reset_index()
grouped_gdf.to_csv("/Users/sergiomarconi/Documents/Data/TOH/csv/TALL_LMA.csv")