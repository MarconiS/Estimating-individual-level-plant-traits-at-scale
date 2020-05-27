#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Apr 29 17:32:00 2020

@author: sergiomarconi
"""


# coding: utf-8
def generate_geom_kernel(solar_az,solar_zn,sensor_az,sensor_zn,li):
    '''Calculate the Li geometric scattering kernel. 
    Parameters
    ----------
    hyObj : HyTools file object with observables loaded.       
    li :        str 
                Geometric scattering kernel type [dense,sparse]
    Returns
    -------
    Volume and geomtric scattering kernels as m x n numpy array
    '''
    relative_az = sensor_az - solar_az 
    #Li kernels
    ############
    # Constants used from Colgan  et al. RS 2012 
    #
    # Eq. 37,52. Wanner et al. JGRA 1995
    solar_zn_ = np.arctan(10* np.tan(solar_zn))
    sensor_zn_ = np.arctan(10* np.tan(sensor_zn))
    # Eq 50. Wanner et al. JGRA 1995
    D = np.sqrt((np.tan(solar_zn_)**2) + (np.tan(sensor_zn_)**2) - 2*np.tan(solar_zn_)*np.tan(sensor_zn_)*np.cos(relative_az))    
    # Eq 49. Wanner et al. JGRA 1995
    t_num = 2. * np.sqrt(D**2 + (np.tan(solar_zn_)*np.tan(sensor_zn_)*np.sin(relative_az))**2) 
    t_denom = (1/np.cos(solar_zn_))  + (1/np.cos(sensor_zn_))
    t = np.arccos(np.clip(t_num/t_denom,-1,1))
    # Eq 33,48. Wanner et al. JGRA 1995
    O = (1/pi) * (t - np.sin(t)*np.cos(t)) * t_denom
    # Eq 51. Wanner et al. JGRA 1995
    cosPhase_ =  np.cos(solar_zn_)*np.cos(sensor_zn_) + np.sin(solar_zn_)*np.sin(sensor_zn_)*np.cos(relative_az)
#
    if li == 'sparse':
        # Eq 32. Wanner et al. JGRA 1995
        k_geom = O - (1/np.cos(solar_zn_)) - (1/np.cos(sensor_zn_)) + .5*(1+ cosPhase_) * (1/np.cos(sensor_zn_))
    elif li == 'dense':
        # Eq 47. Wanner et al. JGRA 1995
        k_geom = (((1+cosPhase_) * (1/np.cos(sensor_zn_)))/ (t_denom - O)) - 2
    #
    return k_geom


def generate_volume_kernel(solar_az,solar_zn,sensor_az,sensor_zn, ross):
    '''Calculate the Ross vlumetric scattering kernel. 
    Parameters
    ----------
    solar_az:   float or np.array
                Solar zenith angle in radians
    solar_zn:   float or np.array 
                Solar zenith angle in radians
    sensor_az:  np.array
                Sensor view azimuth angle in radians
    sensor_zn:  np.array
                Sensor view zenith angle in radians          
    ross:       str 
                Volume scattering kernel type [thick,thin]
    Returns
    -------
    Volume scattering kernel as m x n numpy array
    ''' 
    relative_az = sensor_az - solar_az 
    #Ross kernels 
    ############
    # Eq 2. Schlapfer et al. IEEE-TGARS 2015
    phase = np.arccos(np.cos(solar_zn)*np.cos(sensor_zn) + np.sin(solar_zn)*np.sin(sensor_zn)*  np.cos(relative_az))  
    if ross == 'thick':
        # Eq 13. Wanner et al. JGRA 1995
        k_vol = ((pi/2 - phase)*np.cos(phase) + np.sin(phase))/ (np.cos(sensor_zn)*np.cos(solar_zn)) - pi/4
    elif ross == 'thin':
        # Eq 13. Wanner et al. JGRA 1995
        k_vol = ((pi/2 - phase)*np.cos(phase) + np.sin(phase))/ (np.cos(sensor_zn)*np.cos(solar_zn)) - pi/2
    return k_vol



# In[1]:
def generate_brdf_coeffs_band(bnd, mask, k_vol,k_geom, topo_x): #
    '''Return the BRDF coefficients for the input band.
#    
    Parameters
    ----------
    band:       m x n np.array
                Image band
    mask:       m x n np.array
                Binary image mask
    k_vol:      m x n np.array
                Volume scattering kernel image
    k_geom:     m x n np.array
                Geometric scattering kernel image
    Returns
    -------
    brdf_coeff: list
                BRDF coefficients
    '''
 #   
    # Mask kernels
    k_vol = k_vol[mask]
    k_geom = k_geom[mask]
    topo_x = topo_x[mask]
    # Reshape for regression
    k_vol = np.expand_dims(k_vol,axis=1)
    k_geom = np.expand_dims(k_geom,axis=1)
    topo_x = np.expand_dims(topo_x,axis=1)
    #X = np.concatenate([k_vol,k_geom,np.ones(k_geom.shape)],axis=1)
    X = np.concatenate([k_vol,k_geom,np.ones(k_geom.shape)],axis=1)
    # Mask input band
    mask = np.squeeze(mask)
    y = bnd[mask]/10000
    # Calculate BRDF coefficients
    brdf_coeff = np.linalg.lstsq(X, y)[0].flatten()
    topo_coeff = np.linalg.lstsq(topo_x, y)[0].flatten()
 #   
    return brdf_coeff, topo_coeff
    

def tile_solar_angle(full_path):
    hdf5_file = h5py.File(full_path, 'r')
    file_attrs_string = str(list(hdf5_file.items()))
    file_attrs_string_split = file_attrs_string.split("'")
    sitename = file_attrs_string_split[-1]    
    flight_paths = hdf5_file[sitename]["Reflectance/Metadata/Ancillary_Imagery/Data_Selection_Index"].attrs["Data_Files"]
    flight_paths=str(flight_paths).split(",")
    which_paths = np.unique(hdf5_file[sitename]["Reflectance/Metadata/Ancillary_Imagery/Data_Selection_Index"].value)
    solar_angle = []
    for pt in which_paths:
        #if pt is negative, get any from the available to avoid error(the pixel is blank anyway)
        if pt < 0:
            flight = (flight_paths)[ which_paths[1]].split("_")[5]
        else:
            flight = (flight_paths)[pt].split("_")[5]
          #  
        sol_az = hdf5_file[sitename]["Reflectance/Metadata/Logs/"][str(flight)]["Solar_Azimuth_Angle"].value
        sol_zn = hdf5_file[sitename]["Reflectance/Metadata/Logs/"][str(flight)]["Solar_Zenith_Angle"].value
        solar_angle.append([pt, sol_az, sol_zn])
    return(solar_angle)
    


def print_attrs(name, obj):
    print(name)
    


def h5refl2array(full_path, epsg):
    #refl, refl_md, wavelengths, sol_az, sol_zn, sns_az, sns_zn, slope, aspect = h5refl2array(full_path, epsg = epsg)
    hdf5_file = h5py.File(full_path, 'r')
    file_attrs_string = str(list(hdf5_file.items()))
    file_attrs_string_split = file_attrs_string.split("'")
    sitename = file_attrs_string_split[1]
    epsg = hdf5_file[sitename]["Reflectance/Metadata/Coordinate_System/EPSG Code"].value
    solar_angles = tile_solar_angle(full_path)
    #Extract the reflectance & wavelength datasets
    reflArray = hdf5_file[sitename]['Reflectance']
    refl =reflArray['Reflectance_Data'].value
    wavelengths = reflArray['Metadata']['Spectral_Data']['Wavelength'].value
    # Create dictionary containing relevant metadata information
    refl_md = {}
    refl_md['mapInfo'] = reflArray['Metadata']['Coordinate_System']['Map_Info'].value
    refl_md['wavelength'] = reflArray['Metadata']['Spectral_Data']['Wavelength'].value
    refl_md['shape'] = refl.shape
    #Extract no data value & scale factor
    refl_md['noDataVal'] = float(reflArray['Reflectance_Data'].attrs['Data_Ignore_Value'])
    refl_md['scaleFactor'] = float(reflArray['Reflectance_Data'].attrs['Scale_Factor'])
    #metadata['interleave'] = reflData.attrs['Interleave']
    refl_md['bad_band_window1'] = np.array([1340, 1445])
    refl_md['bad_band_window2'] = np.array([1790, 1955])
    refl_md['epsg'] = str(epsg).split("'")[1]
#    
    #get tiles for BRDF correction
    sns_az = hdf5_file[sitename]['Reflectance/Metadata/to-sensor_azimuth_angle']
    sns_zn = hdf5_file[sitename]['Reflectance/Metadata/to-sensor_zenith_angle']
    slope = hdf5_file[sitename]['Reflectance/Metadata/Ancillary_Imagery/Slope']
    aspect = hdf5_file[sitename]['Reflectance/Metadata/Ancillary_Imagery/Aspect']
#    
    #get solar angles as array to leverage flightpaths mosaic
    flightpaths = hdf5_file[sitename]['Reflectance/Metadata/Ancillary_Imagery/Data_Selection_Index'].value
    sol_zn = hdf5_file[sitename]['Reflectance/Metadata/Ancillary_Imagery/Data_Selection_Index'].value
    sol_az = hdf5_file[sitename]['Reflectance/Metadata/Ancillary_Imagery/Data_Selection_Index'].value
    for pt in range(len(solar_angles)):
        sol_az[flightpaths==solar_angles[pt][0]] = solar_angles[pt][1]
        sol_zn[flightpaths==solar_angles[pt][0]] = solar_angles[pt][2]
#
    mapInfo_string = str(refl_md['mapInfo']);
    mapInfo_split = mapInfo_string.split(",")
    mapInfo_split
#
    # Extract the resolution & convert to floating decimal number
    refl_md['res'] = {}
    refl_md['res']['pixelWidth'] = float(mapInfo_split[5])
    refl_md['res']['pixelHeight'] = float(mapInfo_split[6])
    # Extract the upper left-hand corner coordinates from mapInfo
    xMin = float(mapInfo_split[3])  # convert from string to floating point number
    yMax = float(mapInfo_split[4])
#
    # Calculate the xMax and yMin values from the dimensions
    xMax = xMin + (refl_md['shape'][1] * refl_md['res']['pixelWidth'])  # xMax = left edge + (# of columns * resolution)",
    yMin = yMax - (refl_md['shape'][0] * refl_md['res']['pixelHeight'])  # yMin = top edge - (# of rows * resolution)",
    refl_md['extent'] = (xMin, xMax, yMin, yMax)  # useful format for plotting
    refl_md['ext_dict'] = {}
    refl_md['ext_dict']['xMin'] = xMin
    refl_md['ext_dict']['xMax'] = xMax
    refl_md['ext_dict']['yMin'] = yMin
    refl_md['ext_dict']['yMax'] = yMax
    hdf5_file.close
#
    return refl, refl_md, wavelengths, sol_az, sol_zn, sns_az, sns_zn, slope, aspect


def stack_subset_bands(reflArray, reflArray_metadata, bands, clipIndex):
    subArray_rows = clipIndex['yMax'] - clipIndex['yMin']
    subArray_cols = clipIndex['xMax'] - clipIndex['xMin']
#
    stackedArray = np.zeros((subArray_rows, subArray_cols, len(bands)), dtype=np.int16)
    band_clean_dict = {}
    band_clean_names = []
#
    for i in range(len(bands)):
        band_clean_names.append("b" + str(bands[i]) + "_refl_clean")
        band_clean_dict[band_clean_names[i]] = subset_clean_band(reflArray, reflArray_metadata, clipIndex, bands[i])
        stackedArray[..., i] = band_clean_dict[band_clean_names[i]]
#
    return stackedArray



def subset_clean_band(reflArray, reflArray_metadata, clipIndex, bandIndex):
    bandCleaned = reflArray[clipIndex['yMin']:clipIndex['yMax'], clipIndex['xMin']:clipIndex['xMax'],
                  bandIndex - 1].astype(np.int16)
#
    return bandCleaned



def array2raster(newRaster, reflBandArray, reflArray_metadata, extent, ras_dir, epsg):
    NP2GDAL_CONVERSION = {
        "uint8": 1,
        "int8": 1,
        "uint16": 2,
        "int16": 3,
        "uint32": 4,
        "int32": 5,
        "float32": 6,
        "float64": 7,
        "complex64": 10,
        "complex128": 11,
    }
#
    pwd = os.getcwd()
    os.chdir(ras_dir)
    cols = reflBandArray.shape[1]
    rows = reflBandArray.shape[0]
    bands = reflBandArray.shape[2]
    pixelWidth = float(reflArray_metadata['res']['pixelWidth'])
    pixelHeight = -float(reflArray_metadata['res']['pixelHeight'])
    originX = extent['xMin']
    originY = extent['yMax']
#
    driver = gdal.GetDriverByName('GTiff')
    gdaltype = NP2GDAL_CONVERSION[reflBandArray.dtype.name]
    outRaster = driver.Create(newRaster, cols, rows, bands, gdaltype)
    outRaster.SetGeoTransform((originX, pixelWidth, 0, originY, 0, pixelHeight))
    # outband = outRaster.GetRasterBand(1)
    # outband.WriteArray(reflBandArray[:,:,x])
    for band in range(bands):
        outRaster.GetRasterBand(band + 1).WriteArray(reflBandArray[:, :, band])
#
    outRasterSRS = osr.SpatialReference()
    #outRasterSRS.ImportFromEPSG(reflArray_metadata['epsg'])
    #outRasterSRS.ExportToWkt()
    outRasterSRS.ImportFromEPSG(epsg)
    outRaster.SetProjection(outRasterSRS.ExportToWkt())
    outRaster.FlushCache()
    os.chdir(pwd)


def calc_clip_index(clipExtent, h5Extent, xscale=1, yscale=1):
    h5rows = h5Extent['yMax'] - h5Extent['yMin']
    h5cols = h5Extent['xMax'] - h5Extent['xMin']
#
    ind_ext = {}
    ind_ext['xMin'] = round((clipExtent['xMin'] - h5Extent['xMin']) / xscale)
    ind_ext['xMax'] = round((clipExtent['xMax'] - h5Extent['xMin']) / xscale)
    ind_ext['yMax'] = round(h5rows - (clipExtent['yMin'] - h5Extent['yMin']) / yscale)
    ind_ext['yMin'] = round(h5rows - (clipExtent['yMax'] - h5Extent['yMin']) / yscale)
#
    return ind_ext

def extract_hsi(full_path, itc_id, itc_xmin, itc_xmax, itc_ymin, itc_ymax, epsg, ras_dir = './outdir/plots/hsi/'):
#
    #print(itc_id, itc_xmin, itc_xmax, itc_ymin, itc_ymax, epsg)
    #extract array in h5
    refl_md, refl,wavelengths, solar_angles = h5refl2array(full_path, epsg = epsg)
#    
    #delete water bands
    rgb = np.r_[0:425]
    rgb = np.delete(rgb, np.r_[419:425])
    rgb = np.delete(rgb, np.r_[283:315])
    rgb = np.delete(rgb, np.r_[192:210])
    xmin, xmax, ymin, ymax = refl_md['extent']
    #print(xmin, xmax, ymin, ymax)
    #get extent 
    clipExtent = {}
    clipExtent['xMin'] = itc_xmin
    clipExtent['yMin'] = itc_ymin
    clipExtent['yMax'] = itc_ymax
    clipExtent['xMax'] = itc_xmax
    #print(clipExtent)
    #and then define which cell arrays they belong to
    subInd = calc_clip_index(clipExtent, refl_md['ext_dict'])
    subInd['xMax'] = int(subInd['xMax'])
    subInd['xMin'] = int(subInd['xMin'])
    subInd['yMax'] = int(subInd['yMax'])
    subInd['yMin'] = int(subInd['yMin'])
    #print(subInd)
#    
#    
    refl = refl[(subInd['yMin']):subInd['yMax'], (subInd['xMin']):subInd['xMax'], :]
    refl.shape
    #print(refl.shape)
 #   
    #initialize new raster
    subArray_rows = subInd['yMax'] - subInd['yMin']
    subArray_cols = subInd['xMax'] - subInd['xMin']
    hcp = np.zeros((subArray_rows, subArray_cols, len(rgb)), dtype=np.int16)
 #   
    #load info in multi-layer array
    band_clean_dict = {}
    band_clean_names = []
    for i in range(len(rgb)):
        band_clean_names.append("b" + str(rgb[i]) + "_refl_clean")
        band_clean_dict[band_clean_names[i]] = refl[:, :, rgb[i]].astype(np.int16)
        hcp[..., i] = band_clean_dict[band_clean_names[i]]
#    
    sub_meta = refl_md
    ii = str(itc_id) + '.tif'
    #save array to raster
    array2raster(ii, hcp, sub_meta, clipExtent, ras_dir, int(epsg))
 #   
    return 0
  


def create_tiff_datasets(full_path, wd, epsg, ross, li):
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
    import warnings
    warnings.filterwarnings("ignore")
    try:
        refl, refl_md, wavelengths, sol_az, sol_zn, sns_az, sns_zn, slope, aspect = h5refl2array(full_path, epsg = epsg)
        rgb = np.r_[0:425]
        rgb = np.delete(rgb, np.r_[419:426])
        rgb = np.delete(rgb, np.r_[283:315])
        rgb = np.delete(rgb, np.r_[192:210])
        xmin, xmax, ymin, ymax = refl_md['extent']
     #   
        itc_xmin = xmin
        itc_ymin = ymin
        itc_ymax = ymax
        itc_xmax = xmax
    #
        clipExtent = {}
        clipExtent['xMin'] = itc_xmin
        clipExtent['yMin'] = itc_ymin
        clipExtent['yMax'] = itc_ymax
        clipExtent['xMax'] = itc_xmax
        print(clipExtent)
        subInd = calc_clip_index(clipExtent, refl_md['ext_dict'])
        subInd['xMax'] = int(subInd['xMax'])
        subInd['xMin'] = int(subInd['xMin'])
        subInd['yMax'] = int(subInd['yMax'])
        subInd['yMin'] = int(subInd['yMin'])
        refl = refl[:,(subInd['yMin']):subInd['yMax'], (subInd['xMin']):subInd['xMax']]
        sns_az = sns_az[(subInd['yMin']):subInd['yMax'], (subInd['xMin']):subInd['xMax']]
        sns_zn = sns_zn[(subInd['yMin']):subInd['yMax'], (subInd['xMin']):subInd['xMax']]
        slope = slope[(subInd['yMin']):subInd['yMax'], (subInd['xMin']):subInd['xMax']]
        aspect = aspect[(subInd['yMin']):subInd['yMax'], (subInd['xMin']):subInd['xMax']]
        # mask away bad pixels
        ndvi = (refl[:, :,90] - refl[:,:,58])/(refl[:, :,58] +refl[:, :,90]) > 0.25
        nir860 = (refl[:, :,96] + refl[:, :,97])/20000 > 0.1
        mask = (sns_zn < 10000) * (sns_zn > 0) * (aspect < 10000) * (aspect > 0) * (slope < 10000) * (slope > 0) * (sns_az < 10000) * (sns_az > 0) * ndvi * nir860
        #
        # convert degrees in radiants
        slope = (slope * pi) / 180
        aspect  = (aspect * pi) / 180
        sns_az  = (sns_az * pi) / 180
        sns_zn  = (sns_zn * pi) / 180
        sol_az = (sol_az[0] * pi) / 180
        sol_zn = (sol_zn[0] * pi) / 180
        # topographic correction
        rel_az = aspect - sol_az
        cos_i = np.cos(sol_zn) * np.cos(slope) + np.sin(sol_zn) * np.sin(slope) * np.cos(rel_az)
        c1 =  np.cos(sol_zn) * np.cos(slope)
        # Generate scattering kernels
        k_vol = generate_volume_kernel(sol_az, sol_zn, sns_az,sns_zn, ross = ross)
        k_geom = generate_geom_kernel(sol_az, sol_zn, sns_az,sns_zn,li = li)
        # Generate scattering kernels at Nadir
        k_vol_nadir = generate_volume_kernel(sol_az, sol_zn, sns_az,0, ross = ross)
        k_geom_nadir = generate_geom_kernel(sol_az, sol_zn, sns_az,0,li = li)
        #
        subArray_rows = subInd['yMax'] - subInd['yMin']
        subArray_cols = subInd['xMax'] - subInd['xMin']
        #hcp = np.zeros((subArray_rows, subArray_cols, len(rgb)), dtype=np.int16)
        brd = np.zeros((subArray_rows, subArray_cols, len(rgb)), dtype=np.float32)
        band_clean_dict = {}
        band_clean_names = []
        brdf_coeffs = []
        topo_coeffs = [] 
        #mask = np.squeeze(mask)
        for i in range(len(rgb)):
            band_clean_names.append("b" + str(rgb[i]) + "_refl_clean")
            bnd = refl[:, :, rgb[i]].astype(np.int16)
            band_clean_dict[band_clean_names[i]] = bnd
            #hcp[..., i] = band_clean_dict[band_clean_names[i]]
            #calculate brfd coefficients
            ith_brdf, ith_topo = generate_brdf_coeffs_band(bnd,mask,k_vol,k_geom, cos_i)
            brdf_coeffs.append(ith_brdf)
            topo_coeffs.append(ith_topo)
            #apply brdf and topographic correction
            brdf = ith_brdf[0] * k_vol + ith_brdf[1] * k_geom  + ith_brdf[2]
            brdf_nd = ith_brdf[0] * k_vol_nadir + ith_brdf[1] * k_geom_nadir  + ith_brdf[2]
            bdrf_cor = brdf_nd/brdf
            topo_cor = (c1 + ith_topo) / (cos_i + ith_topo) 
            bnd = bnd/10000
            bnd[~np.squeeze(mask)] = np.nan
            brd[..., i] = bnd * np.squeeze(bdrf_cor*topo_cor)
        #
        #
        brdf_df =  pd.DataFrame(brdf_coeffs,columns=['k_vol','k_geom','k_iso'])
        brdf_df.to_csv(wd+"/test_brdf.csv")
        #save hcp into a tiff file [reflectance]
        sub_meta = refl_md
        #wd = "/orange/ewhite/s.marconi/Chapter1/2015_Campaign/D03/OSBS/L4/"
        itc_id = str(int(itc_xmin/1000)*1000) + "_" + str(int(itc_ymin/1000)*1000)
        ii = itc_id + "_" + '.tif'
        ras_dir = wd+"/HSI/"
        #array2raster(ii, hcp, sub_meta, clipExtent, ras_dir)
        #
        ras_dir = wd+"/corrHSI"
        array2raster(ii, brd, sub_meta, clipExtent, ras_dir = ras_dir, epsg = int(refl_md['epsg']))
    except:
        print("ATTENTION!!")
        print("tile "+str(full_path))
        print("end exception")
