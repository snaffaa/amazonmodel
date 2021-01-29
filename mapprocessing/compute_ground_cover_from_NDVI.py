import numpy as np
import netCDF4 as nc
import pandas as pd
from netCDF4 import Dataset 
import os, sys, datetime, calendar, csv
from time import gmtime, strftime
import netCDF4
from pandas import compat
from scipy.io import netcdf

import pcraster as pcr

"""
compute the fraction cover using NDVI based on W.J. Timmermans, 2007 
"""
workdir = os.getcwd()
mapname = os.path.join(workdir, 'amazon_fine')
inputdir= '/scratch/safaa/ground_ndvi/'
outputdir = '/scratch/safaa/ground_ndvi/output'

def compute_minmax_ndvi(data_array):
    #read the minimum and maximum values of the ndvi variable
    print('the shape of the array is', np.shape(ndvi_data))
    min_value= data_array.min()
    max_value= data_array.max()
    print("Minimum = %.2f, Maximum = %.2f" % (min_value, max_value))
    return [min_value, max_value]
    
def extrat_stuff_ndvi(rootgrp, ndvi_varname):
    # read the data of the time variable 'mt'
    mt = rootgrp.variables['time']
    times = mt[:]
    dtimes = nc.num2date(mt[:], mt.units)
    
    # get the geographical coordinates
    lat_ar = rootgrp.variables['lat'][:]
    long_ar = rootgrp.variables['lon'][:]
    return [lat_ar, long_ar]
    
def update_NetCDF(input_filename,output_filename,var_name,updated_data):
    src = nc.Dataset(input_filename)
    dst = nc.Dataset(output_filename,"w")
    # copy attributes
    for name in src.ncattrs():
        dst.setncattr(name, src.getncattr(name))
        dst.sync()
    #copy dimensions
    for name, dimension in src.dimensions.iteritems():
        dst.createDimension(name, (len(dimension) if not dimension.isunlimited else None))
        dst.sync()
    # copy all file data 
    for name, variable in src.variables.iteritems():
        print('procesing variable: ', name, np.shape(src.variables[name][:]))
        fill_value = None
        if hasattr(variable, "_FillValue"):
            fill_value = variable._FillValue
        if name==var_name:
            x = dst.createVariable(name, 'float32', variable.dimensions, fill_value=fill_value)
        else:
            x = dst.createVariable(name, variable.datatype, variable.dimensions, fill_value=fill_value)
        print('the variable attribute is', variable.ncattrs())
        dst.variables[name].setncatts({k: variable.getncattr(k) for k in variable.ncattrs() if k not in ["_FillValue"]})
        dst.variables[name][:] = src.variables[name][:]
        if name==var_name:
            print('the datatype of the ndvi array is', variable.datatype)
            print('the shape of the original ndvi array is', np.shape(src.variables[name][:]))
            print('the shape of the updated data is', np.shape(updated_data))
            dst.variables[name][:] = updated_data
            print('the shape of the updated ndvi array is', np.shape(dst.variables[name][:]))
            print('first element is (src, new, dst):', src.variables[name][0,0,0], updated_data[0,0,0], dst.variables[name][0,0,0])
        dst.sync()
    # close files
    dst.close()
    src.close()

if __name__ == "__main__":
    ndvi_filename ="ndvi"
    ndvi_varname ="ndvi"
    print ("Extracting data from netCDF: " + ndvi_filename)
    ndvi_pathname = os.path.join(inputdir, ndvi_filename + '.nc');
    rootgrp = Dataset(ndvi_pathname)
    ndvi_data = rootgrp.variables[ndvi_varname][:]
    minmax = compute_minmax_ndvi(ndvi_data)
    ndvi_min = minmax[0]
    ndvi_max = minmax[1]
    cover_data = ndvi_data.copy()
    
    # Iterate over every cell and apply the equation to it.
    with np.nditer(cover_data, flags=['multi_index'], op_flags=['readwrite']) as ndvi_it:
        while not ndvi_it.finished:
            ndvi_it[0] = (ndvi_it[0] - ndvi_min) / (ndvi_max - ndvi_min )
            ndvi_it.iternext()
        #elihw
    #htiw

    print cover_data

    output_file_path = os.path.join(outputdir,"cover_fraction_modis_2.nc")
    update_NetCDF(ndvi_pathname, output_file_path, 'ndvi', cover_data)
    sys.exit('all done')
