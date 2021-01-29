'''
Reads temporal information from a netcdf and puts it into a pcraster
field in memory. In this script we extract the monthly land 
cover fraction per land-use type for the selected year (2001).
'''

# standard modules, imported in full
import os, sys, shutil, datetime, string

# import packages, numpy and pcraster
import numpy as np
import pcraster as pcr
from pcraster.framework import generateNameT
import os

# basic functions to process input mainly
from spatialDataSet2PCR import spatialAttributes, spatialDataSet, setClone
from ncRecipes_fixed import getNCDates
from read_temporal_info_to_pcr import create_date_list, getTimedPCRData
from zonalStatistics import zonal_statistics_pcr


workdir = os.getcwd()
mapdir = os.path.join(workdir, 'mapinput')
clonedir = os.path.join(workdir, 'amazon_fine')
outputdir = os.path.join(workdir, 'amazon_input')
nc_filename = os.path.join(mapdir, 'MODIS-C006_MCD12C1_landcover__LPDAAC__0.05deg_fv0.02.nc4')

# main
def main():

    # initialization
    # set the mask of the Amazon
    amazonmask_filename = os.path.join(mapdir, 'amazon_5min_mask.map')
    
    print ("Getting clone map attrs")
    # set the clone and get its attributes
    clone_attributes = spatialAttributes(amazonmask_filename)
    print ("Setting clone")
    setClone(clone_attributes)
    
    # and read in the map
    print ("Read the map")
    amazonmask = pcr.readmap(amazonmask_filename)
    
    id_list = np.unique(pcr.pcr2numpy(amazonmask, 0))
    id_list = (id_list[id_list != 0]).tolist()
    
    #
    vnames = ["landcover_igbp", "confidence_igbp",\
    "water_igbp", "evergreen_needleleaf_forest_igbp", "evergreen_broadleaf_forest_igbp",\
    "deciduous_needleleaf_forest_igbp","deciduous_broadleaf_forest_igbp",\
    "mixed_forest_igbp", "closed_shrublands_igbp","open_shrublands_igbp", "woody_savannas_igbp", "savannas_igbp",\
    "grasslands_igbp", "permanent_wetlands_igbp", "croplands_igbp", "urban_and_builtup_igbp",\
    "cropland_natural_vegetation_mosaic_igbp", "snowandice_igbp","barren_sparsely_vegetated_igbp"]
    for vn in vnames: 
        nc_dataset = 'NETCDF:"%s":%s' % (nc_filename, vn)
        nc_date_list = getNCDates(nc_filename)
        print (nc_date_list)
        # and create a list of dates we want to extract
        selected_dates, message_str = create_date_list(datetime.datetime(2001, 1, 1),\
            time_increment= 'monthly')
        print ("MSG STR: " + message_str)
        selected_date = selected_dates[0]
        print("Date: " + str(selected_date))
        pcrfield, returned_date, message_str = getTimedPCRData(nc_dataset,\
            nc_date_list, selected_date, dateSelectionMethod= 'interpolation',\
            dataAttributes= clone_attributes)
        print ('   for the date %s the following is read from the actual data:\n%s' % \
            (selected_date, message_str))
        pcrfield = pcr.ifthen(amazonmask, pcrfield)
        fn = vn.split("_igbp")[0] + ".map"
        pcr.report(pcrfield, os.path.join(outputdir, fn))
    #rof
 
if __name__ == "__main__":
    main()
    sys.exit('all done')
