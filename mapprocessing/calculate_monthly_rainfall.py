'''
this function reads temporal information from a netcdf and puts it into a
pcraster field in memory. In this function we'll extract monthly rainfall
for the selected year.
Much of the functionality is coming from the module 
read_temporal_info_to_pcr and was developed for the land cover
parameterization of PCR-GLOBWB in pcrLandCoverParameterization_Main.py in
/home/beek0120/PCRGLOBWB/LandCoverParameterization/dev
'''



# standard modules, imported in full
import os, sys, shutil, datetime, string




# import packages, numpy and pcraster
import numpy as np
import pcraster as pcr
from pcraster.framework import generateNameT
import os


from spatialDataSet2PCR import spatialAttributes, spatialDataSet, setClone
from ncRecipes_fixed import getNCDates
from read_temporal_info_to_pcr import create_date_list, getTimedPCRData
from zonalStatistics import zonal_statistics_pcr


workdir = os.getcwd()
mapdir = os.path.join(workdir, 'mapinput')
outputdir= '/home/naffa002/projects/finalmodel/sediment_transport_model/rainfall_hybam/output/1990'
sys.path.insert (0,"/home/beek0120/PCRGLOBWB/scripts/e2o/model")

from virtualOS import find_match_in_datenum_array

# main
def main():

    # initialization
    # set the mask of the Amazon
    amazonmask_filename = (os.path.join(mapdir ,'amazon_5min_mask.map'))
    # and the netCDF file of interest, here monthly normalized differentianted vegetation index
   
    variablename = 'P'
    #variablename = 'SoilLoss'
    nc_filename =  os.path.join(mapdir, 'HOP_CAI_1deg_monsum.nc')
    #nc_filename =  os.path.join(inputdir, 'SoilLoss.nc')
    nc_dataset = 'NETCDF:"%s":%s' % (nc_filename, variablename)
    # start
    # set the clone and get its attributes
    clone_attributes = spatialAttributes(amazonmask_filename)
    setClone(clone_attributes)
    
    # and read in the map
    amazonmask = pcr.readmap(amazonmask_filename)
    
    id_list = np.unique(pcr.pcr2numpy(amazonmask, 0))
    id_list = (id_list[id_list != 0]).tolist()
    
    # next, create a list of dates in the netCDF file
    nc_date_list = getNCDates(nc_filename , 'Time', units='days since 1980-01-01')
    print (nc_date_list)
    
    # and create a list of dates we want to extract
    selected_dates, message_str = create_date_list(datetime.datetime(1990, 1, 1),\
        time_increment= 'monthly')
    print (message_str)
    poscnt = 0
    pcrfield_avg = pcr.scalar(0)
    print (selected_dates)
    for selected_date in selected_dates:
        
        poscnt =  poscnt + 1
        
        print(' * retrieving info for %s' % selected_date)

        pcrfield, returned_date, message_str= \
            getTimedPCRData(nc_dataset,\
            nc_date_list, selected_date,\
            dateSelectionMethod= 'interpolation',\
            dataAttributes= clone_attributes)
    
        print ('   for the date %s the following is read from the actual data:\n%s' % \
            (selected_date, message_str))
        print (' the monthly values of rainfall are:')
        pcrfield = pcr.ifthen(amazonmask, pcr.cover(pcrfield))
        pcr.report(pcrfield, os.path.join(outputdir, generateNameT('Pre90', poscnt)))

if __name__ == "__main__":
    main()
    sys.exit('all done')
