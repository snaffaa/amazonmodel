import os, sys, types
import glob, subprocess
import pcraster as pcr
from pcraster.framework import generateNameT
from copy import copy
from argparse import ArgumentParser
from spatialDataSet2PCR import * 
import numpy as np
import PCRGLOBST
import pandas as pd 
import csv

# basic functions to process input mainly
from spatialDataSet2PCR import spatialAttributes, spatialDataSet, setClone
from ncRecipes_fixed import getNCDates
from read_temporal_info_to_pcr import create_date_list, getTimedPCRData
from zonalStatistics import zonal_statistics_pcr

parser = ArgumentParser()
parser.add_argument('-o', '--outputdir', help="Output directory name")
parser.add_argument("-v", "--verbose", action="store_true", dest="verbose", default=False, help="Verbose output")
parser.add_argument('filenames',nargs='+')
args = parser.parse_args()

if args.filenames:
    inputfiles = args.filenames
else:
    sys.exit("No input filename for station data was given")    
#fi

if args.outputdir:
    outputdir = args.outputdir
    if (not os.path.isdir(outputdir)):
        print("Output directory %s does not exist yet; will try to create it." % outputdir)
        try:
            os.mkdir(outputdir)
        except OSError:
            sys.exit("Error: creation of the directory %s failed" % outputdir)
        else:
            print("Successfully created the output directory %s " % outputdir)
        #yrt
    #fi
else:
    outputdir = os.getcwd()
#fi

mapdir = '/home/naffa002/projects/finalmodel/mapinput/halfdegree_maps'
ldd_filename = os.path.join(mapdir,'LDD.map')

cellarea_filename = os.path.join(mapdir, 'cellarea_30min.map')
station_ids = [3630050, 3630300, 3625310, 3625220, 3624250,3630200,3624160, \
    3624122, 3628500, 3618721, 3631210, 3627010, 3618710,3627451, \
    3618730, 3625350, 3624120, 3627110, 3625000, 3629152, 3625340,\
    3621000, 3618110, 3629000, 3622401, 3618121,3630120,3617810,\
    3627810, 3627440,3627812, 3627420, 3617800, 3627814, 3627040, \
    3627650, 3625150, 3620000, 3631100, 3623100, 3625350, 3618051,3618951,\
    3621400, 3618715, 3628100]
# station information is organized as follows:
# per station the longitude, latitude, GRDC ID and upstream area in km2 is given
station_info = { }
for inputfile in inputfiles:
    csvfile = open(inputfile, 'r')
    stations = list(csv.reader(csvfile, delimiter=';', quotechar='|'))
    for s in stations:
        try:
            sid = float(s[0])
            if sid in station_ids:
                station_info[s[3]] = (float(s[7]), float(s[6]), sid, float(s[5]))
                print (s[3] + "\t" + str(station_info[s[3]]))
            #fi
        except:
            pass
    #rof
#rof

pcr.setglobaloption("unittrue")
pcr.setglobaloption('large')

# set the clone
pcr.setclone(ldd_filename)

# read in the ldd and cell area [km2]
ldd = pcr.readmap(ldd_filename)
cellarea = 1.0e-6 * pcr.readmap(cellarea_filename)

# get the upstream area
upstream_area = pcr.accuflux(ldd, cellarea)

# compute the latitudes and longitudes
latitudes = pcr.ycoordinate(pcr.boolean(1)) 
longitudes = pcr.xcoordinate(pcr.boolean(1))

# create an empty station map
stations = pcr.ifthen(pcr.defined(ldd), pcr.nominal(0))
pcr_station_info = {}

# next iterate over the stations
for station_name, (longitude, latitude, station_id, catchment_area) in station_info.items():
    print (' * processing %d: %s with an upstream area of %.2f km2' % (station_id, station_name, catchment_area))
    
    # locate the station as a boolean map
    # give the minimum numnber to all cells that has the nearst distance from the cell of x coordinate 
    #and then give all the cells with minimum values a true boolean value
    station_location_x = (pcr.mapminimum(pcr.abs(latitudes - latitude)) == pcr.abs(latitudes - latitude))
    station_location_y = (pcr.mapminimum(pcr.abs(longitudes - longitude)) == pcr.abs(longitudes - longitude))
    station_location = station_location_x & station_location_y
    pcr.report(station_location_x, os.path.join(mapdir, 'station_locationx.map'))
    pcr.report(station_location_y, os.path.join(mapdir, 'station_locationy.map'))
    pcr.report(station_location, os.path.join(mapdir, 'station_location.map'))

    # identify a window around the station
    station_window = pcr.spread(station_location, 0, 1) <= (3 / 12.)
    
    upstream_area_deviation = pcr.abs(upstream_area - catchment_area) / catchment_area
    
    new_station_location = (upstream_area_deviation == pcr.areaminimum(upstream_area_deviation, station_window)) & station_window
    
    # locate the station in the map
    stations = pcr.ifthenelse(new_station_location, pcr.nominal(station_id), stations)
    
    print ('   the station was allocated within an error of %.3f of the original value on a cell with an upstream area of %.2f km2' % \
        ( \
        pcr.cellvalue(pcr.mapmaximum(pcr.ifthen(new_station_location, upstream_area_deviation)), 1)[0], \
        pcr.cellvalue(pcr.mapmaximum(pcr.ifthen(new_station_location, upstream_area)), 1)[0], \
        ))

    new_latitude = pcr.cellvalue(pcr.mapmaximum(pcr.ifthen(new_station_location, latitudes)), 1)[0]
    new_longitude = pcr.cellvalue(pcr.mapmaximum(pcr.ifthen(new_station_location, longitudes)), 1)[0]
    
    pcr_station_info[station_name] = {'longitude': new_longitude, 'latitude': new_latitude}

    print ('   and the updated latitude and longitude are %.4f and %.4f' % \
        ( \
        pcr_station_info[station_name]['latitude'], \
        pcr_station_info[station_name]['longitude'], \
        ))
#rof 
        
df = pd.DataFrame(pcr_station_info)
csvpath = os.path.join(outputdir, 'station_info_30arcmin_new_stations_updated.csv')
df.to_csv(csvpath)
print(df)
    
# create the subcatchments of the stations
station_subcatchment = pcr.subcatchment(ldd, stations)
pcr.report(station_subcatchment, os.path.join(outputdir, 'St_30arcmin_sub_new_stations.map'))
pcr.report(stations, os.path.join(outputdir, 'Sta_30arcmin_new_stations.map'))
#pcr.aguila(stations, station_subcatchment)

