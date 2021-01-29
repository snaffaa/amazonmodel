'''
Extract the model result for the selected stations from the resulted
netcdf file of Amazon Sediment Transport Model run.
'''

import os, sys 
import netCDF4 as nc
from netCDF4 import Dataset 
import numpy as np
import pandas as pd 
import csv

from argparse import ArgumentParser

# find the element in the given array that has value closest to the supplied value.
def find_nearest(array, value):
    array = np.asarray(array)
    idx = (np.abs(array - value)).argmin()
    return idx, array[idx]
#def

# read the netCDF file using Dataset command

def extract_data_from_netCDF(inputfile, outputfile, dataname, locations):
    rootgrp = Dataset(inputfile, format="NETCDF4")
    if args.verbose:
        print('the data of the netcdf includes', rootgrp)
        print ('the variables of the netcdf file are', rootgrp.variables.keys())
    #fi
    data_vars = rootgrp.variables[dataname][:]
    if args.verbose:
        print ('the data variables are', data_vars)
    #fi
    # inspect the variables associated with time dimension
    mt = rootgrp.variables['time']
    # read the data of the time variable 'mt'
    times = mt[:]
    dtimes = nc.num2date(mt[:], mt.units)
    if args.verbose:
        print('times data are', times)
        print('the date times are ', dtimes)
    #fi
    # the start and the stop time/date
    istart = 0
    istop = len(dtimes)
    print("the start time is %s and the stop time is %s" % (istart, istop))
    # get the geographical coordinates
    lat_ar = rootgrp.variables['latitude'][:]
    long_ar = rootgrp.variables['longitude'][:]
    if args.verbose:
        print ('the latitude data for ' + dataname + ' are', lat_ar)
        print ('the longitude data for ' + dataname + ' are ', long_ar)
    #fi
 
    data_per_station = {}
    for index, row in locations.iterrows():
        latidx, lat = find_nearest(lat_ar, row.latitude)
        lonidx, lon = find_nearest(long_ar, row.longitude)
        data_per_station[index] = data_vars[:, latidx, lonidx]
        if args.verbose:
            print(dataname, index, lat, latidx, lon, lonidx)
            print(dataname, index, data_per_station[index])
        #fi
    #rof
    
    df = pd.DataFrame(data_per_station, index = dtimes[istart:istop])
    df.plot.line()
    if args.verbose:
        print('the ' + dataname + ' per station', df)
    #fi
    ofn = outputfile + ".csv"
    print("Writing file: " + ofn)
    df.to_csv(ofn, float_format="%10.0f")
    # get the mean for all years
    grouped_by_year = df.groupby(lambda x: x.year)
    yearly_means = grouped_by_year.mean()
    ofn = outputfile + "_yearly_means.csv"
    if args.verbose:
        print("means: ", yearly_means)
    #fi
    print("Writing file: " + ofn)
    # save in csv file
    yearly_means.to_csv(ofn, float_format="%10.0f")
    # get the sum and the mean per year
    yearly_sums = grouped_by_year.sum()
    yearly_sums.loc['Means'] = yearly_sums.mean()
    ofn = outputfile + "_yearly_totals.csv"
    if args.verbose:
        print("totals: ", yearly_sums)
    #fi
    print("Writing file: " + ofn)
    #save in csv file
    yearly_sums.to_csv(ofn, float_format="%10.0f")
#fed

if __name__ == "__main__": 
    #          
    parser = ArgumentParser()
    parser.add_argument('-o', '--outputdir', help="Output directory name")
    parser.add_argument('-s', '--stations', default="/scratch/safaa/output-validation/station_info_from_netcdf_30arcmin_new_stations_updated.csv",
        help="File with locations of stations")
    parser.add_argument('-vn', '--variablename', dest='varnames', action='append', 
        help="a variable name present in the inputfile that should be extracted")
    parser.add_argument("-v", "--verbose", action="store_true", dest="verbose", default=False,
                    help="print lots of status messages to stdout")

    parser.add_argument('filename')
    args = parser.parse_args()

    if args.filename:
        inputfile = args.filename
    else:
        sys.exit("No input NETCDF4 filename was given")    
    #fi

    inputbase = os.path.basename(inputfile)
    inputfn = os.path.splitext(inputbase)[0]

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

    if args.varnames:
        datanames = args.varnames
    else:
        sys.exit("No data variables to extract were given on the command line")
    #fi

    station_inf = pd.read_csv(args.stations, index_col=0)
    # write the column and row of the dataframe 
    station_locations = station_inf.T

    for dataname in datanames:
        print ("Extracting data %s for each station location" % (dataname))
        outputfile = os.path.join(outputdir, dataname + '_from_' + inputfn)
        extract_data_from_netCDF(inputfile, outputfile, dataname, station_locations)
    #rof
    sys.exit('all done')
#fi


