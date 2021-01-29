# -*- coding: utf-8 -*-
"""
Created on Fri May 08 2015

Amazon Sediment Production Model.

This model developed from PCRGLOB-SET (hoch, 2014).It computes sediment input to the amazon river by using RUSLE equation . 
It employs the RUSLE approach and hence models long-term results with monthly temporal resolution.
This is as well represented by the use of long-term monthly averages.

"""

import os, sys, datetime, csv
import numpy as np
from pcraster import *
from pcraster.framework import *
import calendar

# basic support functions for reading and converting input (maps) etc
from spatialDataSet2PCR import spatialAttributes, spatialDataSet, setClone
from ncRecipes_fixed import getNCDates
from read_temporal_info_to_pcr import create_date_list, getTimedPCRData
from zonalStatistics import zonal_statistics_pcr
from convert2NetCDF import createNetCDF, data2NetCDF

# import functions that are part of our model
from amazon_factors import *

class RUSLE(DynamicModel):
    def __init__(self, clone_pathname, mapdir, txtdir, outdir, startdate):
        # Call the constructor of our parent class.
        DynamicModel.__init__(self)

        # Before do anything else, set the clone map to the map with supplied pathname
        setclone(clone_filename)
        self.clone_mask = readmap(clone_pathname)

        # Input and output directories
        self.mapdir = self.check_dir(mapdir)    # Directory with input maps
        self.txtdir = self.check_dir(txtdir)    # Directory with text files/data
        self.outputdir = self.check_dir(outdir) # Directory for all the output to disk

        # The date where our modelling starts
        self.startdate = startdate
    #fed

    def initial(self):
        print('Initial Section Started')
        self.config_model()  
        self.load_input_initial()
        self.setup_model()  
        print('Initial Section Finished')
    #fed

    def config_model(self):
        print("Start model configuration")
        # Start with some options you can change as you wish
        self.tests = True            
        # Routing choice
        self.accufractionflux = False   
        self.accufraction = 0.9         

        # Model options
        self.apply_Pfactor = False   

        # Reporting/output options
        self.yearly_reports = True   
        self.monthly_reports = True  
        self.netcdf_output = True    

        # Specify the data that we want to write to the netcdf file
        self.dataProducts = {}
        self.dataProducts['SedTrans'] = {\
           'variable': 'SedTransport','unit': 'tonnes per month' }
        self.dataProducts['DelRatio'] = {\
            'variable': 'DeliveryRatio','unit': 'non' }
        self.dataProducts['SoilLoss'] = {\
            'variable': 'SoilLoss','unit': 'tonnes per month' }
        self.netcdf_init()

        print('Configuration is:')
        if self.accufractionflux:
            print('\tRouting Function is Accufractionflux With Fraction: %d' % self.accufraction)
        else:
            print('\tRouting Function is Accuflux')
        #fi
        print('\tP factor is %s' % ("applied" if self.apply_Pfactor else "ignored"))            

       
        print("\tWriting " + ", ".join(map(lambda x: x[0], filter(lambda x: x[1], list(zip( \
            ["yearly maps", "monthly maps", "netcdf file"], \
            [self.yearly_reports, self.monthly_reports, self.netcdf_output])
        )))) + " to: %s" % self.outputdir)

        print("Finished model configuration")
    #fed

    def setup_model(self):
        # add a function to compute the interception fraction
        self.interception_fraction = compute_interception_fraction( \
            self.cover_fraction_info, self.part_interception_info)

        print("Start computing overall yearly (not dynamic) factors")
        
        self.organic_matter_content = scalar(1) 
        self.LS = compute_slope_factor(self.slope_length, self.slope_gradient)
        self.K_year = return_default_soil_erodibility( \
            self.mass_fraction_sand, self.mass_fraction_silt, self.mass_fraction_clay,\
            self.organic_matter_content)
        self.kfact_seasonal = return_soil_erodibility_seasonality(self.mass_fraction_sand,\
            self.mass_fraction_silt,self.mass_fraction_clay)
        self.R_year = compute_annual_erosivity((1.0 - self.interception_fraction) * self.rain_annual)
        mv = -999.9
        self.R_year = cover(self.R_year, mv)
        self.P_month  = compute_conservation_factor(self.cover_fraction_info, self.P_factor_info,self.selected_slope)
        self.Manning_n1 = Manning_n1(self.slope_gradient, self.Manning_n1_TBL)
        self.Manning_n4 = Manning_n4(self.cover_fraction_info, self.manning_n_info)
        print("Finished computing overall yearly (not dynamic) factors")

        # Report static (yearly) factors
        map(self.report_yearly_result, ["LS", "K_year", "R_year"])
        self.check_extent()
    #fed

    def dynamic(self):
        if self.currentTimeStep() == 1:
            print('Dynamic Section Starts')
        #fi
        self.calc_currentdate()
        self.load_input_dynamic()

        # Compute all dynamic factors
        print("Start computing monthly (dynamic) factors for timestep %d (date = %s)" % (self.currentTimeStep(), self.current_date))
        self.R_month = partition_erosivity(self.rain_month, self.rain_annual,self.R_year)
        self.melt_time_fraction = compute_melt_time_fraction(self.tmin, self.tmax)
        self.K_month = compute_kfact_monthly(self.K_year, self.kfact_seasonal, self.melt_time_fraction)
        self.C_month = compute_vegcover_factor2(self.ground_cover_fraction,self.R_month, self.R_year )
        self.SoilLoss = compute_soil_loss(self.R_month, self.K_month, self.LS, self.C_month, self.P_month, self.CellArea)
        self.Manning = compute_Manning(self.Manning_n1, self.Manning_n4, self.ndvi)
        self.HydroCff = compute_hydro_coeff(self.Qsurface, self.rain_month)
        self.DelRatio = compute_delivery_ratio(self.HydroCff, self.slope_gradient, self.Manning, self.slope_length)
        self.SedTrans = \
            compute_erosion(self.SoilLoss, self.DelRatio)
        print('Finished computing for timestep %d (date = %s)' % (self.currentTimeStep(), self.current_date))

        # Report dynamic factors
        map(self.report_monthly_result, ["R_month", "K_month", "C_month", 'P_month','SoilLoss', 'Manning', 'HydroCff','DelRatio', 'SedTrans'])
        self.netcdf_writedata()
        self.check_extent()


    def load_input_initial(self):
        print('Start reading input files')
        self.slope_length = self.rd_map('globalbcat')
        self.LDD = self.rd_map('LDD')
        self.CellArea = self.rd_map('cellarea_05min')
        self.slope_gradient = self.rd_map('slope05min_avgFrom30sec') # slope in ratio (*100 to get percentage)
        self.selected_slope = self.rd_map('sel_slope')
        mv = -999.9
        self.selected_slope= cover(self.selected_slope, mv)
        self.rain_annual = self.rd_map('Pan_90')
        # We replace any missing values in the rainfall map
        mv = -999.9
        self.rain_annual= cover(self.rain_annual, mv)

        self.mass_fraction_sand = self.rd_map('fsand_05deg')
        self.mass_fraction_silt = self.rd_map('fsilt_05deg')
        self.mass_fraction_clay = self.rd_map('fclay_05deg')

        # First read the data sets...
        self.manning_n_info = \
            self.load_landcoverdata(self.check_txtfile("Landcover_n4.txt"))
        self.P_factor_info = \
            self.load_landcoverdata(self.check_txtfile("P_factor.txt"))
        self.part_interception_info = \
            self.load_landcoverdata(self.check_txtfile("part_interception.txt"))

        mapnames = self.manning_n_info.keys()
        

        # ... then read in all the land cover maps.
        self.cover_fraction_info = self.load_landcovermaps(mapnames)
        self.Manning_n1_TBL = self.check_txtfile('ManningN1.txt')
        
        # Adjusting the input map data for our uses
        # The fraction of sand , silt and clay should be in percentage, therfore, divided by 100.
        self.mass_fraction_sand = self.mass_fraction_sand / 100
        self.mass_fraction_silt =  self.mass_fraction_silt / 100
        self.mass_fraction_clay = self.mass_fraction_clay / 100
        print('Finsihed reading input files')
    #fed

    def load_input_dynamic(self):
        print('Reading input files for timestep %d (date = %s)' % (self.currentTimeStep(), self.current_date))
        ndvi = self.rd_map('ndvi0000')
        qsurface = self.rd_map('Qsurface')
        tavg = self.rd_map('tavg0000')
        tmax = self.rd_map('tmax0000')
        tmin = self.rd_map('tmin0000')
        rain_month = self.rd_map('Pre90000')
        self.ground_cover_fraction= self.rd_map('cfr00000')
        
        
        # We replace any missing values in the rainfall map
        mv = -999.9
        self.rain_month = cover(rain_month, mv)

        # The temperature (max and min) data were in °C * 10 but we need 
        # them in the equation in °C... so the temperature divided by 10.
        half_tmprange = 0.5 * (tmax - tmin) / 10 
        # to calculate the maximum and minimum temperature we use the average temperature
        self.tmax = tavg + half_tmprange
        self.tmin = tavg - half_tmprange

        # We replace any missing values in the NDVI map
        fillup_value = areamaximum(ndvi, spreadzone(nominal(uniqueid(defined(ndvi))), 0, 1))
        self.ndvi = cover(ndvi, fillup_value)

        # Convert the discharge from m/month to mm/month
        self.Qsurface = qsurface *1000  
        self.rain_month= self.rain_month 
        print('Finished reading input files for timestep %d (date = %s)' % (self.currentTimeStep(), self.current_date))
    #fed

    """
        This method reads in a textfile that contains a table with two columns. The first
        column contains the names of land cover types and the second contains some value 
        value for it (like manning value or P-factor). The columns are seperated by a tab.
        Here are a few example lines from such a file:

        barren_sparsely_vegetated       0.0113
        closed_shrublands       0.4
        croplands       0.04

        The table in this file is read and returned as a dictionary as result.
    """
    def load_landcoverdata(self, lctxtfn):
        # Read from file the mapping from land cover to roughness value
        landcover_data = {}    
        with open(lctxtfn, "rt") as f:
            reader = csv.reader(f, delimiter='\t')
            try:
                for row in reader:
                    landcover_data[row[0]] = float(row[1])
                #rof
            except csv.Error as e:
                sys.exit('Error: Failed to process file %s, at line %d: %e' % (lctxtfn, reader.line_num, e))
            #yrt
        #htiw
        return landcover_data
    #fed

    """
        For each land cover name in the supplied list, a map with that name is
        read from disk which contains the cover fraction data for the land cover
        type. These maps are returned as a dictionary with the land cover name
        as key and the map as value.

        NB: The fractions of all land covers should add up to one. If not, the
        maps for each type are multiplied by the amount that will make their sum 
        add up to the value one. 
    """
    def load_landcovermaps(self, mapnames):
        # dictiomary with the maps of each land cover
        cover_fraction_info = {}

        # set map with the total fraction: we check here on complete coverage
        # as these fractions may be used throughout the script and not
        # only in the computation of Manning's n
        total_vegetation_cover = scalar(0)

        # iterate over the land cover classes and read each land cover map, and 
        # compute the map with the sum of all land cover fractions.
        for key in mapnames:
            # add the fraction for the land cover denoted by key
            cover_fraction_info[key] = self.rd_map(key)

            # add the cover fraction to the total vegetation cover map
            total_vegetation_cover = total_vegetation_cover \
                + cover_fraction_info[key]
        #rof

        # now, correct by dividing each land cover fraction by the total land cover fraction,
        # so that the total of the vegetation covers becomes unity. We do not process cells
        # for which no significant land cover data is available.
        for key in cover_fraction_info.keys():        
            cover_fraction_info[key] = ifthen(total_vegetation_cover > 1.0e-12, \
                cover_fraction_info[key] / max(1.0e-12, total_vegetation_cover))
        #rof

        return cover_fraction_info
    #fed

    # Calculate date we are modelling in this timestep
    def calc_currentdate(self):
        year = self.startdate.year
        month = self.startdate.month + self.currentTimeStep() - 1
        if month > 12:
            month = month - 12
            year = year + 1
        #fi 
        self.current_date = datetime.datetime(year, month, 1)
    #fed

    def netcdf_init(self):
        if self.netcdf_output:
            self.MV = -999.9 # value to use for missing values
            ncAttributes = {} # add netcdf settings/attributes in here
            latitudes = pcr2numpy(ycoordinate(boolean(1)), 0)[:, 0]
            longitudes = pcr2numpy(xcoordinate(boolean(1)), 0)[0, :]
            for key in self.dataProducts.keys():
                varName = self.dataProducts[key]['variable']
                ncFileName = os.path.join(self.outputdir, '%s.nc' % varName) 
                createNetCDF(ncFileName,\
                    longitudes, latitudes, timedData = True, attributes = ncAttributes)
            #rof
        #fi
    #fed

    def netcdf_writedata(self):
        if self.netcdf_output:
            for key in self.dataProducts.keys():
                variableArray = pcr2numpy(getattr(self, key), self.MV)
                varName = self.dataProducts[key]['variable']
                ncFileName = os.path.join(self.outputdir, '%s.nc' % varName)
                data2NetCDF(ncFileName, varName,{'units': self.dataProducts[key]['unit']},\
                    variableArray, posCnt = self.currentTimeStep() - 1, timeStamp = self.current_date, MV = self.MV)
            #rof
        #fi
    #fed

    def report_yearly_result(self, factorname):
        if self.yearly_reports:
            if (len(factorname) > 8):
                print('Warning: %s has more than 8 characters, cannot use it as output name' % factorname)
            else:
                if type(getattr(self, factorname)) is pcraster.Field:
                    self.report(getattr(self, factorname), os.path.join(self.outputdir, factorname))
                    print('Annual factor %s successfully written' % (factorname))
                else:
                    print('Warning: self.%s not found in model class or is not a pcraster map' % (factorname))
                #fi
            #fi
        #fi
    #fed

    def report_monthly_result(self, factorname):
        if self.monthly_reports:
            ts = self.currentTimeStep()
            if (len(factorname) > 8):
                print('Warning: %s has more than 8 characters, cannot use it as output name' % factorname)
            else:
                if type(getattr(self, factorname)) is pcraster.Field:
                    self.report(getattr(self, factorname), os.path.join(self.outputdir, factorname))
                    print('Factor %s successfully written for timestep %d' % (factorname, ts))
                else:
                    print('Warning: self.%s not found or not a pcraster map for timestep %d' % (factorname, ts))
                #fi
            #fi
        #fi
    #fed

    def check_dir(self, pathname):
        if not os.path.isdir(pathname):
            sys.exit("Error: directory %s does not exist or is not a directory" % pathname)
        #fi
        return pathname
    #fed

    def rd_map(self, mapname):
        pathname = os.path.join(self.mapdir, mapname) 
        if os.path.isfile(pathname + ".001"):
            return self.readmap(pathname)
        else:
            pathname = pathname 
            if os.path.isfile(pathname + ".map"):
                return self.readmap(pathname)
            else:
                sys.exit("Error: map %s in %s does not exist" % (mapname, self.mapdir))
            #fi
        #fi
    #fed

    def check_txtfile(self, filename):
        pathname = os.path.join(self.txtdir, filename)
        if not os.path.isfile(pathname):
            sys.exit("Error: %s in %s does not exist or is not a regular file" % (filename, self.txtdir))
        #fi
        return pathname
    #fed

    # Test of the extent of all the data
    #
    # iterate over all information in maps and check whether areas
    # are missing over the clone
    #       
    # initialize here a counter of the pcraster maps encountered in
    # the datasets read and keep a map where an error occurs;
    # in this case, this is initialized as a nominal map with ID 0 and
    # the counter starts at 1. Values once set, are not overwritten
    def check_extent(self):
        if self.tests:
            var_cnt = 0
            var_mv_id = ifthen(self.clone_mask, nominal(0))

            # iterate over all variables and keep the maps
            for variable_name in dir(self):
                # test on type
                if type(getattr(self, variable_name)) is pcraster.Field:
                    #print(variable_name)
                    # update the counter for the new PCRaster map
                    var_cnt = var_cnt + 1

                    # get a test mask 
                    test_mask = ifthen(self.clone_mask, \
                        defined(getattr(self, variable_name)))
                        
                    # update the map with IDs of maps with missing values
                    var_mv_id = ifthenelse(pcrnot(test_mask) & \
                        (var_mv_id == 0), nominal(var_cnt), var_mv_id)

                    # print message if the minimum value in the test mask is zero
                    if cellvalue(mapminimum(scalar(test_mask)), 1)[0] == 0:
                        print('Warning: missing values were encountered in %d: %s' % \
                            (var_cnt, variable_name))
                    #fi
                #fi
            #rof

            # if the maximum value of the map with identified IDs is not zero,
            # then stop with the error
            if cellvalue(mapmaximum(scalar(var_mv_id)), 1)[0] > 0:

                sys.exit('Error: missing values in input, cannot proceed!')
            #fi
        #fi
    #fed
#ssalc

# Main program
if __name__ == "__main__":
    # First set the directory names where our model should find its input data, and 
    # set the name of the clonemap to be used by our model. 
    workdir = os.getcwd()
    mapdname = os.path.join(workdir, 'amazon_fine')
    txtdname = os.path.join(workdir, 'txtdata')
    mapdir=os.path.join(workdir,' mapinput')
    outdname = '/scratch/safaa/output/1990/30_arc_min_run/'

    if (not os.path.isdir(outdname)):
        print("Output directory %s does not exist yet; will try to create it." % outdname)
        try:
            os.mkdir(outdname)
        except OSError:
            print("Error: creation of the directory %s failed" % outdname)
        else:
            print("Successfully created the directory %s " % outdname)
        #yrt
    #fi

    clonedir = os.path.join(workdir, 'mapinput')
    clone_filename = os.path.join(clonedir, 'amazon_5min_mask.map')
    if (not os.path.isfile(clone_filename)):
        sys.exit('Error: %s is not a file.' % clone_filename) 
    #fi

    # The date where our model starts modelling.
    startdate = datetime.date(1990, 1, 1)
    # the number of timesteps (months) our dynamic model should make.
    nrOfTimeSteps = 12

    # initialize our modeland run it
    myModel = RUSLE(clone_filename, mapdname, txtdname, outdname, startdate)
    dynamicModel = DynamicFramework(myModel, nrOfTimeSteps)
    dynamicModel.run()

    print('Hurray, Model Is Finished!')
#fi
