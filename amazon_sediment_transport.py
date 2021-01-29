'''
Amazon Sediment Transport Model developed from PCRGLOB-WB using 
erosion and sedimentation module.
-PCR-GLOBWB Routing scheme: basic implementation for python version
with the following characteristics:
 * total water storage [m3] of which volume in excess of channel storage is used as flood volume
   and distributed over the surrounding area using smooth floodplain transitions based on Kavetski
   and Kuczera (2007)
 * flood wave propagation can be described by means of the kinematic wave or a simpler traveltime function
 * lakes and reservoirs are included by means of the waterbodies class
'''



import os, sys, datetime, calendar
import pcraster as pcr
import pcraster.framework as pcrm
from time import gmtime, strftime

from routing_wb_dev import pcrglobRoutingSediment

##################
# Model settings #
##################

def set_general_settings(startDate, endDate, currentYear, duration):
    settings = {\
        'testLocalWaterBalance': True,\
        'getSurfaceWaterAttributes': True,\
        'outputPath': os.path.join('/scratch/safaa/sediment_transport_amazon_30arcmin_/without_reservoir/', strftime("%Y%m%d", gmtime())),\
        'duration': duration,\
        'timeSec': 86400,\
        #-area fractions
        'areaFractions': [0.0,0.01,0.05,0.10,0.20,0.30,0.40,\
        0.50,0.60,0.70,0.80,0.90,1.00],\
        # reduction parameter of smoothing interval and error threshold 
        'reductionKK': 0.5,\
        'criterionKK': 40.0,\
        # manning's N
        'channelManN': 0.04,\
        'floodplainManN': 0.10,\
        #-particle diameter [m],
        # particle density [kg/m3] and 
        # sediment uptake factor [-]
        'particleDiameter': 5.0e-6,\
        'particleDensity': 1000.0,\
        'sedimentUptakeFactor': 1.0,\
        #-parameters related to the computation of the transport capacity
        # parameters of unit stream power [m/s] and particle density [kg/m3]
        # applyScalingTC determines whether the unit stream power is scaled
        # below the specified threshold unit stream power
        'criticalUnitStreamPower': 0.004,\
        'thresholdUnitStreamPower': 0.007,\
        'applyScalingTC': True,\
        #-dates
        'currentYear': currentYear,\
        'startDate': startDate,\
        'endDate': endDate
    }
    return settings
#fed

#-map settings
def set_mapSettings(mapdir, fine_maps,inputdir, cloneMapFileName): 
    settings = {\
        'clone': cloneMapFileName,\
        'LDD': os.path.join(mapdir, 'LDD.map'),\
        'cellArea': os.path.join(mapdir, 'cellarea_30min.map'),\
        'channelWidth': os.path.join(mapdir, 'bankfull_width.map'),\
        'channelDepth': os.path.join(mapdir, 'bankfull_depth.map'),\
        'channelGradient': os.path.join(mapdir, 'channel_gradient.map'),\
        #-root of file name with maps of relative elvation above floodplain
        # and associated fractions
        'relZFileName': os.path.join(mapdir, 'dzRel%04d.map'),\
        #-runoff
        #'runoffVarName': 'Runoff',\
        'runoffVarName': 'land_surface_runoff',\
        'runoff': os.path.join(mapdir, 'runoff_dailyTot_reprojected_30arcmin.nc'),\
        #-sediment production [kg/cell/day]
        'sedimentProduction': os.path.join(mapdir, 'sediment_input_daily_new_stations_1990_30arcmin.map')
        
    }
    
    return settings
#fed

#-water bodies settings
def set_waterBodySettings(mapdir, fine_maps): 
    settings = {\
        'waterBodyID': os.path.join(mapdir, 'waterbodyid_2010.map'),\
        'waterBodyOutlet': os.path.join(mapdir, 'waterbodyoutlet_2010.map'),\
        'waterBodyType': os.path.join(mapdir, 'waterbodytype_2010.map'),\
        'bankfullDischarge': os.path.join(mapdir, 'discharge_bankfull.map'),\
        'waterBodyParamTBL': '/home/beek0120/PCRGLOBWB/SedimentTransport/StaticMaps/grand_parameters.tbl',\
        'reservoirDemandTSS': '/home/beek0120/PCRGLOBWB/Amazon/SouthAmerica/maps/reservoirdemand.tss'
     }
    return settings
#fed
    
#-initial settings
def set_initialSettings(mapdir):
    settings= {\
        'longtermAverageDischarge': os.path.join(mapdir,'dis_avg_1990.map'),\
        'Q': os.path.join(mapdir, 'discharge_1990_without_oneyear-12-31.map'),\
        'waterStorage': os.path.join(mapdir,'waterStorage_1990_without_oneyear-12-31.map'),\
        #-sediment stocks of river and floodplain
        'availableSedimentStock': os.path.join(mapdir, 'availableSedimentStock_1990_without_oneyear-12-31.map'),\
        'floodplainSedimentStock': os.path.join(mapdir, 'floodplainSedimentStock_1990_without_oneyear-12-31.map'),\
        'sedimentLoad': os.path.join(mapdir, 'sedimentLoad_1990_without_oneyear-12-31.map'),\
    }
    return settings
#fed

def main():
    #-where to find stuff
    inputdir= None
    mapdir = '/home/naffa002/projects/finalmodel/mapinput/halfdegree_maps'
    fine_maps = None
    cloneMapFileName = os.path.join(mapdir, 'mask_amazon_30arcmin.map')

    #-parameters
    duration    = 1.0
    currentYear = 1990
    startDate   = datetime.datetime(1990,01,1)
    endDate     = datetime.datetime(1999,12,31)
    
    # create dictionaries with settings
    generalSettings = set_general_settings(startDate, endDate, currentYear, duration)
    mapSettings = set_mapSettings(mapdir, fine_maps,inputdir, cloneMapFileName)
    waterBodySettings = set_waterBodySettings(mapdir, fine_maps)
    initialSettings = set_initialSettings(mapdir)

    #-start
    pcr.setclone(cloneMapFileName)
    routingscheme = pcrglobRoutingSediment(startDate, endDate,\
        generalSettings, mapSettings, waterBodySettings, initialSettings)
    dynRouting = pcrm.DynamicFramework(routingscheme, lastTimeStep = int((endDate.toordinal() - startDate.toordinal()+1)/duration),\
        firstTimestep = 1)
    dynRouting.run()
#fed

if __name__ == "__main__":
    main()
#fi
