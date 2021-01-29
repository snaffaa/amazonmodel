import os, sys, types
import glob, subprocess
import pcraster as pcr
from copy import copy

from spatialDataSet2PCR import * 
 
def read_spatialattributes_from_map(inputFileName, maxLength= -1):
    '''Returns the map attributes for the spatial file name specified'''

    #-max length specifies a cutoff to speed up processing of large datasets with multible bands
    #-set local variables: mapInformation holds identifier strings and offset for the variables of interest
    mapInformation={}
    mapInformation['dataFormat']= 'Driver', 0, types.StringTypes
    mapInformation['numberRows']= 'Size is', 1, types.IntType
    mapInformation['numberCols']= 'Size is', 0, types.IntType
    mapInformation['xResolution']= 'Pixel Size =', 0, types.FloatType
    mapInformation['yResolution']= 'Pixel Size =', 1, types.FloatType
    mapInformation['xLL']= 'Lower Left', 0, types.FloatType
    mapInformation['yLL']= 'Lower Left', 1, types.FloatType
    mapInformation['xUR']= 'Upper Right', 0, types.FloatType
    mapInformation['yUR']= 'Upper Right', 1, types.FloatType
    mapInformation['dataType']= 'Type=', 0, types.StringTypes
    mapInformation['minValue']= 'Min=', 0, types.FloatType
    mapInformation['maxValue']= 'Max=', 0, types.FloatType
    mapInformation['noDataValue']= 'NoData Value=', 0, types.FloatType
    #-set dictionary to hold relevant information
    mapAttributes={}
    #-get information with gdalinfo
    command= 'gdalinfo "%s"' % inputFileName
    cOut,err = subprocess.Popen(command, stdout= subprocess.PIPE, stderr= subprocess.PIPE,shell=True).communicate()
    if err != '' and not err[:7].lower() == 'warning':
        sys.exit('Error: no information could be retrieved for the spatial dataset %s' % inputFileName)
    #-get map information
    for mapAttribute, info in mapInformation.iteritems():
        matched = False
        rawList = cOut.split('\n')[:maxLength]
        (key, entryNumber, typeInfo) = info
        #-iterate until the entry is found or the list of entries is empty
        while not matched and len(rawList) > 0:
            #-pop first entry with white space, including \r removed
            entry= rawList.pop(0).strip()
            #-if found, set matched true and process
            if key in entry:
                matched = True
                posCnt = entry.find(key)
                #-if entry found, get actual value
                if posCnt >= 0:
                    posCnt += len(key)
                    rawStr = entry[posCnt:].split(',')[entryNumber]
                    rawStr = rawStr.strip('=:() \t')
                    if typeInfo in [types.IntType,types.FloatType,types.LongType]:
                        rawStr = rawStr.split()[0]
                    #-rawStr of entry returned, process accordingly
                    rawStr= rawStr.strip('=:() \t')
                    if typeInfo == types.IntType or typeInfo ==  types.LongType:
                        try:
                            mapAttributes[mapAttribute]= int(rawStr)
                        except:
                            sys.exit('Error: map attributes could not be processed for %s' % mapAttribute)
                    elif typeInfo == types.FloatType:
                        try:
                            mapAttributes[mapAttribute]= float(rawStr)
                        except:
                            sys.exit('Error: map attributes could not be processed for %s' % mapAttribute)
                    else:
                        mapAttributes[mapAttribute]= rawStr
        #elihw
        if mapAttribute in ['xResolution','yResolution']:
            mapAttributes[mapAttribute]= abs(mapAttributes[mapAttribute])
        #fi    
    #rof
    return mapAttributes
#fed

def remove_tmpfile(tmpfilename):    
    #-remove temporary file if possible
    try:
        os.remove(tmpfilename)
    except:
        pass
      
    
if __name__ == "__main__":
    # just the name we use to construct temporary filenames/maps
    dummy = "dummy"

    workdir = os.getcwd()
    inputdir = os.path.join(workdir, 'amazon_input')
    mapdir=os.path.join(workdir,'mapinput')
    outputdir = os.path.join(workdir, 'amazon_fine')

    # This is the clone map for the amazon at coarse resolution  0.5 deg (30 arc)
    coarse_clone_amazon = os.path.join(mapdir, 'mask_M45.map')
    tmpfile = os.path.join(outputdir, "temp_" + dummy + ".map")

    # retrieve coarse spatial attributes of the amazon clonemap
    cas = read_spatialattributes_from_map(coarse_clone_amazon)
    # and make a copy and change the resolution to fine (by factor 6)
    fas = copy(cas)
    fas["numberCols"] = cas["numberCols"] * 6
    fas["numberRows"] = cas["numberRows"] * 6
    fas["xResolution"] = cas["xResolution"] / 6.0
    fas["yResolution"] = cas["yResolution"] / 6.0
    
    # create a fine amazon clone map using the coarse map as basis
    remove_tmpfile(tmpfile) # remove tmpfile if one is still lying around after a crash
    command= 'mapattr -s -R %d -C %d -x %f -y %f -l %f -P "yb2t" -B %s' % \
        (fas["numberRows"], fas["numberCols"], fas["xLL"], fas["yUR"], fas["xResolution"], tmpfile)
    os.system(command)
    # set the fine map of the amazon as clone, so that our output maps will have location attributes from this map
    pcr.setclone(tmpfile)
    remove_tmpfile(tmpfile)
    
    # We iterate over all the maps (regular and timeseries) in the inputdir
    for imap in glob.glob(os.path.join(inputdir, "*.map")) + \
                glob.glob(os.path.join(inputdir, "*.[0-9][0-9][0-9]")):
	mapname = os.path.split(imap)[1]
	omap = os.path.join(outputdir, mapname)
        # read attributes from each map and check what kind of map it is
        mapinfo = read_spatialattributes_from_map(imap)
        if 'ldd' in mapname or 'LDD' in mapname: 
            # We don't convert this map, we just fetch a finer resolution version of it
            print "Replacing low resolution " + mapname + " with a higher resolution version from PCRGLOBWB"
            lddmap = '/data/hydroworld/PCRGLOBWB20/input5min/routing/lddsound_05min.map'
            new_map = getattr(spatialDataSet(dummy,\
                lddmap, 'Byte', 'LDD', fas["xLL"], fas["xUR"], fas["yLL"], fas["yUR"],\
                fas["xResolution"], fas["yResolution"],\
                pixels= fas["numberCols"], lines= fas["numberRows"], test = True), dummy)
            # Each cell on a local drain direction map must have a pit at the end of its downstream path. 
            # Because we cropped the global ldd map to the amazon area only, this is no longer the case at
            # the edges of the map. We use lddrepair to make sure (again) that all downstream paths will end
            # in a pit cell. (http://pcraster.geo.uu.nl/pcraster/4.2.1/documentation/pcraster_manual/sphinx/op_lddrepair.html)
            fixed_map = pcr.lddrepair(new_map)   
            pcr.report(fixed_map, omap)
        elif(mapinfo["dataType"] == "Byte"): 
            # It is a Boolean map that we convert using gdal_translate
            print "Changing resolution boolean map:" + mapname
            new_map = getattr(spatialDataSet(dummy,\
                imap, 'Byte', 'Boolean', fas["xLL"], fas["xUR"], fas["yLL"], fas["yUR"],\
                fas["xResolution"], fas["yResolution"],\
                pixels= fas["numberCols"], lines= fas["numberRows"], test = True), dummy)
            pcr.report(new_map, omap)
        else:
            # It is a scalar map that we translate using gdal_translate
            print "Translating: " + mapname
            command= 'gdal_translate -ot %s -of PCRaster -mo VALUESCALE=VS_%s -projwin %f %f %f %f -outsize %s %s "%s" "%s" -q' %\
                ('FLOAT32', 'Scalar', fas["xLL"], fas["yUR"], fas["xUR"], fas["yLL"], \
                '%d' % fas["numberCols"],'%d' % fas["numberRows"], imap, tmpfile)
            #print command
            cOut,err = subprocess.Popen(command, stdout= subprocess.PIPE, stderr= subprocess.PIPE,shell=True).communicate()
            #print err
            if err != '' and not err[:7].lower() == 'warning':
                sys.exit('Error: could not read in and convert %s: %s' % (mapname, err))
            #read resulting map
            new_map = pcr.readmap(tmpfile)
            pcr.report(new_map, omap)
        #fi
    #rof
    
    remove_tmpfile(tmpfile)
