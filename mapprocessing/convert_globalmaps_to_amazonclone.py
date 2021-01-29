#
# Reduce the global maps to the amazon area.
#
#   - Limit the maps based on the attributes in the amazon clone map.
#

import os, sys, shutil
from spatialDataSet2PCR import *

#-initialization
workdir = os.getcwd()
inputPath = os.path.join(workdir, './globalmaps')
mapdir = os.path.join(workdir, './mapinput')
outputPath = os.path.join(workdir, './amazon_input/')

cloneFileName = os.path.join(mapdir, 'amazon_5min_mask.map')
sourceCloneFileName = os.path.join(inputPath, 'globalclone.map')
cloneAttributes = spatialAttributes(cloneFileName)
dummyAttributes = spatialAttributes(sourceCloneFileName)
resampleRatio = dummyAttributes.xResolution/cloneAttributes.xResolution
dummyVariableName = 'var'

for root, dirs, files in os.walk(inputPath):
    sourcePath= root
    targetPath= root.replace(inputPath,outputPath)
    print '*\tprocessing %s -> %s' % (sourcePath, targetPath)
    if not os.path.isdir(targetPath):
        os.makedirs(targetPath)
    #fi
    for fileName in os.listdir(sourcePath):
        sourceFileName= os.path.join(sourcePath,fileName)
        targetFileName= os.path.join(targetPath,fileName)
        if os.path.isfile(sourceFileName):
            try:
                dummyAttributes= spatialAttributes(sourceFileName)
                if dummyAttributes.dataType ==  'Byte':
                    valueScale= pcr.Boolean
                if dummyAttributes.dataType ==  'Int32':
                    valueScale= pcr.Nominal
                if dummyAttributes.dataType ==  'Float32':
                    valueScale= pcr.Scalar
                spatialDataSet(dummyVariableName,sourceFileName,\
                    dummyAttributes.dataType, valueScale,\
                    cloneAttributes.xLL, cloneAttributes.xUR,\
                    cloneAttributes.yLL, cloneAttributes.yUR,\
                    cloneAttributes.xResolution, cloneAttributes.yResolution,\
                    lines = cloneAttributes.numberRows,\
                    pixels = cloneAttributes.numberCols, warp= True,\
                    outputFileName = targetFileName, test= True)
                command= 'mapattr -c %s %s' % (cloneFileName,targetFileName)
                os.system(command)
                print '\tfile resampled from %s to %s' % (sourceFileName,targetFileName)
            except:
                shutil.copy(sourceFileName,targetFileName)
                print '\tfile copied from %s to %s' % (sourceFileName,targetFileName)
        #fi
    #rof

#-print all done
print 'all done'
