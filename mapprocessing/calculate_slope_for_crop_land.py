# 1. create with a clone of 5' a map of unique ids
# 2. read in this map using a clone at 30" (creates 10 x 10 30" cells with the same id)
# 3. for the 100 cells we then have the slope angles. with this we compute the areaaverage 
#    but also the percentiles with pcr.areaorder. 
# So, when we have x% agriculture in an area, we can select the corresponding equivalent slope angle.

import os, sys
import pcraster as pcr

sys.path.insert(0, '/home/beek0120/netData/pylib')

from spatialDataSet2PCR import spatialAttributes, spatialDataSet, setClone

def main():

    # set the clone
    clonemap = '/scratch/rens/for_safaa/amazon_input/globalclone.map'
    cell_id_filename = 'cellid.map'
    slope_filename = '/home/beek0120/netData/GlobalDataSets/gtopo30/gtopo30_slope_gaez.map'
    selected_slope_filename = 'sel_slope.map'
    cropfraction_filename = '/scratch/rens/for_safaa/mirca_agriculture/mirca_total_cropland_fraction.map'
    selected_slope_average_filename='sel_slope_ave.map'
   
    

    # clone at 5 arc minute and cell id
    
    cloneattributes = spatialAttributes(clonemap)
    
    cloneattributes.numberCols = cloneattributes.numberCols * 6
    cloneattributes.numberRows = cloneattributes.numberRows * 6
    cloneattributes.xResolution = 0.0833333
    cloneattributes.yResolution = 0.0833333
    
    setClone(cloneattributes)
    
    cell_id = pcr.uniqueid(pcr.boolean(1))
    pcr.report(cell_id, cell_id_filename)
    
    
    
    # clone at 30 arc sec
    cloneattributes = spatialAttributes(clonemap)

    cloneattributes.numberCols = cloneattributes.numberCols * 60
    cloneattributes.numberRows = cloneattributes.numberRows * 60
    cloneattributes.xResolution = 0.00833333
    cloneattributes.yResolution = 0.00833333
    
    setClone(cloneattributes)

    # read in the cell id and the slope
    cell_id = getattr(spatialDataSet( \
                    'dummy', cell_id_filename, \
                    'INT32', 'NOMINAL', \
                    cloneattributes.xLL, cloneattributes.xUR, \
                    cloneattributes.yLL, cloneattributes.yUR, \
                    cloneattributes.xResolution, cloneattributes.yResolution, \
                    pixels= cloneattributes.numberCols, lines= cloneattributes.numberRows, \
                    test= True), 'dummy')
    # read in the cell id and the slope
    slope = getattr(spatialDataSet( \
                    'dummy', slope_filename, \
                    'FLOAT32', 'SCALAR', \
                    cloneattributes.xLL, cloneattributes.xUR, \
                    cloneattributes.yLL, cloneattributes.yUR, \
                    cloneattributes.xResolution, cloneattributes.yResolution, \
                    pixels= cloneattributes.numberCols, lines= cloneattributes.numberRows, \
                    test= True), 'dummy')
    
    slope_order = pcr.areaorder(slope, cell_id)
    slope_order = slope_order / pcr.areamaximum(slope_order, cell_id) * 100.0
 
    desired_threshold = getattr(spatialDataSet( \
                    'dummy', cropfraction_filename, \
                    'FLOAT32', 'SCALAR', \
                    cloneattributes.xLL, cloneattributes.xUR, \
                    cloneattributes.yLL, cloneattributes.yUR, \
                    cloneattributes.xResolution, cloneattributes.yResolution, \
                    pixels= cloneattributes.numberCols, lines= cloneattributes.numberRows, \
                    test= True), 'dummy')
    
    desired_threshold = pcr.min(1.0, desired_threshold) * 100.0
    
    selection = slope_order <= desired_threshold
    selected_slope = pcr.areamaximum(pcr.ifthen(selection, slope), cell_id)
    
    pcr.report(selected_slope, selected_slope_filename)
 
    # clone at 5 arc minute and cell id
    
    cloneattributes = spatialAttributes(clonemap)
    
    cloneattributes.numberCols = cloneattributes.numberCols * 6
    cloneattributes.numberRows = cloneattributes.numberRows * 6
    cloneattributes.xResolution = 0.0833333
    cloneattributes.yResolution = 0.0833333
    
    setClone(cloneattributes)
    # read in the cell id and the slope
    selected_slope = getattr(spatialDataSet( \
                    'dummy', selected_slope_filename, \
                    'FLOAT32', 'SCALAR', \
                    cloneattributes.xLL, cloneattributes.xUR, \
                    cloneattributes.yLL, cloneattributes.yUR, \
                    cloneattributes.xResolution, cloneattributes.yResolution, \
                    pixels= cloneattributes.numberCols, lines= cloneattributes.numberRows, \
                    test= True), 'dummy')
    cell_id = getattr(spatialDataSet( \
                    'dummy', cell_id_filename, \
                    'INT32', 'NOMINAL', \
                    cloneattributes.xLL, cloneattributes.xUR, \
                    cloneattributes.yLL, cloneattributes.yUR, \
                    cloneattributes.xResolution, cloneattributes.yResolution, \
                    pixels= cloneattributes.numberCols, lines= cloneattributes.numberRows, \
                    test= True), 'dummy')
                    
    
    selected_slope_average = pcr.areaaverage(selected_slope,cell_id)
   

    pcr.report(selected_slope, selected_slope_filename)
    pcr.report(selected_slope_average, selected_slope_average_filename)
    
    
if __name__ == "__main__":
    
    main()
