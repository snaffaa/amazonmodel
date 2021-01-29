import pandas as pd 
import csv 
import os
import datetime 
import numpy as np
import sys

inputdir = '/home/naffa002/projects/finalmodel/sediment_flux/new_stations'
outputdir = '/home/naffa002/projects/finalmodel/sediment_flux/output'

def calculate_sediment_flux( ):
    dateparse = lambda x: pd.datetime.strptime(x, '%d-%m-%Y')
    # read the csv file of discharge observations
    df_dm = pd.read_csv(os.path.join(inputdir,'updated_discharge_data_for_all_new_stations.csv'), converters={'id_station': lambda x: str(x)}, \
            dtype={'Vazao': np.float64}, usecols=['id_station', 'Data', 'Vazao'], parse_dates=[1], date_parser=dateparse, decimal = '.',\
             delimiter=';')
    print(df_dm)
    # remove the na (missing) values
    df_dm = df_dm.dropna()
    print(df_dm)
    # select year, month from the date column based on the station id
    df_dm["year"] = pd.DatetimeIndex(df_dm['Data']).year
    print("years in the results are", df_dm["year"])
    df_dm["month"] = pd.DatetimeIndex(df_dm['Data']).month
    print("months in the results are", df_dm["month"])
    filtered_data = df_dm[df_dm["id_station"]=="12230000"]
    print('filtered_data is ', filtered_data)
    print(filtered_data.iloc[0:5, 1:4])
    # read the csv file of sediment observations
    df_sediment= pd.read_csv(os.path.join(inputdir,'all_sediment_stations_observation_data.csv'), converters={'id_station': lambda x: str(x)}, \
            dtype={'ConcentracaoMatSuspensao': np.float64}, usecols=['id_station', 'date', 'ConcentracaoMatSuspensao'], parse_dates=[1], date_parser=dateparse, decimal = ',', \
             delimiter=';')
    print(df_sediment)
    
    
  
    df_sediment = df_sediment.dropna()
    df_sediment["year"] = pd.DatetimeIndex(df_sediment['date']).year
    df_sediment["month"] = pd.DatetimeIndex(df_sediment['date']).month
     # merge discharge and sediment csv files together based on
     # year, month and station_id
    combo = pd.merge(
        df_sediment,
        df_dm,
        on=['id_station','year', 'month'])
    print(combo)
       
    # convert suspended sediment concentration from mg/l to kg/l and 
    #the discharge from m3/s to l/month to calculate the suspended sediment supply in kg/month
    
    
    combo['sed_flux'] = (combo['ConcentracaoMatSuspensao']*1e-6) *(combo['Vazao']*2.628e+9)
    print(combo)
    # calculate the annual flux
    yearly_sums = combo.groupby(['id_station', 'year'])['sed_flux'].sum().\
        reset_index(name='yearly_sums')
    # calculate the annual discharge
    discharge_annual= combo.groupby(['id_station'])['Vazao'].mean().\
        reset_index(name='discharge_annual')
    # add the annual sediment to the merged dataframe
    combo = pd.merge(combo, yearly_sums, on=['id_station', 'year'], how='left')
    
    average = combo.groupby(['id_station'])['yearly_sums'].mean().\
        reset_index(name='annual_average')
        
    print(average.columns)
    #add the avaerage to the merged dataframe 
    combo = pd.merge(combo, average, on=['id_station'], how='left')
    #add the annual discharge to the merged dataframe
    combo = pd.merge(combo, discharge_annual, on=['id_station'], how='left')
    print(combo)
    #convert the dataframe to csv file
    combo.to_csv(os.path.join(outputdir,'sediment_flux_new_stations_updated.csv'), sep = ';')
    return combo
#fed
    

if __name__ == "__main__":
     calculate_sediment_flux()

     sys.exit('all done')
#fi
    







