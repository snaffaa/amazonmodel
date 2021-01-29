
# Full process for running Amazon Sediment Production Model

# Preparing the input maps for the Amazon Sediment Production Model
python mapprocessing/convert_globalmaps_to_amazonclone.py 
python mapprocessing/calculate_monthly_rainfall.py 
python mapprocessing/calculate_annual_rainfall.py
python mapprocessing/compute_ground_cover_from_NDVI.py
python mapprocessing/extract_monthly_ground_cover_fraction.py
python mapprocessing/calculate_slope_for_crop_land.py
python mapprocessing/extract_land_cover.py 
python mapprocessing/increase_resolution_to_5arcmin.py

# Running the Amazon Sediment Production Model
python amazon_sediment

# Use output of the Amazon Sediment Production Model as input to the Amazon sediment transport model to 
# estimate the sediment transport in the Amazon river.
python amazon_sediment_transport.py

# Extract Observational Data for number of stations
python extract_observations_for_stations.py

# Extract from the model output, the estimation data of sediment transport for the (same) stations.
python extract_location_for_stations.py
python extract_estimations_for_stations.py
