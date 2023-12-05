"""
(old - This script is now integrated in "generate_blocks.py")

Load shapefile from china area and convert to CESAR-P input files

"""
import os
import geopandas as gpd

import helper as hp

print("imports loaded")

# ---Define paths
main_path = "C:/Users/eggv/OneDrive - ZHAW/ZBP_shared/1_Research/1015_ZEB_China_Impact_Study/04_Working_phase/geometry_collection/"
path_shape = os.path.join(main_path, "test_buildings.shp")
path_out = os.path.join(main_path, "_results")
path_processed_buildings = os.path.join(path_main_data, "buildings_osm")

# Automated path
path_buildinfo = os.path.join(path_out, "BuildingInformation.csv")
path_shp = os.path.join(path_out, "buildings_to_simulate.shp")
path_csv = os.path.join(path_out, "SiteVertices.csv")

# Read data
data = gpd.read_file(path_shape)

data = hp.calculate_metrics(data)

# Conver to hong kong coordinate system
data = data.to_crs(2326)

# --------------------------------------------------------------------------------
# (a) Assign building specifications (e.g. using rule-based classification)
# --------------------------------------------------------------------------------
data['ORIG_FID'] = range(data.shape[0])
data['levels_f'] = 2
data['height'] = 10
data['floor_level'] = 2
data['c_type'] = 'MFH'
data['c_age'] = 2
data['unique_id'] = data['ORIG_FID']

#TODO: Assing China building archetype based on rule-based classification
#TODO: Assign height
#TODO: Assing c_age

# --------------------------------------------------------------------------------
# (b) Assign building ageclass
# --------------------------------------------------------------------------------

# Create final shapefiles with the building to be simulated and BuildingInformation.csv file
pd_out, gdf_shp = hp.shapefile_as_cesarinput(data)

columns_to_keep = ['ORIG_FID', 'SIA2024BuildingType','BuildingAge']
pd_out = pd_out[columns_to_keep]

pd_out.to_csv(path_buildinfo, index=None, header=True)
gdf_shp.to_file(path_shp)

# Create SiteVertices.csv file
df_geometry = hp.read_sitevertices_from_shp(path_shp,round_digits=2)
df_geometry.to_csv(path_csv, index=False)

print(" -- finished --")
