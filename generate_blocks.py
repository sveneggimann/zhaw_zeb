"""
This script uses the street blocks and the processed building dataset to write out the data to be used for CESAR_P
"""

import os
import numpy as np
import geopandas as gpd
from matplotlib import pyplot
import matplotlib.pyplot as plt
import rasterio
import rasterio.mask
from rasterio.plot import show
from rasterio.merge import merge as rio_merge
from rasterio.plot import show
import shutil
import fiona
from scipy.spatial import KDTree

from progress.bar import Bar

import helper as hp

# Inputs
path_main_data = "C:/Users/eggv/OneDrive - ZHAW/ZBP_shared/1_Research/1015_ZEB_China_Impact_Study/02_Data"

path_blocks = "C:/Users/eggv/OneDrive - ZHAW/ZBP_shared/1_Research/1015_ZEB_China_Impact_Study/02_Data/blocks_data/DT41/3d_urban_form_of_chinese_cities.shp" # Path to street blocks
path_result = "C:/Users/eggv/OneDrive - ZHAW/ZBP_shared/1_Research/1015_ZEB_China_Impact_Study/00_results" # Output_path 

# Blocks to prepare the data (corresponds to the "PARCEL_ID" in the 3d urban form of chinese data layer)
regions = {
    #"Harbin": [],
    #"Beijing": [],
    "Shanghai": [1310120],
    #"Shenzen": []
}

# Read blocks
print("Reading Chinese raw blocks")
blocks = gpd.read_file(path_blocks)

for region in regions:
    print("Region: {}".format(region))
    
    print(" ---reading buildings")
    buildings_region = gpd.read_file(os.path.join(path_main_data, "buildings_osm", "{}_h.geojson".format(region)))

    print(" ---Creating kd search tree for buildings")
    x_coord = buildings_region.geometry.centroid.x.tolist()
    y_coord = buildings_region.geometry.centroid.y.tolist()
    build_kd_tree = KDTree([list(a) for a in zip(x_coord, y_coord)])
    build_rTree = hp.build_rTree(buildings_region)
    
    # Create folder
    path_region_out = os.path.join(path_result, '{}'.format(region))
    if not os.path.exists(path_region_out):
        os.mkdir(path_region_out)
    
    block_ids = regions[region]

    # TODO: if you want to turn, add to turn the buildings

    for block_id in block_ids:

        path_out_block = os.path.join(path_region_out, '{}'.format(block_id))
        if not os.path.exists(path_out_block):
            os.mkdir(path_out_block)

        # Write out block
        block = blocks.loc[blocks['PARCEL_ID'] == block_id].geometry

        # Get all buildings within block and write out
        tree_intersec = build_rTree.intersection(block.geometry.bounds.values[0]) #block.geometry.bounds)
        index_to_take_out = []

        for tree_index in tree_intersec:
            if block.geometry.intersects(buildings_region.iloc[tree_index].geometry):
                index_to_take_out.append(tree_index)

        buildings_block = buildings_region.iloc[[index_to_take_out]]
        buildings_block.to_file(os.path.join(path_out_block, "block.shp"))

        # Write out CESAR-P files

        # Automated path
        path_buildinfo = os.path.join(path_out_block, "BuildingInformation.csv")
        path_shp = os.path.join(path_out_block, "buildings_to_simulate.shp")
        path_csv = os.path.join(path_out_block, "SiteVertices.csv")

        # Read data
        data = gpd.read_file(buildings_block)

        data = hp.calculate_metrics(data)

        # Conver to hong kong coordinate system
        #data = data.to_crs(2326)

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


        


