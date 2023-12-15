"""
This script uses the street blocks and the processed building dataset to write out the data to be used for CESAR_P

Todos that remain to be eimplemented
------
- The function on how the archetypes are classified can be refined
- Turn building orientation (4 times 90°C)
- Add building age category (e.g. from dataset that maybe released from Chinese reserachers next year (see emails))
"""
import os
import numpy as np
import math
import geopandas as gpd
from rasterio.plot import show
from rasterio.merge import merge as rio_merge
from rasterio.plot import show
from scipy.spatial import KDTree
from progress.bar import Bar

import helper as hp

# Inputs
path_main_data = "C:/Users/eggv/OneDrive - ZHAW/ZBP_shared/1_Research/1015_ZEB_China_Impact_Study/02_Data"
path_blocks = "C:/Users/eggv/OneDrive - ZHAW/ZBP_shared/1_Research/1015_ZEB_China_Impact_Study/02_Data/blocks_data/DT41/3d_urban_form_of_chinese_cities_reprojected.shp" # Path to street blocks
path_result = "C:/Users/eggv/OneDrive - ZHAW/ZBP_shared/1_Research/1015_ZEB_China_Impact_Study/00_results" # Output_path 

# Write out CESAR-P
write_out_cesar_p = True

# Assumptions
mean_h_assumption = 3 # [m]

# Blocks to prepare the data (corresponds to the "PARCEL_ID" in the 3d urban form of chinese data layer)
regions = {
    "Harbin": [
        599772, 600174, 600227, 600927, 601278, 601114, 600845, 600426, 600463, 600460, 600325, 602626, 600899, 601094, 601358,
        599655, 598446, 601952, 602439, 603716, 602443, 602965, 603952, 601498, 601312, 601111, 600883, 601254, 601274
        ],
    "Beijing": [
        1190217, 1190324, 1190361, 1193878, 1197489, 1197884, 1198122, 1198127, 1198278, 1199030, 1199041, 1199894, 1200598, 1201109,
        1201354, 1201568, 1201607, 1202433, 1203844, 1204098, 1204099, 1204148, 1205334, 1205515, 1205990, 1205991, 1205998, 1206013,
        1207286, 1208751, 1209006, 1209136, 1209973, 1210095, 1211158, 1211427, 1211699, 1211704, 1211844, 1212426, 1212916, 1213054,
        1213055, 1213308, 1213310, 1213391, 1215247, 1215950, 1216763, 1217108, 1217146, 1218360, 1218628, 1218909, 1219261, 1220907,
        1221655, 1222128, 1222613, 1223007, 1223038, 1230998, 1231496, 1231504, 1231513, 1231558, 1231586, 1231612, 1231617, 1231658,
        1232126, 1239666
        ],
    "Shanghai": [
        1313173, 1313165, 1313132, 1313102, 1313084, 1312643, 1312542, 1312327, 1311966, 1311874, 1310662, 1310504, 1310503,
        1310492, 1310390, 1310306, 1310244, 1310243, 1310217, 1310171, 1310120, 1310108, 1310069, 1310002, 1309861, 1309855,
        1308044, 1308019, 1308014, 1307883, 1303096, 1303059, 1303048, 1303025, 1301877, 1298109, 1297821, 1297805, 1289381,
        1288829, 1285367, 1285279, 1285256, 528903, 528853, 528747, 528642
        ],
    "Shenzen": [
        65846, 65859, 65872, 65875, 65888, 65900, 65945, 65965, 66502, 66616, 66703, 66900, 66914, 66917, 66953, 67264, 67267,
        67756, 67814, 67970, 68048, 68357, 68499, 68774, 68829, 68843, 68870, 69127, 69937, 70353, 70364, 70419, 70467, 70516,
        70581, 70628, 70860, 72233, 73059, 724506],
}

# Read blocks
print("Reading Chinese raw blocks")
blocks = gpd.read_file(path_blocks)

for region in regions.keys():
    print("--------------------------")
    print("Region: {}".format(region))
    print("--------------------------")
    print(" ---reading buildings")
    buildings_region = gpd.read_file(os.path.join(path_main_data, "buildings_osm", "{}_h.geojson".format(region)))

    # -----------------------
    print(" ---classify building archetypes")
    # -----------------------
    progress_par = Bar('     Classingying archetpyes: ', max=(len(buildings_region)))
    archetypes = []
    for index in buildings_region.index:
        try:
            build_geom = buildings_region.loc[index].geometry.buffer(0)
            
            # Mean height
            mean_h_round = buildings_region.loc[index]['tiff_h']
            
            # MAx height
            max_h = buildings_region.loc[index]['max_h']


            # Function how buildings are chlasified
            archetpye = hp.classify_archetype(
                build_geom,
                max_h, #mean_h_round,
                assumed_height_floor=mean_h_assumption)
            archetypes.append(archetpye)
        except:
            print("something went wrong with the building classificaiton")
            archetypes.append("error")
        progress_par.next()
    progress_par.finish()
    buildings_region['archetype'] = archetypes

    print("    reating kd search tree for buildings")
    x_coord = buildings_region.geometry.centroid.x.tolist()
    y_coord = buildings_region.geometry.centroid.y.tolist()
    build_kd_tree = KDTree([list(a) for a in zip(x_coord, y_coord)])
    build_rTree = hp.build_rTree(buildings_region)
    
    # Create folder
    path_region_out = os.path.join(path_result, '{}'.format(region))
    if not os.path.exists(path_region_out):
        os.mkdir(path_region_out)
    
    block_ids = regions[region]

    # ITerate blocks to write out files (buildings, block, cesar-p etc)
    for block_id in block_ids:

        path_out_block = os.path.join(path_region_out, '{}'.format(block_id))
        if not os.path.exists(path_out_block):
            os.mkdir(path_out_block)

        # Check that the same crs
        assert blocks.crs == buildings_region.crs

        # Write out block
        block = blocks.loc[blocks['PARCEL_ID'] == block_id]
        block.to_file(os.path.join(path_out_block, "block.shp"))

        # Get all buildings within block and write out
        tree_intersec = build_rTree.intersection(block.geometry.bounds.values[0]) #block.geometry.bounds)
        index_to_take_out = []

        for tree_index in tree_intersec:
            block_geometry = block.geometry.tolist()[0]
            # check whether building centroid is within block
            if block_geometry.contains(buildings_region.iloc[tree_index].geometry.centroid):
                index_to_take_out.append(tree_index)
        buildings_block = buildings_region.iloc[index_to_take_out]
        buildings_block.to_file(os.path.join(path_out_block, "buildings_block.shp"))

        # -----------------------
        # Write out CESAR-P files
        # -----------------------
        if write_out_cesar_p:
            path_buildinfo = os.path.join(path_out_block, "BuildingInformation.csv")
            path_shp = os.path.join(path_out_block, "buildings_to_simulate.shp")
            path_csv = os.path.join(path_out_block, "SiteVertices.csv")

            # Read data
            data = buildings_block #gpd.read_file(os.path.join(path_out_block, "buildings_block.shp"))

            data = hp.calculate_metrics(data)

            # --------------------------------------------------------------------------------
            # Assign building specifications for CEAR-P
            # --------------------------------------------------------------------------------
            data['ORIG_FID'] = range(data.shape[0])
            data['height'] = data['tiff_h']                              # Height from tiff files

            data['levels_f'] = data['tiff_h'] / mean_h_assumption        # Calculate number of floors
            data['levels_f'] = data['levels_f'].apply(np.ceil)
            data['floor_level'] = data['tiff_h'] / mean_h_assumption     # Calculate number of floors
            data['floor_level'] = data['floor_level'].apply(np.ceil)
            
            data['c_type'] = 'MFH' # TODO: replace by chinese archetypes data['archetype']
            data['c_age'] = 2      # Age class. Note: Needs to be derived from somewhere
            data['unique_id'] = data['ORIG_FID']                                    
    

            # Create final shapefiles with the building to be simulated and BuildingInformation.csv file
            pd_out, gdf_shp = hp.shapefile_as_cesarinput(data)

            columns_to_keep = ['ORIG_FID', 'SIA2024BuildingType','BuildingAge']
            pd_out = pd_out[columns_to_keep]

            pd_out.to_csv(path_buildinfo, index=None, header=True)
            gdf_shp.to_file(path_shp)

            # Create SiteVertices.csv file
            df_geometry = hp.read_sitevertices_from_shp(path_shp,round_digits=2)
            df_geometry.to_csv(path_csv, index=False)


    

print(" ---- finished scripts ----")

