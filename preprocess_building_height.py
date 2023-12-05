"""
Load Chinese regions, corresponding openstreetmap data and assign building heights through spatial overlay.


Used code-snipplets / online-help:
- merge tif: https://trac.osgeo.org/gdal/wiki/CatalogueForQIS
- masking tifs: https://rasterio.readthedocs.io/en/latest/topics/masking-by-shapefile.html 
- https://automating-gis-processes.github.io/CSC18/lessons/L6/raster-mosaic.html
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
from progress.bar import Bar

import helper as hp

# This is a required preprocessing step if a city is covered by several tifs that need to be merged

pre_processing_merge = False # Only once to merge tifs

# ----------------
# Paths
# ----------------
path_main_data = "C:/Users/eggv/OneDrive - ZHAW/ZBP_shared/1_Research/1015_ZEB_China_Impact_Study/02_Data"

path_china_regions_1 = os.path.join(path_main_data, "boundaries", "geoBoundaries-CHN-ADM1.geojson")
path_china_regions_2 = os.path.join(path_main_data, "boundaries", "geoBoundaries-CHN-ADM2.geojson")
path_osm_data = os.path.join(path_main_data, "buildings_osm")
path_building_height = os.path.join(path_main_data, "CH10")

path_result = "C:/Users/eggv/OneDrive - ZHAW/ZBP_shared/1_Research/1015_ZEB_China_Impact_Study/00_results"

# Chinese regions (Harbin, Beijing, Shanghai, Shenzen)
regions = [
    "Harbin",
    "Beijing",
    "Shanghai",
    "Shenzen"]

# Load regions
regions_1 = gpd.read_file(path_china_regions_1)
regions_2 = gpd.read_file(path_china_regions_2)

# Lookup regions for geometry. You need to explore whether to use regions_1 or regions_2 in the 
# geometry files for China
regions_selection = {
    "Harbin": regions_2.loc[regions_2['shapeName'] == 'Haerbinshi'],
    "Beijing": regions_1.loc[regions_1['shapeName'] == 'Beijing Municipality'],
    "Shanghai": regions_1.loc[regions_1['shapeName'] == 'Shanghai Municipality'],
    "Shenzen": regions_2.loc[regions_2['shapeName'].isin(['Shenzhenxian','Shenzhenshi'])]}


# ----------------
# (0) Preprocessing step
# In this step, different merging tifs are merged into one big tif. The reason for this is as 
# as some regions/cities contains several tifs. You need to check which tifs are to be merged and assign them to the
# Lookup-_tifs: C:\Users\eggv\OneDrive - ZHAW\ZBP_shared\1_Research\1015_ZEB_China_Impact_Study\02_Data\CH10
# ----------------
if pre_processing_merge:
    print("preproceeing merging tifs")
    
    # Look-up TIF building heigh with region
    # Note: If no the same CRS, need to be transformed first.
    lookup_tifs = {
        "Harbin": ["CNBH10m_X127Y45.tif"],
        "Beijing": ['CNBH10m_X117Y39.tif', 'CNBH10m_X115Y41.tif', 'CNBH10m_X115Y39.tif', 'CNBH10m_X117Y41.tif'],
        "Shanghai": ['CNBH10m_X121Y31.tif'],
        "Shenzen": ['CNBH10m_X115Y23.tif', 'CNBH10m-X113Y23_32650.tif', ]} # Not same CRS initioally, so transformed to 32650 both (32650, 32649 )

    for region in regions:
        print(" ... {}".format(region))
        tiff_names = lookup_tifs[region]
    
        # https://medium.com/spatial-data-science/how-to-mosaic-merge-raster-data-in-python-fb18e44f3c8
        if len(tiff_names) > 1:
            input_tiff_paths = [os.path.join(path_building_height, p) for p in tiff_names]
            output_tiff_path = os.path.join(path_building_height, '_merge_{}.tif'.format(region))

            hp.merge_tiffs(
                input_tiff_paths,
                output_tiff_path,
                nodata_value=np.nan)
        else:
            src = os.path.join(path_building_height, tiff_names[0])
            dst = os.path.join(path_building_height, '_merge_{}.tif'.format(region))
            shutil.copyfile(src, dst)
    raise Exception("finished preprocessing merging tifs")

# ----------------
# (1) Assign raster height
# ----------------

lookup_buildings = {
    "Harbin": os.path.join(path_osm_data, 'harbin.geojson'),
    "Beijing": os.path.join(path_osm_data, 'beijing.geojson'),
    "Shanghai": os.path.join(path_osm_data, 'shanghai.geojson'),
    "Shenzen": os.path.join(path_osm_data, 'shenzen.geojson')}

# Iterate regions to assing heights
for region in regions:
    print("assigning heights for region: {}".format(region))
    
    # Geometry of region
    #region_geometry = regions_selection[region] 

    # Get all buildings within region
    path_buildings = lookup_buildings[region]
    all_osm_region = gpd.read_file(path_buildings)

    # Get height projection
    path_heights = os.path.join(path_building_height, '_merge_{}.tif'.format(region))   
    with rasterio.open(path_heights) as src:
        crs_height = int(str(src.meta['crs']).split(":")[1])

    # Make same crs projection
    all_osm_region_new_crs = all_osm_region.to_crs(crs_height) #32652)

    shapes = hp.getFeatures(all_osm_region_new_crs)

    tiff_h = []
    max_hs = []
    archetypes = []

    with rasterio.open(path_heights) as src:
       
        # Iterate every building and get average height
        progress_par = Bar('reading shape: ', max=(len(shapes)))

        for counter, build_geom in enumerate(shapes):

            # Masking the tif
            try:
                out_image, out_transform = rasterio.mask.mask(src, [build_geom], crop=True)
                out_image[np.isnan(out_image)] = 0

                # Get all overlapping tif pixels
                out_image = out_image[out_image>0]  #Remove all zeros
                all_pixel_h_list = list(out_image.flatten())

                # Calculate the mediam of all tiff pixels
                mean_h = np.median(all_pixel_h_list)

                # Round building height
                mean_h_round = np.round(mean_h, 0)

                max_h = np.max(all_pixel_h_list)
    
                if np.isnan(mean_h_round):
                    mean_h_round = -9999
                
                if np.isnan(max_h):
                    max_h = -9999
            except:
                print("Building doesn't intersect")
                mean_h_round = -9999
                max_h = -9999           # Note: Dummy height values are asssigned

            #archetypes.append(archetpye)
            tiff_h.append(mean_h_round)
            max_hs.append(max_h)

            progress_par.next()
        progress_par.finish()
    all_osm_region['tiff_h'] = tiff_h
    all_osm_region['max_h'] = max_hs

    all_osm_region.to_file(path_buildings.replace('.geojson', '_h.geojson'))
    print("finished region: {}".format(region))

print("---- finished script ------")
