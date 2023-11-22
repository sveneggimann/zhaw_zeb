"""
Load chinese regions, corresponding openstreetmap data and assign building heights.

Info:
- merge tif: https://trac.osgeo.org/gdal/wiki/CatalogueForQIS

- masking tifs: https://rasterio.readthedocs.io/en/latest/topics/masking-by-shapefile.html 
"""

import os
import numpy as np
import geopandas as gpd
from matplotlib import pyplot
import matplotlib.pyplot as plt
import rasterio
from rasterio.plot import show
from rasterio.merge import merge as rio_merge
from rasterio.plot import show
import shutil
import fiona
from progress.bar import Bar

import helper as hp


test_path = "C:/Users/eggv/OneDrive - ZHAW/ZBP_shared/1_Research/1015_ZEB_China_Impact_Study/02_Data/CH10/test/small.tif"


'''src1 = rasterio.open("C:/Users/eggv/OneDrive - ZHAW/ZBP_shared/1_Research/1015_ZEB_China_Impact_Study/02_Data/CH10/test/test.tif")
src2 = rasterio.open("C:/Users/eggv/OneDrive - ZHAW/ZBP_shared/1_Research/1015_ZEB_China_Impact_Study/02_Data/CH10/test/test1.tif")
srcs_to_mosaic = [src1, src2]
arr, out_trans = rio_merge(srcs_to_mosaic)
arr[np.isnan(arr)] = 0'''

pre_processing_merge = False # Only once

# Paths
path_main_data = "C:/Users/eggv/OneDrive - ZHAW/ZBP_shared/1_Research/1015_ZEB_China_Impact_Study/02_Data"

path_china_regions_1 = os.path.join(path_main_data, "boundaries", "geoBoundaries-CHN-ADM1.geojson")
path_china_regions_2 = os.path.join(path_main_data, "boundaries", "geoBoundaries-CHN-ADM2.geojson")
path_osm_data = os.path.join(path_main_data, "buildings_osm") # , "gis_osm_buildings_a_free_1.shp")
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

# Lookup regions for geometry
regions_selection = {
    "Harbin": regions_2.loc[regions_2['shapeName'] == 'Haerbinshi'],
    "Beijing": regions_1.loc[regions_1['shapeName'] == 'Beijing Municipality'],
    "Shanghai": regions_1.loc[regions_1['shapeName'] == 'Shanghai Municipality'],
    "Shenzen": regions_2.loc[regions_2['shapeName'] == 'Shenzhenxian']}

# Look-up TIF building heigh with region (TODO: Merge TIFS)
lookup_tifs = {
    "Harbin": ["CNBH10m_X127Y45.tif"],
    "Beijing": ['CNBH10m_X117Y39.tif', 'CNBH10m_X115Y41.tif', 'CNBH10m_X115Y39.tif', 'CNBH10m_X117Y41.tif'],
    "Shanghai": ['CNBH10m_X121Y31.tif'],
    "Shenzen": ['CNBH10m_X113Y23.tif', 'CNBH10m_X115Y23.tif'],
}

# ----------------
# (0) Preprocessing step, merging tifs
# https://automating-gis-processes.github.io/CSC18/lessons/L6/raster-mosaic.html
# ----------------
if pre_processing_merge:
    print("preproceeing merging tifs")
    for region in regions:
        print(" ... {}".format(region))

        tiff_names = lookup_tifs[region]
    
        # https://medium.com/spatial-data-science/how-to-mosaic-merge-raster-data-in-python-fb18e44f3c8
        if len(tiff_names) > 1:
            raster_to_mosiac = []
            for p in tiff_names:
                path_tif = os.path.join(path_building_height, p)
                raster = rasterio.open(path_tif)
                raster_to_mosiac.append(raster)
            
            mosaic, output = rio_merge(
                raster_to_mosiac,
                nodata=0)
            
            # Replace all nan values with 0
            mosaic[np.isnan(mosaic)] = 0

            output_path = os.path.join(path_building_height, '_merge_{}.tif'.format(region))
            
            output_meta = raster.meta.copy()
            output_meta.update(
                {"driver": "GTiff",
                    "height": mosaic.shape[1],
                    "width": mosaic.shape[2],
                    "transform": output,
                })

            with rasterio.open(output_path, 'w', **output_meta) as m:
                m.write(mosaic)
        else:
            src = os.path.join(path_building_height, tiff_names[0])
            dst = os.path.join(path_building_height, '_merge_{}'.format(region))
            shutil.copyfile(src, dst)
        print("finished preprocessing merging tifs")


lookup_buildings = {
    "Harbin": os.path.join(path_osm_data, 'harbin.geojson'),
    "Beijing": os.path.join(path_osm_data, 'beijing.geojson'),
    "Shanghai": os.path.join(path_osm_data, 'shangahi.geojson'),
    "Shenzen": os.path.join(path_osm_data, 'shenzen.geojson')}

# Iterate regions to assing heights
for region in regions:
    print("assigning heights for region: {}".format(region))
    
    # Geometry of region
    region_geometry = regions_selection[region] 

    # Get all buildings within region
    path_buildings = lookup_buildings[region]
    all_osm_region = gpd.read_file(path_buildings)

    # Make same crs projection
    all_osm_region_new_crs = all_osm_region.to_crs(32650)

    shapes = hp.getFeatures(all_osm_region_new_crs)

    #with fiona.open(path_buildings, "r") as shapefile:
    #    shapes = [feature["geometry"] for feature in shapefile]

    # Heights
    path_heights = os.path.join(path_building_height, '_merge_{}.tif'.format(region))   
    
    tiff_h = []

    with rasterio.open(path_heights) as src:
       
        # Iterate every building and get average height
        #for i in all_osm_region.index:
        #    build_geom = all_osm_region.loc[i].geometry
        progress_par = Bar('reading shape: ', max=(len(shapes)))

        for build_geom in shapes:

            # Masking the tif
            out_image, out_transform = rasterio.mask.mask(src, [build_geom], crop=True)
            out_image[np.isnan(out_image)] = 0

            # Get all overlapping tif pixels
            all_pixel_h = out_image.remove(0)

            # Calculate average of all tiff pixels
            mean_h = np.mean(all_pixel_h)

            tiff_h.append(mean_h)
            progress_par.next()
        progress_par.finish()


    all_osm_region['tiff_h'] = tiff_h

    all_osm_region.to_file(all_osm_region.replace('.geojson', '_h.geojson'))

    



    
    '''
from osgeo import gdal
filepath = test_path

# Open the file: https://automating-gis-processes.github.io/2016/Lesson7-read-raster.html
raster = gdal.Open(filepath)

# Check type of the variable 'raster'
type(raster)

print(raster.GetProjection())
print(raster.RasterXSize)
print(raster.RasterYSize)
print(raster.RasterCount)
print(raster.GetMetadata())
band = raster.GetRasterBand(1)
print(type(band))
print(gdal.GetDataTypeName(band.DataType))

# https://automating-gis-processes.github.io/2016/Lesson7-read-raster-array.html
rasterArray = raster.ReadAsArray()
nodata = band.GetNoDataValue()
rasterArray = np.ma.masked_equal(rasterArray, nodata)
type(rasterArray)

# Check again array statistics
rasterArray.min()
'''

##from osgeo import gdal
#dataset = gdal.Open(test_path, gdal.GA_ReadOnly)
#for x in range(1, dataset.RasterCount + 1):
#    band = dataset.GetRasterBand(x)
#    array = band.ReadAsArray()
#with rasterio.open(test_path, mode="r") as src:
#    data_masked = src.read(masked=True)
#    nodata_value = src.nodata
#data = data_masked.filled(nodata_value) 

#src = rasterio.open(test_path)
#msk = src.read_masks(1)