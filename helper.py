"""Helper functions
"""
import os
import random
import pandas as pd
import geopandas as gpd
import json

import haversine as hs   
from haversine import Unit

from shapely.geometry import Polygon

def age_classes_lookup():
    """Defined GWR age classes
    GBAUP
    """
    age_cat_classified = {
        0: [8011, 8012],                    # < 1945
        1: [8013],                          # 1945 - 1960
        2: [8014, 8015, 8016],              # 1961 - 1985
        3: [8017, 8018, 8019, 8020, 8021],  # 1986 - 2010
        4: [8022, 8023],                    # > 2010
        }
    timespans = {
        0: [1900, 1945],
        1: [1945, 1960],
        2: [1961, 1985],
        3: [1986, 2010],
        4: [2011, 2022]}

    age_classes_lookup = {}
    for key, values in age_cat_classified.items():
        for value in values:
            age_classes_lookup[value] = key

    return age_classes_lookup, timespans


def shapefile_as_cesarinput(
        gdf_shp,
        attr_def={
            'floor_level': 'levels_f',
            'height': 'height',
            'c_age': 'c_age',
            'c_type': 'c_type',
            'unique_id': 'unique_id'}):
    """Create shapefile as cesar input
    
    - Assign random year to building within the defined building age class
    """
  
    # -----------------------------------
    # Create CESAR building information
    # -----------------------------------
    info_dict_buildingtype = {
        'SFH':              'SFH',
        'MFH':              'MFH',
        'SHOP':             'SHOP',
        'OFFICE':           'OFFICE',
        'SCHOOL':           'SCHOOL',
        'RESTAURANT':       'RESTAURANT',
        'HOSPITAL':         'HOSPITAL',
        'OTHER':            'OTHER',
        'INDUSTRIAL':       'INDUSTRIAL'
        }

    # Middle year of range
    _, info_dict_age = age_classes_lookup()

    rows = []
    for index in gdf_shp.index:
        if gdf_shp.loc[index][attr_def['c_type']] in ['OTHER', 'INDUSTRIAL', 'none']:
            continue
        
        unique_id = int(gdf_shp.loc[index][attr_def['unique_id']])
        gdf_shp.at[index, 'TARGET_FID'] = unique_id
        building_type = gdf_shp.at[index, 'c_type']
        age_class = gdf_shp.at[index, 'c_age']
        BuildingType = info_dict_buildingtype[building_type]
        GroundFloorArea = round(gdf_shp.loc[index].geometry.area, 2)

        # Get random year for building age
        timespan = info_dict_age[age_class]
        BuildingAge = random.randint(timespan[0], timespan[1])

        nr_floors = int(gdf_shp.loc[index][attr_def['floor_level']])

        height = gdf_shp.loc[index][attr_def['height']]
        rows.append([
            unique_id,
            #EGID,
            BuildingType,
            BuildingAge,
            #LastRetrofit,
            GroundFloorArea,
            #ECarrierHeating,
            #ECarrierDHW,
            #GlazingRatio,
            nr_floors,
            height,
            height,
            unique_id,
            ])
    
    pd_out = pd.DataFrame(
        rows, 
        columns=[
            "ORIG_FID",
            #"EGID",
            "SIA2024BuildingType",
            "BuildingAge",
            #"LastRetrofit",
            "GroundFloorArea",
            #"ECarrierHeating",
            #"ECarrierDHW",
            #"GlazingRatio",
            "nr_floors",
            'height',
            'h_final',
            "unique_id"
            ])

    assert pd_out.ORIG_FID.is_unique

    return pd_out, gdf_shp

def read_sitevertices_from_shp(path_gdp, round_digits=2):
    """
    Read building shape information from shp and aggregate to a DataFrame.
    To each building a unique bld_id is assigned.

    Expected entries per row in csv, each representing one vertex of a building
    'gis_fid': fid identifying building in the external gis tool
    'height': height of building in meter
    'x': x coordinate of vertex, meter
    'y': y coordinate of vertex, meter

    :param file_path: full path to shp file
    :return: pandas DataFrame with one row for each building, columns being 'gis_fid', 'height', 'footprint_shape' and 'bld_id' as index.
             'footprint_shape' is a pandas DataFrame[columns=[x,y]] holding all building vertices
    """
    gdf_shp = gpd.read_file(path_gdp)

    required_keys = ["ORIG_FID", "height"]
    gdf_columns = gdf_shp.columns.tolist()
    for required_key in required_keys:
        assert required_key in gdf_columns, "Attribute: '{}' is missing in shapefile".format(required_key)

    container_list = []

    for building_index in gdf_shp.index:
        building_geometry = gdf_shp.loc[building_index].geometry
        target_fid = gdf_shp.loc[building_index]["ORIG_FID"]
        height = gdf_shp.loc[building_index]["height"]

        # Check if closed polygon
        if not building_geometry.boundary.is_ring:
            if building_geometry.interiors:

                # If inner ring, ignore it.
                building_geometry = Polygon(list(building_geometry.exterior.coords))
            if not building_geometry.boundary.is_ring:
                print(f"Polygon with target_fid {target_fid} is still not closed. Filling did not help. The building will now be skipped.")
                continue

        # Get boundary coordinates
        clockwise_coordinates = list(building_geometry.boundary.coords)

        # Iterate building vertices anti-clockwise
        for vertex in clockwise_coordinates:
            

            container_list.append([target_fid, round(vertex[0],round_digits ), round(vertex[1], 2), height])

    df_geometry = pd.DataFrame(container_list, columns=["TARGET_FID", "POINT_X", "POINT_Y", "HEIGHT"], dtype=float)

    return df_geometry

def calculate_metrics(gdf):
    """CAlculate width (m), and length (m) of polygon, as well as area (m2)
    """
    # Convert to metric coordinate system https://stackoverflow.com/questions/19412462/getting-distance-between-two-points-based-on-latitude-longitude
    lengths = []
    widths = []
    areas = []
    for i in gdf.index:
        geom = gdf.loc[i].geometry
        maxx = geom.bounds[2]
        minx = geom.bounds[0]
        maxy = geom.bounds[3]
        miny = geom.bounds[1]
        #widht = maxx - minx # maxx - minx
        #length = maxy - miny # maxy - miny

        length = hs.haversine((miny, maxx), (miny, minx), unit=Unit.METERS) # (lat, lon)
        width = hs.haversine((maxy, minx), (miny, minx), unit=Unit.METERS)

        area = length * width
        lengths.append(length)
        widths.append(width)
        areas.append(area)
    
    gdf['length_m'] = lengths
    gdf['width_m'] = widths
    gdf['footprint_m2'] = areas

    return gdf


def getFeatures(gdf):
    """Function to parse features from GeoDataFrame in such a manner that rasterio wants them"""
    return [json.loads(gdf.to_json())['features'][0]['geometry']