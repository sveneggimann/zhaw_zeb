"""Download and extrct amtliche Data vermessung for City of St Gallen and city of Zürich

Info
-----
AV files to download: https://data.geo.admin.ch/ch.swisstopo-vd.amtliche-vermessung/meta.txt

Converting interlis to shapfiles is done with GDAL (ogr2ogr). In order
to run the interlis part, the file needs to be executed with the 
OSGeo4W64 shell. Set "mode to 'OSGeo4W64

set CONDA_DLL_SEARCH_MODIFICATION_ENABLE=1
e.g. python3 ../../scripts/create_parzellen.py

works: ogr2ogr -f "ESRI Shapefile" "C:/users/eggv/_temp/_temp/extract" "C:/Users/eggv/_temp/_temp/de/0100.itf,C:/Users/eggv/_temp/_temp/de/DM01AVCH24LV95D.ili"
not  : ogr2ogr -f "ESRI Shapefile" "C:/users/eggv/_temp/_temp/extract" "C:/Users/eggv/_temp/_temp/de/0100.itf,C:/Users/eggv/_temp/_temp/de/DM01AVCH24LV95D.ili"
"""
mode = "normal" #'normal' #'OSGeo4W64'
import os, sys
import subprocess
import requests
import urllib
from urllib.request import urlretrieve
import zipfile
from shutil import copyfile
import configparser

subprocess.call('dir', shell=True)

#path_staVerdi = os.path.abspath(os.path.join(os.path.dirname(__file__), '..', '..', ))
#sys.path.append(#path_staVerdi)
   # 
if mode == 'normal':
    from shapely import wkt
    from shapely.ops import linemerge, unary_union, polygonize
    import geopandas as gpd
    import pandas as pd
    from shapely.geometry import LineString, Polygon, Point
    import helper as hp
    #import staVerdi.helper.helper_functions as hp




# -------------------------
# Download raw data
# -------------------------
# Path to config file
path_ini_file = 'C:/Users/eggv/OneDrive - ZHAW/Sven/00_code/zhaw_zeb/stess/config_open_av.ini'

# Read configuration
config = configparser.ConfigParser()
config.read(path_ini_file)

target_crs = '21781'

download_raw = False
extract_fitting = True
polygonzie = False

merge_and_street_clip = False

cantons_to_pick = ['SG', 'ZH']

if download_raw:

    hp.create_folder(config['PATHS']['path_out_folder'])

    #data = urllib2.urlopen(target_url) # it's a file like object and works just like a file
    response = requests.get(config['PATHS']['target_url']) 
    data = response.text.split("\n")
    for canton_to_pick in cantons_to_pick:

        for line in data: # files are iterable

            if canton_to_pick in line:
                line_split = line.split(" ")
                url = line_split[0]
                file_name = url.split("/")[-1]

                path_save_file = os.path.join(config['PATHS']['path_out_folder'], file_name)

                if not os.path.exists(path_save_file):
                    print('Downloading file {}'.format(file_name))
                    # Copy a network object to a local file
                    urlretrieve(url, path_save_file)
                else:
                    print(" .... {} already downloaded".format(file_name))

    print("... finished downloading file")

# -------------------------
# Extract fitting data
# -------------------------
if extract_fitting:
    extract_interlis = True
    hp.create_folder(config['PATHS']['path_zipped'])

    for file_name in os.listdir(config['PATHS']['path_out_folder']):

        print('Extracting file {}'.format(file_name))
        path_zip = os.path.join(config['PATHS']['path_out_folder'], file_name)
        path_temp = os.path.join(config['PATHS']['path_zipped'], '_temp')
        path_out = os.path.join(config['PATHS']['path_zipped'], file_name[:-4])

        # Files which could not be extracted. Take old datafiles
        #files_with_problems = [6004, 6010, 6076, 6110, 6119, 6142, 6202, 6235, 6248]

        #if int(file_name[:-4]) in files_with_problems:
        #    continue

        if not os.path.exists(path_out):
            hp.delete_folder(path_temp) #DELETE FOLDER
            hp.create_folder(path_temp)
            hp.create_folder(path_out)

            with zipfile.ZipFile(path_zip, 'r') as zip_ref:
                zip_ref.extractall(path_temp)

            # Extract interlis file
            if extract_interlis:
                path_temp_interlis = os.path.join(path_temp, "extract")
                hp.create_folder(path_temp_interlis)


                # Get itf file
                sub_files = os.listdir(os.path.join(path_temp, "de"))
                
                for f in sub_files:
                    if f.endswith('itf'):
                        path_itf_file = os.path.join(os.path.join(path_temp, "de", f))
                #path_itf_file = os.path.join(path_temp, "de", "000{}.itf".format(file_name[:-4]))
                    if f.endswith("ili"):
                        path_ili_file = os.path.join(path_temp, "de", f)
                #path_ili_file = os.path.join(path_temp, "de", "DM01AVCH24LV95D.ili")
                print("----")
                print(path_temp_interlis)
                print(path_itf_file)
                print(path_ili_file)

                path_temp_interlis = path_temp_interlis.replace("\\","/")
                path_itf_file = path_itf_file.replace("\\","/")
                path_ili_file = path_ili_file.replace("\\","/")
                
                path_temp_interlis = path_temp_interlis.replace("/","\\")
                path_itf_file = path_itf_file.replace("/","\\")
                path_ili_file = path_ili_file.replace("/","\\")

                #command = [
                #    'ogr2ogr',
                #    "-f",
                #    '"ESRI Shapefile"',
                #    path_temp_interlis,
                #    path_itf_file,
                #    path_ili_file]
                #Info: https://giswiki.hsr.ch/HowTo_OGR2OGR#INTERLIS_2-Reader_und_-Writer --Note: See interlis info 1
                args = [
                    'ogr2ogr',
                    '-f',
                    'ESRI Shapefile',
                    'C:\\Users\\eggv\\_temp\\_temp\\extract',
                    '{},{}'.format(
                        path_itf_file,
                        path_ili_file)]
                
                #args2 = [
                #    'ogr2ogr',
                #    '-f',
                #    '"ESRI Shapefile"',
                #    'C:/Users/eggv/_temp/_temp/extract',
                #    "C:/Users/eggv/_temp/_temp/de/0100.itf,C:/Users/eggv/_temp/_temp/de/DM01AVCH24LV95D.ili"]

                #t = 'ogr2ogr -f "ESRI Shapefile" "C:/Users/eggv/_temp/_temp/extract" "C:/Users/eggv/_temp/_temp/de/0100.itf,C:/Users/eggv/_temp/_temp/de/DM01AVCH24LV95D.ili"'

                #befehl = 'ogr2ogr -f "ESRI Shapefile" C:\\Users\\eggv\\_temp\\_temp\\extract {},{}'.format(path_itf_file, path_ili_file)
                #subprocess.run(command, shell=True, capture_output=False, text=True)
                #subprocess.check_call(command2)
                #args = ["{}/polygonize.sh".format(config['PATHS']['path_script']), INPUT, OUTPUT, config['PATHS']['temp_path']]
                #print("----1")
                #exec(befehl)
                print("----2")
                subprocess.run(args)
                #print("----3")
                #subprocess.call(args)
                prnt("fff")

            # Files to keep Liegenschaften__Liegenschaft_Geometrie.shp
            files_to_keep = [
                os.path.join("Liegenschaften__Liegenschaft_Geometrie.dbf"),
                os.path.join("Liegenschaften__Liegenschaft_Geometrie.shp"),
                os.path.join("Liegenschaften__Liegenschaft_Geometrie.shx")]

            for file_path in files_to_keep:
                split_names = file_path.split("\\")
                old_path = os.path.join(path_temp_interlis, file_path)
                new_path = os.path.join(path_out, split_names[-1])
                copyfile(old_path, new_path)
        else:
            print("... already extracted {}".format(path_out))
  
        #hp.delete_folder(path_temp)

    print("... finished unzipping fitting files")

# --------------------
# Conver to polygons
# --------------------
if polygonzie:

    hp.create_folder(config['PATHS']['path_polygonized'])
    for folder in os.listdir(config['PATHS']['path_zipped']):

        print("Folder to polygonize: " + str(folder))
        out_path = os.path.join(config['PATHS']['path_polygonized'], "{}.shp".format(folder))
        if os.path.exists(out_path):
            continue
        else:
            file_name = 'Liegenschaften__Liegenschaft_Geometrie.shp'

            gdf = gpd.read_file(os.path.join(config['PATHS']['path_zipped'], folder, file_name))
            geometries = []
            gdf_intersecting = unary_union(gdf.geometry)

            if gdf_intersecting.type == 'MultiLineString':
                gdf_intersecting = linemerge(gdf_intersecting)
                
                if gdf_intersecting.type == 'LineString':
                    geometries.append(gdf_intersecting)
                else:
                    for line in gdf_intersecting:
                        geometries.append(wkt.loads(line.wkt))

            if gdf_intersecting.type == 'LineString':
                geometries.append(gdf_intersecting)

            merged = linemerge(geometries)
            borders = unary_union(merged)
            all_polygons = list(polygonize(borders))

            # Write to shapefile
            parzellen = pd.DataFrame(all_polygons, columns=['geometry'])
            
            parzellen = gpd.GeoDataFrame(
                parzellen, crs={'init': 'epsg:{}'.format('2056')})

            parzellen = parzellen.to_crs({'init': 'epsg:{}'.format(target_crs)})
            
            parzellen.to_file(out_path)

    print("... finished polygonize")


# --------------------
# Merge and clip streets
# --------------------
if merge_and_street_clip:
    try_clipping_strees = False
                
    hp.create_folder(config['PATHS']['parzellen_path_merged'])

    # Sliver Cirterium
    sliver_crit = 5 # Area / perimeter

    # Get only major street
    ##street_network_full = street_network_full[street_network_full['OBJEKTART'].isin()]
    if try_clipping_strees:      

        print("Load street network")
        street_network_full = gpd.read_file(config['PATHS']['street_path'])

        print("Create retree for street network")
        street_rtree = hp.build_rTree(street_network_full)

    # Merge gdf
    merged_parzellen = gpd.GeoDataFrame()

    for file_name in os.listdir(config['PATHS']['path_polygonized']):
        if file_name.endswith('.shp'):
            print("Merging bfs folder: {} length: {}".format(file_name, merged_parzellen.shape[0]))
            gdf = gpd.read_file(os.path.join(config['PATHS']['path_polygonized'], file_name))

            # ------------------------------
            # Calculate sliver polygon criteria to remove only all intersecting polygons which are
            # streets (long and narrow)
            # ------------------------------
            gdf['sliver_r'] = gdf.geometry.area / gdf.geometry.length

            # ------------------------------
            # Get intersecting polygons only
            # ------------------------------
            if try_clipping_strees:
                index_to_delete = set()
                for poly_index in gdf.index:

                    streets_in_tree = street_rtree.intersection(
                        gdf.loc[poly_index].geometry.bounds)
                    
                    for street_index in streets_in_tree:
                        line = street_network_full.loc[street_index].geometry
                        poly = gdf.loc[poly_index].geometry

                        if poly.intersects(line):

                            # Only delete if sliver crit is fulfilled
                            if gdf.loc[poly_index].sliver_r <= sliver_crit:
                                index_to_delete.add(poly_index)
                
                index_to_delete = list(index_to_delete)

                # Drop intersecting polygons
                gdf = gdf.drop(index=index_to_delete)

            # Merge
            merged_parzellen = merged_parzellen.append(gdf)

    # Reset index
    merged_parzellen = merged_parzellen.reset_index(drop=True)

    merged_parzellen.to_file(os.path.join(config['PATHS']['parzellen_path_merged'], "parzellen_merged.shp"))