'''

Create Statistics of library of neighbourhoods

'''

# Read all blocks and buildings

import sys
import os
import pyproj
import geopandas as gpd
import shapely
import math
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

path_result = "C:/Users/eggv/OneDrive - ZHAW/ZBP_shared/1_Research/1015_ZEB_China_Impact_Study/00_results"

case_studies = os.listdir(path_result)

assumed_height_floor = 3 # [m]

# Region, Block, Building_numbers, FAR, 
result_stats = []
result_columns = [
    'region', 'block_id', 'far', 'floor_area', 'build_nrs', 'mean_height', 'nr_terrace_house', 'nr_low_rise_apartment',
    'nr_high_rise_slab', 'nr_high_rise_tower', 'nr_non_residential']

for case_study in case_studies:
    print("case_study: {}".format(case_study))
    path_case_study = os.path.join(path_result, case_study)

    blocks = os.listdir(path_case_study)

    for block in blocks:
        path_block = os.path.join(path_case_study, block)

        path_block_geom = os.path.join(path_block, 'block.shp')
        path_build_blocks = os.path.join(path_block, 'buildings_block.shp')

        block_geom = gpd.read_file(path_block_geom)
        build_blocks = gpd.read_file(path_build_blocks)

        # Replace all -9999 values with one minimum height of 3 m
        index_faulty_height = build_blocks.loc[build_blocks['tiff_h'] == -9999].index
        if len(index_faulty_height) > 0:
            build_blocks.loc[index_faulty_height, 'tiff_h'] = assumed_height_floor


        # Reproject to metric coodrinate sytsem
        block_geom_metric_reprojection = block_geom.to_crs(7564)  # metric crs for China
        build_blocks_reprojection = build_blocks.to_crs(7564)  # metric crs for China
        
        # Calculate floor level (round-up)
        build_blocks_reprojection['floor_level'] = build_blocks_reprojection['tiff_h'] / assumed_height_floor
        build_blocks_reprojection['floor_level'] = build_blocks_reprojection['floor_level'].apply(np.ceil)
    
        # Stats
        build_nrs = build_blocks_reprojection.shape[0]     # nr of buildings
        block_geom_metric_reprojection = block_geom_metric_reprojection.geometry.area[0]     # area block
        floor_area = np.sum(build_blocks_reprojection.geometry.area  * build_blocks_reprojection['floor_level'])  #[m2]
        far = floor_area / block_geom_metric_reprojection
        mean_floor_level = np.mean(build_blocks_reprojection['floor_level'])
        mean_height = np.mean(build_blocks_reprojection['tiff_h'])

        nr_archetype1 = build_blocks.loc[build_blocks['archetype'] == 'terrace_house'].shape[0]
        nr_archetype2 = build_blocks.loc[build_blocks['archetype'] == 'low_rise_apartment'].shape[0]
        nr_archetyp3 = build_blocks.loc[build_blocks['archetype'] == 'high_rise_slab'].shape[0]
        nr_archetype4 = build_blocks.loc[build_blocks['archetype'] == 'high_rise_tower'].shape[0]
        nr_archetype5 = build_blocks.loc[build_blocks['archetype'] == 'non_residential'].shape[0]

        stats_block = [
            case_study,
            block,
            far,
            floor_area,
            build_nrs,
            mean_height,
            nr_archetype1,
            nr_archetype2,
            nr_archetyp3,
            nr_archetype4,
            nr_archetype5
            ]
        result_stats.append(stats_block)


pd_results = pd.DataFrame(result_stats, columns=result_columns)

print("Pandas resutls")
print(pd_results)


# --Plot share of archetypes
stats = pd_results.groupby('region').sum()
stats = stats[['build_nrs', 'nr_terrace_house', 'nr_low_rise_apartment', 'nr_high_rise_slab', 'nr_high_rise_tower', 'nr_non_residential']]

# Calculate percentage
for index in stats.index:
    stats.loc[index] = (100 / stats.loc[index]['build_nrs']) * stats.loc[index]
stats['region'] = stats.index
stats = stats.drop(columns=['build_nrs'])
ax = stats.plot( 
    x = 'region', 
    kind = 'barh', 
    stacked = True, 
    title = 'Distribution of archetypes', 
    mark_right = True)

plt.ylabel("Region")
plt.xlabel("Archetype (as percentage of all buildings)")

plt.show()


# --PLot individual number of buildings

# --PLot number of blocks
block_nr = pd_results.groupby('region').count()  # Count number of blocks
block_nr['region'] = block_nr.index
ax = block_nr.plot.bar(x='region', y='block_id', rot=0)
plt.xlabel("Region")
plt.ylabel("Number of sampled blocks")
plt.legend('',frameon=False)
plt.show()

# --Plot FAR as boxplots
pd_results.boxplot("far", by="region")
plt.ylim(0,)
plt.xlabel("Region")
plt.ylabel("Floor Area Ratio (FAR)")
plt.legend('',frameon=False)
plt.show()



print("--- finished script ----")


        