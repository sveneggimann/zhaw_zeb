"""

Not used

"""

# Google scholar = 3857

import os
import pprint
import geopandas as gpd
print("-load blocks--")
google_scholar_crs = 3857



blocks = gpd.read_file("C:/Users/eggv/OneDrive - ZHAW/ZBP_shared/1_Research/1015_ZEB_China_Impact_Study/02_Data/blocks_data/DT41/3d_urban_form_of_chinese_cities.shp")


#blocks_new = blocks.to_crs(google_scholar_crs)
#blocks_new.to_file("C:/Users/eggv/OneDrive - ZHAW/ZBP_shared/1_Research/1015_ZEB_China_Impact_Study/02_Data/blocks_data/DT41/transforrmed.shp")


# Filter out all with less than five buildings, cleaning steps
blocks = blocks.loc[blocks['BLD_Count'] >= 5]
blocks = blocks.loc[blocks['BLD_FAR'] >= 0.2]
blocks.to_file("C:/Users/eggv/OneDrive - ZHAW/ZBP_shared/1_Research/1015_ZEB_China_Impact_Study/02_Data/blocks_data/DT41/blocks_selection.shp")

# Inward & outwards buffer
# Remove sliver polygon criteria

# Plots
stats = blocks[['FORM_TYPE', 'BLD_FLOORN', 'BLD_FAR']].groupby(["FORM_TYPE"]).mean()

pprint.pprint(blocks.crs)
blocks.drop(columns=['geometry']).to_csv("C:/Users/eggv/OneDrive - ZHAW/ZBP_shared/1_Research/1015_ZEB_China_Impact_Study/02_Data/blocks_data/DT41/data.csv")

print("-finished-")