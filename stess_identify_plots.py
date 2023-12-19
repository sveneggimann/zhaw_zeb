"""

Quick and dirty approach to identify suitable land plots


"""
import geopandas as gpd
from progress.bar import Bar
from scipy.spatial import KDTree
import helper as hp

# TODO: Integrate Bauzonen criteria (not yet done)


# ---------
# Criteria
# ---------
#Note: Criteria are not applied yet to filter, but oculd be easily aded at the end (the different plots can also be visualized in QIGS by showing the selection (see data_overview.qgz file))
crit_min_area = 3000        # Minimum area [m2]
crit_max_nr_build = 1       # Maximum number of buildings on plot [nr]
sliver_assumption = 2.5     # REatio to remove odd long and thin polygons, i.e. streets [-]
footprint_coverage_f = 50   # Max covered by building footprint area [%]
# TODO: add other data lyeras and other criteria

# Paths
path_plots = "C:/Users/eggv/OneDrive - ZHAW/ZBP_shared/1_Research/1016_SwissSTES/00_sven_geometries/av_data/av_data_clipped.shp"
path_bauzonen_buffer = "C:/Users/eggv/OneDrive - ZHAW/ZBP_shared/1_Research/1016_SwissSTES/00_sven_geometries/bauzonen_buffer_200m.shp"
path_buildings = "C:/Users/eggv/OneDrive - ZHAW/ZBP_shared/1_Research/1016_SwissSTES/00_sven_geometries/street/swissTLM3d_buildings_project21781.shp"
path_building_assigned_plot = "C:/Users/eggv/OneDrive - ZHAW/ZBP_shared/1_Research/1016_SwissSTES/00_sven_geometries/street/swissTLM3d_buildings_with_plot.shp"
path_result = "C:/Users/eggv/OneDrive - ZHAW/ZBP_shared/1_Research/1016_SwissSTES/00_sven_geometries/plot_selection.shp"

print("Load buildings and convert to dots")
plots = gpd.read_file(path_plots)

plots['unique_id'] = range(plots.shape[0])

# Remove obvious slivers
plots = plots.loc[plots['sliver_r'] > sliver_assumption]

plots['area'] = plots.geometry.area
plots['build_area_ftp'] = 0
plots['build_nrs'] = 0
plots['build_200m'] = 0

# Reduce selection of plots
plots_min_size = plots.loc[plots['area'] > crit_min_area]


print("Read buildings")
buildings = gpd.read_file(path_buildings)
buildings['area_ftp'] = buildings.area          # Individual buildling area [m2]


# Search tree for buildings
print("    reating kd search tree for buildings")
x_coord = buildings.geometry.centroid.x.tolist()
y_coord = buildings.geometry.centroid.y.tolist()
build_kd_tree = KDTree([list(a) for a in zip(x_coord, y_coord)])
build_rTree = hp.build_rTree(buildings)


cont_plot_id = []
progress_par = Bar('Assigning buildings to plots:', max=(plots.shape[0]))

for plot_index in plots.index:
    plot_id = plots.loc[plot_index]['unique_id']
    plot_geom = plots.loc[plot_index].geometry

    # Get all buildings within block and write out
    tree_intersec = build_rTree.intersection(plot_geom.bounds)
    index_to_take_out = []
    sum_area = 0
    for tree_index in tree_intersec:
        if plot_geom.contains(buildings.iloc[tree_index].geometry.centroid):
            index_to_take_out.append(tree_index)
            sum_area += buildings.iloc[tree_index]['area_ftp']
    
    # TODO: Get number of buildings within X distance
            


    plots.loc[plot_index, 'build_area_ftp'] = sum_area
    plots.loc[plot_index, 'build_nrs'] = len(index_to_take_out)

    progress_par.next()
progress_par.finish()

# Calculate ftp coverage factor
plots['ftp_coverage_f'] = (plots['build_area_ftp'] / plots['area']) * 100

#plots = plots.loc[plots['footprint_coverage_f'] > footprint_coverage_f]

plots.to_file(path_result)


print(" --- finished ---")

