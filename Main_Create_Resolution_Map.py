# To add a new cell, type '# %%'
# To add a new markdown cell, type '# %% [markdown]'
# %% [markdown]
# # Example on how to create resolution maps
# 
# @author: Max Felius
# 
# This Notebook will show how to combine the different MRSS datasets in order to create a complete point density map.

# %%
#imports
from decouple import config
import CreateResolutionMap
import os

# %% [markdown]
# ## Select the datasets to be combined
#
# Make sure to configure the .env file with the correct path.

# %%
absolute_path = config('abs_path_downloads')

#datasets used to combine
datasets = [
    'full-pixel_mrss_rsat2_asc_t109_v3_b728bb51760435c769a9d45b3e8e74d6325b90fa.csv',
    'full-pixel_mrss_rsat2_dsc_t302_v3_6cb798ad333fdf7ad686b31443834e676ea5a0bc.csv'#,
    #'full-pixel_mrss_s1_asc_t88_v4_080a1cbf7de1b6d42b3465772d9065fe7115d4bf.csv',
    #'full-pixel_mrss_s1_dsc_t37_v5_6b2ebff998208d769573c02548be69075447f465.csv',
    #'full-pixel_mrss_s1_dsc_t139_v5_47cf4b28795f4dcd36525ceb1e4741a4e9809ac0.csv',
    #'limburg_xf_asc_v8_noqc_psds_hl.csv',
    #'limburg_xf_dsc_v6_psds_hl.csv'
]

#add the absolute path to the datasets
for idx,data in enumerate(datasets):
    datasets[idx] = os.path.join(absolute_path,data)

# %% [markdown]
# ## Create the resolution map for looking at the point density
# 
# Map is made using the following procedure.
# - First the data is read with the use of pandas and appended to eachother in order to make one big data file
# - The maximum and minimum values of the rijksdriehoek columns (pnt_rdx, pnt_rdy) are taken from the data file
# - A grid is created using the spacing of 2 times the 'cell_radius' in the x and y direction
# - Using the rijksdriehoek coordinates (pnt_rdx, pnt_rdy), the KD-Tree is initialized
# - Per each grid point a subset is created using the KD-Tree. The KD-Tree can only work with search <b>radius</b>, thus here the search radius is defined as $\sqrt{2} \cdot cell\_radius \cdot 1.01$. This way all points that fall within the search square will be used.
# - A simple if-statement defines a cell and all the points that fall within this cell are counted.


# %%
#cell radius is the distance between the center of the square to the edge. Resolution is 2*cell_radius
cell_radius = 100

#creating the data object
dataset_obj = CreateResolutionMap.create_resolution_map(datasets,cell_radius,'combined_dataset_rsat.csv')

#start creating the map
data = dataset_obj.start()

#save the dataset
data.to_csv('Resolution_Map_Combined_Dataset_RSAT_200x200.csv')


# %% [markdown]
# ## Sentinel 1 track map

datasets = [
	'full-pixel_mrss_s1_asc_t88_v4_080a1cbf7de1b6d42b3465772d9065fe7115d4bf.csv',
    'full-pixel_mrss_s1_dsc_t37_v5_6b2ebff998208d769573c02548be69075447f465.csv',
    'full-pixel_mrss_s1_dsc_t139_v5_47cf4b28795f4dcd36525ceb1e4741a4e9809ac0.csv'
]

#cell radius is the distance between the center of the square to the edge. Resolution is 2*cell_radius
cell_radius = 100

#creating the data object
dataset_obj = CreateResolutionMap.create_resolution_map(datasets,cell_radius,'combined_dataset_s1.csv')

#start creating the map
data = dataset_obj.start()

#save the dataset
data.to_csv('Resolution_Map_Combined_Dataset_s1_200x200.csv')
