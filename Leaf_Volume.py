import numpy as np
import matplotlib.pyplot as plt

import open3d as o3d
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

from mpl_toolkits.mplot3d import Axes3D
from scipy.spatial import ConvexHull
from skspatial.plotting import plot_3d

def compute_height(voxel_df):
    """
    Compute the height of a tree
    
    Parameters
    ----------------
    voxel_df: 
        dataframe containing x, y, and z voxel coordinates. Must be in 
        this exact order (x,y,z)
        
    Returns
    ----------------
    tree height 
    """
    return voxel_df[2].max() - voxel_df[2].min()

def calc_convex_hull(sliced_vox, voxel_size, visualize=False):
    """
    Perform convex hull to find contact frequency
    
    Parameters
    ----------------
    sliced_vox:
        voxelized leaf dataframe slices
    voxel_size:
        size used to voxelize tree (0.0325 m is default)
    visualize:
        True or False
        visualize the convex hull 
        
    Returns
    ----------------
    sum of leaf volumes
    
    """ 
 
    leaf_vols = []

    for vox_con in sliced_vox:
        vols = ConvexHull(vox_con).volume # units are in m3 

        if visualize:
            fig = plt.figure()
            ax = fig.add_subplot(111, projection = "3d")

            ax.plot(vox_con.T[0], vox_con.T[1], vox_con.T[2], "gs")

            hull = ConvexHull(vox_con)
            for s in hull.simplices:
                s = np.append(s, s[0])
                ax.plot(vox_con[s, 0], vox_con[s, 1], vox_con[s, 2], "r-")

            for i in ["x", "y", "z"]:
                eval("ax.set_{:s}label('{:s}')".format(i,i))

            plt.show()

        leaf_vols.append(vols)
        
    return leaf_vols

def LeafVolume(file_path, visualize=False, voxel_size = 0.0325):
    """
    voxelize and slice trees to estimate their leaf volumes 
    
    Parameters
    ----------------
    file_path:
        the location of the point cloud data files
    voxel_size:
        size used to voxelize tree (0.0325 m is default)
    visualize:
        True or False
        visualize slicing and convex hull
        
    Returns
    ----------------
    sum of leaf volumes
    """  
   
    # read point cloud
    pcd = o3d.io.read_point_cloud(file_path) 
    
    # voxelize
    downpcd = pcd.voxel_down_sample(voxel_size = voxel_size)
    
    # generate an array of coordinates for voxelized point cloud
    df_vox = pd.DataFrame(downpcd.points)
    
    # normal calculation
    downpcd.estimate_normals(search_param = o3d.geometry.KDTreeSearchParamHybrid(radius = 0.1, max_nn =30))
    df_norm = pd.DataFrame(downpcd.normals)

    # consolidate voxels and normals into a single dataframe
    sorted_df = pd.concat([df_vox, df_norm], axis = 1)
    sorted_df.columns = ["x", "y", "z", "norm_x", "norm_y", "norm_z"]
    
    # and sort by elevation
    sorted_df = sorted_df.sort_values(by = ['z'])
    
    #separate the now sorted voxels and normals 
    voxs = sorted_df.iloc[:,:-3]
    norms = sorted_df.iloc[:,-3:]
    
    # compute tree height
    slice_width = 0.5
    h = compute_height(df_vox)

    # how many half meter slices, round up  
    num_slices = np.ceil(h/slice_width)
#     print(f"how many half meter slices: {num_slices}")
    
    # slice the arrays 
    sliced_norms = np.array_split(norms.to_numpy(), num_slices) 
    sliced_vox = np.array_split(voxs.to_numpy(), num_slices)

    if visualize:
        fig = plt.figure()
        ax = fig.add_subplot(111, projection = "3d")

        for vox_slice in sliced_vox:
            ax.scatter(vox_slice[:,0], vox_slice[:,1], vox_slice[:,2])
        plt.show()
        
    # convex hull
    leaf_vols = calc_convex_hull(sliced_vox,voxel_size,visualize)
    
    return leaf_vols

def generate_df(file_name_list, value_list):
    """
    voxelize and slice trees to estimate their leaf volumes 
    
    Parameters
    ----------------
    file_path:
        the location of the point cloud data files
    voxel_size:
        size used to voxelize tree (0.0325 m is default)
    visualize:
        True or False
        visualize slicing and convex hull
        
    Returns
    ----------------
    file containing leaf volumes for plots containing all tree species or just one species
    """ 
    values_df = pd.DataFrame(value_list, columns = ['leaf volumes'])
    
    file_name_df = pd.DataFrame(
        {
            'plot_name':file_name_list,
        }
    )
    file_name_df = file_name_df['plot_name'].str.split('_',n=3,expand=True)
    file_name_df = file_name_df.drop(3, axis=1)
                                    
    file_name_df.loc[file_name_df[2] == "leaves.pcd",2] = "all"
    file_name_df.columns = ["plot","year","species"]
    
    return pd.concat([file_name_df,values_df],axis=1)