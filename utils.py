import math
import numpy as np
import matplotlib.pyplot as plt

import open3d as o3d
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

from mpl_toolkits.mplot3d import Axes3D
from scipy.spatial import ConvexHull
from skspatial.plotting import plot_3d

def compute_projected_area(voxel_df):
    """
    Compute the projected area of a tree
    
    Parameters
    ----------------
    voxel_df: 
        dataframe containing x, y, and z voxel coordinates. Must be in 
        this exact order (x,y,z)
    
    Returns
    ----------------
    projected area
    """
    width = voxel_df[0].max() - voxel_df[0].min()
    radius = width/2
    
    return math.pi * radius ** 2

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

def calculate_leaf_angles(slices_voxels, slices_normals, visualize=False): 
    """
    Compute the needle angles
    
    Parameters
    ----------------
    slices_voxels: 
        voxel coordinates datafame
    slices_normals:
        normal vectors dataframe
    visualize:
        True or False
        visualize the histogram of individual needle angles 
        
    Returns
    ----------------
    leaf angles in degrees
    """   
    all_leafangles = []

    for index,(vox, norm) in enumerate(zip(slices_voxels, slices_normals)): # iterate over slices
        slice_values = [] # each iteration of the slice open a new list 

        for i, j in zip(vox, norm): # iterate over values in slices 
            dot_product = (np.dot(i, j)) # numerator
            vox_mag = math.sqrt(np.sum(np.power(i, 2))) # denominator
            norm_mag = math.sqrt(np.sum(np.power(j, 2))) # denominator
            denomi_prod = vox_mag * norm_mag
            angle = np.arccos(dot_product/denomi_prod)
            angle_deg = np.degrees(angle)
            if angle_deg >= 90:
                angle_deg = 180 - angle_deg
            slice_values.append(angle_deg) # must append outside of if statement to get all angles, so all_leafangles.append(angle_deg) combines all slice angles

        all_leafangles.append(slice_values)

        if visualize:
            plt.figure # initialize figure so the information is separate
            plt.hist(slice_values)
            plt.title(f'slice {index+1} inclination angles')
            plt.xlabel("leaf angle (degrees)")
            plt.ylabel("density")

            plt.show()
        
    return all_leafangles

def calc_convex_hull(slices_voxels_conv, voxs, voxel_size, visualize=False): 
    """
    Perform convex hull to find contact frequency
    
    Parameters
    ----------------
    slices_voxels_conv: 
        total possible voxel volume
    voxs:
        actual leaf voxels
    voxel_size:
        size used to voxelize tree (0.0325 m is default)
    visualize:
        True or False
        visualize the convex hull 
        
    Returns
    ----------------
    all_actual_voxs:
        leaf occupied space 
    all_null_vox:
        gaps
    
    """       
    vox_size = np.power(voxel_size, 3)

    all_actual_voxs = []
    all_null_vox = []

    for vox_con in slices_voxels_conv:
        vols = ConvexHull(vox_con).volume
        poss_voxs = vols/vox_size 
        actual_vox = len(vox_con)
        null_vox = poss_voxs - actual_vox
#         print('Possible voxels:', poss_voxs)
#         print('Actual voxels:', actual_vox)
#         print('Null voxels:', null_vox)

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

        all_actual_voxs.append(poss_voxs)
        all_null_vox.append(null_vox)
    return all_actual_voxs, all_null_vox

def whole_pipeline(file_path_1, file_path_2, visualize=False, voxel_size = 0.0325):
    """
    Import all functions into one pipeline 
    
    Parameters
    ----------------
    slices_voxels_conv: 
        total possible voxel volume
    voxs:
        actual leaf voxels
    voxel_size:
        size used to voxelize tree (0.0325 m is default)
    visualize:
        True or False
        visualize the convex hull 
        
    Returns
    ----------------
    LAI:
        leaf area index 
    LAD:
        leaf area density
    allangles:
        average leaf angle for each slice
    LA_correction_factors:
        leaf angle-based correction factor 
    """
    # read point cloud
    pcd = o3d.io.read_point_cloud(file_path_1) # load point cloud, must be in .pcd format or .ply format

    # voxelize
    #voxel_size = 0.0325
    downpcd = pcd.voxel_down_sample(voxel_size = voxel_size)
    # extract voxel centroids
    df_vox = pd.DataFrame(downpcd.points) # makes the voxelized point cloud an array of the coordinates that make it up

    # normal calculation
    downpcd.estimate_normals(search_param = o3d.geometry.KDTreeSearchParamHybrid(radius = 0.1, max_nn =30))
    df_norm = pd.DataFrame(downpcd.normals)

    # compute projected area
    PA = compute_projected_area(df_vox)  

    # consolidate voxels and normals into a single dataframe
    # and sort by elevation
    sorted_df = pd.concat([df_vox, df_norm], axis = 1)
    sorted_df.columns = ["x", "y", "z", "norm_x", "norm_y", "norm_z"]
    sorted_df = sorted_df.sort_values(by = ['z'])

    #separate the now sorted voxels and normals 
    voxs = sorted_df.iloc[:,:-3]
    norms = sorted_df.iloc[:,-3:]
    
    # compute tree height
    slice_width = 0.5
    h = compute_height(df_vox)

    # how many half meter slices, round up  
    num_slices = np.ceil(h/slice_width)
    print(f"how many half meter slices: {num_slices}") # tree L2, 5 slices

    # slice the arrays 
    sliced_norms = np.array_split(norms.to_numpy(), num_slices) # put these in the function 
    sliced_vox = np.array_split(voxs.to_numpy(), num_slices)

    if visualize:
        fig = plt.figure()
        ax = fig.add_subplot(111, projection = "3d")

        for vox_slice in sliced_vox:
            ax.scatter(vox_slice[:,0], vox_slice[:,1], vox_slice[:,2])
        plt.show()
    
    # leaf angles
    slices_voxs = sliced_vox
    slices_norms = sliced_norms
    all_leafangles = []
    dot_products_all = []
    denominator_prods = []
    for index,(vox, norm) in enumerate(zip(slices_voxs, slices_norms)): # iterate over slices
        slice_values = [] # each iteration of the slice we want to open a new list 
        dot_prods_slices = []
        denom_prods_slices = []

        for i, j in zip(vox, norm): # iterate over values in slices 
            dot_product = (np.dot(i, j)) # numerator
            vox_mag = np.sqrt(np.sum(np.power(i, 2))) # denominator
            norm_mag = np.sqrt(np.sum(np.power(j, 2))) # denominator
            denomi_prod = vox_mag * norm_mag
            angle = np.arccos(dot_product/denomi_prod)
            angle_deg = np.degrees(angle)
            if angle_deg >= 90:
                angle_deg = 180 - angle_deg
            slice_values.append(angle_deg) # must append outside of if statement to get all angles, so all_leafangles.append(angle_deg) combines all slice angles
            dot_prods_slices.append(dot_product)
            denom_prods_slices.append(denomi_prod)

        all_leafangles.append(slice_values)
        dot_products_all.append(dot_prods_slices)
        denominator_prods.append(denom_prods_slices)


        if visualize:
            plt.figure # initialize figure so the information is separate
            plt.hist(slice_values)
            plt.title(f'slice {index+1} inclination angles')
            plt.xlabel("leaf angle (degrees)")
            plt.ylabel("density")

            plt.show()

    conifer_leaf_angles = calculate_leaf_angles(sliced_vox, sliced_norms,visualize)
    
    # prints mean leaf angle for each slice
    allangles = []
    for x in all_leafangles:
        allangles.append(np.mean(x))
#         allangles.append(x)

        
    ## Correction Factor
    # define zenith angles, is there a better way of generating these? I want to able to perform all of the following on all the slices at once 
    conif_zeniths = np.array_split(np.arange (90, 0, -1), num_slices) # put these in the function 
    
    # mikel help 
    # def correction_fact_calc(conif_angles, conif_slices, conif_dtprods, conif_dmprods, conif_zeniths)
    conif_angles = conifer_leaf_angles
    conif_slices = slices_voxs
    conif_dtprods = dot_products_all
    conif_dmprods = denominator_prods 
    LA_correction_factors = []
    cf_stdvs = []
    G_function = []
    
    for index in range(len(conifer_leaf_angles)):
        which_incang = np.deg2rad(conif_angles[index])
        which_slice = conif_slices[index]
        which_dtp = conif_dtprods[index] 
        which_dnmp = conif_dmprods[index]
        zenan = np.deg2rad(conif_zeniths[index])

        # spherical coords
        x_TLS_ar = which_slice[:,0]
        y_TLS_ar = which_slice[:,1]
        z_TLS_ar = which_slice[:,2]

        scanner_radius = []
        for i, j, k in zip(x_TLS_ar, y_TLS_ar, z_TLS_ar):
            i_sq = np.power(i, 2)
            j_sq = np.power(j, 2)
            k_sq = np.power(k, 2)
            scanner = i_sq + j_sq + k_sq
            radius = np.sqrt(scanner)
            scanner_radius.append(radius)
        avg_scanner_radius = np.mean(scanner_radius) # average must be outside the second for loop

        # find phi or azimuth of scanner arctan(y/x)
        phi_laser = []
        for i, j in zip(y_TLS_ar, x_TLS_ar):
            phi_laser.append(np.arctan(i/j)) 

        # find phi or azimuth of the leaves 
        # divide the dot product by the denominator and take inverse tan
        phi_leaves = []
        for i, j in zip(which_dtp, which_dnmp):
            phi_leaves.append(np.arctan(i/j) )

        rand_phi_laser = np.random.choice(phi_laser, size=len(zenan),replace=False)

        # nB from Hosoi 2007 nB = sin(thta)cos(phi), sin(thta), cos(thta), unit vector corresponding to the direction of the laser beam
        nB = []
        for i, j in zip(zenan, rand_phi_laser):
            y = np.sin(i)*np.cos(j), np.sin(i)*np.sin(j), np.cos(i)
            nB.append(y)
        nB = np.asarray(nB)

        # nL from Hosoi 2007 nL = sin(thta)cos(phi), sin(thta), cos(thta), unit vector corresponding to the direction of the normal to the leaf surface
        nL = []
        for i, j in zip(which_incang, phi_leaves): #CHANGE INCLIN_ANG ACCORDING TO WHAT SLICE YOU ARE ON
            y = np.sin(i)*np.cos(j), np.sin(i)*np.sin(j), np.cos(i)
            nL.append(y)
        nL = np.asarray(nL)

        # find G and S
        # where is zenith laser <= pi/2 - zenith leaves(leaf inclination angle) (hosoi, 2006)
        pi_inc_value = [] # pi/2 - inclination angle values 
        for i in which_incang:
            y = (np.pi/2) - i
            x = y
            pi_inc_value.append(x)

        rand_inc_leaves = np.random.choice(pi_inc_value, size=len(zenan),replace=False)

        S = []
        for i, j in zip(zenan, rand_inc_leaves):
            y = (np.cos(i)) * (np.cos(j))
            if i > (np.pi/2) - j:    
                tan = (np.tan(i)*np.tan(j)) # cotangent needs to be calculated as 1/tan
                cotan = (1/tan)
                x = np.arccos(cotan)
                bracket = 1 + (2*((np.tan(x)-x)))/np.pi
                y = y * bracket
            S.append(y)

        # correction factor 
        corr_fac = []
        for i, s in zip(zenan, S):
            y = ((np.cos(i))/s)/3
            corr_fac.append(y)

        LA_correction_factors.append(np.mean(corr_fac))

        stdvs = np.std(corr_fac)
        cf_stdvs.append(stdvs)
        G_function.append(S)

    print('Slice Correction Factors:', LA_correction_factors)
    print('CF Standard Deviations:', cf_stdvs)
   
    # Convex Hull
    all_actual_voxs, all_null_vox = calc_convex_hull(slices_voxs,voxs,voxel_size,visualize)    
    
    # WHOLE PLOT
    pcd_whole = o3d.io.read_point_cloud(file_path_2)
    downpcd_whole = pcd_whole.voxel_down_sample(voxel_size = voxel_size)
    
    df_vox_whole = pd.DataFrame(downpcd_whole.points)
    
    # compute tree height
    h_whole = compute_height(df_vox_whole)
    
    df_vox_whole.columns = ["x", "y", "z"]
    sorted_vox_whole = df_vox_whole.sort_values(by = ['z'])
    vox_array_whole = sorted_vox_whole.to_numpy()
    
    # how many half meter slices, round up  
    num_slices_whole = np.ceil(h_whole/slice_width)
    print(f"how many half meter slices for the second file: {num_slices_whole}")
    sliced_vox_whole = np.array_split(vox_array_whole, num_slices_whole)
    
    ## Convex Hull whole
    _, whole_plot_nulls = calc_convex_hull(sliced_vox_whole,vox_array_whole,voxel_size,visualize)   
    
    all_actual_voxs = np.asarray(all_actual_voxs)
    # whole plot, use for species
    all_null_vox = np.asarray(whole_plot_nulls)
    
    LAD = []
    for j, k, n in zip(all_actual_voxs, all_null_vox, LA_correction_factors): 
#         y = n * (1/slice_width) * (j/(j+k))
        y = n * (1/slice_width) * (j/(j+(3000000 - j))) # specific to plot footprint, uncomment if you want to just use convex hull
        LAD.append(y)

    print('The LAD is:', LAD)
    LAI = np.sum(LAD)
    print('The LAI is:', LAI)
    
    return LAI,LAD,allangles,LA_correction_factors
    
    
