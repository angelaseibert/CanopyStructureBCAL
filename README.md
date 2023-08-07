Estimating Leaf Area Index (LAI), leaf area density, and leaf inclination angles using Terrestrial Laser Scanning (Boise Center Aerospace Lab and SPRUCE)

**Welcome!**

Goal
------------

The main goal of this project is to estimate annual canopy structural metrics for SPRUCE plots from August 2015 - August 2022. To do so, we use a voxel-based contact frequency model. In this model, the  occupied voxels are quantified for each canopy layer, providing the leaf area density of each layer. Summing these layers results in the LAI. The vegetation structural data were measured from twelve plots (12 m x 12 m) in the Marcell Experimental Forest in northern Minnesota, USA. The lidar data comes from a Riegl VZ-1000 TLS with a 1550 nm laser. The original data is publically available:

* [Link to the field data](https://mnspruce.ornl.gov/datasets/spruce-terrestrial-laser-scanning-of-experimental-plots-beginning-in-2015)

The scripts provided in this project perform the following steps:

1. Separate leaf and wood points.
2. Voxelize the point cloud data.
3. Compute normal vectors and estimate leaf angles.
4. Slice the tree plots and estimate contact frequencies for each slice to get leaf area density.
5. Sum leaf area densities to estimate the LAI of each plot.


**Below you can find the project organization.**

Project Organization 
------------
    ├── LICENSE
    ├── README.md                <- The top-level README file overview of the project.
    ├── Canopy_Structure.ipynb   <- Runs full workflow to estimate all structural metrics using functions from utils.py.
    ├── Leaf_Volume.py           <- Contains functions to estimate the voxel volume of plots.
    ├── SPRUCE_RF.ipynb          <- Random Forest Classifier to separate leaves and wood.
    ├── destructive_trees.zip    <- Destructively harvested tree files used as model validation data.
    ├── leaf_volume_output.ipynb <- Notebook containing Leaf_Volume.py functions applied.
    ├── utils_final.py                 <- Contains functions to estimate all structural metrics. 

  ----------
