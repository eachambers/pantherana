#!/usr/bin/python
# The following code is modified from that provided by Erik Enbody
# His repo is available here: https://github.com/erikenbody/ccgp_feems/

# base
import numpy as np
import pkg_resources
from sklearn.impute import SimpleImputer
from pandas_plink import read_plink
from os.path import basename, dirname, join
import pandas as pd
import sys
import pickle
import statsmodels.api as sm

# viz
import matplotlib.pyplot as plt
from matplotlib import gridspec
import cartopy.crs as ccrs

# feems
from feems.utils import prepare_graph_inputs
from feems import SpatialGraph, Viz, Objective, query_node_attributes
from feems.cross_validation import run_cv, comp_mats

# change matplotlib fonts
plt.rcParams["font.family"] = "Arial"
plt.rcParams["font.sans-serif"] = "Arial"

# set a path to where the data will live
data_path = "Anne/run_feems/data"
output_path = "Anne/run_feems/output"

# ==== READ IN INPUT FILES
# read the genotype data and mean impute missing data
plink_path = "{}/forreri_FILT".format(data_path)
(bim, fam, G) = read_plink(plink_path)
imp = SimpleImputer(missing_values=np.nan, strategy="mean")
genotypes = imp.fit_transform((np.array(G)).T)
n, p = genotypes.shape
print("n_samples={}, n_snps={}".format(genotypes.shape[0], genotypes.shape[1]))

# read in coords, retaining only x and y
coord = np.loadtxt("../data/forreri_coords.txt", usecols = (1,2))

# ==== SET UP GRAPH
# coords but added 0.01 buffer around points
#outer = np.array([[-105.7535, 9.949587],
#				[-105.7535, 22.54783],
#				[-82.96517, 22.56783],
#				[-82.96517, 9.949587]])

outer = np.loadtxt("../data/forreri_outer.txt")
grid_path = "../data/forr_grid.shp"  # path to discrete global grid

# graph input files
outer, edges, grid, _ = prepare_graph_inputs(coord=coord, 
                                             ggrid=grid_path,
                                             translated=True, 
                                             buffer=0,
                                             outer=outer)

# Set scale_snps to false when there is pervasive missing data and/or when invariant data 
# is present. But it must be true in order to run the lambda optimizer
sp_graph = SpatialGraph(genotypes, coord, grid, edges, scale_snps=True)

# ==== RUN CV ANALYSIS

LoadCVFromDisk = False

# define grid. lamba values of 0.001-20
lamb_grid = np.geomspace(1e-2, 2e1, 20)[::-1]
pd.DataFrame(lamb_grid).to_csv('../output/lamb_grid.csv', header=False, index=False)

# run cross-validation
if not LoadCVFromDisk:
    cv_err = run_cv(sp_graph, lamb_grid, n_folds=sp_graph.n_observed_nodes, factr=1e10)
    pickle.dump(cv_err,open("cv_err.pkl","wb"))

# average over folds
mean_cv_err = np.mean(cv_err, axis=0)
pd.DataFrame(mean_cv_err).to_csv('../output/mean_cv_err.csv', header=False, index=False)

# Select lambda based on min cv error
lamb_cv = float(lamb_grid[np.argmin(mean_cv_err)]) # 20

# ==== EXPORT FILES
# https://github.com/karolisr/pitcairnia-dr-nrv/blob/93534b6906d2a59a821c19deaef649eb1a3fb283/31-feems/31-feems-run.py#L196

feems_nodes = sp_graph.nodes
pd.DataFrame(feems_nodes).to_csv('../output/feems_nodes.csv', header=False, index=False)

feems_node_pos = sp_graph.node_pos
pd.DataFrame(feems_node_pos).to_csv('../output/feems_node_pos.csv', header=False, index=False)

#feems_node_pos_T = sp_graph.node_pos.T
#pd.DataFrame(feems_node_pos_T).to_csv('../output/feems_node_pos_T.csv', header=False, index=False)

feems_edges = sp_graph.edges
pd.DataFrame(feems_edges).to_csv('../output/feems_edges.csv', header=False, index=False)

feems_w = sp_graph.w
pd.DataFrame(feems_w).to_csv('../output/feems_w.csv', header=False, index=False)
