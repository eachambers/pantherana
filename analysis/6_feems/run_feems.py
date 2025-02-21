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
data_path = "/media/wanglab/798f0e01-89d1-4d0f-8ed8-ef323be70ab9/Anne/run_feems/data"
output_path = "/media/wanglab/798f0e01-89d1-4d0f-8ed8-ef323be70ab9/Anne/run_feems/output"

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

# ==== GENERATE CV PLOT

# Plot CV error
fig, ax = plt.subplots(dpi=300)
ax.plot(np.log10(lamb_grid), mean_cv_err, ".");
ax.set_xlabel("log10(lambda)");
ax.set_ylabel("L2 CV Error");
ax.axvline(np.log10(lamb_cv), color = "orange")

plt.savefig("forreri_cv.png")

# ==== GENERATE MAP

# Figure setup
projection = ccrs.EquidistantConic(central_longitude=-94.35934, central_latitude=16.25871)
fig = plt.figure(dpi=300)
ax = fig.add_subplot(1, 1, 1, projection=projection)
v = Viz(ax, sp_graph, projection=projection, edge_width=0.5, 
        edge_alpha=1, edge_zorder=100, sample_pt_size=10, 
        obs_node_size=5, sample_pt_color="black", 
        cbar_font_size=10)

v.draw_map()
v.draw_samples()
v.draw_edges(use_weights=False) # if True, will colorize edges according to migration rate weights
v.draw_obs_nodes(use_ids=False) # will number points if True

# plt.savefig("test_weights_ids.png")

# Add fit to the map
sp_graph.fit(float(lamb_cv)) # may take a little while to run

fig = plt.figure(dpi=300)
ax = fig.add_subplot(1, 1, 1, projection=projection)  
v = Viz(ax, sp_graph, projection=projection, edge_width=0.5, 
  edge_alpha=1, edge_zorder=100, sample_pt_size=20, 
  obs_node_size=7.5, sample_pt_color="black", 
  cbar_font_size=10)
v.draw_map()
v.draw_edges(use_weights=True)
v.draw_obs_nodes(use_ids=False) 
v.draw_edge_colorbar()
ax.text(.2, .85, "lambda={:.5f}\ncv l2 error={:.5f}".format(lamb_cv, mean_cv_err[np.argmin(mean_cv_err), 0]), 
           fontdict={"fontsize": 4}, transform = ax.transAxes)

# Save the figure
plt.savefig("forreri_feems.png")

# ==== LOOK AT THE FIT

def cov_to_dist(S):
    s2 = np.diag(S).reshape(-1, 1)
    ones = np.ones((s2.shape[0], 1))
    D = s2 @ ones.T + ones @ s2.T - 2 * S
    return(D)
    
from cov_to_dist import *
# from misc_functions import cov_to_dist
# import math

# Initialize figure
fig = plt.figure(constrained_layout=True, dpi=300, figsize=(6, 6))
spec = gridspec.GridSpec(ncols=2, nrows=2, figure=fig)

# from scipy.spatial.distance import pdist, squareform
# from sklearn.metrics.pairwise import haversine_distances

# (A) Genetic distance vs geographic distance
D_geno = squareform(pdist(genotypes, metric="sqeuclidean")) / p
coord_rad = coord[:,::-1] * math.pi / 180.0
D_geo = haversine_distances(coord_rad) * 6371000/1000
tril_idx = np.tril_indices(n, k=-1)
x = D_geo[tril_idx]
y = D_geno[tril_idx]
X = sm.add_constant(x)
mod = sm.OLS(y, X)
res = mod.fit()
muhat, betahat = res.params

ax_00 = fig.add_subplot(spec[0, 0])
ax_00.set_title("A", loc='left')
ax_00.scatter(x, 
              y, 
              marker=".", 
              alpha=1, 
              zorder=0, 
              color="black",
              s=3)

x_ = np.linspace(np.min(x), np.max(x), 20)
ax_00.plot(x_, muhat + betahat * x_, zorder=2, color="orange", linestyle='--', linewidth=1)
ax_00.text(3500, .6, "R²={:.4f}".format(res.rsquared))
ax_00.set_xlabel("great circle distance (km)")
ax_00.set_ylabel("genetic distance")

plt.savefig("forreri_gendistgeodist.png")

# (B) Genetic distance vs fitted distance for constant w model 
tril_idx = np.tril_indices(sp_graph.n_observed_nodes, k=-1)
ax_01 = fig.add_subplot(spec[0, 1])
ax_01.set_title("B", loc='left')
sp_graph.fit_null_model() # takes a few seconds to run; constant-w/variance fit, converged in 108 iterations, train_loss=489512.7140462
sp_graph.comp_graph_laplacian(sp_graph.w)

obj = Objective(sp_graph)
fit_cov, _, emp_cov = comp_mats(obj)
fit_dist = cov_to_dist(fit_cov)[tril_idx] # name np is not defined issue
emp_dist = cov_to_dist(emp_cov)[tril_idx]
X = sm.add_constant(fit_dist)
mod = sm.OLS(emp_dist, X)
res = mod.fit()
muhat, betahat = res.params
ax_01.scatter(fit_dist, 
              emp_dist, 
              marker=".", 
              alpha=1, 
              zorder=0, 
              color="black",
              s=3)

x_ = np.linspace(np.min(fit_dist), np.max(fit_dist), 20)
ax_01.plot(x_, muhat + betahat * x_, zorder=2, color="orange", linestyle='--', linewidth=1)
ax_01.text(4.2, 2, "R²={:.4f}".format(res.rsquared))
ax_01.set_xlabel("fitted distance (constant w)")

# (C) Genetic distance vs fitted distance for lambda = lambda_cv
tril_idx = np.tril_indices(sp_graph.n_observed_nodes, k=-1)
ax_10 = fig.add_subplot(spec[1, 0])
ax_10.set_title("C", loc='left')
lamb = lamb_cv
sp_graph.fit(lamb=lamb,
             lb=math.log(1e-6), 
             ub=math.log(1e+6))
sp_graph.comp_graph_laplacian(sp_graph.w)

obj = Objective(sp_graph)
fit_cov, _, emp_cov = comp_mats(obj)
fit_dist = cov_to_dist(fit_cov)[tril_idx]
emp_dist = cov_to_dist(emp_cov)[tril_idx]
X = sm.add_constant(fit_dist)
mod = sm.OLS(emp_dist, X)
res = mod.fit()
muhat, betahat = res.params
ax_10.scatter(fit_dist,
              emp_dist,
              marker=".", 
              alpha=1, 
              zorder=0, 
              color="black",
              s=3)

x_ = np.linspace(np.min(fit_dist), np.max(fit_dist), 20)
ax_10.plot(x_, muhat + betahat * x_, zorder=2, color="orange", linestyle='--', linewidth=1)
ax_10.text(4.2, 2, "R²={:.4f}".format(res.rsquared))
ax_10.set_xlabel("fitted distance ($\lambda = \lambda_{CV}$)")
ax_10.set_ylabel("genetic distance")

# (D) Genetic distance vs fitted distance for lambda = 1e-3*lambda_cv
tril_idx = np.tril_indices(sp_graph.n_observed_nodes, k=-1)
ax_11 = fig.add_subplot(spec[1, 1])
ax_11.set_title("D", loc='left')
lamb = 0.001*lamb_cv
sp_graph.fit(lamb=lamb,
             lb=math.log(1e-6), 
             ub=math.log(1e+6))
sp_graph.comp_graph_laplacian(sp_graph.w)

obj = Objective(sp_graph)
fit_cov, _, emp_cov = comp_mats(obj)
fit_dist = cov_to_dist(fit_cov)[tril_idx]
emp_dist = cov_to_dist(emp_cov)[tril_idx]
X = sm.add_constant(fit_dist)
mod = sm.OLS(emp_dist, X)
res = mod.fit()
muhat, betahat = res.params
ax_11.scatter(fit_dist,
              emp_dist,
              marker=".", 
              alpha=1, 
              zorder=0, 
              color="black",
              s=3)

x_ = np.linspace(np.min(fit_dist), np.max(fit_dist), 20)
ax_11.plot(x_, muhat + betahat * x_, zorder=2, color="orange", linestyle='--', linewidth=1)
ax_11.text(4.2, 2, "R²={:.4f}".format(res.rsquared))
ax_11.set_xlabel("fitted distance ($\lambda = \lambda_{CV}\cdot 10^{-3}$)")
ax_11.set_ylabel("genetic distance")

# axis 00 
ax_00 = fig.add_subplot(spec[0, 0], projection=projection)
ax_00.set_title("A", loc=title_loc, pad=title_pad, fontdict={"fontsize": title_fontsize})
sp_graph.fit(float(lamb_grid[0]))
v = Viz(ax_00, sp_graph, projection=projection, edge_width=edge_width, 
        edge_alpha=1, edge_zorder=100, sample_pt_size=20, 
        obs_node_size=obs_node_size, sample_pt_color="black", 
        cbar_font_size=10)
v.draw_map()
v.draw_edges(use_weights=True)
v.draw_obs_nodes(use_ids=False)
v.draw_edge_colorbar()
 
ax_00.text(.2, .85, "lambda={:.5f}\ncv l2 error={:.5f}".format(lamb_grid[0], mean_cv_err[0, 0]), 
           fontdict={"fontsize": 4}, transform = ax_00.transAxes)

# axis 10
ax_10 = fig.add_subplot(spec[1, 0], projection=projection)
ax_10.set_title("B", loc=title_loc, pad=title_pad, fontdict={"fontsize": title_fontsize})
sp_graph.fit(float(lamb_grid[3]))
v = Viz(ax_10, sp_graph, projection=projection, edge_width=edge_width, 
        edge_alpha=1, edge_zorder=100, sample_pt_size=20,
        obs_node_size=obs_node_size, sample_pt_color="black", 
        cbar_font_size=10)
v.draw_map()
v.draw_edges(use_weights=True)
v.draw_obs_nodes(use_ids=False) 
ax_10.text(.2, .85, "lambda={:.5f}\ncv l2 error={:.5f}".format(lamb_grid[3], mean_cv_err[3, 0]), 
           fontdict={"fontsize": 4}, transform = ax_10.transAxes)

# axis 01
ax_01 = fig.add_subplot(spec[0, 1], projection=projection)
ax_01.set_title("C", loc=title_loc, pad=title_pad, fontdict={"fontsize": title_fontsize})
sp_graph.fit(float(lamb_cv))
v = Viz(ax_01, sp_graph, projection=projection, edge_width=edge_width, 
        edge_alpha=1, edge_zorder=100, sample_pt_size=20, 
        obs_node_size=obs_node_size, sample_pt_color="black", 
        cbar_font_size=10)
v.draw_map()
v.draw_edges(use_weights=True)
v.draw_obs_nodes(use_ids=False) 
ax_01.text(.2, .85, "lambda={:.5f}\ncv l2 error={:.5f}".format(lamb_cv, mean_cv_err[np.argmin(mean_cv_err), 0]), 
           fontdict={"fontsize": 4}, transform = ax_01.transAxes)

# axis 11
ax_11 = fig.add_subplot(spec[1, 1], projection=projection)
ax_11.set_title("D", loc=title_loc, pad=title_pad, fontdict={"fontsize": title_fontsize})
sp_graph.fit(float(lamb_grid[10]))
v = Viz(ax_11, sp_graph, projection=projection, edge_width=edge_width, 
        edge_alpha=1, edge_zorder=100, sample_pt_size=20, 
        obs_node_size=obs_node_size, sample_pt_color="black", 
        cbar_font_size=10)
v.draw_map()
v.draw_edges(use_weights=True)
v.draw_obs_nodes(use_ids=False)
v.cbar_font_size = cbar_font_size
v.cbar_orientation = cbar_orientation
v.cbar_ticklabelsize = cbar_ticklabelsize
v.draw_edge_colorbar()
ax_11.text(.2, .85, "lambda={:.5f}\ncv l2 error={:.5f}".format(lamb_grid[10], mean_cv_err[10, 0]), 
           fontdict={"fontsize": 4}, transform = ax_11.transAxes)

plt.savefig("{}/{prefix}_cv_maps.png".format(output_path, prefix=prefix))

# ==== ERRORS

# If you get an error indicating PROJECTION issue, be sure to remove the .prj file from wherever
#your input files are.

# spatial_graph.py:81: RuntimeWarning: divide by zero encountered in true_divide
#   self.frequencies = self.frequencies / np.sqrt(self.mu * (1 - self.mu))
# ^ Requires verification that invariant sites aren't included in genotype matrix;
# also be sure to stringently filter for missing data

