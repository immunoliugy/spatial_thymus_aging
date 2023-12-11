import ast
import bisect
from collections import Counter
import copy
import csv
from datetime import datetime
import importlib
import itertools
import math
from matplotlib import pyplot as plt
from matplotlib.collections import PatchCollection
from mpl_toolkits.axes_grid1 import make_axes_locatable
from mpl_toolkits.axes_grid1.anchored_artists import AnchoredSizeBar
from matplotlib.colors import LinearSegmentedColormap
import matplotlib.pyplot as plt
from matplotlib.ticker import ScalarFormatter
from matplotlib_venn import venn2, venn2_circles
from matplotlib_scalebar.scalebar import ScaleBar
import multiprocess
import numpy as np
import openpyxl
import os
from os import path
import pandas as pd
import pickle5 as pkl
import random as rand
import random
import scanpy as sc
import scipy.sparse
import scipy
from scipy import spatial
import scanpy.external as sce
from scipy.spatial import distance
from scipy.stats import kstest, pearsonr, gaussian_kde as kde
from scipy.ndimage import gaussian_filter
import seaborn as sns
from sklearn.decomposition import NMF
from sklearn.neighbors import KNeighborsClassifier
from sklearn.cluster import KMeans
import squidpy as sq
import statsmodels
from statistics import mean, stdev
from statsmodels import stats
import statsmodels.stats.multitest
import sys
from tqdm.notebook import tqdm
from scipy.stats import ttest_ind
from scipy.spatial import cKDTree


def report_ct_counter(locs, loc_to_ct):
    cts = []
    for loc in locs:
        cts.append(loc_to_ct[loc])
    return (Counter(cts))


def merge_in_cter(new_cter, all_cter):
    for ct in all_cter:
        if ct not in new_cter:
            new_cter[ct] = 0
        all_cter[ct].append(new_cter[ct] / sum(new_cter.values()))
    #all_cter[ct] = [np.average(all_cter[ct])]
    return all_cter


def record_averages(measured, recorder):
    for ct in measured:
        if ct not in recorder:
            recorder[ct] = []
        recorder[ct].append(np.average(measured[ct]))
    return recorder


def run_prox_test_on_ct(ct,
                        region_locs,
                        ct_to_loc,
                        loc_to_ct,
                        rand_iter_num=100):  # medulla_locs or cortex_locs
    print(ct)

    # Cortex/Medulla test
    point_tree = cKDTree(region_locs)
    dist = 50  # microns
    dist = dist / 0.65

    nearby_ct_counter_all_means = {}
    rand_ct_counter_all_means = {}

    nearby_ct_counter_all = {i: [] for i in set(ct_to_loc)}

    region_locs_ct = set(set(ct_to_loc[ct]) & set(region_locs))

    if len(region_locs_ct) < 10:  # If fewer than 10 of this cell type
        print(ct + ' | Not enough cells')
        return 'Not enough cells'

    loc_to_num_nearby_locs = {}
    for loc in region_locs_ct:
        nearby_locs_idx = point_tree.query_ball_point(loc, dist)
        nearby_locs = [region_locs[i] for i in nearby_locs_idx]
        nearby_locs.remove(loc)
        if len(nearby_locs) < 5:  # If fewer than 5 cells nearby
            continue
        loc_to_num_nearby_locs[loc] = len(nearby_locs)

        cter = report_ct_counter(nearby_locs, loc_to_ct)
        nearby_ct_counter_all = merge_in_cter(cter, nearby_ct_counter_all)

    nearby_ct_counter_all_means = record_averages(nearby_ct_counter_all,
                                                  nearby_ct_counter_all_means)

    # Record averages

    # For each loc, sample randomly within a certain distance

    for it in tqdm(range(rand_iter_num)):
        rand_ct_counter_all = {i: [] for i in set(ct_to_loc)}

        for loc in region_locs_ct:
            rand_locs_idx = point_tree.query_ball_point(loc, dist * 3)
            rand_locs = [region_locs[i] for i in rand_locs_idx]
            rand_locs.remove(loc)
            if loc not in loc_to_num_nearby_locs:
                continue
            rand_locs = random.sample(rand_locs, loc_to_num_nearby_locs[loc])
            #print(nearby_locs)

            rand_cter = report_ct_counter(rand_locs, loc_to_ct)

            rand_ct_counter_all = merge_in_cter(rand_cter, rand_ct_counter_all)
        rand_ct_counter_all_means = record_averages(rand_ct_counter_all,
                                                    rand_ct_counter_all_means)

    return nearby_ct_counter_all_means, rand_ct_counter_all_means


def calc_p_vals(observed, rand, show_plot=False, show_sig_plots_only=False):
    p_val_all = {}

    for ct in observed:
        sig = False
        if len(rand[ct]) == 0:
            p_val = 1
        else:
            p_val = len([i for i in rand[ct] if i > observed[ct]]) / len(
                rand[ct])
        if p_val < 0.05:
            sig = True
        if p_val > 0.95:
            sig = True
        p_val_all[ct] = p_val
        if show_plot:
            if (show_sig_plots_only & ~sig):
                continue
            else:
                plt.figure()
                plt.hist(rand[ct], bins=5, alpha=0.5)
                plt.title('Mac close to ' + ct)
                plt.axvline(x=observed[ct], color='r', linestyle='--')
                plt.title(ct + f' | pval={p_val}')

    return p_val_all


def perm_test_to_final_plot(input_info, load_previous=False): # input_info = [slideseq_name, region_type]
    #puck = load_puck(slideseq_name)
    slideseq_name = input_info[0]
    region_type = input_info[1]
    directory = input_info[2]
    # load adata
   
    adata = sc.read(f'{directory}/filt_annotated_h5ad/{slideseq_name}_cortex_medulla_xy.h5ad')

    bcs = adata.obs.index
    locs = [(x, y) for x, y in zip(adata.obs['x'], adata.obs['y'])]
    bc_loc_dict_s1 = {i: j for i, j in zip(bcs, locs)}
    #locs = [puck.bc_loc_dict_s1[puck.matched_bead_barcodes[i]] for i in bcs]
    x, y = zip(*locs)
    adata.obs['x'] = x
    adata.obs['y'] = y
    adata.obs['UMI'] = adata.obs['total_counts']
    adata.obs

    rctd_cell_annot_fn = f'{directory}/RCTD_outputs/{slideseq_name}_RCTD_results.csv'
    rctd_df = pd.read_csv(rctd_cell_annot_fn)
    rctd_df

    index = [i[:-2] + '-1' for i in rctd_df['Unnamed: 0']]
    rctd_df.index = index
    rctd_df.drop(['Unnamed: 0'], inplace=True, axis=1)
    rctd_df

    bc_to_ct = {i: j for i, j in zip(rctd_df.index, rctd_df.first_type)}

    all_locs = locs
    medulla_bcs = adata.obs[adata.obs['cortex_medulla'] == 'medulla'].index
    cortex_bcs = adata.obs[adata.obs['cortex_medulla'] == 'cortex'].index

    all_locs = [i for i in all_locs if i in rctd_df.index]
    medulla_bcs = [i for i in medulla_bcs if i in rctd_df.index]
    cortex_bcs = [i for i in cortex_bcs if i in rctd_df.index]

    medulla_locs = [bc_loc_dict_s1[i] for i in medulla_bcs]
    cortex_locs = [bc_loc_dict_s1[i] for i in cortex_bcs]

    loc_to_ct = {bc_loc_dict_s1[i]: bc_to_ct[i] for i in bc_to_ct}
    ct_to_loc = {}

    for loc in set(loc_to_ct):
        ct = loc_to_ct[loc]
        if ct not in ct_to_loc:
            ct_to_loc[ct] = []
        ct_to_loc[ct].append(loc)

    observed_all = {}
    rand_all = {}

    ct1 = set(list(ct_to_loc.keys()))
    new_ct1 = copy.deepcopy(ct1)

    if region_type == 'cortex':
        region_locs_to_use = cortex_locs
    elif region_type == 'medulla':
        region_locs_to_use = medulla_locs
    elif region_type == 'both':
        region_locs_to_use = cortex_locs + medulla_locs
    for ct in ct1:
        #     if ct not in ['B','gdT']: # for testing on fewer cases
        #         continue

        result = run_prox_test_on_ct(ct,
                                     region_locs_to_use,
                                     ct_to_loc,
                                     loc_to_ct,
                                     rand_iter_num=500)
        if result == 'Not enough cells':
            new_ct1.remove(ct)
            continue
        else:
            observed = result[0]
            rand = result[1]
        observed_all[ct] = observed
        rand_all[ct] = rand
    ct1 = copy.deepcopy(new_ct1)

    z_score_all = {}
    p_val_all = {}
    for ct in ct1:
        z_score_all[ct] = {}
        #     if ct not in ['B','gdT']: # for testing on fewer cases
        #         continue
        obs = observed_all[ct]
        rand = rand_all[ct]
        p_val = calc_p_vals(
            obs, rand)  #,show_plot = True,show_sig_plots_only = True
        p_val_all[ct] = p_val

        for ctct2 in obs:
            z_score_all[ct][ctct2] = scipy.stats.zscore(obs[ctct2] +
                                                        rand[ctct2])[0]

    p = [0, 0.05, 0.05, 1]
    f = lambda x: np.interp(x, p, [1, 0, 0, 0])

    cmap = LinearSegmentedColormap.from_list(
        'map_white',
        list(
            zip(np.linspace(0, 1), plt.cm.magma(f(np.linspace(min(p),
                                                              max(p)))))))

    ct1_arr = [i + '-1' for i in ct1]
    ct2_arr = [i[:-2] for i in ct1_arr]

    arr = np.ones([len(ct1_arr), len(ct2_arr)])

    for i in range(len(ct1_arr)):
        for j in range(len(ct2_arr)):
            arr[i, j] = p_val_all[ct1_arr[i][:-2]][ct2_arr[j]]

    plot_df = pd.DataFrame(arr, columns=ct2_arr, index=ct1_arr)
    print(plot_df)

    # plt.figure(figsize=(1,1))
    p = [0, 0.05, 0.05, 1]
    f = lambda x: np.interp(x, p, [1, 0, 0, 0])

    cmap = LinearSegmentedColormap.from_list(
        'map_white',
        list(
            zip(np.linspace(0, 1), plt.cm.magma(f(np.linspace(min(p),
                                                              max(p)))))))

    cmap_to_use = cmap  #sns.color_palette(220,20,as_cmap=True)
    cm = sns.clustermap(plot_df,
                        cmap=cmap_to_use,
                        figsize=(10, 10),
                        yticklabels=True,
                        xticklabels=True,
                        row_cluster=False,
                        col_cluster=False,
                        linewidths=0.5,
                        linecolor='gray')  #
    hm = cm.ax_heatmap.get_position()

    fn = open(f'prox_arr_{slideseq_name}_{region_type}.pkl', 'wb')
    pkl.dump(arr, fn)

    fn = open(f'prox_arr_{slideseq_name}_pd_{region_type}.pkl', 'wb')
    pkl.dump(plot_df, fn)

    ##########

    ct1_arr = [i + '-1' for i in ct1]
    ct2_arr = [i[:-2] for i in ct1_arr]

    arr = np.ones([len(ct1_arr), len(ct2_arr)])

    for i in range(len(ct1_arr)):
        for j in range(len(ct2_arr)):
            arr[i, j] = z_score_all[ct1_arr[i][:-2]][ct2_arr[j]]

    z_score_df = pd.DataFrame(arr, columns=ct2_arr, index=ct1_arr)

    fn = open(f'zscore_{slideseq_name}_pd_{region_type}.pkl', 'wb')
    pkl.dump(z_score_df, fn)
