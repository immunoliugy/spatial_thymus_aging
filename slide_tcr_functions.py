import ast
import bisect
from collections import Counter
import copy
import csv
from datetime import datetime
import editdistance
from google.cloud import storage
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
import matplotlib.font_manager as fm
from matplotlib_scalebar.scalebar import ScaleBar
import multiprocess
from multiprocess import Pool
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
import skbio
from sklearn.decomposition import NMF
from sklearn.neighbors import KNeighborsClassifier
from sklearn.cluster import KMeans
import statsmodels
from statistics import mean, stdev
from statsmodels import stats
import statsmodels.stats.multitest
import sys
from tqdm import tqdm
from scipy.stats import ttest_ind
from scipy.spatial import cKDTree

# Set up functions
def download_data(slideseq_name):
    data_direc = '/PHShome/sx931/scratch/Thymus_data/move_to_mgh'
    slideseq_name_abbrev = '_'.join(slideseq_name.split('_')[1:])
    clonotype_name = '' #/with_failed
    clonotype_fastq_name = ''
    clonotype_mixcr_name = ''
    file_type = 'slide-seq-tools'
    bucket_name = 'gs://fc-secure-2423a594-f5fe-40c6-b998-f28ac17ea012'

    directory_clonotype = './{}'.format(clonotype_name)

    os.system('/apps/released/built-by-outside-authors-x86linuxtarget/gcloud/gcloud-420.0.0/google-cloud-sdk/bin/gsutil cp $bucket_name/$clonotype_name/$clonotype_mixcr_name\_paired_clones.txt $data_direc/$clonotype_mixcr_name\_paired_clones.txt')
    os.system('/apps/released/built-by-outside-authors-x86linuxtarget/gcloud/gcloud-420.0.0/google-cloud-sdk/bin/gsutil cp $bucket_name/$clonotype_name/$clonotype_mixcr_name\_paired_cloneID.txt  $data_direc/$clonotype_mixcr_name\_paired_cloneID.txt')
    os.system('/apps/released/built-by-outside-authors-x86linuxtarget/gcloud/gcloud-420.0.0/google-cloud-sdk/bin/gsutil cp $bucket_name/$clonotype_name/$clonotype_fastq_name\_R1_001.fastq.gz  $data_direc/$clonotype_fastq_name\_R1_001.fastq.gz')
    os.system('gunzip -f ./$clonotype_fastq_name\_R1_001.fastq.gz')
    os.system('/apps/released/built-by-outside-authors-x86linuxtarget/gcloud/gcloud-420.0.0/google-cloud-sdk/bin/gsutil cp $bucket_name/trac_trbc2_constant_umis_$slideseq_name\.pickle  $data_direc/trac_trbc2_constant_umis_$slideseq_name\.pickle')
    os.system('/apps/released/built-by-outside-authors-x86linuxtarget/gcloud/gcloud-420.0.0/google-cloud-sdk/bin/gsutil cp $bucket_name/$slideseq_name/myRCTD_$slideseq_name_abbrev\.rds $data_direc/')
        
    if file_type == 'jilong':
        os.system('/apps/released/built-by-outside-authors-x86linuxtarget/gcloud/gcloud-420.0.0/google-cloud-sdk/bin/gsutil cp $bucket_name/$slideseq_name/hg19.exonic+intronic/barcode_matching/$slideseq_name_abbrev\_matched_bead_barcodes.txt $data_direc/$slideseq_name_abbrev\_matched_bead_barcodes.txt')
        os.system('/apps/released/built-by-outside-authors-x86linuxtarget/gcloud/gcloud-420.0.0/google-cloud-sdk/bin/gsutil cp $bucket_name/$slideseq_name/hg19.exonic+intronic/barcode_matching/$slideseq_name_abbrev\_matched_bead_locations.txt  $data_direc/$slideseq_name_abbrev\_matched_bead_locations.txt')
        os.system('/apps/released/built-by-outside-authors-x86linuxtarget/gcloud/gcloud-420.0.0/google-cloud-sdk/bin/gsutil cp $bucket_name/$slideseq_name/hg19.exonic+intronic/alignment/$slideseq_name_abbrev\.digital_expression.txt.gz  $data_direc/$slideseq_name_abbrev\.digital_expression.txt.gz')
        os.system('gunzip -f ./$slideseq_name_abbrev\.digital_expression.txt.gz')
        os.system('/apps/released/built-by-outside-authors-x86linuxtarget/gcloud/gcloud-420.0.0/google-cloud-sdk/bin/gsutil cp $bucket_name/$slideseq_name/hg19.exonic+intronic/alignment/$slideseq_name_abbrev\.digital_expression.txt  $data_direc')

    else: # change below
        os.system('/apps/released/built-by-outside-authors-x86linuxtarget/gcloud/gcloud-420.0.0/google-cloud-sdk/bin/gsutil cp $bucket_name/$slideseq_name/barcode_matching/$slideseq_name_abbrev\_barcode_matching.txt.gz  $data_direc/$slideseq_name_abbrev\_barcode_matching.txt.gz')
        os.system('/apps/released/built-by-outside-authors-x86linuxtarget/gcloud/gcloud-420.0.0/google-cloud-sdk/bin/gsutil cp $bucket_name/$slideseq_name/barcode_matching/$slideseq_name_abbrev\_barcode_matching.txt  $data_direc/$slideseq_name_abbrev\_barcode_matching.txt')
        os.system('gunzip -f ./$slideseq_name_abbrev\_barcode_matching.txt.gz')
        os.system('/apps/released/built-by-outside-authors-x86linuxtarget/gcloud/gcloud-420.0.0/google-cloud-sdk/bin/gsutil cp $bucket_name/$slideseq_name/$slideseq_name_abbrev\.matched.digital_expression.txt.gz  $data_direc/$slideseq_name_abbrev\.matched.digital_expression.txt.gz')
        os.system('gunzip -f ./$slideseq_name_abbrev\.matched.digital_expression.txt.gz')
        os.system('!/apps/released/built-by-outside-authors-x86linuxtarget/gcloud/gcloud-420.0.0/google-cloud-sdk/bin/gsutil cp $bucket_name/$slideseq_name/$slideseq_name_abbrev\.matched.digital_expression.txt  $data_direc/$slideseq_name_abbrev\.matched.digital_expression.txt')

def make_puck(slideseq_name, directory_slideseq):
    puck = PuckReplicate(
        puck_name=slideseq_name,  # Name of the Slide-seq puck
        genome="hg19",  # Genome"
        clonotype_sample_name="",
        clonotype_mixcr_name="",  # Name of the MiXCR outputs
        clonotype_fastq_name="",  # Name of the rhPCR sequencing results
        umi_filter_number=
        1,  # We consider clonotypes on beads with > umi_filter_number UMIs.
        rctd_fn="",directory_slideseq = directory_slideseq,
        resave=True,  # See note below
        resave_sparse=True,
        file_type = 'slide-seq-tools',
        pub=True,
        skip_clonotype = True)  # Set note below
    return puck



class PuckReplicate:
    '''A class for each puck'''

    def __init__(self,
                 puck_name,
                 genome,
                 clonotype_sample_name,
                 clonotype_fastq_name,
                 umi_filter_number,
                 resave,
                 rctd_fn,
                 clonotype_mixcr_name,
                 directory_slideseq,
                 cluster_num_to_name_list=False,
                 cell_type_list=False,
                 abbrev_rep_name=False,
                 rcc=False,
                 resave_sparse=False,
                 skip_clonotype=False,
                 file_type='jilong',
                 pub=False):
        """Sets up variables"""
        self.puck_name = puck_name
        self.abbrev_sample_reference = "_".join(puck_name.split("_")[1:])
        self.genome = genome
        self.clonotype_sample_name = clonotype_sample_name
        self.clonotype_fastq_name = clonotype_fastq_name
        self.resave_state = resave
        self.abbrev_rep_name = abbrev_rep_name
        self.cluster_num_to_name_list = cluster_num_to_name_list
        self.clonotype_mixcr_name = clonotype_mixcr_name
        self.cell_type_list = cell_type_list
        self.umi_filter_number = umi_filter_number
        self.resave_sparse = resave_sparse
        self.rctd_fn = rctd_fn  ##### CHANGE HERE
        self.file_type = file_type
        if self.file_type == 'jilong':
            self.dge_filename = "{}/{}.digital_expression.txt".format(
                directory_slideseq, self.abbrev_sample_reference)

        elif self.file_type == 'slide-seq-tools':
            self.dge_filename = "{}/{}.matched.digital_expression.txt".format(
                directory_slideseq, self.abbrev_sample_reference)
        if self.file_type == 'jilong':
            self.spatial_barcodes_s1 = "{}/{}_matched_bead_barcodes.txt".format(
                directory_slideseq, self.abbrev_sample_reference)
            self.spatial_locations_s1 = "{}/{}_matched_bead_locations.txt".format(
                directory_slideseq, self.abbrev_sample_reference)
            self.bc_loc_dict_s1 = generate_barcode_loc_dictionary(
                self.spatial_barcodes_s1,
                self.spatial_locations_s1,
                old_status=False)
        elif self.file_type == 'slide-seq-tools':
            self.bc_loc_df = pd.read_csv('{}/{}_barcode_matching.txt'.format(
                directory_slideseq, self.abbrev_sample_reference),
                                         sep='\t',
                                         header=None)
            self.spatial_barcodes_s1 = list(self.bc_loc_df[0])
            self.spatial_locations_s1 = [
                (i, j) for i, j in zip(self.bc_loc_df[2], self.bc_loc_df[3])
            ]
            self.bc_loc_dict_s1 = {
                i: j
                for i, j in zip(self.spatial_barcodes_s1,
                                self.spatial_locations_s1)
            }
            self.matched_bead_barcodes = {
                i: j
                for i, j in zip(self.bc_loc_df[1], self.bc_loc_df[0])
            }

        self.loc_to_bc_s1 = {y: x for x, y in self.bc_loc_dict_s1.items()}

        if skip_clonotype is False:
            self.identify_clonotypes()
        self.open_dge()
        if rcc:
            self.run_rcc()

    def identify_clonotypes(self):
        """Identifies clonotypes from rhTCR sequencing"""
        if (path.exists("{}/cloneids_readids_{}.pickle".format(
                directory_slideseq, self.puck_name))
                is False) or (self.resave_state):
            self.cloneids, self.readids = read_in_clonotypes(
                sample_name=self.clonotype_mixcr_name)
            with open(
                    "{}/cloneids_readids_{}.pickle".format(
                        directory_slideseq, self.puck_name), "wb") as handle:
                pkl.dump([self.cloneids, self.readids],
                         handle,
                         protocol=pkl.HIGHEST_PROTOCOL)

        else:
            with open(
                    "{}/cloneids_readids_{}.pickle".format(
                        directory_slideseq, self.puck_name), "rb") as p_f:
                self.cloneids, self.readids = pkl.load(p_f)

        if path.exists("{}/readIDtobarcode_dict_{}.pickle".format(
                directory_slideseq,
                self.puck_name)) is False or (self.resave_state):
            self.readIDtobarcode_dict = readIDtobarcode(
                "{}/{}_R1_001.fastq".format(directory_slideseq,
                                            self.clonotype_fastq_name),
                up_filter=False,
            )
            with open(
                    "{}/readIDtobarcode_dict_{}.pickle".format(
                        directory_slideseq, self.puck_name), "wb") as handle:
                pkl.dump(self.readIDtobarcode_dict,
                         handle,
                         protocol=pkl.HIGHEST_PROTOCOL)
        else:
            with open(
                    "{}/readIDtobarcode_dict_{}.pickle".format(
                        directory_slideseq, self.puck_name), "rb") as p_f:
                self.readIDtobarcode_dict = pkl.load(p_f)

        self.tcr_bcumi_dict = tcr_bcumi_dict_f(self.readIDtobarcode_dict,
                                               self.cloneids, self.readids)

        self.clonotypes = {}
        self.clonotypes["1a"] = list(
            self.cloneids[self.cloneids.topChains == "TRA"].aaSeqImputedCDR3)
        self.clonotypes["1b"] = list(
            self.cloneids[self.cloneids.topChains == "TRB"].aaSeqImputedCDR3)

        self.tcr_bcumi_dict = merge_hamming_clonotype(self.tcr_bcumi_dict)

        filename = "{}/tcr_loc_dict_s1_{}.pickle".format(
            directory_slideseq, self.puck_name)

        if path.exists(filename) is False or (self.resave_state):
            self.tcr_loc_dict_s1, _ = clonotype_to_spatial(
                self.tcr_bcumi_dict,
                self.bc_loc_dict_s1,
                display_missed_barcodes=False,
                hamming_correction=True,
            )
            with open(filename, "wb") as handle:
                pkl.dump(self.tcr_loc_dict_s1,
                         handle,
                         protocol=pkl.HIGHEST_PROTOCOL)
        else:
            with open(filename, "rb") as p_f:
                self.tcr_loc_dict_s1 = pkl.load(p_f)

        self.tcr_loc_dict_s1 = merge_hamming_clonotype(self.tcr_loc_dict_s1)

        # Add back BC and UMIs of constant reads
        with open(
                "{}/trac_trbc2_constant_umis_{}.pickle".format(
                    directory_slideseq, self.puck_name), "rb") as p_f:
            self.trac_trbc2_constant = pkl.load(p_f)

        ### add back constant UMIs
        if self.genome == "hg19":
            tcra_gene_name = "TRAC"
        else:
            tcra_gene_name = "Trac"
        if self.genome == "hg19":
            tcrb_gene_name = "TRBC2"
        else:
            tcrb_gene_name = "Trbc2"
        all_constant_bc_umis = list(
            tuple(i) for i in self.trac_trbc2_constant[tcra_gene_name]) + list(
                tuple(i) for i in self.trac_trbc2_constant[tcrb_gene_name])
        all_constant_bc_umis = set(all_constant_bc_umis)

        self.tcr_bcumi_dict_filt4 = {}
        self.tcr_loc_dict_s1_filtered4 = {}

        for tcr in self.tcr_bcumi_dict:
            umis_to_include = set()
            cter = Counter(self.tcr_bcumi_dict[tcr])
            # Include if at least 2 reads
            umis_to_include.update(set(i for i in cter if cter[i] >= 2))
            #print(f'2 reads added to: {len(umis_to_include)}')
            # Include if at least 2 UMIs/bead
            deduped_bcs = [i[0] for i in cter]
            cter2 = Counter(deduped_bcs)
            bcs_to_include = set([i for i in cter2 if cter2[i] >= 2])
            umi_to_bc_dict = {i: i[0] for i in cter}
            umis_to_include.update(
                set(i for i in cter if umi_to_bc_dict[i] in bcs_to_include))
            #print(f'2 umis added to: {len(umis_to_include)}')
            # TODO: Include if in constant
            umis_to_include.update(
                set(i for i in cter if i in all_constant_bc_umis))

            self.tcr_bcumi_dict_filt4[tcr] = list(set(umis_to_include))

        self.tcr_loc_dict_s1_filtered4, errors = clonotype_to_spatial(
            self.tcr_bcumi_dict_filt4,
            self.bc_loc_dict_s1,
            display_missed_barcodes=False,
            hamming_correction=True,
        )
        self.tcr_loc_dict_s1_filtered4 = merge_hamming_clonotype(
            self.tcr_loc_dict_s1_filtered4)

        if path.exists("{}/dfs1_filtered_{}.csv".format(
                directory_out,
                self.puck_name)) is False or (self.resave_state):
            self.df_s1_filtered4 = save_to_csv(
                "{}/dfs1_filtered_{}.csv".format(directory_out,
                                                 self.puck_name),
                self.tcr_loc_dict_s1_filtered4,
                self.clonotypes,
                self.bc_loc_dict_s1,
                save=True,
            )
        else:
            self.df_s1_filtered4 = pd.read_csv(
                "{}/dfs1_filtered_{}.csv".format(
                    directory_out,
                    self.puck_name))  # s1, s08, s9 & _filtered mods

    def open_dge(self):
        """Opens DGE and saves a sparse matrix if needed"""
        if self.resave_sparse:
            save_as_sparse(
                self.dge_filename,
                self.abbrev_sample_reference,
                skip_lines=0,
            )
        self.s1_dge_fmtd = open_sparse_matrix(
            "{}/{}_sparse_matrix.npz".format(directory_slideseq,
                                             self.abbrev_sample_reference),
            self.dge_filename,
            sample_abbrev=self.abbrev_sample_reference)

        self.s1_dge_fmtd_norm = self.s1_dge_fmtd.div(
            self.s1_dge_fmtd.sum(axis=1), axis=0)

        if self.file_type == 'jilong':
            self.s1_dge_fmtd, x = dge_formatting(self.s1_dge_fmtd,
                                                 self.spatial_barcodes_s1,
                                                 self.spatial_locations_s1)
        if self.file_type == 'slide-seq-tools':
            self.s1_dge_fmtd = dge_formatting_sst(self.s1_dge_fmtd,
                                                  self.spatial_barcodes_s1,
                                                  self.spatial_locations_s1)
        self.s1_total_counts = self.s1_dge_fmtd.sum(axis=1)

    def run_rcc(self):
        """Analysis for RCC samples"""
        self.cluster_labels = make_knn_assignments(
            "{}/{}".format(directory_slideseq, self.rctd_fn),
            self.abbrev_rep_name,
            self,
            self.cluster_num_to_name_list,
            self.cell_type_list,
        )

        # Identify tumor boundary
        d = 50
        distarray = cdist(
            self.tumor_beads,
            [i for i in self.tll_loc_cluster_dict])  # tumor beads vs all beads

        all_beads = [i for i in self.tll_loc_cluster_dict]
        boundary_points = []

        skip_these_beads = set(
            []
        )  # for the tumor beads that are just too far away from the boundary

        tumor_beads_loc_to_index = {}
        for i in range(len(self.tumor_beads)):
            tumor_beads_loc_to_index[self.tumor_beads[i]] = i

        for each_loc in tqdm(self.tumor_beads):

            if each_loc in skip_these_beads:
                continue

            nearby_cell_types = [0,
                                 0]  # First is tumor, second is anything else
            # Get nearest points within distance

            all_beads_indices = [
                i for i, x in enumerate(distarray[
                    self.tumor_beads.index(each_loc),
                ]) if x < d
            ]
            cell_types = [
                self.tll_loc_cluster_dict[all_beads[idx]]
                for idx in all_beads_indices
            ]
            for i in cell_types:
                if i == "Tumor":
                    nearby_cell_types[0] += 1
                else:
                    nearby_cell_types[1] += 1
            if nearby_cell_types[0] / sum(nearby_cell_types) <= 0.9:
                boundary_points.append(each_loc)

            else:
                skip_these_beads |= set([
                    all_beads[i] for i in all_beads_indices
                    if all_beads[i] in self.tumor_beads
                ])
        self.boundary_points = boundary_points
        self.bc_to_cluster_label_dict = dict(
            zip(self.cluster_labels.barcode, self.cluster_labels.cluster))

        ##### Randomness between lung, tumor, and TIL

        self.locations_reference = switch_dict(
            copy.deepcopy(self.tll_loc_cluster_dict))

        self.ltt_dist_dict = {}
        for tcr in tqdm(self.tcr_loc_dict_s1_filtered4):
            if tcr not in self.clonotypes["1b"]:
                continue
            ltt_dist = get_locations_distribution(
                set(self.tcr_loc_dict_s1_filtered4[tcr]),
                self.locations_reference)
            if sum(ltt_dist) < 10:
                continue

            self.ltt_dist_dict[tcr] = ltt_dist

        ##### EXPECTED VALUES

        # Sampling expected only from tcr locations
        all_tcr_locs = []
        for clonotype in self.tcr_loc_dict_s1_filtered4:
            all_tcr_locs = all_tcr_locs + list(
                set(self.tcr_loc_dict_s1_filtered4[clonotype]))
        expected = get_locations_distribution(all_tcr_locs,
                                              self.locations_reference)

        #### SET UP CLONO PLOT DF FOR Clonotype Analysis Plot
        cl_all = [
            i for i in self.tcr_loc_dict_s1_filtered4
            if ((len(set(self.tcr_loc_dict_s1_filtered4[i])) >= 10) and (
                i in self.clonotypes["1b"]))
        ]
        num_locs = [len(list(set(self.tcr_loc_dict_s1[i]))) for i in cl_all]
        self.x_ltt = [
            self.ltt_dist_dict[clonotype][0] /
            sum(self.ltt_dist_dict[clonotype]) / expected[0] * sum(expected) -
            1 for clonotype in cl_all
        ]  # lung
        self.y_ltt = [
            self.ltt_dist_dict[clonotype][1] /
            sum(self.ltt_dist_dict[clonotype]) / expected[1] * sum(expected) -
            1 for clonotype in cl_all
        ]  # til
        self.z_ltt = [
            self.ltt_dist_dict[clonotype][2] /
            sum(self.ltt_dist_dict[clonotype]) / expected[2] * sum(expected) -
            1 for clonotype in cl_all
        ]  # tumor

        clono_plot_df = pd.DataFrame.from_dict({
            "Lung": self.x_ltt,
            "TIL": self.y_ltt,
            "Tumor": self.z_ltt,
            "clonotype": cl_all,
            "num_points": num_locs,
        })
        self.clono_plot_df = clono_plot_df
        
##


def cohens_d(i, j):
    """Calculates the Cohen's d effect size between two groups"""
    return (mean(i) - mean(j)) / (math.sqrt((stdev(i)**2 + stdev(j)**2) / 2))


def abline(slope, intercept, ax_name=None):
    """Plot a line from slope and intercept"""
    if ax_name is None:
        axes = plt.gca()
        x_vals = np.array(axes.get_xlim())
        y_vals = intercept + slope * x_vals
        plt.plot(x_vals, y_vals, "--")
    else:
        x_vals = np.array(ax_name.get_xlim())
        y_vals = intercept + slope * x_vals
        ax_name.plot(x_vals, y_vals, "--")


def closest_node(node, nodes):
    """Finds the closest point in a list of points to a given point"""
    closest_index = distance.cdist([node], nodes).argmin()
    return nodes[closest_index]


def distance_between_points(pt1, pt2):
    """Distance between two points"""
    return math.sqrt((pt1[0] - pt2[0])**2 + (pt1[1] - pt2[1])**2)


def switch_dict(orig_dict):
    """Switches the keys and values of a dictionary and merges to a list"""
    new_dict = {}
    for i in orig_dict:
        if orig_dict[i] in new_dict:
            new_dict[orig_dict[i]].append(i)
        else:
            new_dict[orig_dict[i]] = [i]
    return new_dict


def generate_barcode_loc_dictionary(barcode_file,
                                    locations_file,
                                    old_status=False):
    """Creates barcode and location dictionary from Slide-seq data"""
    barcode_loc_dict = {}

    if old_status:
        with open(barcode_file, "r") as file1, open(locations_file,
                                                    "r") as file2:
            for line1, line2 in zip(file1, file2):
                line1 = line1.rstrip("\n")
                line2 = line2.rstrip("\n")
                line2_split = list(line2.split("\t"))
                barcode_loc_dict[line1] = (line2_split[1], line2_split[2])

    else:
        barcodes = []
        with open(barcode_file, "r") as file1, open(locations_file,
                                                    "r") as file2:
            for line in file1:
                line = line.rstrip("\n")
                line = line.replace(",", "")
                barcodes.append(line)

        with open(locations_file, "r") as file1:
            locations = []
            for line in file1:
                line = line.rstrip("\n")
                line_split = list(line.split("\t"))
                locations.append((line_split[1], line_split[2]))

    file1.close()
    file2.close()

    for i, j in enumerate(barcodes):
        barcode_loc_dict[j] = (float(locations[i][0]), float(locations[i][1]))

    return barcode_loc_dict


def hamming_distance(s1, s2):
    """Reports if hamming distance between two strings is 1"""
    if len(s1) != len(s2):
        return 'Fail'
    return sum(ch1 != ch2 for ch1, ch2 in zip(s1, s2))


def clonotype_to_spatial(
    clonotype_to_barcode,
    barcode_to_spatial,
    display_missed_barcodes=True,
    hamming_correction=True,
    verbose=False,
):
    hamming_dict = {}
    """Returns dictionary of clonotypes to all spatial locations"""
    cl_to_spatial = {}
    # For each clonotype
    print(len(clonotype_to_barcode), "total to analyze")

    errors = []

    all_barcodes = set(barcode_to_spatial)
    all_barcodes_firsthalf = {}
    all_barcodes_secondhalf = {}
    for bc in all_barcodes:
        if bc[0:7] not in all_barcodes_firsthalf:
            all_barcodes_firsthalf[bc[0:7]] = set()
        all_barcodes_firsthalf[bc[0:7]].add(bc)
        if bc[7:] not in all_barcodes_secondhalf:
            all_barcodes_secondhalf[bc[7:]] = set()
        all_barcodes_secondhalf[bc[7:]].add(bc)

    for tcr_clonotype in tqdm(clonotype_to_barcode):
        missed_barcodes = []
        bc_umis_per_cl = set(clonotype_to_barcode[tcr_clonotype])

        if tcr_clonotype not in cl_to_spatial:
            cl_to_spatial[tcr_clonotype] = []

        for read in bc_umis_per_cl:
            bc = read[0][0:14]
            if bc in barcode_to_spatial:
                cl_to_spatial[tcr_clonotype].append(barcode_to_spatial[bc])
                continue
            else:
                if display_missed_barcodes:
                    if bc not in missed_barcodes:
                        print(bc,
                              " was not found in in situ barcode sequencing")
                        missed_barcodes.append(bc)

                if hamming_correction:
                    all_barcodes_searchspace = set()
                    if bc[7:] in all_barcodes_secondhalf:
                        all_barcodes_searchspace.update(
                            all_barcodes_secondhalf[bc[7:]])
                    if bc[0:7] in all_barcodes_firsthalf:
                        all_barcodes_searchspace.update(
                            all_barcodes_firsthalf[bc[0:7]])
                    all_barcodes_pass_hamming = [
                        bc2 for bc2 in all_barcodes_searchspace
                        if hamming_distance(bc, bc2) == 1
                    ]

                    if len(all_barcodes_pass_hamming) > 1:
                        if verbose:
                            print("bc:", bc, "match?",
                                  all_barcodes_pass_hamming)
                            print("error not 1 d apart")
                        errors.append(all_barcodes_pass_hamming)
                        continue
                    if len(all_barcodes_pass_hamming) == 1:
                        cl_to_spatial[tcr_clonotype].append(
                            barcode_to_spatial[all_barcodes_pass_hamming[0]])

                continue
    return cl_to_spatial, errors


def cl_loc_convert_loc_cl(tcr_loc_dict):
    """Converts clonotype-location dictionary to location-clonotype"""
    loc_to_tcr = {}  # Dictionary of locations to clonotypes
    for tcr in tcr_loc_dict:
        loc_all = tcr_loc_dict[tcr]
        for loc in loc_all:
            if loc in loc_to_tcr:
                loc_to_tcr[loc].append(tcr)
            else:
                loc_to_tcr[loc] = [tcr]

    # Remove duplicates
    loc_to_tcr_dedup = {}
    for loc in loc_to_tcr:
        loc_to_tcr_dedup[loc] = set(loc_to_tcr[loc])

    return loc_to_tcr_dedup


def bc_loc_convert_loc_bc(bc_loc_dict):
    """Converts bc-location dictionary to location-bc"""
    loc_bc_dict = {}
    for bc in bc_loc_dict:
        loc_all = bc_loc_dict[bc]

        if loc_all not in loc_bc_dict:
            loc_bc_dict[loc_all] = []

        loc_bc_dict[loc_all].append(bc)

    return loc_bc_dict


def check_clonotype(cl_string, a_clonotypes_list, b_clonotypes_list):
    """Checks if a clonotype is in alpha or beta list"""
    if cl_string in a_clonotypes_list:
        return "a"
    if cl_string in b_clonotypes_list:
        return "b"


def save_to_csv(filename,
                tcr_loc_dict,
                clonotypes,
                bc_umi_dict,
                output_dataframe=True,
                save=True):
    """Output a csv file with bead barcodes, locations, and clonotypes
    of tcr_a and tcr_b"""
    loc_to_tcr = cl_loc_convert_loc_cl(
        tcr_loc_dict)  # Convert tcr:locations to locations:tcr

    # Make list of all locations
    locations = loc_to_tcr.keys()
    loclist = list(locations)

    # Make list of all clonotypes at each location
    clonotypes_all = []

    for loc in loclist:
        clonotypes_all.append(loc_to_tcr[loc])

    # Separate list of all clonotypes into tcr_a and tcr_b
    clonotypes_a = []
    clonotypes_b = []
    for i in clonotypes_all:
        insert_a = []
        insert_b = []
        for j in i:
            outcome = check_clonotype(j, clonotypes["1a"], clonotypes["1b"])
            if outcome == "a":
                insert_a.append(j)
                insert_b.append("")
            elif outcome == "b":
                insert_b.append(j)
                insert_a.append("")

        clonotypes_a.append(insert_a)
        clonotypes_b.append(insert_b)

    # Make list of bead_barcodes at each location
    bead_barcodes = []
    for i in loclist:
        loc_bc = bc_loc_convert_loc_bc(bc_umi_dict)
        bead_barcodes.append(loc_bc[i])

    # Save as csv
    if save:
        with open(filename, "w") as csvfile:
            fieldnames = ["x", "y", "tcr_a", "tcr_b", "bead_barcodes"]

            writer = csv.DictWriter(csvfile, fieldnames=fieldnames)
            writer.writeheader()
            for i, j in enumerate(loclist):
                writer.writerow({
                    "x": j[0],
                    "y": j[1],
                    "tcr_a": clonotypes_a[i],
                    "tcr_b": clonotypes_b[i],
                    "bead_barcodes": bead_barcodes[i][0],
                })

    # Output result as dataframe
    if output_dataframe:
        d = {
            "x": [i[0] for i in loclist],
            "y": [i[1] for i in loclist],
            "tcr_a": clonotypes_a,
            "tcr_b": clonotypes_b,
            "bead_barcodes": [i[0] for i in bead_barcodes],
        }

        df = pd.DataFrame(data=d)
        return df


def import_dge(dge_filename, skip_lines=0, print_progress=True):
    """Imports a dge as a dictionary and barcodes"""
    csv.field_size_limit(sys.maxsize)
    cter = 0

    reader = csv.reader(open(dge_filename), delimiter='\t')

    result = {}
    cter = 0
    for row in reader:
        cter += 1
        key = row[0]
        if key in result:
            # implement your duplicate row handling here
            pass
        result[key] = row[1:]
        if cter % 1000 == 0 and print_progress:
            print(cter)

    barcodes = list(result['GENE'])
    result.pop('GENE')
    for item in result:
        result[item] = [int(i) for i in result[item]]

    return result, barcodes


def dge_formatting(dge, barcodes_fn, locations_fn):
    """Formats a dge as a dataframe and returns coordinates for locations"""
    coords_barcodes = pd.read_csv(barcodes_fn, sep="\t", header=None)
    coords_barcodes = [i.replace(",", "")
                       for i in list(coords_barcodes[0])]  # Remove commas
    coords_locations = pd.read_csv(locations_fn, sep="\t", header=None)
    coords = pd.DataFrame(coords_locations)
    coords.set_axis(["id", "xcoord", "ycoord"], axis=1, inplace=True)
    coords["barcodes"] = coords_barcodes
    coords.set_index("barcodes", inplace=True)

    barcodes_in_dge = dge.index
    barcodes_in_coords = coords.index
    common_barcodes = [i for i in barcodes_in_dge if i in barcodes_in_coords]
    dge = dge.loc[common_barcodes]
    coords = coords.loc[common_barcodes]

    return dge, coords


def dge_formatting_sst(dge, barcodes, locations):
    """Formats a dge as a dataframe and returns coordinates for locations"""
    coords = pd.DataFrame({'barcodes': barcodes})
    coords['xcoord'] = [i[0] for i in locations]
    coords['ycoord'] = [i[1] for i in locations]
    coords.set_index("barcodes", inplace=True)

    barcodes_in_dge = dge.index
    barcodes_in_coords = coords.index
    #common_barcodes = [i for i in barcodes_in_dge if i in barcodes_in_coords]
    #dge = dge.loc[common_barcodes]
    #     coords = coords.loc[common_barcodes]

    return dge


def find_specific_barcodes(which_tcr, df_scb):
    """Finds the barcodes of a certain clonotype"""

    tcr_spec = df_scb[which_tcr]

    tcr_spec_counter = []
    ct = 0
    for i in tcr_spec:

        if isinstance(i, list):
            cl = i
        else:
            cl = ast.literal_eval(i)
        cl = [n.strip() for n in cl]
        filt_cl = [j for j in cl if j != ""]
        if len(filt_cl) == 0:
            tcr_spec_counter.append(-1)
        else:
            tcr_spec_counter.append(ct)
        ct += 1

    tcr_specific_barcodes = [
        df_scb["bead_barcodes"][b] for b in tcr_spec_counter if b != -1
    ]

    print(len(tcr_specific_barcodes))
    return tcr_specific_barcodes


def open_sparse_matrix(filename, dge_filename, sample_abbrev):
    """Opens a previously saved sparse matrix as a dataframe, used for DGE"""
    sparse_matrix = scipy.sparse.load_npz(filename)
    sparse_matrix_dense = sparse_matrix.todense()
    sparse_matrix_dense = sparse_matrix_dense.astype(float)

    if (path.exists("{}/{}_genenames_barcodes.csv".format(
            directory_slideseq, sample_abbrev)) is False):
        new_dict, barcodes = import_dge(dge_filename, print_progress=True)
        df_dge = pd.DataFrame(data=new_dict, index=barcodes)
        df_dge.index.name = "barcodes"
        dge_fmtd = df_dge
        cols = [list(dge_fmtd.columns)]
        rows = [list(dge_fmtd.index)]
        f = pd.DataFrame({"Genes": cols, "Barcodes": rows})
        f.to_csv(
            "{}/{}_genenames_barcodes.csv".format(directory_slideseq,
                                                  sample_abbrev),
            header=True,
            index=False,
        )

    rows_cols = pd.read_csv("{}/{}_genenames_barcodes.csv".format(
        directory_slideseq, sample_abbrev))
    cols = list(rows_cols["Genes"])
    cols = [ast.literal_eval(i) for i in cols][0]
    rows = list(rows_cols["Barcodes"])
    rows = [ast.literal_eval(i) for i in rows][0]
    new_df = pd.DataFrame(data=sparse_matrix_dense, columns=cols, index=rows)
    return new_df


def save_as_sparse(filename, abbrev_sample_reference, skip_lines=0):
    """Saves a dataframe as a sparse matrix, used for DGE"""
    new_dict, barcodes = import_dge(filename, skip_lines, print_progress=True)
    df_dge = pd.DataFrame(data=new_dict, index=barcodes)
    df_dge.index.name = "barcodes"
    dge_fmtd = df_dge
    cols = [list(dge_fmtd.columns)]
    rows = [list(dge_fmtd.index)]
    f = pd.DataFrame({"Genes": cols, "Barcodes": rows})
    f.to_csv(
        "{}/{}_genenames_barcodes.csv".format(directory_slideseq,
                                              abbrev_sample_reference),
        header=True,
        index=False,
    )
    array = dge_fmtd.values
    sparse_matrix = scipy.sparse.csc_matrix(array)
    scipy.sparse.save_npz(
        "{}/{}_sparse_matrix.npz".format(directory_slideseq,
                                         abbrev_sample_reference),
        sparse_matrix,
    )


def read_in_clonotypes(sample_name, full_directory=False):
    """Reads in clonotypes from MiXCR output"""
    if full_directory is False:
        fn1 = "{}{}_paired_clones.txt".format(directory, sample_name)
        fn2 = "{}{}_paired_cloneID.txt".format(directory, sample_name)
    if full_directory:
        fn1 = "{}{}_paired_clones.txt".format(directory, sample_name)
        fn2 = "{}{}_paired_cloneID.txt".format(directory, sample_name)
    cloneids = pd.read_csv(fn1, sep="\t")
    readids = pd.read_csv(fn2, sep="\t", error_bad_lines=False)
    readids = readids[(readids.topChains == "TRB") |
                      (readids.topChains == "TRA")]
    cloneids = cloneids[(cloneids.topChains == "TRB") |
                        (cloneids.topChains == "TRA")]
    return cloneids, readids


def readIDtobarcode(fastq_file_r1, up_filter):
    """Converts read IDs from MiXCR output to barcodes"""
    readIDtobarcode = {}
    line_count = 0
    all_read_count = 0
    up_present_ct = 0
    up_not_present_ct = 0
    up_shifted_ct = 0
    with open(fastq_file_r1, "r") as file1:
        for line1 in file1:
            if line_count % 4 == 0:  # readID name
                line1 = line1.strip()
                read_name = line1[1:].split(' ')[0]
            if line_count % 4 == 1:  # Check if fastq sequence
                all_read_count += 1
                barcode = line1[0:8] + line1[26:32]
                if line1[8:26] == "TCTTCAGCGTTCCCGAGA":
                    up_present_ct += 1

                else:

                    if "TCTTCAGCGTTCCCGAGA" in line1:
                        up_shifted_ct += 1
                    else:
                        up_not_present_ct += 1
                    if up_filter:
                        line_count += 1
                        continue

                umi = line1[32:41]

                readIDtobarcode[read_name] = [barcode, umi]
            line_count += 1
            if line_count % 1000000 == 0:
                print(line_count)
                print(datetime.now())
    file1.close()
    print("Total reads:", all_read_count)
    print(
        "UP present reads:",
        up_present_ct,
        "UP not present reads:",
        up_not_present_ct,
        ". UP shifted: ",
        up_shifted_ct,
    )
    return readIDtobarcode


def tcr_bcumi_dict_f(readbarcode, cloneids, readids):
    """Makes a dictionary of TCR clonotypes to barcode-UMIs"""
    aa_to_clid = {}
    aa_to_tcr = {}
    cloneids_dict = cloneids.to_dict()
    for i in cloneids_dict['cloneId'].keys():
        clid = cloneids_dict['cloneId'][i]
        aa_seq = cloneids_dict['aaSeqImputedCDR3'][i]
        tcr_seq = cloneids_dict['targetSequences'][i]
        if aa_seq not in aa_to_clid:
            aa_to_clid[aa_seq] = []
        aa_to_clid[aa_seq].append(clid)

        if aa_seq not in aa_to_tcr:
            aa_to_tcr[aa_seq] = []
        aa_to_tcr[aa_seq].append(tcr_seq)

    clid_to_r1name = {}
    r1name_all = [i.split(' ')[0] for i in list(readids["descrsR1"])]
    clid_all = list(readids["cloneId"])
    for key, value in tqdm(zip(clid_all, r1name_all)):
        if key in clid_to_r1name:
            clid_to_r1name[key].append(value)
        else:
            clid_to_r1name[key] = [value]
    cl_aa_to_barcode_umis = {}
    for aa in tqdm(aa_to_clid):
        if aa not in cl_aa_to_barcode_umis:
            cl_aa_to_barcode_umis[aa] = []
        for clid in aa_to_clid[aa]:
            if clid not in clid_to_r1name:
                continue
            for r1name in clid_to_r1name[clid]:
                if r1name in readbarcode:
                    cl_aa_to_barcode_umis[aa].append(tuple(
                        readbarcode[r1name]))

    return cl_aa_to_barcode_umis


def merge_hamming_clonotype(clonotype_dict):
    """Merges clonotypes with a hamming distance of 1 into
    the clonotype with more UMIs"""
    hamming_sets = []
    found = False
    clonotype_dict_hamming = {}
    for i in clonotype_dict:
        found = False
        for hset in hamming_sets:
            for cl in hset:
                if hamming_distance(i, cl) == 1:
                    hset.append(i)
                    found = True
                    break
            if found:
                break
        if found:
            continue
        hamming_sets.append([i])

    for hset in hamming_sets:
        hset_lengths = [len(clonotype_dict[i]) for i in hset]
        pref_clonotype = hset[hset_lengths.index(max(hset_lengths))]

        clonotype_dict_hamming[pref_clonotype] = []
        for cl in hset:
            clonotype_dict_hamming[pref_clonotype] = (
                clonotype_dict_hamming[pref_clonotype] + clonotype_dict[cl])
    return clonotype_dict_hamming


def plot_variable_vs_constant(dge_output,
                              tcr_type,
                              tcr_loc_dict_s1,
                              loc_to_bc_s1,
                              human=False,
                              cutoff=15,
                              plot_histogram=True,
                              save=False):
    """Plots the variable recovery"""
    if tcr_type == "tcr_b":
        gene_name = "Trbc2"
    if tcr_type == "tcr_a":
        gene_name = "Trac"
    if human:
        if tcr_type == "tcr_b":
            gene_name = "TRBC2"
        if tcr_type == "tcr_a":
            gene_name = "TRAC"

    bc_umis = (
        {}
    )  # This converts dictionary to barcodes:UMIs instead of clonotypes:[bc,umi]
    for cl in tcr_loc_dict_s1:
        for loc in tcr_loc_dict_s1[cl]:
            bc = loc_to_bc_s1[loc]
            if bc not in bc_umis:
                bc_umis[bc] = []
            bc_umis[bc].append(loc)
    # tcr_bc = list(df_output['bead_barcodes'])

    tcr_barcodes = list(dge_output[dge_output[gene_name] >=
                                   1].index)  # all constant tcr barcodes

    # for every barcode record constant counts
    cons_counts = []
    for j in tcr_barcodes:
        cons_counts.append(dge_output.loc[j][gene_name])

    var_counts = []  # record umi counts
    for i in tcr_barcodes:
        if i in bc_umis:
            var_counts.append(len(bc_umis[i]))
        if i not in bc_umis:
            var_counts.append(0)

    frac_count = [j / i for i, j in zip(cons_counts, var_counts)]

    frac_count_adj = []

    for i in frac_count:
        if i > 0:
            frac_count_adj.append(1)
        else:
            frac_count_adj.append(0)

    ##### groups counts based on constant counts

    t = int(max(cons_counts))  # maximum number of constant counts

    output = [0] * t  # Summed beads passing test for this bin
    counts = [0] * t  # Number of bead barcodes for each bin

    for i, j in enumerate(cons_counts):
        idxcount = int(j) - 1
        output[idxcount] += frac_count_adj[i]
        counts[idxcount] += 1

    frac_fin = []
    cutoff_sum = [[], []]  # output, beads
    for i, j in enumerate(output):
        if i < cutoff:
            if counts[i] == 0:  # If there are no beads, append 0
                frac_fin.append(0)
            else:
                frac_fin.append(
                    float(j) / float(counts[i])
                )  # number of beads passing test divided by total constant beads
        else:
            cutoff_sum[0].append(j)
            cutoff_sum[1].append(counts[i])
    frac_fin.append([sum(cutoff_sum[0]) / sum(cutoff_sum[1])])
    _, ax = plt.subplots(figsize=(cutoff, 5))
    print(cutoff)
    if tcr_type == "tcr_a":
        color_to_use = "#9FB7CD"
    if tcr_type == "tcr_b":
        color_to_use = "#0F4C81"

    ax.bar(np.array([i + 1 for i in range(cutoff + 1)]),
           frac_fin,
           color=color_to_use)
    labels = [str(i) for i in counts]
    updated_labels = []
    for v, i in enumerate(frac_fin):
        if v < cutoff:
            ax.text(v + 1,
                    1.03,
                    labels[v],
                    ha="center",
                    va="bottom",
                    rotation="vertical")
            updated_labels.append(int(labels[v]))
        if v == cutoff:
            ax.text(
                v + 1,
                1.03,
                str(sum(cutoff_sum[1])),
                ha="center",
                va="bottom",
                rotation="vertical",
            )
            updated_labels.append(sum(cutoff_sum[1]))
    plt.ylim([0, 1])
    if tcr_type == "tcr_b":
        #         plt.yticks([])
        plt.ylabel("")
        plt.xlabel("Number of constant reads")

    else:
        plt.ylabel("Fraction")
        plt.ylabel("")
        plt.yticks([])
        plt.xlabel("")
    plt.xticks([i + 1 for i in range(cutoff + 1)])
    print(frac_fin)
    if plot_histogram:
        plt.figure(figsize=(cutoff, 1))
        print(updated_labels)
        plt.plot(range(len(updated_labels)), updated_labels)
        plt.axis("off")
    if save == True:
        plt.savefig(f'./plots/{tcr_type}.pdf')


def cdist(i, j):
    """Calculates the distance between every pair of points in two lists"""
    return scipy.spatial.distance.cdist(i, j)


def find_le(a, x):
    """Finds index of rightmost value less than or equal to x"""
    i = bisect.bisect_right(a, x)
    if i:
        return i - 1
    return 0


def find_ge(a, x):
    """Finds index of leftmost value greater than or equal to x"""
    i = bisect.bisect_left(a, x)
    if i <= len(a):
        return i
    raise ValueError


def bc_list_to_clusters_vector(barcode_list,
                               cell_type_list,
                               bc_to_cluster_label_dict,
                               verbose=False):
    """Makes barcodes to clusters dictionary"""
    clusters = [0 for i in range(len(cell_type_list))]
    for barcode in barcode_list:
        if barcode not in bc_to_cluster_label_dict:
            if verbose:
                print(barcode, "not in cluster labels")
            continue

        cl = bc_to_cluster_label_dict[barcode]
        clusters[cl] += 1
    return clusters


def make_knn_assignments(rctd_filename, puck_abbreviation, puck,
                         cluster_num_to_name_list, cell_type_list):
    """Generates K-nearest neighbors assignments for compartments on puck"""
    puck.all_locs = list(puck.loc_to_bc_s1.keys())
    aggregated_clusters = pd.read_csv(rctd_filename)
    rep = [i.split("_")[0] for i in list(aggregated_clusters.barcode)]
    barcode = [i.split("_")[1] for i in list(aggregated_clusters.barcode)]
    aggregated_clusters["rep"] = rep
    aggregated_clusters["barcode"] = barcode

    cluster_labels = aggregated_ctcrlusters[aggregated_clusters.rep ==
                                            puck_abbreviation]

    cluster_labels["x"] = [
        puck.bc_loc_dict_s1[bc][0] for bc in cluster_labels.barcode
    ]
    cluster_labels["y"] = [
        puck.bc_loc_dict_s1[bc][1] for bc in cluster_labels.barcode
    ]

    bc_to_cluster_label_dict = dict(
        zip(cluster_labels.barcode, cluster_labels.cluster))
    loc_to_cluster_dict = {
        puck.bc_loc_dict_s1[i]: bc_to_cluster_label_dict[i]
        for i in list(bc_to_cluster_label_dict.keys())
    }

    cluster_labels["cluster_name"] = [
        cluster_num_to_name_list[i] for i in cluster_labels.cluster
    ]
    # Load in cell type assignments for posttreatment
    # Perform KNN to assign tumor, lung, and TIL regions

    ###### KNN Classifier
    loc_cluster_dict = {}
    loc_cluster_dict_all = {}

    # cell_type_list = ['Tumor','Normal']
    # cell_type_list = ['TumorCell_CA9','NonImmuneCell','PlasmaCell']

    for index, row in cluster_labels.iterrows():
        loc_cluster_dict_all[(float(row["x"]),
                              float(row["y"]))] = row["cluster_name"]
        if row["cluster_name"] in cell_type_list:
            loc_cluster_dict[(float(row["x"]),
                              float(row["y"]))] = row["cluster_name"]

    unassigned_points = set(
        [i for i in puck.all_locs if i not in set(loc_cluster_dict.keys())])

    neigh = KNeighborsClassifier(
        n_neighbors=500
    )  # nearest neighbors was 500 for post treatment. trying 10
    print("Begin classifier process")
    X = np.array(list(loc_cluster_dict.keys()))
    y = np.array(list(loc_cluster_dict.values()))
    neigh.fit(X, y)
    cl_assignment = {}
    for pt in tqdm(unassigned_points):
        cl_assignment[pt] = neigh.predict([pt])[0]
    for pt in tqdm(loc_cluster_dict):
        cl_assignment[pt] = neigh.predict([pt])[0]

    # Plots results of classifier
    testl = []
    testtil = []
    testtu = []
    for loc in cl_assignment:
        if cl_assignment[loc] == cell_type_list[0]:
            testl.append(loc)
        if cl_assignment[loc] == cell_type_list[1]:
            testtu.append(loc)
        if cl_assignment[loc] == cell_type_list[2]:
            testtil.append(loc)

    loc_cluster_dict = copy.deepcopy(cl_assignment)
    puck.tll_loc_cluster_dict = loc_cluster_dict
    puck.tumor_beads = [
        i for i in loc_cluster_dict if loc_cluster_dict[i] == "Tumor"
    ]
    puck.loc_to_cluster_dict = loc_to_cluster_dict
    return cluster_labels


def get_locations_distribution(all_locations, locations_reference):
    """Gets the cell types for a list of locations"""
    lung_locs = set(locations_reference["Lung"])
    til_locs = set(locations_reference["TIL Chemokines"])
    tumor_locs = set(locations_reference["Tumor"])
    ltt_dist = [0, 0, 0]

    for location in all_locations:
        if location in lung_locs:
            ltt_dist[0] += 1

        elif location in til_locs:
            ltt_dist[1] += 1

        elif location in tumor_locs:
            ltt_dist[2] += 1
    return ltt_dist


def normalize_by_total_reads(list_name, total_reads):
    ''' Normalize a list by the total number of reads'''
    return [i / total_reads for i in list_name]


def clono_vd(pre_list, post_list, labels, filename):
    ''' Makes a venn diagram for comparing pre and post-treatment clonotypes'''
    plt.rcParams.update({"font.size": 12})
    set1 = set(pre_list)
    set2 = set(post_list)

    v = venn2([set1, set2], labels)
    c = venn2_circles([set1, set2])
    c[0].set_color("white")
    c[0].set_ls("solid")
    c[0].set_edgecolor("blue")

    c[1].set_color("white")
    c[1].set_ls("solid")
    c[1].set_edgecolor("red")
    plt.savefig(filename)
    plt.show()


def add_cp_puck(puck, cl_agg_enrichment):
    '''Adds a clono_plot from a puck to an aggregated enrichment'''
    for cl in puck.clono_plot_df.clonotype:
        for c in compartments:
            if cl not in cl_agg_enrichment[c]:
                cl_agg_enrichment[c][cl] = []
            cl_agg_enrichment[c][cl].append(
                float(
                    puck.clono_plot_df[puck.clono_plot_df.clonotype == cl][c]))
    return cl_agg_enrichment


def make_plotting_df():
    '''Makes plotting df for violin plots below and above median distance'''
    pval_recorder = {}
    plt.rcParams.update({"font.size": 22})
    plotting_df = {}
    plotting_df["clonotype"] = []
    plotting_df["greater_or_less"] = []
    plotting_df["expression"] = []

    for cl in list(gene_group_cl_below_median.keys()
                   ):  ### changed to different group to test
        for i in range(len(gene_group_cl_below_median[cl])):
            plotting_df["clonotype"].append(cl)
            plotting_df["greater_or_less"].append("less")
            plotting_df["expression"].append(gene_group_cl_below_median[cl][i])

        for i in range(len(gene_group_cl_above_median[cl])):
            plotting_df["clonotype"].append(cl)
            plotting_df["greater_or_less"].append("greater")
            plotting_df["expression"].append(gene_group_cl_above_median[cl][i])

        pval = kstest(gene_group_cl_below_median[cl],
                      gene_group_cl_above_median[cl])[1]
        pval_recorder[cl] = (pval)
        if pval < 0.05:
            print("KS", cl, group, pval)  #pretty_clonotype_name[cl],

    plotting_df = pd.DataFrame(plotting_df)
    return plotting_df, pval_recorder


def find_highest_expressing_clonotype(barcode,
                                      list_of_clonotypes,
                                      verbose=False):
    '''Finds highest expressing clonotype'''
    errors = []
    counts = 0
    cl_keep = ''
    for cl in list_of_clonotypes:
        if (cl, barcode) not in tcr_bc_counter_dict:
            errors.append((cl, barcode))
        else:
            bc_count = tcr_bc_counter_dict[(cl, barcode)]
            if bc_count > counts:
                counts = bc_count
                cl_keep = cl
    if (errors != [] and verbose):
        print(errors)
    return cl_keep


def save_puck(slideseq_name, puck,direc_name = None):
    if direc_name == None:
        with open("Puck_object_{}.pickle".format(slideseq_name), "wb") as handle:
            pkl.dump(puck, handle, protocol=pkl.HIGHEST_PROTOCOL)
    else:
         with open(f"{direc_name}Puck_object_{slideseq_name}.pickle", "wb") as handle:
            pkl.dump(puck, handle, protocol=pkl.HIGHEST_PROTOCOL)
    return puck


def load_puck(slideseq_name,direc_name = None):
    if direc_name == None:
        with open("Puck_object_{}.pickle".format(slideseq_name), "rb") as handle:
            puck = pkl.load(handle)
    else:
        with open(f"{direc_name}Puck_object_{slideseq_name}.pickle", "rb") as handle:
            puck = pkl.load(handle)
    return puck


def make_adata(slideseq_name, puck, save=True):
    adata = sc.AnnData(puck.s1_dge_fmtd)
    adata.var_names_make_unique()
    sc.pl.highest_expr_genes(
        adata,
        n_top=20,
    )
    sc.pp.filter_cells(adata, min_genes=50)
    sc.pp.filter_genes(adata, min_cells=3)
    adata.var['mt'] = adata.var_names.str.startswith(
        'MT-')  # annotate the group of mitochondrial genes as 'mt'
    sc.pp.calculate_qc_metrics(adata,
                               qc_vars=['mt'],
                               percent_top=None,
                               log1p=False,
                               inplace=True)
    sc.pl.violin(adata, ['n_genes_by_counts', 'total_counts', 'pct_counts_mt'],
                 jitter=0.4,
                 multi_panel=True)
    sc.pl.scatter(adata, x='total_counts', y='pct_counts_mt')
    sc.pl.scatter(adata, x='total_counts', y='n_genes_by_counts')
    #adata = adata[adata.obs.n_genes_by_counts < 2500, :]
    adata = adata[adata.obs.pct_counts_mt < 5, :]
    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)
    sc.pp.highly_variable_genes(adata,
                                min_mean=0.0125,
                                max_mean=3,
                                min_disp=0.5)
    sc.pl.highly_variable_genes(adata)
    adata.raw = adata
    adata = adata[:, adata.var.highly_variable]
    sc.pp.regress_out(adata, ['total_counts', 'pct_counts_mt'])
    sc.pp.scale(adata, max_value=10)
    sc.tl.pca(adata, svd_solver='arpack')
#     sc.pl.pca(adata, color='TRBC2')
    sc.pl.pca_variance_ratio(adata, log=True)
    sc.pp.neighbors(adata, n_neighbors=10, n_pcs=40)
    sc.tl.umap(adata)
#     sc.pl.umap(adata, color=['TRBC2'], cmap='viridis')
    sc.tl.leiden(adata)
    sc.pl.umap(adata, color=['leiden'])
    sc.tl.rank_genes_groups(adata, 'leiden', method='wilcoxon')
    sc.pl.rank_genes_groups(adata, n_genes=25, sharey=False)
    # Show top 10 ranked genes in clusters in a dataframe
    pd.DataFrame(adata.uns['rank_genes_groups']['names']).head(5)
    if save == True:
        results_file = './{}/{}.h5ad'.format(directory_slideseq, slideseq_name)
        adata.write(results_file)
    return adata


def load_adata(slideseq_name,directory_slideseq):
    with open('./{}/{}.h5ad'.format(directory_slideseq, slideseq_name),
              "rb") as handle:
        adata = sc.read_h5ad(handle)
    return adata


def pkl_dump(obj,file_name):
    file = open(f'{file_name}.pkl', 'wb')
    # dump information to that file
    pkl.dump(obj, file)
    # close the file
    file.close()  
    
def pkl_load(file_name):
    file = open(f'{file_name}.pkl','rb')
    obj = pkl.load(file)
    file.close()
    return obj