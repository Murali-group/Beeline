# scSGL - a python package for fene regulatory network inference using graph signal processing based
# signed graph learning
# Copyright (C) 2021 Abdullah Karaaslanli <evdilak@gmail.com>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <https://www.gnu.org/licenses/>.

from itertools import permutations

import numpy as np
import pandas as pd

from scipy.spatial.distance import squareform

from . import unsigned
from . import signed
from ..associations import correlation, dotprod, proprho, zikendall

def _find_bs_upper_bound(k, d, density):
    for i in range(1, 100):
        w = unsigned.learn(k, d, i)
        densities_est = np.count_nonzero(w)/len(w)
        if densities_est > density:
            return i

def _binary_search(k, d, density_pos, density_neg):
    apos_done = False
    aneg_done = False

    apos_min = 0
    aneg_min = 0

    # TODO: Find better upper bound
    if density_pos > 0:
        apos_max = _find_bs_upper_bound(k, d, density_pos)
    else:
        apos_max = 0
        apos_done = True
        wpos = np.zeros((len(k, 1)))

    if density_neg > 0:
        aneg_max = _find_bs_upper_bound(k, d, density_neg)
    else:
        aneg_max = 0
        aneg_done = True
        wneg = np.zeros((len(k, 1)))

    densities_pos = np.zeros(50)
    densities_neg = np.zeros(50)
    for i in range(50):
        apos = (apos_min + apos_max)/2
        aneg = (aneg_min + aneg_max)/2

        if aneg == 0: # Learn only positive edges
            wpos = unsigned.learn(k, d, apos)
            densities_pos[i] = np.count_nonzero(wpos)/len(wpos)
        elif apos == 0: # Learn only negative edges
            wneg = unsigned.learn(-k, d, aneg)
            densities_neg[i] = np.count_nonzero(wneg)/len(wneg)
        else: 
            wpos, wneg= signed.learn(k, d, apos, aneg, lpos_init="zeros", 
                lneg_init="zeros")
            densities_pos[i] = np.count_nonzero(wpos)/len(wpos)
            densities_neg[i] = np.count_nonzero(wneg)/len(wneg)

        # print("Current densities: {:.2f}    {:.2f}".format(densities_pos[i], densities_neg[i]))

        # Check if desired density is obtained for positive part
        if not apos_done:
            if np.abs(density_pos - densities_pos[i]) < 1e-2:
                apos_done = True
            elif density_pos < densities_pos[i]:
                apos_max = apos
            elif density_pos > densities_pos[i]:
                apos_min = apos

        # Check if desired density is obtained for negative part
        if not aneg_done:
            if np.abs(density_neg - densities_neg[i]) < 1e-2:
                aneg_done = True
            elif density_neg < densities_neg[i]:
                aneg_max = aneg
            elif density_neg > densities_neg[i]:
                aneg_min = aneg

        # If desired densities are obtained, break
        if (apos_done and aneg_done):
            break

        # If binary search stuck, break
        if i>2:
            if np.abs(densities_pos[i] - densities_pos[i-1]) < 1e-3 and \
            np.abs(densities_pos[i] - densities_pos[i-2]) < 1e-3 and \
            np.abs(densities_neg[i] - densities_neg[i-1]) < 1e-3 and \
            np.abs(densities_neg[i] - densities_neg[i-2]) < 1e-3:
                break

    return wpos, -wneg   

def learn_signed_graph(X, pos_density, neg_density, assoc="dotprod", gene_names = None, 
                       verbose=False):
    # TODO: Docstring
    # TODO: Input check

    assocs = {"dotprod": dotprod.calc,
              "correlation": correlation.calc,
              "proprho": proprho.calc,
              "zikendall": zikendall.calc}

    if gene_names is None:
        gene_names = np.arange(1, X.shape[0]+1)

    # Check if there is any genes that has no expression at all
    nnzeros = np.count_nonzero(X, axis=1) != 0
    X_nnzeros = X[nnzeros, :]

    # Calculate association matrix
    K = assocs[assoc](X_nnzeros)
    k = K[np.triu_indices_from(K, k=1)]
    k /= np.max(np.abs(k))
    d = np.diag(K)/np.max(np.abs(k))

    # Learn graph with desired density
    if verbose:
        print("Estimating a graph whose positive and negative edges densities are",
              "{:.3f} and {:.3f}...".format(pos_density, neg_density))

    wpos, wneg = _binary_search(k, d, pos_density, neg_density)

    pos_density_est = np.count_nonzero(wpos)/len(wpos)
    neg_density_est = np.count_nonzero(wneg)/len(wneg)

    if verbose:
        print("Graph is found. Its positive and negative edge densities are {:.3f} and {:.3f}"\
            .format(pos_density_est, neg_density_est))

    return convert_df(gene_names[nnzeros], wpos, wneg)

def convert_df(gene_names, lpos, lneg):
    gene1 = [i for i, _ in permutations(gene_names, r=2)]
    gene2 = [j for _, j in permutations(gene_names, r=2)]
    L = squareform(np.squeeze(lpos + lneg))
    edge_weights = [L[i, j] for i, j in permutations(range(len(gene_names)), r=2) if i != j]

    grn_df = pd.DataFrame({"Gene1": gene1, "Gene2": gene2, "EdgeWeight": edge_weights})
    grn_df = grn_df[grn_df.EdgeWeight != 0]
    
    return grn_df

