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

from itertools import combinations

import numba
from numba.np.ufunc import parallel
import numpy as np
from numpy.core.numeric import count_nonzero

from scipy import stats

from . import utils

from rpy2 import robjects
from rpy2.robjects.packages import importr
from rpy2.robjects import numpy2ri
numpy2ri.activate()
# dismay = importr("dismay")
pcapp = importr("pcaPP")

@numba.njit
def _calc_kendall(x, y):
    P = 0
    Q = 0
    T = 0
    U = 0
    for i in range(len(x)):
        for j in range(i+1, len(x)):
            P += int(((x[i] > x[j]) and (y[i] > y[j])) or ((x[i] < x[j]) and (y[i] < y[j])))
            Q += int(((x[i] > x[j]) and (y[i] < y[j])) or ((x[i] < x[j]) and (y[i] > y[j])))
            T += int((x[i] == x[j]) and (y[i] != y[j]))
            U += int((x[i] != x[j]) and (y[i] == y[j]))

    denominator = np.sqrt((P+Q+T)*(P+Q+U)) 
    return (P-Q)/denominator if denominator > 0 else 0

def _calc_nonzero_kendall(x, y):
    n_samples = len(x)
    if n_samples == 1:
        return 0
    elif len(np.unique(x)) == 1 or len(np.unique(y)) == 1:
        return 0
    else:
        tau, _ = stats.kendalltau(x, y)
        return 0 if np.isnan(tau) else tau

@numba.njit
def _calc_nonzero_kendall_mat(counts):
    n_vars = counts.shape[0] 
    results = np.eye(n_vars)
    for i in range(n_vars):
        for j in range(i+1, n_vars):
            x = counts[i, :]
            y = counts[j, :]
            nz = (x != 0) & (y != 0)
            results[i, j] = _calc_kendall(x[nz], y[nz])
            results[j, i] = results[i, j]
    return results

@numba.njit
def _calc_p(counts):
    n_vars = counts.shape[0] 
    results = np.zeros((n_vars, n_vars))
    for i in range(n_vars):
        for j in range(n_vars):
            if i==j:
                continue
            x = counts[i, :]
            y = counts[j, :]

            y = y[x!=0]
            x = x[x!=0]

            x10 = x[y==0]

            x11 = np.expand_dims(x[y!=0], axis=1)

            if len(x11) > 0 and len(x10) > 0:
                results[i, j] = np.count_nonzero(x10 > x11)/(len(x11)*len(x10))
            elif len(x11) == 0:
                results[i, j] = 1
            elif len(x10) == 0:
                results[i, j] = 0

    return results

def _nonzero_kendall(X):
    n_samples = X.shape[0]
    result = np.eye(n_samples)
    for i in range(n_samples):
        x = X[i, :]
        sorting_indices = np.argsort(x)
        x = x[sorting_indices]
        for j in range(i+1, n_samples):
            y = X[j, :]
            y = y[sorting_indices]
            nnzeros = (x != 0) & (y != 0)

            result[i, j] = pcapp.cor_fk(x[nnzeros], y[nnzeros])
            result[j, i] = result[i, j]

    return result

def calc(counts):
    # tau11 = dismay.cor_fk_nz(counts.T)
    tau11 = _nonzero_kendall(counts)

    n_samples = counts.shape[1]
    
    nz = counts != 0
    p11 = (nz.astype(np.int)@(nz.T).astype(np.int))/n_samples
    p00 = (np.logical_not(nz).astype(np.int)@np.logical_not(nz.T).astype(np.int))/n_samples
    p01 = (np.logical_not(nz).astype(np.int)@(nz.T).astype(np.int))/n_samples
    p10 = (nz.astype(np.int)@np.logical_not(nz.T).astype(np.int))/n_samples

    p1 = _calc_p(counts)

    if np.any(np.isnan(tau11)):
        tau11[np.isnan(tau11)] = 0

    return p11**2 * tau11 + 2*(p00*p11 - p01*p10) + 2*p11*(p10*(1 - 2*p1) + p01*(1 - 2*p1.T))


def permutations(counts, k, tau_neg, tau_pos):
    # TODO: Docstring

    return utils._permutations(counts, calc, k, tau_neg, tau_pos)

def associations(counts, k):
    # TODO: Docstring
    
    return utils._associations(counts, calc, k)