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

import numpy as np

from . import utils

def calc(counts):
    # TODO: Docstring

    # Input validation
    if np.isnan(counts).any():
        raise ValueError("Argument 'counts' cannot have NANs.")

    if np.any(counts < 0):
        raise ValueError("Argument 'counts' cannot have negative numbers.")

    # Replace zeros in counts with minimum non-zero entry if alpha is not provided
    lr = counts.copy()
    lr[counts == 0] = np.min(counts[counts != 0])

    # Log Transformation
    lr = np.log(lr)
    ref = np.mean(lr, axis=0) # mean over samples
    lr -= ref
    
    # Calculate the proportionality metric rho
    covlr = np.cov(lr) # Covariance matrix with rows being variables
    vars = np.diag(covlr) # Variance of each row
    rho = 2*covlr/(vars + vars[..., None])
    rho[np.diag_indices_from(rho)] = 1

    return rho

def permutations(counts, k, tau_neg, tau_pos):
    # TODO: Docstring

    return utils._permutations(counts, calc, k, tau_neg, tau_pos)

def associations(counts, k):
    # TODO: Docstring
    
    return utils._associations(counts, calc, k)