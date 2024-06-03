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
import scipy as sp

from ..utils import rowsum_matrix

def _project_hyperplane(v, n):
    """ Project v onto the hyperplane defined by np.sum(v) = -n
    """
    return v - (n + np.sum(v))/(len(v))

def _qp_admm(b, inv_mat_part1, inv_mat_part2, rho=1, max_iter=1000):
    m = len(b) # number of node pairs
    n = (1 + np.sqrt(8*m + 1))//2 # number of nodes

    # Initialization
    w = np.zeros((m, 1)) # slack variable
    lambda_ = np.zeros((m, 1)) # Lagrange multiplier

    for iter in range(max_iter):
        # Update l
        l_temp = b + rho*w + lambda_
        l = inv_mat_part1@l_temp + inv_mat_part2(l_temp)
        l = _project_hyperplane(l, n)

        # Update slack variable
        w = l - lambda_/rho
        w[w>0] = 0
        
        # Update Lagrange multiplier
        lambda_ += rho*(w - l)

        residual = np.linalg.norm(w-l)
        if residual < 1e-4:
            break

    w[w>-1e-4] = 0 # Remove very small edges

    # returns adjacency matrix
    return np.abs(w)

def learn(k, d, alpha):
    # TODO: Docstring

    n = len(d) # number of nodes
    m = len(k) # number of node pairs

    S = rowsum_matrix.rowsum_matrix(n)

    # ADMM Parameter
    rho = 1 

    # Inverse matrix for li subproblems
    a = 4*alpha + rho
    b = 2*alpha
    c1 = 1/a
    c2 = b/(a*(a+n*b-2*b))
    c3 = 4*b**2/(a*(a+2*n*b-2*b)*(a+n*b-2*b))
    inv_mat_part1 = c1*sp.eye(m) - c2*S.T@S
    inv_mat_part2 = lambda x : c3*np.sum(x)*np.ones((m, 1))

    b = S.T@d - 2*k
    if np.ndim(b) == 1:
        b = b[..., None]
    
    return _qp_admm(b, inv_mat_part1, inv_mat_part2, rho)