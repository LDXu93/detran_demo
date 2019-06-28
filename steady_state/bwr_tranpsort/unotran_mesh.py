#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import matplotlib.pyplot as plt
from pydmd import DMD
import scipy as sp
from scipy.linalg import norm
import numpy as np
import pydgm
import numpy as np
np.set_printoptions(precision=16, linewidth=132)


class solver:
    """
    Initialize Unotran setup, include geometry, quadrature, and boundary conditon. 
    """

    def __init__(self):

        # set up mesh geometry
        pydgm.control.fine_mesh_x = [6, 18, 18, 18, 18, 6, 6, 18, 18, 18, 18, 6, 6, 18, 18, 18, 18,
                                     6, 6, 18, 18, 18, 18, 6, 6, 18, 18, 18, 18, 6, 6, 18, 18, 18, 18, 6, 6, 18, 18, 18, 18, 6]
        pydgm.control.coarse_mesh_x = [0.,      1.1176,  4.3688,  7.62,    10.8712,  14.1224,  15.24,    16.3576,
                                       19.6088,  22.86, 26.1112,  29.3624,  30.48,    31.5976,  34.8488,  38.1,
                                       41.3512,  44.6024, 45.72,   46.8376,  50.0888,  53.34,    56.5912,  59.8424,
                                       60.96,  62.0776, 65.3288,  68.58,    71.8312,  75.0824,  76.2,     77.3176,
                                       80.5688,  83.82,  87.0712,  90.3224,  91.44,    92.5576,  95.8088,  99.06,
                                       102.3112, 105.5624, 106.68]
        pydgm.control.material_map = [5, 1, 2, 2, 1, 5, 5, 1, 3, 3, 1, 5, 5, 1, 2, 2,
                                      1, 5, 5, 1, 3, 3, 1, 5, 5, 1, 2, 2, 1, 5, 5, 1, 3, 3, 1, 5, 5, 1, 2, 2, 1, 5]

        pydgm.control.angle_order = 4  # set to s_4 quadrature
        pydgm.control.angle_option = pydgm.angle.gl  # gauss legendre Quadrature
        pydgm.control.boundary_east = 0.0  # set east boundary vacuum
        pydgm.control.boundary_west = 0.0  # set west boundary vacuum
        pydgm.control.spatial_dimension = 1  # set to a 1d problem
        pydgm.control.allow_fission = True  # compute the fission source
        pydgm.control.eigen_print = 1  # print information of eigen iteration
        pydgm.control.outer_print = 0  # Dont print information of eigen iteration
        pydgm.control.eigen_tolerance = 1e-8  # residual of eigen iteration
        pydgm.control.outer_tolerance = 1e-8  # residual of source iteration
        pydgm.control.max_eigen_iters = 1  # maxinum number of eigen iteration
        pydgm.control.max_outer_iters = 1  # maxinum number of outer iteration
        pydgm.control.store_psi = False  # don't store the angular flux
        pydgm.control.solver_type = 'eigen'.ljust(256)  # eigen solver
        pydgm.control.source_value = 0.0  # no external source
        pydgm.control.equation_type = 'DD'  # diamond difference
        pydgm.control.scatter_leg_order = 0  # only 0th scattering legendre order
        pydgm.control.ignore_warnings = True  # silence warning of outer iteration
        pydgm.control.xs_name = '2gXS.anlxs'.ljust(256)

        # Initialize the dependancies
        pydgm.solver.initialize_solver()
        self.init_flux = pydgm.state.phi



