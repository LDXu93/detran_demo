from detran import *
import cPickle as pickle
import matplotlib.pyplot as plt
import math
import time
import sys
import os
import numpy as np

fluxes = []
precursors = []
times = []
number_groups = 1
number_material = 1
order_legendre = 1


mat = Material.Create(number_material, number_groups)

# set lambda for each precurser
mat.set_lambda(0, 0.0127)
mat.set_lambda(1, 0.317)
mat.set_lambda(2, 0.115)
mat.set_lambda(3, 0.311)
mat.set_lambda(4, 1.4)
mat.set_lambda(5, 3.87)

# set neutron velocity for each energy group
mat.set_velocity(0, 4.48e9)
m = 0
s_t = np.array([[2.134]])
s_s = np.array([[[0.334]]])
# Set the cross sections in detran
for m in range(number_material):
    for g in range(number_groups):
        mat.set_sigma_t(m, g, s_t[m, g])
        for gp in range(number_groups):
            mat.set_sigma_s(m, g, gp, s_s[m, g, gp])
            
mat.set_diff_coef(m, vec_dbl([1.4389]))
mat.set_sigma_f(m, vec_dbl([0.0074527]))
mat.set_chi(m, 0, 1.0)
# set beta for each precurser
mat.set_beta(m, 0, 2.470e-04)
mat.set_beta(m, 1, 1.3845e-03)
mat.set_beta(m, 2, 1.222e-03)
mat.set_beta(m, 3, 2.6455e-03)
mat.set_beta(m, 4, 8.320e-04)
mat.set_beta(m, 5, 1.690e-04)
mat.compute_diff_coef()
mat.compute_sigma_a()
