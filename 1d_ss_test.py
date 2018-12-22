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

def get_mesh() :
    """ | mod | f1 | f1 | f1 | f2 | f1 | mod | """
        
    fm = [10]
    cm = [0.0, 20.0]
    mt = [0]
    
    mesh = Mesh1D.Create(fm, cm, mt)

    return mesh
    
def get_material() :

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
    return mat

def get_db() :

    """ Provides the base nodal db (or for use in the reference)
    """
    db = InputDB.Create()
    db.put_dbl("inner_tolerance",                1.0e-16)
    db.put_dbl("outer_tolerance",                1.0e-16)
    db.put_int("number_groups",                  number_groups)
    db.put_int("dimension",                      1)
    db.put_int("quad_number_polar_octant",       4)
    db.put_str("equation",                       "dd")
    db.put_str("bc_west",                        "vacuum")
    db.put_str("bc_east",                        "vacuum")
    db.put_int("inner_max_iters",                50000)
    db.put_int("outer_max_iters",                50000)
    return db
    
    

def run() :
    """ Run the reference eigenvalue problem using Detran.  Also, produce
        the reference modes for producing an orthogonal basis.
    """
    db   = get_db()

    mat  = get_material()

    mesh = get_mesh()


    solver = Fixed1D(db, mat, mesh, False)
    solver.setup()
    S = ConstantSource.Create(number_groups, mesh, 1)
    solver.set_source(S)
    solver.set_solver()

    solver.solve()


    ng = solver.material().number_groups()
    nc = solver.mesh().number_cells()

    phi = np.zeros((nc, ng))
    for i in range(nc) :
        for g in range(ng) :
            phi[i, g] = solver.state().phi(g)[i]

    #print 'phi = '
    print phi
    np.savetxt('phi', phi.flatten())
    

    

if __name__ == "__main__" :
    Manager.initialize(sys.argv)
    run()
