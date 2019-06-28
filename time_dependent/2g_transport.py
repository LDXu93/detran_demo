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


def monitor(ts, step, t, dt, it, conv) :
    # Extract the state vector and precursor concentrations
    # from the time stepper and save the first element of
    # the flux and precursors along with the time.  Since this
    # is a purely reflective, homogeneous problem, the flux and
    # precursors are constant.  Hence, this is really just a point
    # kinetics problem, i.e., kinetics with no space/angle dependence.
    c = np.asarray(ts.precursor().C(0)[0])
    f = ts.state().phi(0)[0]
    fluxes.append(f)
    precursors.append(c)
    times.append(t)
    if step == 0 :
        print "       step     time          phi             C"
        print "   -----------------------------------------------------"
    print "     %4i   %12.6e  %12.6e   %12.6e " % (step, t, f, c)

def get_mesh() :
    """ | mod | f1 | f1 | f1 | f2 | f1 | mod | """
        
    fm = [10, 10, 10, 10, 10, 10, 10]
    cm = [0.0, 20.0, 40.0, 60.0, 80.0, 100.0, 120.0, 140.0]
    mt = [0, 0, 0, 0, 0, 0, 0]
    
    mesh = Mesh1D.Create(fm, cm, mt)

    return mesh


def get_db() :
    
    basis_inp = InputDB.Create("basis_data")
    """ Provides the base nodal db (or for use in the reference)
    """
    db = InputDB.Create()
    db.put_int("number_groups",                  2)
    db.put_int("dimension",                      1)
    db.put_str("equation",                       "dd")
    db.put_str("quad_type",                      "gl")
    db.put_int("quad_number_polar_octant",       2)
    db.put_str("inner_solver",                   "SI")
    db.put_int("inner_print_level",              0)
    db.put_int("inner_max_iters",                50000)
    db.put_dbl("inner_tolerance",                1.0e-16)
    db.put_dbl("eigen_tolerance",                1.0e-16)
    db.put_int("eigen_max_iters",                1000)
    db.put_str("outer_solver",                   "GS")
    db.put_int("outer_max_iters",                10000)
    db.put_int("outer_print_level",              0)
    db.put_dbl("outer_tolerance",                1.0e-16)
    db.put_int("outer_krylov_group_cutoff",      0)
    db.put_str("basis_p_type",                   "jacobi")
    db.put_str("bc_west",                        "vacuum")
    db.put_str("bc_east",                        "vacuum")
    db.put_int("erme_angular_expansion",         1)
    db.put_int('adjoint',                        0)
    db.put_int("store_angular_flux",             1)
    db.put_int("ts_scheme",                  Time1D.BDF1)
    db.put_dbl("ts_step_size",               0.1)
    db.put_dbl("ts_final_time",              3.0)

    return db


def get_steady_state() :
    """ Run the reference eigenvalue problem using Detran.  Also, produce
        the reference modes for producing an orthogonal basis.
    """
    inp = get_db()
    mat = TDMaterial()
    mesh = get_mesh()
    mat.update(0.0, 0.0, 1, False)
    
    solver_type = 'fixed'
    # Flag controls the presence of fission or not
    if solver_type == 'fixed':
        solver = Fixed1D(inp, downcast(mat.material), mesh, True)
        solver.setup()
        
        S = ConstantSource.Create(2, mesh, 1)
        solver.set_source(S)
        solver.set_solver()
    elif solver_type == 'eigen':
        solver = Eigen1D(db, mat, mesh)
    
    solver.solve()
    print '#####'
    state = solver.state()
    print(dir(state))


    ng = solver.material().number_groups()
    nc = solver.mesh().number_cells()
    
    phi = np.zeros((nc, ng))
    precursor = np.zeros((6, ng))
    for i in range(nc) :
        for g in range(ng) :
            phi[i, g] = solver.state().phi(g)[i]
    np.savetxt('flux_ss.txt',phi)
    for i in range(6) :
        for g in range(ng) :
            precursor[i, g] = solver.state().precursor(g)[i]

            
    
    np.savetxt('precursor_ss.txt',precursor_ss)
    #--------------------------------------------------------------------------#
    # PLOTS AND DATA
    #--------------------------------------------------------------------------#

    data = {}
    data['phi_full'] = np.mean(phi, axis=1)
    data['keff'] = keff
    pickle.dump(data, open('ss_reference_data.p', 'wb'))

    return state


def run() :
    inp = get_db()
    mat = TDMaterial()
    mesh = get_mesh()
    # Set up the time stepper.  Note that mat.material is passed and
    # not just mat.  That's because the material class defined above
    # based on the UserMaterial is more of a wrapper class around
    # the actual time-dependent material.  It's sort of an annoying
    # structure that will be simplified with the new smart pointer
    # design under development.
    ts = Time1D(inp, mat.material, mesh, True)
    ts.set_monitor(monitor)
    # Compute the steady-state solution.  By definition, we should
    # be able to step in time without any change in the flux if
    # the materials are not changing.  We see that is the case
    # here for the first 0.1 seconds before the step-insertion
    # of reactivity.
    print 'start'
    ic = get_steady_state()
    print'test'
    # Make sure the material is updated with the eigenvalue.  Basically,
    # the fission cross section is divided by the eigenvalue.
    mat.material.set_eigenvalue(ic.eigenvalue())
    mat.material.update(0, 0, 1, False)

    # Solve the transient problem.
    ts.solve(ic)


class TDMaterial(UserMaterial) :
    """ To solve transient problems in Python, there is a UserMaterial
        class that must be used.  The update_material function must be
        present, and it assigns all the material properties at any 
        point in time.  Because the class constructor is arbitrary,
        the user can define it and assign to it the TimeStepper 
        (or state, or whatever) so that materials can be functions
        of the flux, etc.
    """

    def __init__(self) :
        # nine materials, seven energy groups, and six precursor groups
        super(TDMaterial, self).__init__(3, 2, 6, "biblis")
        self.update_material()

    def update_material(self) :

        mat = self.material

        # set lambda for each precurser
        mat.set_lambda(0, 0.0127)
        mat.set_lambda(1, 0.317)
        mat.set_lambda(2, 0.115)
        mat.set_lambda(3, 0.311)
        mat.set_lambda(4, 1.4)
        mat.set_lambda(5, 3.87)

        # set neutron velocity for each energy group
        mat.set_velocity(0, 4.48e9)
        mat.set_velocity(1, 8.33e5)

        # --------------------------------------------
        # Material 0: constant fuel
        # --------------------------------------------
        m = 0
        mat.set_sigma_t(m, vec_dbl([0.0274640, 0.0914080]))
        mat.set_sigma_s(m, 1, 0, 0.017101)
        mat.set_diff_coef(m, vec_dbl([1.4389, 0.3638]))
        mat.set_sigma_f(m, vec_dbl([0.0074527, 0.1323600]))
        mat.set_chi(m, 0, 1.0)
        # set beta for each precurser
        mat.set_beta(m, 0, 2.470e-04)
        mat.set_beta(m, 1, 1.3845e-03)
        mat.set_beta(m, 2, 1.222e-03)
        mat.set_beta(m, 3, 2.6455e-03)
        mat.set_beta(m, 4, 8.320e-04)
        mat.set_beta(m, 5, 1.690e-04)

        # --------------------------------------------
        # Material 1: varying fuel
        # --------------------------------------------
        m = 1
        mat.set_sigma_t(m, vec_dbl([0.0274640, 0.0914080]))
        if self.material.time() < 2:
            mat.set_sigma_s(m, 0, 1, 0.017101)
        else:
            f = 2 - self.material.time() / 2.0
            mat.set_sigma_s(m, 0, 1, 0.0914080 * (1 - f) + 0.017101 * f)
        mat.set_diff_coef(m, vec_dbl([1.4389, 0.3638]))
        mat.set_sigma_f(m, vec_dbl([0.0074527, 0.1323600]))
        mat.set_chi(m, 0, 1.0)
        # set beta for each precurser
        mat.set_beta(m, 0, 2.470e-04)
        mat.set_beta(m, 1, 1.3845e-03)
        mat.set_beta(m, 2, 1.222e-03)
        mat.set_beta(m, 3, 2.6455e-03)
        mat.set_beta(m, 4, 8.320e-04)
        mat.set_beta(m, 5, 1.690e-04)

        # --------------------------------------------
        # Material 2: water
        # --------------------------------------------
        m = 2
        mat.set_sigma_t(m, vec_dbl([0.0257622, 0.0715960]))
        mat.set_sigma_s(m, 1, 0, 0.023106)
        mat.set_diff_coef(m, vec_dbl([1.3200, 0.2772]))

        mat.compute_diff_coef()
        mat.compute_sigma_a()
        print mat.sigma_a(1), mat.sigma_s(1)

def plot():
    plt.semilogy(times, fluxes)
    plt.show()


if __name__ == "__main__" :
    Manager.initialize(sys.argv)
    run()
    data = np.column_stack((times, fluxes, precursors))
    np.savetxt('data', data, header='{:<23s} {:<25s} {:<25s}'.format('times', 'fluxes', 'precursors'))
    plot()
