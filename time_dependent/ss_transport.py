from detran import *

np.set_printoptions(precision=16)

number_groups = 1

def get_material() :
    """ Return the three material, 238 group library generated from SCALE 6.1.
        Materials 0 and 2 are UO2 and MOX while 1 is moderator.
    """
    # three materials, two energy groups
    mat = Material.Create(3, number_groups)

    # Define the cross sections

    s_t = np.array([[2.134],
                    [10.],
                    [10.]])
      
    s_s = np.array([[[0.334]],
                     [[0.334]],
                     [[0.334]]])
    # Set the cross sections in detran
    for m in range(3):
        for g in range(number_groups):
            mat.set_sigma_t(m, g, s_t[m, g])
            for gp in range(number_groups):
                mat.set_sigma_s(m, g, gp, s_s[m, g, gp])

    # Recompute the diffusion coefficient and absorption XS
    mat.compute_diff_coef()
    mat.compute_sigma_a()

    # Finalize the material
    mat.finalize()

    return mat

def get_mesh() :
    """ Get a fuel, moderator, or core mesh
    """
    
    fm = [36]
    cm = [0.0, 1.3]
    mt = [0]
    mesh = Mesh1D.Create(fm, cm, mt)

    return mesh

def get_db() :

    """ Provides the base nodal db (or for use in the reference)
    """
    db = InputDB.Create()
    db.put_dbl("inner_tolerance",                1.0e-16)
    db.put_dbl("outer_tolerance",                1.0e-16)
    db.put_int("number_groups",                  number_groups)
    db.put_int("dimension",                      1)
    db.put_str("equation",                       "dd")
    db.put_str("bc_west",                        "vacuum")
    db.put_str("bc_east",                        "vacuum")
    db.put_int("inner_max_iters",                50000)
    db.put_int("outer_max_iters",                50000)
    return db

def run_reference() :
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

if __name__ == '__main__':
    Manager.initialize(sys.argv)
    run_reference()
