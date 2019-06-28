# pyexamples/slab_reactor/slab_reactor_materials.py
#
# This implements the materials for the test cases published 
# in Scott Mosher's Ph.D. thesis, "A Variational Coarse Mesh 
# Transport Method". All data and reference values are from 
# Appendix A of that work.

from detran import *

def get_materials() :
    # Two-group data from 1-d coarse mesh benchmarks (Mosher, Ilas etc.)
    
    # Create the Materials object.
    mat = Material.Create(5, 2, "slabreactor");

    # ---------------------------
    # Material 0: UO2-1           
    # ---------------------------

    # Total
    mat.set_sigma_t(0, 0, 5.373960000000e-01);       # (obj, matid, g, value);
    mat.set_sigma_t(0, 1, 1.044040000000e+00);       

    # Fission 
    mat.set_sigma_f(0, 0, 6.445320000000e-03);         # Note, default is zero
    mat.set_sigma_f(0, 1, 1.385100000000e-01);   
    mat.set_chi(0, 0, 1.0); 
    mat.set_chi(0, 1, 0.0);        

    # Scattering
    mat.set_sigma_s(0, 0, 0, 5.038510000000e-01);    # 1 <- 1
    mat.set_sigma_s(0, 0, 1, 2.672070000000e-04);    # 1 <- 2
    mat.set_sigma_s(0, 1, 0, 1.911150000000e-02);    # 2 <- 1
    mat.set_sigma_s(0, 1, 1, 9.583650000000e-01);    # 2 <- 2

    # ---------------------------
    # Material 1: UO2-2          
    # ---------------------------

    # Total
    mat.set_sigma_t(1, 0, 5.377550000000e-01);       # (obj, matid, g, value);
    mat.set_sigma_t(1, 1, 1.076520000000e+00);       

    # Fission 
    mat.set_sigma_f(1, 0, 8.788020000000e-03);
    mat.set_sigma_f(1, 1, 2.208680000000e-01);   
    mat.set_chi(1, 0, 1.0); 
    mat.set_chi(1, 1, 0.0);        

    # Scattering
    mat.set_sigma_s(1, 0, 0, 5.031960000000e-01);    # 1 <- 1
    mat.set_sigma_s(1, 0, 1, 3.664320000000e-04);    # 1 <- 2
    mat.set_sigma_s(1, 1, 0, 1.882600000000e-02 );    # 2 <- 1
    mat.set_sigma_s(1, 1, 1, 9.515900000000e-01);    # 2 <- 2

    # ---------------------------
    # Material 2: UO2-Gd   
    # ---------------------------

    # Total
    mat.set_sigma_t(2, 0, 5.405050000000e-01);       # (obj, matid, g, value);
    mat.set_sigma_t(2, 1, 2.093440000000e+00);       

    # Fission 
    mat.set_sigma_f(2, 0, 8.789080000000e-03);
    mat.set_sigma_f(2, 1, 1.153600000000e-01);   
    mat.set_chi(2, 0, 1.0); 
    mat.set_chi(2, 1, 0.0);        

    # Scattering
    mat.set_sigma_s(2, 0, 0, 5.025980000000e-01);    # 1 <- 1
    mat.set_sigma_s(2, 0, 1, 3.128760000000e-03);    # 1 <- 2
    mat.set_sigma_s(2, 1, 0, 1.777900000000e-02 );    # 2 <- 1
    mat.set_sigma_s(2, 1, 1, 8.288520000000e-01);    # 2 <- 2

    # ---------------------------
    # Material 3: MOX        
    # ---------------------------

    # Total
    mat.set_sigma_t(3, 0, 5.446680000000e-01);       # (obj, matid, g, value);
    mat.set_sigma_t(3, 1, 1.470290000000e+00);       

    # Fission 
    mat.set_sigma_f(3, 0, 1.990060000000e-02);
    mat.set_sigma_f(3, 1, 9.624560000000e-01);   
    mat.set_chi(3, 0, 1.0); 
    mat.set_chi(3, 1, 0.0);        

    # Scattering
    mat.set_sigma_s(3, 0, 0, 4.990620000000e-01);    # 1 <- 1
    mat.set_sigma_s(3, 0, 1, 1.282830000000e-03);    # 1 <- 2
    mat.set_sigma_s(3, 1, 0, 1.437630000000e-02 );    # 2 <- 1
    mat.set_sigma_s(3, 1, 1, 9.185740000000e-01 );    # 2 <- 2

    # ---------------------------
    # Material 4: Water          
    # ---------------------------
    
    # Total
    mat.set_sigma_t(4, 0, 5.953610000000e-01);       # (obj, matid, g, value);
    mat.set_sigma_t(4, 1, 1.705580000000e+00);       

    # Fission 
    mat.set_sigma_f(4, 0, 0.0);
    mat.set_sigma_f(4, 1, 0.00);   
    mat.set_chi(4, 0, 1.0); 
    mat.set_chi(4, 1, 0.0);        

    # Scattering
    mat.set_sigma_s(4, 0, 0, 5.674490000000e-01);	   # 1 <- 1
    mat.set_sigma_s(4, 0, 1, 4.503090000000e-04);    # 1 <- 2
    mat.set_sigma_s(4, 1, 0, 2.739510000000e-02 );    # 2 <- 1
    mat.set_sigma_s(4, 1, 1, 1.695970000000e+00);    # 2 <- 2

    # ---------------------------
    # FINALIZE     
    # ---------------------------

    mat.finalize();

    return mat

