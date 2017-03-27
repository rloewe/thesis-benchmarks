import climate
from climate.pyom import PyOM
from climate.pyom.core import numerics, external, idemix
from climate.setup.acc2 import ACC2
import numpy as np
import sys
np.set_printoptions(threshold=np.nan)

yt_start = -39.0
yt_end   = 43
yu_start = -40.0
yu_end   = 42

baroTimer = climate.Timer("barotropics")
if len(sys.argv) == 2:
    b = ACC2(fortran=sys.argv[1])
    b.fortran.my_mpi_init(0)
    b.set_parameter()
    b.set_legacy_parameter()
    b.fortran.pe_decomposition()
    b.fortran.allocate_main_module()
    b.fortran.allocate_isoneutral_module()
    b.fortran.allocate_tke_module()
    b.fortran.allocate_eke_module()
    b.fortran.allocate_idemix_module()
    b.set_grid()
    b.fortran.calc_grid()
    b.set_coriolis()
    b.fortran.calc_beta()
    b.set_topography()
    b.fortran.calc_topo()
    b.set_initial_conditions()
    b.fortran.calc_initial_conditions()
    b.fortran.streamfunction_init()
    with baroTimer:
        for i in xrange(1000):
            b.fortran.solve_streamfunction()
else:
    a = ACC2()
    """
     allocate everything
    """
    a.set_parameter()
    a._allocate()

    """
      Grid
    """
    a.set_grid()
    numerics.calc_grid(a)

    """
     Coriolis
    """
    a.set_coriolis()
    numerics.calc_beta(a)

    """
     topography
    """
    a.set_topography()

    numerics.calc_topo(a)
    idemix.calc_spectral_topo(a)

    a.set_initial_conditions()
    numerics.calc_initial_conditions(a)

    a.set_forcing()

    external.streamfunction_init(a)
    with baroTimer:
        for i in xrange(1000):
            external.solve_streamfunction(a)
baroTimer.printTime()
