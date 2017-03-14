import climate
from climate.pyom import PyOM, external, idemix
from climate.pyom import numerics
from climate.setup.acc2 import ACC2
import numpy as np
import sys
np.set_printoptions(threshold=np.nan)

yt_start = -39.0
yt_end   = 43
yu_start = -40.0
yu_end   = 42

baroTimer = climate.Timer("barotropics")
a = ACC2()
"""
 allocate everything
"""
a.set_parameter()
a.allocate()

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
