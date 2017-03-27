from climate.setup.acc2 import ACC2
from climate import Timer
from climate.pyom.core import momentum, numerics, external, idemix, eke, isoneutral
import sys

timer = Timer("momentum")
if len(sys.argv) == 2:
    with timer:
        a = ACC2(fortran=sys.argv[1])
        a.setup()
        for i in xrange(100):
            a.fortran.momentum()
else:
    with timer:
        a = ACC2()
        a.setup()
        #"""
        #allocate everything
        #"""
        #a.set_parameter()
        #a.allocate()

        #"""
        #Grid
        #"""
        #a.set_grid()
        #numerics.calc_grid(a)

        #"""
        #Coriolis
        #"""
        #a.set_coriolis()
        #numerics.calc_beta(a)

        #"""
        #topography
        #"""
        #a.set_topography()
        #numerics.calc_topo(a)
        #idemix.calc_spectral_topo(a)

        #"""
        #initial condition and forcing
        #"""
        #a.set_initial_conditions()
        #numerics.calc_initial_conditions(a)

        #a.set_forcing()
        #if a.enable_streamfunction:
        #    external.streamfunction_init(a)

        #"""
        #initialize EKE module
        #"""
        #eke.init_eke(a)

        #"""
        #initialize isoneutral module
        #"""
        #isoneutral.check_isoneutral_slope_crit(a)
        for i in xrange(100):
            momentum.momentum(a)
timer.printTime()
