from sympy import Function, Eq
from sympy.solvers import solve
from sympy.abc import x, y, d
import matplotlib.pyplot as plt
import numpy
import math
numpy.set_printoptions(threshold=numpy.nan)

def cartesian(x, y):
    return numpy.transpose([numpy.tile(x, len(y)), numpy.repeat(y, len(x))])

#  Stommel's 1961 model of convection in coupled boxes.
#
#  Compute the temperature and salinity difference, here represented
#  by y and x, between two well-stirred boxes that are both
#  in contact with a reservoir that has (y, x) = (1, 1).  Both
#  boxes conduct y and x at rates 1 and delta (delta < 1), and the
#  density difference between the boxes is given by d = -y + R*x, where
#  we are assuming R > 1, typically.  There can be advection (or
#  flushing) between the boxes at a rate d*lambda, where the
#  flushing does not depend upon the sign of the density anomaly.
#
#  This code also has 1) flushing with inertia (not especially
#  interesting), 2) a flickering or random temperature
#  perturbation (slightly interesting), and 3) an oscillating
#  reservoir temperature (hoping for but not yet finding chaos;
#  may be present in parameter regimes not checked).
#
#  This code has not been fully tested, but seems to reproduce
#  S61's two cases fairly well when all three of the 'new' things
#  are turned off, of course.
#
#  Written by Jim Price, April 28, 1999.  Public domain.

nn = 0

R      = 2.     #  abs of the ratio of the expansion coefficients, x/y
delta  = 1./6.  #  conduction rate of salinity wrt temperature
lambdA = 0.2    #  inverse non-d flushing rate; = inf for no flushing
q      = 0.     #  initial flushing rate (0 to 1) 'new'
qdelta = 100.   #  time constant (inertia) for flushing; 'new'
                #    set = 1/dtau for equilibrium flushing as in S61
                #    set = 0.2 for slowly responding flushing
yres   = 1.     #  steady reservoir y, = 1 for S61 case  'new'
resosc = 0.     #  amplitude of reservoir y oscillation  'new'
dtau   = 0.01   #  the time step of non-d time; 0.01 seems OK
nstep  = 1500   #  number of time steps; 1500 is usually
                #    enough to insure convergence to a steady state

#  This model version is set up to do integration over a range
#  of initial T,S or y,x from 0 to 1.  The increment of T and S 
#  are set by delT and delS = 1/ni,  where ni is the number of
#  integrations (n1 = 1 to 20 is reasonable).
yres0 = yres

ni = 6.
delT = 1./ni  #  make ni = 1 or 2 to reduce the number of
              #     integrations to be done
delS = delT

tau = numpy.append(numpy.zeros(1), numpy.arange(2, nstep+1) * dtau) # the non-d time


xVal = numpy.zeros(nstep)
yVal = numpy.zeros(nstep)
tau = numpy.zeros(nstep)
for n1 in numpy.linspace(0, 1, ni+1):
    for n2 in numpy.linspace(0, 1, ni+1):
        if n1 == 0 or n1 == 1 or n2 == 0 or n2 == 1:
            xVal[0] = n1  # set the initial temperature
            yVal[0] = n2  # set the initial salinity

            for t in xrange(1, nstep):
                xVal[t] = 1 + (xVal[0] - 1) * math.exp(delta * tau[t])
                yVal[t] = 1 + (yVal[0] - 1) * math.exp(tau[t])

            d = R*xVal - yVal

            if nn == 0:
                nn = 1

                plt.figure(1)
                plt.clf()
                plt.subplot(2,1,1)
                line1, line2 = plt.plot(tau, xVal, '-', tau, yVal, '--')
                plt.legend((line1, line2), ('salinity', 'temperature'))
                plt.ylabel('T, S diff, non-d')
                plt.title('Experiment 1,1')

                plt.subplot(2,1,2)
                plt.plot(tau, d)
                plt.xlabel('time, non-d')
                plt.ylabel('density diff')

plt.show()
