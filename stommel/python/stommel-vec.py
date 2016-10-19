import numpy
import math
import matplotlib.pyplot as plt
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
delT = ni+1  #  make ni = 1 or 2 to reduce the number of
              #     integrations to be done
delS = delT

ns = numpy.append(cartesian([0, 1], numpy.linspace(0, 1, delS)), cartesian(numpy.linspace(0, 1, delT), [0, 1]), axis=0)
x = numpy.zeros((len(ns), nstep))
y = numpy.zeros((len(ns), nstep))
tau = numpy.append(numpy.zeros(1), numpy.arange(2, nstep+1) * dtau) # the non-d time
x[:, 0] = ns[:, 0]
y[:, 0] = ns[:, 1]

# evaluate the reservoir temperature (y); note that
#   this temperatre is steady if resosc = 0. (the S61 case)

yres = yres0 + resosc * numpy.sin(tau[:-1]*math.pi)

for m in xrange(1, nstep):

    # the first part of a second order R-K method; time step forward
    #    by half a time step to make a first guess at the new times

    dr = numpy.abs(R*x[:, m-1] - y[:, m-1])        # the density anomaly

    qequil = dr/lambdA                 # the equilibrium flushing

    yh = y[:, m-1] + dtau*(yres[m-1] - y[:, m-1])/2 - dtau*y[:, m-1]*q/2

    xh = x[:, m-1] + dtau*delta*(1 - x[:, m-1])/2 - dtau*x[:, m-1]*q/2

    qh = q + dtau*qdelta*(qequil - q)/2

    # the second part; use the half time step values to make a
    #    full step

    dr = numpy.abs(R*xh - yh)

    qequil = dr / lambdA

    y[:, m] = y[:, m-1] + dtau*(yres[m-1] - yh) - dtau*qh*yh

    x[:, m] = x[:, m-1] + dtau*delta*(1 - xh) - dtau*qh*xh

    q = q + dtau*qdelta*(qequil - qh)

    # now add on a flickering temperature if you want to (or comment out)

    # tflickamp = 0.01;   %  set the amplitude here
    # tflick = tflickamp*unifrnd(-1., 1.);
    # y[m] = y[m] + tflick;

d = R*x - y

plt.figure(1)
plt.clf()
plt.subplot(2,1,1)
line1, line2 = plt.plot(tau, x[0], '-', tau, y[0], '--')
plt.legend((line1, line2), ('salinity', 'temperature'))
plt.ylabel('T, S diff, non-d')
plt.title('Experiment 1,1')

plt.subplot(2,1,2)
plt.plot(tau, d[0])
plt.xlabel('time, non-d')
plt.ylabel('density diff')

# contour the density (or flushing rate), and add the (T, S)
# trajectories on top of the contours
step = 0.1
zm = numpy.arange(0, 1+step, step)
dm = (1./lambdA)*(R*zm - zm[numpy.newaxis, :].T)

plt.figure(2)
plt.clf()
plt.axis([0,1,0,1])
step = 2
dc = numpy.arange(-10, 20+step, step)
c = plt.contour(zm, zm, dm, dc, colors='k')
plt.clabel(c)
plt.xlabel('salinity diff, non-d')
plt.ylabel('temp diff, non-d')

(m1, m2) = x.shape
m2 -= 1

#  plot the individual trajectories. 
#  color code according to which equilibrium point
#  the trajectory ends up on. this will likely have be 
#  reset if the model parameters (R, delta, lambda) are changed.
for i in xrange(len(ns)):
    if d[i, m2] >= 0:
        plt.plot(x[i], y[i], 'r--')
        plt.plot(x[i, m2], y[i, m2], '*r')
    else:
        plt.plot(x[i], y[i], 'g')
        plt.plot(x[i, m2], y[i, m2], '*g')

f = numpy.arange(-29, 30) * 0.1
lhs = lambdA * f
rhs = (R/(1 + numpy.abs(f)/delta)) - 1/(1 + numpy.abs(f))

plt.figure(3)
plt.axis([-3, 3, 0., 1.])
plt.clf()
plt.plot(f, rhs, f, lhs)
plt.xlabel('f, flow rate')
plt.ylabel('lhs(f), rhs(f)')
plt.title('roots of S61 model')
plt.grid()


plt.show()
