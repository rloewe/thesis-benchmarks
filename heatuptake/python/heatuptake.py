#import bohrium as np
import numpy as np
import math
import time

# Constants
nz=100
nt=3
nout=1001
nmax=100000001
nf=140

#pi = 3.14159265358979
pi = math.pi

# run with GCM forcing data?
nmodel = False

# Southern ocean switch (0 for no Drake Passage, 1 for Drake passage)
socean = 1

# eddy diffusivities (placed here to scale dt with kappav)
#       kappav0: surface diapycnal diffusivity
#       Bryan-Lewis:
#         nb1: amplification of mixing in abyss
#         c0, c1, c2: Bryan-Lewis parameters
#         kappav = c0 + c1 atan[ c2 ( d - c3 ) ]
#       Simmons:
#         set nsimmons=True to use (ignores other values)
#         kappas: upper ocean diffusivity in Simmons (if used)
#       fractional change in eddy transport (0.1 = +10%)
kappav0  = 1.e-5
nbl      = 10.
c0       = kappav0*(nbl+1.)/2.
c1       = kappav0*(nbl-1.)/pi
c2       = 4.5e-3
c3       = 2500.
nsimmons = False
kappas   = 3.e-5
kappagm0 = 1.e3*socean
deddy    = 0.

# if using GCM forcing data, read in now
if nmodel:
    a = open("qekin.dat", "r")
    b = open("qnin.dat", "r")
    c = open("sstin.dat", "r")
    qekin = []
    qnin = []
    sstin = []
    for i in xrange(nf):
        #TODO: fix this shit
        qekin.append()
        qnin.append()
        sstin.append()
    a.close()
    b.close()
    c.close()

# Southern Ocean Ekman transport, Ekman depth, Drake Passage depth,
#        fractional change in Ekman transport (0.1 = +10%)
qek0 = 30.e6*socean

# Translated from the original code. What have they been smoking?
if nmodel:
#NOTE: What the fuck?
#    qek0 = 0
#    for i in xrange(10):
#        qek0 += qekin[i]/10
#    qek0 = qekin[0]
#    qek0 = 33e6
#    qek0 = 33.744e6
    qek0 = 36.82e6

depth  = 5.e3
ddrake = 0.8*depth
dwind  = 0.

# NADW formation, fractional reduction in NADW formation (0.5 = -50%),
#        temperature range of NADW formation, NADW warming
qn0=20.e6

if nmodel:
#    qn0 = 0
#    for i in xrange(1, 11):
#        qn0 += qnin[i]/10
#    qn0 = qnin[1]
#    qn0 = 15e6
#    qn0 = 14.4255e6
    qn0 = 14.77e6


dqn     = 0.
# NB: 1.5 degrees gets added to the final figures (in matlab script)
tnadw10 = 4.5
tnadw20 = 0.5
if nmodel:
    tnadw10 = 6.
    tnadw20 = 2.
dtnadw = 0.

# set ndeboer = True for dynamic NADW formation following de Boer (2011)
ndeboer = False

# surface temperature, surface interface number (fractional)
#     (theta0 hard-wired to keep levels the same across each scenario;
#      small number added to ensure data output times do not coincide 
#      with a new layer appearing at the surface => noisy diagnostics!)
#   NB: 1.5 degrees gets added to the final figures in non-GCM runs
#      (in matlab script)
thetas0 = 19.5
# NB: do not modify the next line; change dthetas at end of this block
dthetas = 4.
theta0 = thetas0 + dthetas + 1.117e-4
thetas = thetas0
dthetas = 0.

if nmodel:
#NOTE: Still wtf?
    """
    thethas0 = 0
    for i in xrange(10)
        thetas0 += sstin[i] / 10
    """
    thetas0 = sstin[0]
#    thetas0 = 21.5
    theta0 = 26.5
    thetas = thetas0
    dthetas = 0.

surface = nz * (theta0-thetas) / theta0

# length of year, integration time, time step, no steps, output times
year = 31557600. # NOTE: this is not correct as it is calculated 365.25*24*60*60, should have been 365.2425*24*60*60
tmax = 1.e4*year
tmax = 0.5e4*year
dt   = 0.25e-2*year

if (kappav0 * nbl > 1e-5):
    dt = dt * 1.e-5 / (kappav0 * nbl)
if (qek0 < 15e-6):
    dt = 0.1e-2 * year
if (socean < 0.1):
    dt = 0.025e-2 * year
if (nmodel):
    dt = 0.01e-2 * year
nstop = int(tmax/dt)

# output counter
ndump = 0

# define anthropogenic forcing time series:
#   (a) anthropogenic forcing; (b) AMOC collapse
t1 = tmax - 1000.*year
t2 = t1 + 200.*year
t3 = t1 + 100.*year
n1 = int(t1/dt)
n2 = int(t2/dt)
n3 = int(t3/dt)
fanth = np.zeros(nmax)
famoc = np.zeros(nmax)
fsst  = np.zeros(nmax)
if not nmodel:
    fanth[n1:n2-1] = np.sin(0.5 * pi * (np.arange(n1+1, n2)*dt - t1) / (t2 - t1))
    fanth[n2-1:nstop] = 1
    famoc[n3+1:n2-1] = np.sin(0.5 * pi * (np.arange(n3+2, n2)*dt - t3) / (t2 - t3)) ** 2
    famoc[n2-1:nstop] = 1
else:
    n1 = int((tmax-nf*year)/dt)
    for n in xrange(1, nstop + 1):
        t2 = n * dt
        t3 = (t2 - tmax)/year + nf
        n3 = 1
        if t3+0.5 > nf:
            fanth[n-1] = 1.
            famoc[n-1] = 1.
            fsst[n-1] = 1.
        elif t3+0.5 > nf:
            fanth[n-1] = qekin[nf - 1] / qek0
            famoc[n-1] = qnin[nf - 1] / qn0
            fsst[n-1] = sstin[nf - 1] / sstin[0]
        else:
            n3 = int(t3+0.5)
            t1 = t3 - n3 + 0.5
            fanth[n-1] = (t1*qekin[n3] + (1-t1)*qekin[n3 - 1])/qek0
            famoc[n-1] = (t1*qnin[n3] + (1-t1)*qnin[n3 - 1])/qn0
            fsst[n-1] = (t1*sstin[n3] + (1-t1)*sstin[n3 - 1])/sstin[0]

# constants required for heat content
rho0 = 1027.
cp   = 3992.

# surface area north of ACC, length and width of ACC
a  = 2.e14
lx = 2.e7
ly = 1.5e6

# adams bashforth parameters
ab = np.array([(23./12.)*dt/a, (16./12.)*dt/a, (5./12.)*dt/a])
qtot = np.zeros((nz, nt))

# temperature and initial depths fo interfaces
theta = theta0 * np.arange(nz, -1, -1) / nz
dtheta = theta0 / nz

index = np.arange(nz + 1)
mask = theta[index] < thetas
d = np.zeros(nz + 1)
d[mask] = depth*(index[mask] - surface)/(nz - surface)


# minimum upper layer thickness (numerical parameter)
hmin = 1.e-3

# main loop
qn          = np.zeros(nz)
qu          = np.zeros(nz)
qeddy       = np.zeros(nz)
qek         = np.zeros(nz)
kappav      = np.zeros(nz)
heatout     = np.zeros(nout)
heat07out   = np.zeros(nout)
heat72out   = np.zeros(nout)
heat25out   = np.zeros(nout)
thetasout   = np.zeros(nout)
heatupout   = np.zeros(nout+2)
heatekout   = np.zeros(nout+2)
heatuout    = np.zeros(nout+2)
heatnout    = np.zeros(nout+2)
heateddyout = np.zeros(nout+2)
qekout      = np.zeros((nz+1, nout))
qeddyout    = np.zeros((nz+1, nout))
qnout       = np.zeros((nz+1, nout))
quout       = np.zeros((nz+1, nout))
qtotout     = np.zeros((nz+1, nout))
dout        = np.zeros((nz+1, nout))
heatuptake  = 0.

for n in xrange(nstop):
    #if n % 100000 == 0:
    #    print n
    thetas = thetas0 + dthetas * fanth[n]
    if nmodel:
        thetas = sttin[0] * fsst[n]
    surface = nz * (theta0 - thetas) / theta0

    # find surface layer and fraction not outcropped
    kmin = int(surface + 1)
    delta = kmin - surface

    # calculate layer thicknesses
    h = np.maximum(d[1:] - d[:-1], hmin)
    dmid = 0.5 * (d[:-1] + d[1:])

    qn00 = qn0 * (1. - dqn*famoc[n])
    tnadw1 = tnadw10 + dtnadw*famoc[n]
    tnadw2 = tnadw20 + dtnadw*famoc[n]
    if ndeboer:
        nadwfactor = d[80] ** 2
        if n <= n1:
            nadwfactor0 = nadwfactor
        qn00 = qn0 * nadwfactor / nadwfactor0
    if nmodel:
        qn00 = qn0 * famoc[n]
        tnadw1 = tnadw10 + dtnadw*famoc[n]
        tnadw2 = tnadw20 + dtnadw*famoc[n]

    # qn
    qn[kmin-1:nz-1] = 0. # NOTE: this is not optimal
    qnmask = np.arange(kmin-1, nz-1)[(theta[kmin:nz] <= tnadw1) & (theta[kmin:nz] > tnadw2)]
    thetamask = qnmask + 1
    qn[qnmask] = qn00 * np.cos(0.5 * pi * ((tnadw1 - theta[thetamask]) / (tnadw1 - tnadw2)))**2
    qnmask = np.arange(kmin-1, nz-1)[(theta[kmin:nz] <= thetas) & (theta[kmin:nz] > tnadw1)]
    thetamask = qnmask + 1
    qn[qnmask] = qn00 * np.sin(0.5 * pi * ((thetas - theta[thetamask]) / (thetas - tnadw1)))

    # qeddy
    kappagm = kappagm0 * (1. + deddy*fanth[n])
    qeddy[kmin-1:nz-1] = kappagm*d[kmin:nz]*(lx/ly)/np.maximum((thetas-theta[kmin:nz])/thetas, 1.e-10)
    qeddymask = np.arange(kmin-1, nz-1)[d[kmin:nz] > ddrake]
    dmask = qeddymask + 1
    qeddy[qeddymask] *= (depth-d[dmask])/(depth-ddrake)

    # qek
    if nmodel:
        qek00 = qek0 * fanth[n] * socean
    else:
        qek00 = qek0*(1 + dwind*fanth[n])
    qek[kmin-1:nz-1] = qek00 # NOTE: this is not optimal
    qekmask = np.arange(kmin-1, nz-1)[(d[kmin:nz] > ddrake)]
    qek[qekmask] = qek00*(depth - d[qekmask + 1])/(depth - ddrake)
    qekmask = np.arange(kmin-1, nz-1)[theta[kmin:nz] > thetas - 10.]
    qek[qekmask] = qek00*(thetas - theta[qekmask + 1])/10.

    if nsimmons:
        kappav[kmin-1:nz] = np.maximum(kappas, 3.85e-7*(dmid[kmin-1:nz] - 3000.))
    else:
        kappav[kmin-1:nz] = c0 + c1 * np.arctan(c2*(dmid[kmin-1:nz] - c3))

    qu[kmin-1] = a * (delta*kappav[kmin-1]/h[kmin-1] - kappav[kmin]/h[kmin])
    qu[kmin:nz-1] = a * (kappav[kmin:nz-1]/h[kmin:nz-1] - kappav[kmin+1:nz]/h[kmin+1:nz])

    qtot[kmin-1:nz-1, 0] = qek[kmin-1:nz-1] + qu[kmin-1:nz-1] - qn[kmin-1:nz-1] - qeddy[kmin-1:nz-1]

    if ((n + 1) <= n1 and (n + 1) % (nstop/10) == 0) or ((n + 1) >= n1 and (n + 1) % (nstop / 50) == 0):
        for k in xrange(kmin, nz):

            if k == kmin:
                print "time:", (n+1.)/nstop, "n:", (n+1)
            print theta[k], d[k], 1.e-6*qek[k-1],-1.e-6*qeddy[k-1], -1.e-6*qn[k-1], 1.e-6*qu[k-1], 1.e-6*qtot[k-1, 0]
            if k == nz - 1:
                print "heat flux:", heatuptake/a


    # heat content and uptake (including decomposition=
    #   subscripts 1-3 refer to >10, 5-10, <5 degrees
    #   (NB: temperature classes hard-wired here for now)

    heatek     = rho0 * cp * dtheta * np.sum(qek[kmin:nz-1])
    heatn      = rho0 * cp * dtheta * np.sum(qn[kmin:nz-1])
    heatu      = rho0 * cp * dtheta * np.sum(qu[kmin:nz-1])
    heateddy   = rho0 * cp * dtheta * np.sum(qeddy[kmin:nz-1])
    heatuptake = rho0 * cp * dtheta * np.sum(qtot[kmin:nz-1, 0])
    heatstart  = rho0 * cp * dtheta * a * d[kmin]
    heat       = heatstart + rho0 * cp * dtheta * a * np.sum(d[kmin+1:nz])
    heat07     = heatstart + rho0 * cp * dtheta * a * np.sum(np.minimum(d[kmin+1:nz], 700.))
    heat72     = heatstart + rho0 * cp * dtheta * a * np.sum(np.minimum(d[kmin+1:nz], 2000.))
    heat25 = heat - heat72
    heat72 -= heat07

    d[kmin:nz] += ab[0]*qtot[kmin-1:nz-1, 0] + ab[1]*qtot[kmin-1:nz-1, 1] + ab[2]*qtot[kmin-1:nz-1, 2]
    d[kmin] = np.max(d[kmin], hmin)

    if n >= n1 and (n - n1) % ((nstop - n1 - 1) / (nout - 1)) == 0:
        ndump += 1
        nindex = ndump - 1
        for k in xrange(nz+1):
            if k < kmin or k == nz:
                qekout[k, nindex]   = 0.
                qeddyout[k, nindex] = 0.
                qnout[k, nindex]    = 0.
                quout[k, nindex]    = 0.
                qtotout[k, nindex]  = 0.
                dout[k, nindex]     = 0.
            else:
                qekout[k, nindex]   = qek[k-1]*1.e-6
                qeddyout[k, nindex] = qeddy[k-1]*1.e-6
                qnout[k, nindex]    = qn[k-1]*1.e-6
                quout[k, nindex]    = qu[k-1]*1.e-6
                qtotout[k, nindex]  = qtot[k-1, 0]*1.e-6
                dout[k, nindex]     = d[k]
            dout[nz, nindex]    = d[nz]
            heatout[nindex]     = heat
            heat07out[nindex]   = heat07
            heat72out[nindex]   = heat72
            heat25out[nindex]   = heat25
            heatupout[nindex]   = heatuptake
            heatekout[nindex]   = heatek
            heatuout[nindex]    = heatu
            heatnout[nindex]    = heatn
            heateddyout[nindex] = heateddy
            thetasout[nindex]   = thetas

    qtot[kmin-1:nz-1, 2] = qtot[kmin-1:nz-1, 1]
    qtot[kmin-1:nz-1, 1] = qtot[kmin-1:nz-1, 0]

np.savetxt("qek.dat", qekout, delimiter="\n")
np.savetxt("qeddy.dat", qeddyout, delimiter="\n")
np.savetxt("qn.dat", qnout, delimiter="\n")
np.savetxt("qu.dat", quout, delimiter="\n")
np.savetxt("qtot.dat", qtotout, delimiter="\n")
np.savetxt("d.dat", dout, delimiter="\n")
np.savetxt("heat.dat", heatout, delimiter="\n")
np.savetxt("heat07.dat", heat07out, delimiter="\n")
np.savetxt("heat72.dat", heat72out, delimiter="\n")
np.savetxt("heat25.dat", heat25out, delimiter="\n")

heatupout[0] = heatupout[1]
heatupout[nout+1] = heatupout[nout]
heatupwrite = 0.25*(2.*heatupout[1:ndump] + heatupout[2:ndump+1] + heatupout[:ndump-1])/a
np.savetxt("heatup.dat", heatupwrite, delimiter="\n")

heatekout[0] = heatekout[1]
heatekout[nout+1] = heatekout[nout]
heatekwrite = 0.25*(2.*heatekout[1:ndump+1] + heatekout[2:ndump+2] + heatekout[:ndump])/a
np.savetxt("heatek.dat", heatekwrite, delimiter="\n")

heatuout[0] = heatuout[1]
heatuout[nout+1] = heatuout[nout]
heatuwrite = 0.25*(2.*heatuout[1:ndump+1] + heatuout[2:ndump+2] + heatuout[:ndump])/a
np.savetxt("heatu.dat", heatuwrite, delimiter="\n")

heateddyout[0] = heateddyout[1]
heateddyout[nout+1] = heateddyout[nout]
heateddywrite = 0.25*(2.*heateddyout[1:ndump+1] + heateddyout[2:ndump+2] + heateddyout[:ndump])/a
np.savetxt("heateddy.dat", heateddywrite, delimiter="\n")

heatnout[0] = heatnout[1]
heatnout[nout+1] = heatnout[nout]
heatnwrite = 0.25*(2.*heatnout[1:ndump+1] + heatnout[2:ndump+2] + heatnout[:ndump])/a
np.savetxt("heatn.dat", heatnwrite, delimiter="\n")

np.savetxt("thetas.dat", thetasout, delimiter="\n")
