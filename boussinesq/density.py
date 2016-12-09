from climate.boussinesq.density import gsw, linear_eq, nonlinear_eq1, nonlinear_eq2, nonlinear_eq3
import numpy

try:     # try to load module with netcdf bindings
    from netCDF4 import Dataset as NF
except ImportError:
    from Scientific.IO.NetCDF import NetCDFFile as NF

nx=200
ny=200
nz=200

T = numpy.zeros(nx)
S = numpy.zeros(ny)
P = numpy.zeros(nz)
rho1 = numpy.zeros((nx, ny, nz))
rho2 = numpy.zeros((nx, ny, nz))
rho3 = numpy.zeros((nx, ny, nz))
drho1dT = numpy.zeros((nx, ny, nz))
drho2dT = numpy.zeros((nx, ny, nz))
drho3dT = numpy.zeros((nx, ny, nz))
drho1dS = numpy.zeros((nx, ny, nz))
drho2dS = numpy.zeros((nx, ny, nz))
drho3dS = numpy.zeros((nx, ny, nz))
drho1dP = numpy.zeros((nx, ny, nz))
drho2dP = numpy.zeros((nx, ny, nz))
drho3dP = numpy.zeros((nx, ny, nz))
Hd1 = numpy.zeros((nx, ny, nz))
Hd2 = numpy.zeros((nx, ny, nz))
Hd3 = numpy.zeros((nx, ny, nz))
dHd3dT = numpy.zeros((nx, ny, nz))
dHd2dT = numpy.zeros((nx, ny, nz))
dHd3dS = numpy.zeros((nx, ny, nz))
dHd2dS = numpy.zeros((nx, ny, nz))

for i in xrange(nx):
    for j in xrange(ny):
        for k in xrange(nz):
            T[i] = -2.+27.*float(i)/nx
            S[j] = 33+4*float(j)/ny
            #S[j] = 35*(j-1.0)/ny
            P[k] = 5000*float(k)/nz
            rho3[i,j,k] = gsw.gsw_rho(S[j],T[i],P[k])
            rho2[i,j,k] = nonlinear_eq2.nonlin2_eq_of_state_rho(S[j],T[i],P[k])
            drho3dT[i,j,k] = gsw.gsw_drhodT(S[j],T[i],P[k])
            drho2dT[i,j,k] = nonlinear_eq2.nonlin2_eq_of_state_drhodT(T[i],P[k])
            drho3dS[i,j,k] = gsw.gsw_drhodS(S[j],T[i],P[k])
            drho2dS[i,j,k] = nonlinear_eq2.nonlin2_eq_of_state_drhodS()
            drho3dP[i,j,k] = gsw.gsw_drhodP(S[j],T[i],P[k])
            drho2dP[i,j,k] = nonlinear_eq2.nonlin2_eq_of_state_drhodP(T[i])
            Hd3[i,j,k] = gsw.gsw_dyn_enthalpy(S[j],T[i],P[k])
            Hd2[i,j,k] = nonlinear_eq2.nonlin2_eq_of_state_dyn_enthalpy(S[j],T[i],P[k])
            dHd3dT[i,j,k] = gsw.gsw_dHdT(S[j],T[i],P[k])
            dHd3dS[i,j,k] = gsw.gsw_dHdS(S[j],T[i],P[k])
            dHd2dT[i,j,k] = nonlinear_eq2.nonlin2_eq_of_state_int_drhodT(T[i],P[k])*(-9.81/1024.)
            dHd2dS[i,j,k] = nonlinear_eq2.nonlin2_eq_of_state_int_drhodS(P[k])*(-9.81/1024.)

f = NF("check_density.cdf", "r")

if not (numpy.abs(f.variables["T"] - T) < 1e-5).all():
    print "T"
if not (numpy.abs(f.variables["S"] - S) < 1e-5).all():
    print "S"
if not (numpy.abs(f.variables["P"] - P) < 1e-5).all():
    print "P"
if not (numpy.abs(f.variables["rho3"][:].T - rho3) < 1e-4).all():
    print "rho3"
if not (numpy.abs(f.variables["rho2"][:].T - rho2) < 1e-4).all():
    print "rho2"
if not (numpy.abs(f.variables["drho3dT"][:].T - drho3dT) < 1e-4).all():
    print "drho3dT"
if not (numpy.abs(f.variables["drho2dT"][:].T - drho2dT) < 1e-4).all():
    print "drho2dT"
if not (numpy.abs(f.variables["drho3dS"][:].T - drho3dS) < 1e-4).all():
    print "drho3dS"
if not (numpy.abs(f.variables["drho2dS"][:].T - drho2dS) < 1e-5).all():
    print "drho2dS"
if not (numpy.abs(f.variables["drho3dP"][:].T - drho3dP) < 1e-5).all():
    print "drho3dP"
if not (numpy.abs(f.variables["drho2dP"][:].T - drho2dP) < 1e-5).all():
    print "drho2dP"
if not (numpy.abs(f.variables["Hd3"][:].T - Hd3) < 1e-4).all():
    print "Hd3"
if not (numpy.abs(f.variables["Hd2"][:].T - Hd2) < 1e-3).all():
    print "Hd2"
if not (numpy.abs(f.variables["dHd3dT"][:].T - dHd3dT) < 1e-4).all():
    print "dHd3dT"
if not (numpy.abs(f.variables["dHd3dS"][:].T - dHd3dS) < 1e-4).all():
    print "dHd3dS"
if not (numpy.abs(f.variables["dHd2dT"][:].T - dHd2dT) < 1e-4).all():
    print "dHd2dT"
if not (numpy.abs(f.variables["dHd2dS"][:].T - dHd2dS) < 1e-4).all():
    print "dHd2dS"
