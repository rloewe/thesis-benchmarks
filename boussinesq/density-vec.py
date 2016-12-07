from climate.boussinesq.density import gsw, linear_eq, nonlinear_eq1, nonlinear_eq2, nonlinear_eq3
from climate import Timer
import numpy as np

try:     # try to load module with netcdf bindings
    from netCDF4 import Dataset as NF
except ImportError:
    from Scientific.IO.NetCDF import NetCDFFile as NF

nx=200
ny=200
nz=200

T = -2.+27.*np.arange(nx)/nx
S = 33+4.*np.arange(ny)/ny
P = 5000*np.arange(nz)/nz
drho2dS = np.empty((nx, ny, nz))

Sj = np.ones(ny)*S[np.newaxis, :].T*np.ones(ny)[np.newaxis, np.newaxis, :].T
Ti = np.ones(ny)*np.ones(ny)[np.newaxis, :].T*T[np.newaxis, np.newaxis, :].T
Pk = P*np.ones(ny)[np.newaxis, :].T*np.ones(ny)[np.newaxis, np.newaxis, :].T

timer = Timer("total")
with timer:
    rho3 = gsw.gsw_rho(Sj, Ti, P)
    np.flush()
    rho2 = nonlinear_eq2.nonlin2_eq_of_state_rho(Sj, Ti, P)
    np.flush()
    drho3dT = gsw.gsw_drhodT(Sj,Ti,P)
    np.flush()
    drho3dS = gsw.gsw_drhodS(Sj,Ti,P)
    np.flush()
    drho3dP = gsw.gsw_drhodP(Sj,Ti,P)
    np.flush()
    drho2dT = nonlinear_eq2.nonlin2_eq_of_state_drhodT(Ti,P)
    np.flush()
    #drho2dS.fill(nonlinear_eq2.nonlin2_eq_of_state_drhodS())
    np.flush()
    drho2dP = nonlinear_eq2.nonlin2_eq_of_state_drhodP(Ti)
    np.flush()
    Hd3 = gsw.gsw_dyn_enthalpy(Sj,Ti,P)
    np.flush()
    Hd2 = nonlinear_eq2.nonlin2_eq_of_state_dyn_enthalpy(Sj,Ti,P)
    np.flush()
    dHd3dT = gsw.gsw_dHdT(Sj,Ti,P)
    np.flush()
    dHd3dS = gsw.gsw_dHdS(Sj,Ti,P)
    np.flush()
    dHd2dT = nonlinear_eq2.nonlin2_eq_of_state_int_drhodT(Ti,P)*(-9.81/1024.)
    np.flush()
    dHd2dS = nonlinear_eq2.nonlin2_eq_of_state_int_drhodS(Pk)*(-9.81/1024.)
    np.flush()

timer.printTime()

f = NF("check_density.cdf", "r")

if not (np.abs(f.variables["T"] - T) < 1e-5).all() or not f.variables["T"][:].shape == T.shape:
    print "T"
if not (np.abs(f.variables["S"] - S) < 1e-5).all() or not f.variables["S"][:].shape == S.shape:
    print "S"
if not (np.abs(f.variables["P"] - P) < 1e-5).all() or not f.variables["P"][:].shape == P.shape:
    print "P"
if not (np.abs(f.variables["rho3"][:].T - rho3) < 1e-3).all() or not f.variables["rho3"][:].shape == rho3.shape:
    print "rho3"
if not (np.abs(f.variables["rho2"][:].T - rho2) < 1e-4).all() or not f.variables["rho2"][:].shape == rho2.shape:
    print "rho2"
if not (np.abs(f.variables["drho3dT"][:].T - drho3dT) < 1e-4).all() or not f.variables["drho3dT"][:].shape == drho3dT.shape:
    print "drho3dT"
if not (np.abs(f.variables["drho2dT"][:].T - drho2dT) < 1e-4).all() or not f.variables["drho2dT"][:].shape == drho2dT.shape:
    print "drho2dT"
if not (np.abs(f.variables["drho3dS"][:].T - drho3dS) < 1e-4).all() or not f.variables["drho3dS"][:].shape == drho3dS.shape:
    print "drho3dS"
if not (np.abs(f.variables["drho2dS"][:].T - drho2dS) < 1e-5).all() or not f.variables["drho2dS"][:].shape == drho2dS.shape:
    print "drho2dS"
if not (np.abs(f.variables["drho3dP"][:].T - drho3dP) < 1e-5).all() or not f.variables["drho3dP"][:].shape == drho3dP.shape:
    print "drho3dP"
if not (np.abs(f.variables["drho2dP"][:].T - drho2dP) < 1e-5).all() or not f.variables["drho2dP"][:].shape == drho2dP.shape:
    print "drho2dP"
if not (np.abs(f.variables["Hd3"][:].T - Hd3) < 1e-3).all() or not f.variables["Hd3"][:].shape == Hd3.shape:
    print "Hd3"
if not (np.abs(f.variables["Hd2"][:].T - Hd2) < 1e-3).all() or not f.variables["Hd2"][:].shape == Hd2.shape:
    print "Hd2"
if not (np.abs(f.variables["dHd3dT"][:].T - dHd3dT) < 1e-4).all() or not f.variables["dHd3dT"][:].shape == dHd3dT.shape:
    print "dHd3dT"
if not (np.abs(f.variables["dHd3dS"][:].T - dHd3dS) < 1e-4).all() or not f.variables["dHd3dS"][:].shape == dHd3dS.shape:
    print "dHd3dS"
if not (np.abs(f.variables["dHd2dT"][:].T - dHd2dT) < 1e-4).all() or not f.variables["dHd2dT"][:].shape == dHd2dT.shape:
    print "dHd2dT"
if not (np.abs(f.variables["dHd2dS"][:].T - dHd2dS) < 1e-4).all() or not f.variables["dHd2dS"][:].shape == dHd2dS.shape:
    print "dHd2dS"
