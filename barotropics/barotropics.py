import climate
from climate.pyom import PyOM
from climate.pyom.external import solve_stream
from climate.pyom import numerics
import numpy as np
import sys
np.set_printoptions(threshold=np.nan)

yt_start = -39.0
yt_end   = 43
yu_start = -40.0
yu_end   = 42

class simulation(PyOM):

   """ a simple global model with a Southern Ocean and Atlantic part
   """
   def set_parameter(self):
     """set main parameter
     """

     (self.nx,self.ny,self.nz)    = (30,42,15)
     self.dt_mom    = 4800  
     self.dt_tracer = 86400/2.0

     self.coord_degree     = 1
     self.enable_cyclic_x  = 1

     self.congr_epsilon = 1e-12
     self.congr_max_iterations = 5000
     self.enable_streamfunction = 1

     self.enable_neutral_diffusion =  1
     self.k_iso_0 = 1000.0
     self.k_iso_steep = 500.0
     self.iso_dslope=0.005
     self.iso_slopec=0.01
     self.enable_skew_diffusion = 1

     self.enable_hor_friction = 1
     self.a_h = (2*self.degtom)**3*2e-11    
     self.enable_hor_friction_cos_scaling = 1
     self.hor_friction_cospower=1

     self.enable_bottom_friction = 1
     self.r_bot = 1e-5

     self.enable_implicit_vert_friction = 1
     self.enable_tke = 1
     self.c_k = 0.1
     self.c_eps = 0.7
     self.alpha_tke = 30.0
     self.mxl_min = 1e-8
     self.tke_mxl_choice = 2

     self.k_gm_0 = 1000.0
     self.enable_eke = 1
     self.eke_k_max  = 1e4
     self.eke_c_k    = 0.4
     self.eke_c_eps  = 0.5
     self.eke_cross  = 2.
     self.eke_crhin  = 1.0
     self.eke_lmin   = 100.0
     self.enable_eke_superbee_advection = 1
     self.enable_eke_isopycnal_diffusion = 1

     self.enable_idemix = 1
     self.enable_idemix_hor_diffusion = 1
     self.enable_eke_diss_surfbot = 1
     self.eke_diss_surfbot_frac = 0.2
     self.enable_idemix_superbee_advection = 1

     self.eq_of_state_type = 3 
     return


   def set_grid(self):
     ddz  = np.array([50.,70.,100.,140.,190.,240.,290.,340.,390.,440.,490.,540.,590.,640.,690.])
     self.dxt[:] = 2.0
     self.dyt[:] = 2.0
     self.x_origin=  0.0
     self.y_origin= -40.0
     self.dzt[:] = ddz[::-1]/2.5
     return


   def set_coriolis(self):
     #print self.omega, np.pi
     self.coriolis_t[:,:] = 2*self.omega*np.sin(self.yt/180.*np.pi)
     return

   def set_topography(self):
     """ setup topography
     """
     (X,Y)= np.meshgrid(self.xt,self.yt); X=X.transpose(); Y=Y.transpose()
     self.kbot[:]=0

     if climate.is_bohrium:
         kbot = self.kbot.copy2numpy()
     else:
         kbot = self.kbot

     kbot[ X >1.0]=1
     kbot[ Y <-20]=1
     if climate.is_bohrium:
         self.kbot = np.array(kbot)
     return
   
   def set_initial_conditions(self):
     """ setup initial conditions
     """
     (XT,YT)= np.meshgrid(self.xt,self.yt); XT=XT.transpose(); YT=YT.transpose()
     (XU,YU)= np.meshgrid(self.xu,self.yu); XU=XU.transpose(); YU=YU.transpose()

     # initial conditions
     self.temp[:,:,:,self.tau] = (1-self.zt/self.zw[0])*15*self.maskt
     self.temp[:,:,:,self.taum1] = (1-self.zt/self.zw[0])*15*self.maskt
     #for k in range(self.nz):
     #   self.temp[:,:,k,self.tau]   =  (1-self.zt[k]/self.zw[0])*15*self.maskt[:,:,k]
     #   self.temp[:,:,k,self.taum1] =  (1-self.zt[k]/self.zw[0])*15*self.maskt[:,:,k]
     self.salt[:,:,:,self.tau]   = 35.0*self.maskt[:]
     self.salt[:,:,:,self.taum1] = 35.0*self.maskt[:]

     # wind stress forcing
     for j in range(self.js_pe,self.je_pe+1):
          jj = self.jf2py(j)
          taux=0.0
          if  self.yt[jj]<-20 : taux =  .1e-3*np.sin(np.pi*(self.yu[jj]-yu_start)/(-20.0-yt_start))
          if  self.yt[jj]>10  : taux =  .1e-3*(1-np.cos(2*np.pi*(self.yu[jj]-10.0)/(yu_end-10.0)))
          self.surface_taux[:,jj] = taux*self.masku[:,jj,-1]

     # surface heatflux forcing      
     self.t_rest = np.zeros( self.u[:,:,1,0].shape)
     self.t_star = np.zeros( self.u[:,:,1,0].shape)

     self.t_rest = self.dzt[-1]/(30.*86400.)*self.maskT[:,:,-1]
     for j in range(self.js_pe,self.je_pe+1):
           jj = self.jf2py(j) 
           t_star=15
           if self.yt[jj]<-20.0: t_star=15*(self.yt[jj]-yt_start )/(-20.0-yt_start)
           if self.yt[jj]> 20.0: t_star=15*(1-(self.yt[jj]-20)/(yt_end-20) )
           self.t_star[:,jj] = t_star


     T=self.fortran.tke_module   
     if T.enable_tke:
        T.forc_tke_surface[2:self.nx+2, 2:self.ny+2] = np.sqrt((0.5*(self.surface_taux[2:self.nx+2,2:self.ny+2]+self.surface_taux[1:self.nx+1,2:self.ny+2]))**2  \
                +(0.5*(self.surface_tauy[2:self.nx+2,2:self.ny+2]+self.surface_tauy[2:self.nx+2,1:self.ny+1]))**2 )**(3./2.)
        #for j in range(self.js_pe,self.je_pe+1):
        #   for i in range(self.is_pe,self.ie_pe+1):
        #     (ii,jj) = (self.if2py(i), self.jf2py(j) )
        #     T.forc_tke_surface[ii,jj] = np.sqrt( (0.5*(self.surface_taux[ii,jj]+self.surface_taux[ii-1,jj]))**2  \
        #                                      +(0.5*(self.surface_tauy[ii,jj]+self.surface_tauy[ii,jj-1]))**2 )**(3./2.)


     I=self.fortran.idemix_module   
     if I.enable_idemix:
       I.forc_iw_bottom[:] =  1.0e-6*self.maskw[:,:,-1]
       I.forc_iw_surface[:] = 0.1e-6*self.maskw[:,:,-1]
     return


   def set_forcing(self):
     self.forc_temp_surface[:]=self.t_rest*(self.t_star-self.temp[:,:,-1,self.tau-1])
     return

   def set_diagnostics(self):
     self.register_average(name='temp',long_name='Temperature',         units = 'deg C' , grid = 'TTT', var = self.temp)
     self.register_average(name='salt',long_name='Salinity',            units = 'g/kg' ,  grid = 'TTT', var = self.salt)
     self.register_average(name='u',   long_name='Zonal velocity',      units = 'm/s' ,   grid = 'UTT', var = self.u)
     self.register_average(name='v',   long_name='self.ridional velocity', units = 'm/s' ,   grid = 'TUT', var = self.v)
     self.register_average(name='w',   long_name='Vertical velocity',   units = 'm/s' ,   grid = 'TTU', var = self.w)
     self.register_average(name='taux',long_name='wind stress',         units = 'm^2/s' , grid = 'UT',  var = self.surface_taux)
     self.register_average(name='tauy',long_name='wind stress',         units = 'm^2/s' , grid = 'TU',  var = self.surface_tauy)
     self.register_average(name='psi' ,long_name='Streamfunction',      units = 'm^3/s' , grid = 'UU',  var = self.psi)
     return

   def user_defined_signal(self):
       """ this routine must be called by all processors
       """
       a = zeros( (self.nx,self.ny), 'd', order = 'F')
       a[self.is_pe-1:self.ie_pe,0] = self.xt[2:-2]
       self.fortran.pe0_recv_2d(a)
       self.xt_gl = a[:,0].copy()
       
       a[0,self.js_pe-1:self.je_pe] = self.yt[2:-2]
       self.fortran.pe0_recv_2d(a)
       self.yt_gl = a[0,:].copy()
       
       self.psi_gl = zeros( (self.nx,self.ny), 'd', order = 'F')
       self.psi_gl[self.is_pe-1:self.ie_pe,self.js_pe-1:self.je_pe] = where( self.maskz[2:-2,2:-2,-1] >0,  self.psi[2:-2,2:-2,self.tau-1] , NaN) 
       self.fortran.pe0_recv_2d(self.psi_gl)
       
       self.temp_gl = zeros( (self.nx,self.ny,self.nz), 'd', order = 'F')
       for k in range(self.nz):
         a[self.is_pe-1:self.ie_pe,self.js_pe-1:self.je_pe] = where( self.maskt[2:-2,2:-2,k] >0,  self.temp[2:-2,2:-2,k,self.tau-1] , NaN) 
         self.fortran.pe0_recv_2d(a)
         self.temp_gl[:,:,k]=a.copy()

       self.kappa_gl = zeros( (self.nx,self.ny,self.nz), 'd', order = 'F')
       for k in range(self.nz):
         a[self.is_pe-1:self.ie_pe,self.js_pe-1:self.je_pe] = where( self.maskw[2:-2,2:-2,k] >0,  self.kappah[2:-2,2:-2,k] , NaN) 
         self.fortran.pe0_recv_2d(a)
         self.kappa_gl[:,:,k]=a.copy()

       return
   
   def make_plot(self):
       
       self.set_signal('user_defined') # following routine is called by all PEs
       self.user_defined_signal()
       
       self.figure.clf()
       ax=self.figure.add_subplot(221)
       
       co=ax.contourf(self.yt_gl,self.zt,self.temp_gl[self.nx/2-1,:,:].transpose())
       self.figure.colorbar(co)
       ax.set_title('temperature')
       ax.set_ylabel('z [m]')
       ax.axis('tight')

       ax=self.figure.add_subplot(223)
       try:
        co=ax.contourf(self.yt_gl,self.zw,log10(self.kappa_gl[self.nx/2-1,:,:].transpose()) )
       except:
        pass
       self.figure.colorbar(co)
       ax.set_title('Diffusivity')
       ax.set_xlabel('Latitude [deg N]')
       ax.set_ylabel('z [m]')
       ax.axis('tight')
       
       ax=self.figure.add_subplot(122)
       co=ax.contourf(self.xt_gl,self.yt_gl,self.psi_gl.transpose()*1e-6)
       self.figure.colorbar(co)
       ax.set_title('Streamfunction [Sv]')
       ax.set_xlabel('Longitude [deg E]')
       ax.axis('tight')
       
       return

baroTimer = climate.Timer("barotropics")
with baroTimer:
    a = simulation()
    """
     allocate everything
    """
    #print "set parameter"
    a.set_parameter()
    #print "set parameter end"
    #if climate.is_bohrium:
    #    np.flush()
    #print "allocate"
    a.allocate()
    #print "allocate end"
    #if climate.is_bohrium:
    #    np.flush()

    """
      Grid
    """
    #print "set grid"
    a.set_grid()
    #print "done"
    #if climate.is_bohrium:
    #    np.flush()
    #print "calc grid"
    numerics.calc_grid(a)
    #print "done"
    #if climate.is_bohrium:
    #    np.flush()

    """
     Coriolis
    """
    #print "set coriolis"
    a.set_coriolis()
    #print "done"
    #if climate.is_bohrium:
    #    np.flush()
    #print "calc beta"
    numerics.calc_beta(a)
    #print "done"
    #if climate.is_bohrium:
    #    np.flush()

    """
     topography
    """
    #print "set topo"
    a.set_topography()
    #print "done"

    #if climate.is_bohrium:
    #    np.flush()
    #print "calc topo"
    numerics.calc_topo(a)
    #print "done"
    #if climate.is_bohrium:
    #    np.flush()

    a.set_initial_conditions()
    numerics.calc_initial_conditions(a)

    #print "solve init"
    solve_stream.streamfunction_init(a)
    #print "done solve init"
    #if climate.is_bohrium:
    #    np.flush()
    #print "solve stream"
    solve_stream.solve_streamfunction(a)
    #print "done solve stream"
    #if climate.is_bohrium:
    #    np.flush()
