      program gnanadesikan
c     
c     function: multilayer Gnanadesikan model
c     May 2014 - corrected for variable mixing bug
c
c     this version includes Bryan-Lewis profile of diapycnal mixing
c     also incoporate GCM forcing data as an option
c
      implicit none
c
      integer nz,nt,nout,nmax
      integer nf,nmodel,ndeboer,nsimmons
c
      parameter(nz=100,nt=3,nout=1001,nmax=100000001)
      parameter(nf=140)
c
      real*8 h(nz),d(0:nz),theta(0:nz),kappav(nz),dmid(nz)
      real*8 qek(nz),qeddy(nz),qn(nz),qu(nz)
      real*8 qtot(nz,nt),ab(nt),fanth(nmax),famoc(nmax)
      real*8 dout(0:nz,nout)
      real*8 qekout(0:nz,nout),qeddyout(0:nz,nout)
      real*8 qnout(0:nz,nout),quout(0:nz,nout),qtotout(0:nz,nout)
      real*8 heatout(nout)
      real*8 heat07out(nout),heat72out(nout),heat25out(nout)
      real*8 heatupout(0:nout+1),thetasout(nout)
      real*8 heatekout(0:nout+1),heatuout(0:nout+1)
      real*8 heateddyout(0:nout+1),heatnout(0:nout+1)
      real*8 pi,year,tmax,dt,t1,t2,t3,depth,qek0,ddrake,socean
      real*8 theta0,thetas,delta,surface,tnadw10,tnadw20,a
      real*8 nbl,c0,c1,c2,c3
      real*8 dwind,qek00,kappagm0,deddy,dtnadw,tnadw1,tnadw2
      real*8 thetas0,dthetas,dtheta,rho0,cp,heat,heatuptake
      real*8 heat07,heat72,heat25
      real*8 hmin,kappav0,kappagm,lx,ly,qn0,dqn,qn00
      real*8 heatn,heatu,heateddy,heatek
      real*8 nadwfactor,nadwfactor0,kappas
      real*8 qekin(nf),qnin(nf),sstin(nf),fsst(nmax)

c
c     nb: surface area of earth - 5.10072e14 m^2 = included in matlab script
c
	  integer nstop,n,k,kmin,n1,n2,n3,ndump,i
c
      pi=3.14159265358979
c
c     run with GCM forcing data? (1 for yes)
        nmodel=0
c
c     Southern ocean switch (0. for no Drake Passage, 1. for Drake Passage)
        socean=1.
c
c     eddy diffusivities (placed here to scale dt with kappav)
c            kappav0: surface diapycnal diffusivity
c            Bryan-Lewis:
c              nb1: amplification of mixing in abyss
c              c0, c1, c2: Bryan-Lewis parameters
c              kappav = c0 + c1 atan[ c2 ( d - c3 ) ]
c            Simmons:
c              set nsimmons=1 to use (ignores other values)
c              kappas: upper ocean diffusivity in Simmons (if used)
c            fractional change in eddy transport (0.1 = +10%)
        kappav0=1.e-5
        nbl=10.
        c0=kappav0*(nbl+1)/2.
        c1=kappav0*(nbl-1.)/pi
        c2=4.5e-3
        c3=2500.
        nsimmons=0
        kappas=3.e-5
c
        kappagm0=1.e3*socean
        deddy=0.
c
c     if using GCM forcing data, read in now
        if (nmodel.eq.1) then
          open(unit=9,file='qekin.dat',status='old')
          open(unit=10,file='qnin.dat',status='old')
          open(unit=11,file='sstin.dat',status='old')
          do i=1,nf
            read(9,*) qekin(i)
            qekin(i)=qekin(i)/1027.
            read(10,*) qnin(i)
            qnin(i)=qnin(i)/1027.
            read(11,*) sstin(i)
c           NB: 1.5 degrees not added to the final figures in this case
            sstin(i)=sstin(i)-296.4+20.5
          end do
          close(9)
          close(10)
          close(11)
        endif
c
c     Southern Ocean Ekman transport, Ekman depth, Drake Passage depth,
c            fractional change in Ekman transport (0.1 = +10%)
        qek0=30.e6*socean
c
        if (nmodel.eq.1) then
          qek0=0.
          do i=1,10
            qek0=qek0+qekin(i)/10.
          end do
          qek0=qekin(1)
c          qek0=33.e6
c          qek0=37.755e6
          qek0=36.82e6
        endif
c
        depth=5.e3
        ddrake=0.8*depth
        dwind=0.0
c
c     NADW formation, fractional reduction in NADW formation (0.5 = -50%),
c            temperature range of NADW formation, NADW warming 
        qn0=20.e6
c
        if (nmodel.eq.1) then
          qn0=0.
          do i=1,10
            qn0=qn0+qnin(i)/10.
          end do
c          qn0=qnin(1)
c          qn0=15.e6
c          qn0=14.4255e6
          qn0=14.77e6
        endif
c
        dqn=0.
c       NB: 1.5 degrees gets added to the final figures (in matlab script)
        tnadw10=4.5
        tnadw20=0.5
        if (nmodel.eq.1) then
          tnadw10=6.0
          tnadw20=2.0
       endif
        dtnadw=0.
c
c       set ndeboer = 1 for dynamic NADW formation following de Boer (2011)
        ndeboer=0
c                
c     surface temperature, surface interface number (fractional)
c         (theta0 hard-wired to keep levels the same across each scenario;
c          small number added to ensure data output times do not coincide 
c          with a new layer appearing at the surface => noisy diagnostics!)
c       NB: 1.5 degrees gets added to the final figures in non-GCM runs
c          (in matlab script)
        thetas0=19.5
c       NB: do not modify the next line; change dthetas at end of this block
        dthetas=4.
        theta0=thetas0+dthetas+1.117e-4
        thetas=thetas0
        dthetas=0.
c
        if (nmodel.eq.1) then
c          thetas0=0.
c          do i=1,10
c            thetas0=thetas0+sstin(i)/10.
c          end do
          thetas0=sstin(1)
c          thetas0=21.5
          theta0=26.5
          thetas=thetas0
          dthetas=0.
        endif
c
        surface=real(nz)*(theta0-thetas)/theta0
c
c     length of year, integration time, time step, no steps, output times
	    year=31557600.
	    tmax=1.e4*year
	    tmax=0.5e4*year
        dt=0.25e-2*year
c
        if (kappav0*nbl.gt.1.e-5) dt=dt*1.e-5/(kappav0*nbl)
        if (qek0.lt.15.e6) dt=0.1e-2*year
        if (socean.lt.0.1) dt=0.025e-2*year
        if (nmodel.eq.1) dt=0.01e-2*year
        nstop=int(tmax/dt)
c
c     output counter
        ndump=0
c
c     define anthropogenic forcing time series: 
c       (a) antropgoenic forcing; (b) AMOC collapse
        t1=tmax-1000.*year
        t2=t1+200.*year
        t3=t1+100.*year
        n1=int(t1/dt)
        n2=int(t2/dt)
        n3=int(t3/dt)
        do n=1,nstop
          if (n.le.n1) then
            fanth(n)=0.
          elseif (n.gt.n1.and.n.lt.n2) then 
            fanth(n)=(n*dt-t1)/(t2-t1)
            fanth(n)=(sin(0.5*pi*fanth(n)))**2
          else 
            fanth(n)=1.
          endif
          if (n.le.n3) then
            famoc(n)=0.
          elseif (n.gt.n3.and.n.lt.n2) then 
            famoc(n)=(n*dt-t3)/(t2-t3)
            famoc(n)=(sin(0.5*pi*famoc(n)))**2
          else 
            famoc(n)=1.
          endif
        enddo
c       overwrite if using GCM forcing data
        if (nmodel.eq.1) then
          n1=int((tmax-real(nf)*year)/dt)
          do n=1,nstop
            t2=real(n)*dt
            t3=(t2-tmax)/year+real(nf)
            n3=1
            if (t3+0.5.le.1.) then
              fanth(n)=1.
              famoc(n)=1.
              fsst(n)=1.
            elseif (t3+0.5.ge.real(nf)) then
              fanth(n)=qekin(nf)/qek0
              famoc(n)=qnin(nf)/qn0
              fsst(n)=sstin(nf)/sstin(1)
            else
              n3=int(t3+0.5)
              t1=t3-real(n3)+0.5
              fanth(n)=(t1*qekin(n3+1)+(1.-t1)*qekin(n3))/qek0
              famoc(n)=(t1*qnin(n3+1)+(1.-t1)*qnin(n3))/qn0
              fsst(n)=(t1*sstin(n3+1)+(1.-t1)*sstin(n3))/sstin(1)
            endif
          end do
        endif
c
c     constants required for heat content
        rho0=1027.
        cp=3992.
c
c     surface area north of ACC, length and width of ACC
        a=2.e14	  
        lx=2.e7
        ly=1.5e6
c
c     adams bashforth parameters
  	    ab(1)=(23./12.)*dt/a
	    ab(2)=-(16./12.)*dt/a
        ab(3)=(5./12.)*dt/a
        do k=1,nz
          qtot(k,2)=0.
          qtot(k,3)=0.
        enddo
c
c     temperature and initial depths of interfaces
        do k=0,nz
          theta(k)=theta0*(real(nz-k)/real(nz))
          dtheta=theta0/real(nz)
          if (theta(k).ge.thetas) then
            d(k)=0.0
          else
            d(k)=depth*(real(k)-surface)/(real(nz)-surface)
          endif
        enddo
c 
c     minimum upper layer thickness (numerical parameter)
        hmin=1.e-3
c
c     main loop
c
      do n=1,nstop
c
c       set thetas
          thetas=thetas0+dthetas*fanth(n)
          if (nmodel.eq.1) then
            thetas=sstin(1)*fsst(n)
          endif
          surface=real(nz)*(theta0-thetas)/theta0
c
c       find surface layer and fraction not outcropped
          kmin=int(surface+1)
    	  delta=real(kmin)-surface
c          
c       calculate layer thicknesses
          do k=1,nz
            h(k)=max(d(k)-d(k-1),hmin)
            dmid(k)=0.5*(d(k-1)+d(k))
          end do
c
        do k=kmin,nz-1
c
          qn00=qn0*(1.0-dqn*famoc(n))
          tnadw1=tnadw10+dtnadw*famoc(n)
          tnadw2=tnadw20+dtnadw*famoc(n)
          if (ndeboer.eq.1) then
            nadwfactor=d(80)**2
            if (n.le.n1) nadwfactor0=nadwfactor
            qn00=qn0*nadwfactor/nadwfactor0
          endif
          if (nmodel.eq.1) then
            qn00=qn0*famoc(n)
            tnadw1=tnadw10+dtnadw*famoc(n)
            tnadw2=tnadw20+dtnadw*famoc(n)
          endif
c
          if (theta(k).le.thetas.and.theta(k).gt.tnadw1) then
            qn(k)=(thetas-theta(k))/(thetas-tnadw1)
            qn(k)=qn00*sin(0.5*pi*qn(k))
          elseif (theta(k).le.tnadw1.and.theta(k).gt.tnadw2) then
            qn(k)=(tnadw1-theta(k))/(tnadw1-tnadw2)
            qn(k)=qn00*(cos(0.5*pi*qn(k)))**2
          else
            qn(k)=0.0
          endif 
c
          kappagm=kappagm0*(1.+deddy*fanth(n))
          qeddy(k)=kappagm*d(k)*(lx/ly)/
     -             max((thetas-theta(k))/thetas,1.e-10)
          if (d(k).gt.ddrake) then
            qeddy(k)=qeddy(k)*(depth-d(k))/(depth-ddrake)
          endif
c     
          qek00=qek0*(1.+dwind*fanth(n))
          if (nmodel.eq.1) then
            qek00=qek0*fanth(n)*socean
          endif
c
          if (theta(k).gt.thetas-10.) then
            qek(k)=qek00*(thetas-theta(k))/10.
          elseif (d(k).gt.ddrake) then
            qek(k)=qek00*(depth-d(k))/(depth-ddrake)
          else
            qek(k)=qek00
          endif
c
        end do
c
        do k=kmin,nz
          kappav(k)=c0+c1*atan(c2*(dmid(k)-c3))
c          print *,k,kappav(k)
          if (nsimmons.eq.1) then
c            kappav(k)=max(3.e-5,3.85e-7*(dmid(k)-3000.))
            kappav(k)=max(kappas,3.85e-7*(dmid(k)-3000.))
          endif
c
        end do
c
        do k=kmin,nz-1
c
          if (k.eq.kmin) then
            qu(k)=a*(delta*kappav(k)/h(k)-kappav(k+1)/h(k+1))
          else
            qu(k)=a*(kappav(k)/h(k)-kappav(k+1)/h(k+1))
          endif
c
          qtot(k,1)=qek(k)+qu(k)-qn(k)-qeddy(k)
c
          if ((n.le.n1.and.mod(n,nstop/10).eq.0).or.
     -        (n.ge.n1.and.mod(n,nstop/50).eq.0)) then
            if (k.eq.kmin) print *,'time:',real(n)/real(nstop),'n:',n
            print *,theta(k),d(k),1.e-6*qek(k),-1.e-6*qeddy(k),
     -            -1.e-6*qn(k),1.e-6*qu(k),1.e-6*qtot(k,1)

            if (k.eq.nz-1) print *,'heat flux:',heatuptake/a
            if (k.eq.nz-1) print *,' '
          endif
c
        enddo
c
c       heat content and uptake (including decomposition)
c         subscripts 1-3 refer to >10, 5-10, <5 degrees
c         (NB: temperature classes hard-wired here for now) 
c
          heatek=0.
          heatn=0.
          heatu=0.
          heateddy=0.
          heatuptake=0.
          heat=rho0*cp*dtheta*a*d(kmin)
          heat07=rho0*cp*dtheta*a*d(kmin)
          heat72=rho0*cp*dtheta*a*d(kmin)
c
          do k=kmin+1,nz-1
            heatek=heatek+rho0*cp*dtheta*qek(k)
            heatn=heatn+rho0*cp*dtheta*qn(k)
            heatu=heatu+rho0*cp*dtheta*qu(k)
            heateddy=heateddy+rho0*cp*dtheta*qeddy(k)
            heatuptake=heatuptake+rho0*cp*dtheta*qtot(k,1)
            heat=heat+rho0*cp*dtheta*a*d(k)
            heat07=heat07+rho0*cp*dtheta*a*min(d(k),700.)
            heat72=heat72+rho0*cp*dtheta*a*min(d(k),2000.)
          enddo
          heat25=heat-heat72
          heat72=heat72-heat07
c
c     step forward
        do k=kmin,nz-1
          d(k)=d(k)+ab(1)*qtot(k,1)+ab(2)*qtot(k,2)+ab(3)*qtot(k,3)
               if (k.eq.kmin) d(k)=max(d(k),hmin)
        enddo
c   
        if (n.ge.n1.and.mod(n-n1,(nstop-n1)/(nout-1)).eq.0) then
          print *, n
c
c         save data 
            ndump=ndump+1 
c            print *,n,n1,n-n1,nstop-n1,nout-1,(nstop-n1)/(nout-1)
c            print *,ndump,nout,n,nstop
            do k=0,nz
              if (k.lt.kmin.or.k.eq.nz) then
                qekout(k,ndump)=0. 
                qeddyout(k,ndump)=0. 
                qnout(k,ndump)=0. 
                quout(k,ndump)=0.
                qtotout(k,ndump)=0. 
                dout(k,ndump)=0.
              else 
                qekout(k,ndump)=qek(k)*1.e-6
                qeddyout(k,ndump)=qeddy(k)*1.e-6
                qnout(k,ndump)=qn(k)*1.e-6
                quout(k,ndump)=qu(k)*1.e-6
                qtotout(k,ndump)=qtot(k,1)*1.e-6
                dout(k,ndump)=d(k)       
              endif
              dout(nz,ndump)=d(nz)
              heatout(ndump)=heat
              heat07out(ndump)=heat07
              heat72out(ndump)=heat72
              heat25out(ndump)=heat25
              heatupout(ndump)=heatuptake
              heatekout(ndump)=heatek
              heatuout(ndump)=heatu
              heatnout(ndump)=heatn
              heateddyout(ndump)=heateddy
              thetasout(ndump)=thetas
            enddo
        endif
c
c       shuffle forcing	
        do k=kmin,nz-1
          qtot(k,3)=qtot(k,2)
          qtot(k,2)=qtot(k,1)
        end do
c
      end do	  
c
c     end of main loop
c
c     write data to files
c        
        open(unit=9,file='qek.dat',status='unknown') 
        do k=0,nz
         do n=1,ndump
          write(9,*) qekout(k,n)
         end do
        end do
        close(9)
c
        open(unit=9,file='qeddy.dat',status='unknown') 
        do k=0,nz
         do n=1,ndump
          write(9,*) qeddyout(k,n)
         end do
        end do
        close(9)
c
        open(unit=9,file='qn.dat',status='unknown') 
        do k=0,nz
         do n=1,ndump
          write(9,*) qnout(k,n)
         end do
        end do
        close(9)
c
        open(unit=9,file='qu.dat',status='unknown') 
        do k=0,nz
         do n=1,ndump
          write(9,*) quout(k,n)
         end do
        end do
        close(9)
c
        open(unit=9,file='qtot.dat',status='unknown') 
        do k=0,nz
         do n=1,ndump
          write(9,*) qtotout(k,n)
         end do
        end do
        close(9)
c
        open(unit=9,file='d.dat',status='unknown') 
        do k=0,nz
         do n=1,ndump
          write(9,*) dout(k,n)
         end do
        end do
        close(9)
c
        open(unit=9,file='heat.dat',status='unknown') 
        do n=1,ndump
          write(9,*) heatout(n)
        end do
        close(9)
c
        open(unit=9,file='heat07.dat',status='unknown')
        do n=1,ndump
          write(9,*) heat07out(n)
        end do
        close(9)
c
        open(unit=9,file='heat72.dat',status='unknown')
        do n=1,ndump
          write(9,*) heat72out(n)
        end do
        close(9)
c
        open(unit=9,file='heat25.dat',status='unknown')
        do n=1,ndump
          write(9,*) heat25out(n)
        end do
        close(9)
c
        heatupout(0)=heatupout(1)
        heatupout(nout+1)=heatupout(nout)
        open(unit=9,file='heatup.dat',status='unknown') 
        do n=1,ndump
          write(9,*) 0.25*(2.*heatupout(n)
     -                   +heatupout(n+1)+heatupout(n-1))/a
        end do
        close(9)
c
        heatekout(0)=heatekout(1)
        heatekout(nout+1)=heatekout(nout)
        open(unit=9,file='heatek.dat',status='unknown') 
        do n=1,ndump
          write(9,*) 0.25*(2.*heatekout(n)
     -                   +heatekout(n+1)+heatekout(n-1))/a
        end do
        close(9)
c
        heatuout(0)=heatuout(1)
        heatuout(nout+1)=heatuout(nout)
        open(unit=9,file='heatu.dat',status='unknown') 
        do n=1,ndump
          write(9,*) 0.25*(2.*heatuout(n)
     -                   +heatuout(n+1)+heatuout(n-1))/a
        end do
        close(9)
c
        heateddyout(0)=heateddyout(1)
        heateddyout(nout+1)=heateddyout(nout)
        open(unit=9,file='heateddy.dat',status='unknown') 
        do n=1,ndump
          write(9,*) 0.25*(2.*heateddyout(n)
     -                   +heateddyout(n+1)+heateddyout(n-1))/a
        end do
        close(9)
c
        heatnout(0)=heatnout(1)
        heatnout(nout+1)=heatnout(nout)
        open(unit=9,file='heatn.dat',status='unknown') 
        do n=1,ndump
          write(9,*) 0.25*(2.*heatnout(n)
     -                   +heatnout(n+1)+heatnout(n-1))/a
        end do
        close(9)
c
        open(unit=9,file='thetas.dat',status='unknown') 
        do n=1,ndump
          write(9,*) thetasout(n)
        end do
        close(9)
c
      stop 
      end

