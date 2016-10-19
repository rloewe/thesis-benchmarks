

%  Stommel's 1961 model of convection in coupled boxes. 
%
%  Compute the temperature and salinity difference, here represented
%  by y and x, between two well-stirred boxes that are both
%  in contact with a reservoir that has (y, x) = (1, 1).  Both
%  boxes conduct y and x at rates 1 and delta (delta < 1), and the 
%  density difference between the boxes is given by d = -y + R*x, where
%  we are assuming R > 1, typically.  There can be advection (or 
%  flushing) between the boxes at a rate d*lambda, where the 
%  flushing does not depend upon the sign of the density anomaly.  
%  
%  This code also has 1) flushing with inertia (not especially
%  interesting), 2) a flickering or random temperature
%  perturbation (slightly interesting), and 3) an oscillating 
%  reservoir temperature (hoping for but not yet finding chaos;
%  may be present in parameter regimes not checked).  
%
%  This code has not been fully tested, but seems to reproduce 
%  S61's two cases fairly well when all three of the 'new' things 
%  are turned off, of course.
%
%  Written by Jim Price, April 28, 1999.  Public domain.
%

set(0,'DefaultLineLineWidth',1.2)
set(0,'DefaultTextFontSize',14)
set(0,'DefaultAxesLineWidth',1.6)
set(0,'DefaultAxesFontSize',14)
 
clear; format compact; 
nn = 0;

%  set the model parameters (follows S61, aside from 'new')

R = 2.0;       %  abs of the ratio of the expansion coefficients, x/y
delta = 1/6;   %  conduction rate of salinity wrt temperature
lambda = 0.2;  %  inverse non-d flushing rate; = inf for no flushing
q = 0.;        %  initial flushing rate (0 to 1) 'new'
qdelta = 100.; %  time constant (inertia) for flushing; 'new'
               %    set = 1/dtau for equilibrium flushing as in S61
               %    set = 0.2 for slowly responding flushing 
yres = 1.;     %  steady reservoir y, = 1 for S61 case  'new'
resosc = 0.;   %  amplitude of reservoir y oscillation  'new'            
dtau = 0.01;   %  the time step of non-d time; 0.01 seems OK
nstep = 1500;  %  number of time steps; 1500 is usually 
               %    enough to insure convergence to a steady state
%                 
%  This model version is set up to do integration over a range
%  of initial T,S or y,x from 0 to 1.  The increment of T and S 
%  are set by delT and delS = 1/ni,  where ni is the number of
%  integrations (n1 = 1 to 20 is reasonable).

yres0 = yres;

ni = 6;
delT = 1/ni;   %  make ni = 1 or 2 to reduce the number of
               %     integrations to be done
delS = delT;

for n1=0:delT:1  %  n1 and n2 are used to set the initial T,S
    for n2=0:delS:1

        if n1==0 | n1==1 | n2==0 | n2==1  %  skip all non-boundary initial points

            x(1) = n1;   %  set the initial temperature
            y(1) = n2;   %  set the initial salinity

            for m = 2:nstep

                tau(m) = m*dtau;  %  the non-d time

                %  evaluate the reservoir temperature (y); note that
                %   this temperature is steady if resosc = 0. (the S61 case)

                yres = yres0 + resosc*sin(tau(m)*pi);  

                % the first part of a second order R-K method; time step forward
                %    by half a time step to make a first guess at the new times

                dr = abs(R*x(m-1)   - y(m-1));        %  the density anomaly

                qequil = dr/lambda;                   %  the equilibrium flushing

                yh = y(m-1) + dtau*(yres - y(m-1))/2 - ...
                   dtau*y(m-1)*q/2;                   %  time step the temperature 

                xh = x(m-1) + dtau*delta*(1 - x(m-1))/2 - ...
                   dtau*x(m-1)*q/2;                   %  time step the salinity

                qh = q + dtau*qdelta*(qequil - q)/2;  %  time step the flushing

                % the second part; use the half time step values to make a full step

                dr = abs(R*xh  - yh);

                qequil = dr/lambda;

                y(m) = y(m-1) + dtau*(yres - yh) - ...
                   dtau*qh*yh;  

                x(m) = x(m-1) + dtau*delta*(1 - xh) - ...
                   dtau*qh*xh; 

                q = q + dtau*qdelta*(qequil - qh);

                % now add on a flickering temperature if you want to (or comment out)

                % tflickamp = 0.01;   %  set the amplitude here
                % tflick = tflickamp*unifrnd(-1., 1.);
                % y(m) = y(m) + tflick;

            end  %  end of time step loop

            d = R*x  - y;   %  evaluate the density

            if nn == 0;
                nn = 1;


                %  make a time series plot of the first case only 

                figure(1)
                clf reset

                subplot(2,1,1)
                plot(tau, x,'-', tau, y,'--')
                legend('salinity', 'temperature')
                ylabel('T, S diff, non-d')
                title('Experiment 1,1')

                subplot(2,1,2)
                plot(tau, d)
                xlabel('time, non-d')
                ylabel('density diff')


                %  contour the density (or flushing rate), and add the (T,S)
                %  trajectories on top of the contours

                ym = 0:0.1:1;
                xm = 0:0.1:1;
                for k1 = 1:11
                for k2 = 1:11  
                  dm(k1,k2) = (1/lambda)*(R*xm(k2)  - ym(k1));
                end
                end

                figure(2)
                clf reset

                dc = -10:2:20;
                c = contour(xm, ym, dm, dc,'k');
                clabel(c);
                xlabel('salinity diff, non-d')
                ylabel('temp diff, non-d')

                hold on

            end   %  the if on nn = 0 

            [m1 m2] = size(x);

            %  plot the individual trajectories. 
            %  color code according to which equilibrium point
            %  the trajectory ends up on. this will likely have be 
            %  reset if the model parameters (R, delta, lambda) are changed.

            if d(m2) >= 0
            plot(x, y, 'r--')
            plot(x(m2), y(m2), '*r')

            else 
            plot(x, y, 'g')
            plot(x(m2), y(m2), '*g')
            end

            clear x y d

        end %  or to gate out non-boundary values of n1, n2
    end %  loop on n1
end %  loop on n2


%  make some plots to show where roots (equil. points) are

for k=1:60
f(k) = (k-30)*0.1;
lhs(k) = lambda*f(k);
rhs(k) = (R/(1 + abs(f(k))/delta)) - 1/(1 + abs(f(k)));
end

figure(3)
clf reset
plot(f, rhs, f, lhs)
xlabel('f, flow rate')
ylabel('lhs(f), rhs(f)')
title('roots of S61 model')
grid

%   end of the script
