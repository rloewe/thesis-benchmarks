import numpy
import scipy.integrate
import matplotlib.pyplot as plt
numpy.set_printoptions(threshold=numpy.nan)

def q12(X, p):
    DeltaT = p["Temperature"][1] - p["Temperature"][0]
    DeltaS = X[1] - X[0]
    flow = -p["k"]*(p["alpha"]*DeltaT - p["beta"]*DeltaS)
    return flow

def rhs(time, X, p):
    deltaS12 = X[1] - X[0]
    Fs = FW(time, p)
    r = Fs *p["area"]
    r[0] += abs(q12(X,p))*deltaS12
    r[1] -= abs(q12(X,p))*deltaS12
    r = r / p["volume"]

    if p["optimization"]:
        r = r*1.e12

    return r

def FW(time, p):
    ramp = 1.
    amp = p["FW_amp0"]

    Fs = numpy.array([
        ramp * p["S0"] * amp * p["meter"] / p["year"],
        -ramp * p["S0"] * amp * p["meter"] / p["year"]
        ])

    return Fs

def set_parameters():
    p = {}
    p["meter"] = 1.
    p["sec"] = 1.
    p["km"] = 1000.* p["meter"]
    p["area"] = (4000. * p["km"]) ** 2
    p["depth"] = 4.*p["km"]
    p["volume"] = p["area"] * p["depth"]
    p["year"] = 365*24*3600.
    p["time_max"] = 500000. * p["year"]
    p["alpha"] = 4.e-4 # K^-1
    p["beta"] = 4.e-3  # ppt^-1
    p["Sverdrup"] = 10 ** 6 * p["meter"] ** 3 / p["sec"]
    p["k"] = 10 * p["Sverdrup"] / 0.005 # Sv/(delta sigma)
    p["S0"] = 35. # Salinity reference
    p["FW_amp0"] = 0.2 # meter/year for calculating steady state that
                       # initializes the hysteresis.
    p["Nrand"] = 50
    p["optimization"] = False
    return p

p = set_parameters()

# find steady state solutions
tspan = numpy.arange(0,1550,50) * p["year"]
for ii in xrange(3):
    if ii == 0:
        X0 = numpy.array([36, 35])
        p["Temperature"] = [25, 5]
    elif ii == 1:
        X0 = numpy.array([36, 34.255])
        p["Temperature"] = [25, 5]
    else:
        X0 = numpy.array([36, 34.5])
        p["Temperature"] = [10, 5]
    ode = scipy.integrate.ode((lambda t, X: rhs(t, X, p))).set_integrator('dopri5', rtol=1.e-7,atol=1.e-8)
    ode.set_initial_value(X0, tspan[0])
    X = numpy.array([X0])
    for t in tspan[1:]:
        res = ode.integrate(t)
        X = numpy.append(X, numpy.array([res]), axis = 0)

    plt.figure()
    plt.subplot(2,1,1)
    line1 = plt.plot(tspan/p["year"], X[:, 0], 'r')
    line2 = plt.plot(tspan/p["year"], X[:, 1], 'b')
    plt.ylabel('S')
    plt.legend((line1[0], line2[0]), ('S box 1','S box 2'))
    if ii == 0:
        plt.title('Thermally driven solution')
    elif ii == 1:
        plt.title('Thermally driven solution starting near the unstable state')
    else:
        plt.title('Salinity driven solution')

    # Initial DeltaS
    print X[0, 1] - X[0, 0]

    # Final DeltaS
    print X[-1, 1] - X[-1, 0]

    plt.subplot(2,1,2)
    q12plot = numpy.zeros(len(tspan))
    for i in xrange(len(tspan)):
        q12plot[i] = q12(X[i, :], p)/p["Sverdrup"]

    plt.plot(tspan/p["year"], q12plot, 'r')
    plt.xlabel('time (yr)')
    plt.ylabel('Transport [Sv]')

    FW_tmp = p["S0"]*p["FW_amp0"]*p["meter"]/p["year"]
    # Thermally Driven Flow
    SteadyState1=(0.5*p["alpha"]*(p["Temperature"][1]-p["Temperature"][0])+0.5*numpy.sqrt((p["alpha"]*(p["Temperature"][1]-p["Temperature"][0])) ** 2-(4*p["beta"]/p["k"])*FW_tmp*p["area"]))/p["beta"]
    #pltfigure(1); plot(FW_amp0,SteadyState1/p["beta"],'b+');
    SteadyState2=(0.5*p["alpha"]*(p["Temperature"][1]-p["Temperature"][0])- 0.5*numpy.sqrt((p["alpha"]*(p["Temperature"][1]-p["Temperature"][0]))**2-(4*p["beta"]/p["k"])*FW_tmp*p["area"]))/p["beta"]
    #plt.plot(FW_amp0,SteadyState2/p["beta"],'ro');
    # Salt driven flow
    SteadyState3=(0.5*p["alpha"]*(p["Temperature"][1]-p["Temperature"][0])+ 0.5*numpy.sqrt((p["alpha"]*(p["Temperature"][1]-p["Temperature"][0]))**2+(4*p["beta"]/p["k"])*FW_tmp*p["area"]))/p["beta"]
    #plt.plot(FW_amp0,SteadyState3/p["beta"],'k*');
    SteadyState4=(0.5*p["alpha"]*(p["Temperature"][1]-p["Temperature"][0])- 0.5*numpy.sqrt((p["alpha"]*(p["Temperature"][1]-p["Temperature"][0]))**2+(4*p["beta"]/p["k"])*FW_tmp*p["area"]))/p["beta"]
    #ply.plot(FW_amp0,SteadyState4/p["beta"],'g.');
    print "SteadyState1", SteadyState1
    print "SteadyState2", SteadyState2
    print "SteadyState3", SteadyState3
    print "SteadyState4", SteadyState4

plt.show()
