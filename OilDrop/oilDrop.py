import numpy as np
import os
import matplotlib.pyplot as plt
from scipy import interpolate
import qLabMods
from pint import UnitRegistry

ureg = UnitRegistry()

def drop1():
    d = 2 # boxes
    d = d*0.1e-3 # m
    dist = 7.6e-3 # m
    downT = np.array([9.75,7.22,7.49,9.11,8.67,9.73,8.52,7.09,8.64])
    upT = np.array([2.67,2.51,3.16,2.81,2.69,2.77,2.87,2.60,2.64])
    dTime = 0.05 # sec
    V = 302 # V
    dV = 2
    tR = 2.196e6 # Ohms
    dtR = 0.002e6

    return dist,d,downT,upT,dTime,V,dV,tR,dtR

def drop2():
    d = 2
    d = d*0.1e-3
    dist = 7.6e-3
    downT = np.array([20.55,26.38,24.6,33.33,22.47,21.27,20.34,21.89])
    upT = np.array([3.11,2.85,2.94,2.59,2.67,1.29,2.92,2.93])
    dTime = 0.05
    V = 300
    dV = 2
    tR = 2.064e6
    dtR = 0.002e6

    return dist,d,downT,upT,dTime,V,dV,tR,dtR

def drop3():
    d = 2
    d = d*0.1e-3
    dist = 7.6e-3
    downT = np.array([47.32,43.41])
    upT = np.array([0.73,0.95])
    dTime = 0.05
    V = 300
    dV = 2
    tR = 2.0206e6
    dtR = 0.002e6

    return dist,d,downT,upT,dTime,V,dV,tR,dtR

def drop4():
    d = 2
    d = d*0.1e-3
    dist = 7.6e-3
    downT = np.array([8.89,11.9,11.86,9.52,11.75,11.53,10.71,11.68])
    upT = np.array([1.52,1.37,0.98,0.94,1.03,0.87,1.1,1.13])
    dTime = 0.05
    V = 299
    dV = 2
    tR = 2.012e6
    dtR = 0.002e6

    return dist,d,downT,upT,dTime,V,dV,tR,dtR

def drop5():
    d = 2
    d = d*0.1e-3
    dist = 7.6e-3
    downT = np.array([20.33,19.69,18.72,26.2,26.08,22.92,20.87,15.89])
    upT = np.array([0.9,1.0,1.04,0.92,0.85,0.96,0.86,0.94])
    dTime = 0.05
    V = 299
    dV = 2
    tR = 1.99e6
    dtR = 0.002e6

    return dist,d,downT,upT,dTime,V,dV,tR,dtR

def drop6():
    d = 2
    d = d*0.1e-3
    dist = 7.6e-3
    downT = np.array([4.46,5.18,5.66,4.48,5.43,5.09,4.64,5.31])
    upT = np.array([3.75,4.06,3.66,3.27,3.39,3.52,3.7,3.12])
    dTime = 0.05
    V = 299
    dV = 2
    tR = 1.985e6
    dtR = 0.002e6

    return dist,d,downT,upT,dTime,V,dV,tR,dtR

def drop7():
    d = 2
    d = d*0.1e-3
    dist = 7.6e-3
    downT = np.array([6.43,5.58,6.41,5.88])
    upT = np.array([2.55,2.9,2.84,2.58])
    dTime = 0.05
    V = 299
    dV = 2
    tR = 1.98e6
    dtR = 0.002e6
    return dist,d,downT,upT,dTime,V,dV,tR,dtR

def drop8():
    d = 2
    d = d*0.1e-3
    dist = 7.6e-3
    downT = np.array([7.2,7.32,8.86,9.81,6.61,9.85,6.87,8.1])
    upT = np.array([1.15,2.01,1.99,1.52,1.73,1.67,1.75,2.01])
    dTime = 0.05
    V = 299
    dV = 2
    tR = 1.977e6
    dtR = 0.002e6
    return dist,d,downT,upT,dTime,V,dV,tR,dtR

def drop9():
    d = 2
    d = d*0.1e-3
    dist = 7.6e-3
    downT = np.array([13.04,18.7,17.28,16.5,17.46,17.97,16.89,13.9])
    upT = np.array([1.61,0.46,0.65,0.82,0.5,0.56,0.54,0.46])
    dTime = 0.05
    V = 299
    dV = 2
    tR = 1.975e6
    dtR = 0.002e6
    return dist,d,downT,upT,dTime,V,dV,tR,dtR


def calculate_velocities(d,downT,upT,dTime):
    vDown = d/downT
    vDownBar = np.mean(vDown)
    dvDown = vDown*(dTime/downT)
    vUp = d/upT
    vUpBar = np.mean(vUp)
    dvUp = vUp*(dTime/upT)
    # print(vFall,vRise)
    # print(vFallBar,vRiseBar)

    return(vDown,vDownBar,dvDown,vUp,vUpBar,dvUp)

def interp_T(tR,dtR):

    C = np.array([17,18,19,20,21,22,23,24,25,26,27,28,29])
    R = np.array([2.526,2.446,2.371,2.3,2.233,2.169,2.11,2.025,2.00,1.950,1.902,
                1.857,1.815])*1e6

    f = interpolate.interp1d(R,C)
    temp = f(tR)

    return temp

def calculate_q(dist,vDown,vDownBar,dvDown,vUp,vUpBar,dvUp,volts,dTime,T):
    g = 9.8
    rho = 885 # kg/m^3
    # T = 22.5 # Celsius, change later
    dTemp = .3
    b = 0.82e-7

    E = volts/dist
    e = 1.6021766208e-19

    eta = (1.800 + ((T-15)/209))*10**-5
    dEtadT = (10**-5)/209
    etaUnc = dEtadT*dTemp

    r = np.sqrt((9*eta*vDownBar)/(2*g*rho))
    drdEta = (9*vDownBar)*(etaUnc**2)/(8*g*rho*eta)
    drdvDown = (9*eta)*(dvDown**2)/(8*g*rho*vDownBar)
    rUnc = ((drdEta**2)*(etaUnc**2) + (drdvDown**2)*(dvDown**2))**.5

    m = (4/3)*(np.pi*r**3*rho)
    mUnc = np.abs(3*m*(rUnc/r))

    gamma = (1/(1+(b/r)))**(-3/2)
    dGammadr = (3/2)*(-b/(r**2))*(1 + (b/r))**.5
    gammaUnc = np.abs(dGammadr*rUnc)

    qU = gamma*(m*9.8/E)*(vUp/vDown + 1)
    qU = np.mean(qU)
    dqdGamma = (((m*g)/E)*(vUpBar/vDownBar + 1)*gammaUnc)**2
    dqdm = (((gamma*g)/E)*(vUpBar/vDownBar + 1)*mUnc)**2
    dqdvUp = (((gamma*m*g)/(E*vDownBar))*dvUp)**2
    dqdvDown = (((gamma*m*g*vUp)/(E*(vDown**2)))*dvDown)**2
    dq = np.sqrt(dqdGamma + dqdm + dqdvUp + dqdvDown)

    ratio = qU/e
    ratioUnc = np.mean(dq)*(1/e)

    # print('\nq',qU)
    # print('dq/e:',ratioUnc)
    # print('q/e',ratio)

    return qU,ratio,dq,ratioUnc,r,rUnc,eta,etaUnc,m,mUnc,gamma,gammaUnc

def plots(ratio,uncs):

    samples = np.arange(len(ratio))
    plt.figure()
    plt.errorbar(samples,ratio,yerr=uncs,fmt = 'b.',capsize = 3)
    plt.xlabel('Sample Number',fontsize = 12)
    plt.ylabel('Integer Multiple of Electron Charge [# of e]',fontsize = 12)
    plt.grid()
    plt.savefig('oilDrop.pdf')
    # plt.show()


def main():
    dropList = [drop1,drop2,drop3,drop4,drop5,drop6,drop7,drop8,drop9]
    drops = len(dropList)

    etaList = []
    etaUncs = []
    mList = []
    mUncs = []
    gammaList = []
    gammaUncs = []
    rList = []
    rUncs = []
    qList = []
    dqList = []
    ratioList = []
    ratioUncs = []
    for ind in np.arange(drops):
        # print('\nDrop {}'.format(ind+1))
        (dist,d,downT,upT,dTime,volts,dV,tR,dtR) = dropList[ind]()

        temp = interp_T(tR,dtR)

        (vDown,vDownBar,dvDown,vUp,vUpBar,dvUp) = calculate_velocities(d,downT,
                                                                    upT,dTime)
        (qU,ratio,dq,ratioUnc,r,rUnc,
        eta,etaUnc,m,mUnc,gamma,gammaUnc) = calculate_q(dist,vDown,vDownBar,
                                                        dvDown,vUp,vUpBar,dvUp,
                                                        volts,dTime,temp)

        etaList.append(eta)
        etaUncs.append(etaUnc)
        mList.append(m)
        mUncs.append(np.mean(mUnc))
        gammaList.append(gamma)
        gammaUncs.append(gammaUnc)
        rList.append(r)
        rUncs.append(np.mean(rUnc))
        qList.append(qU)
        dqList.append(np.mean(dq))
        ratioList.append(ratio)
        ratioUncs.append(ratioUnc)

        dTime = dTime*np.ones(len(upT))
        varray = [volts]
        dVarray = [dV]
        tRarray = [tR]
        dtRarray = [dtR]
        # qArray = [qU]
        # dqArray = [np.mean(dq)]
        # ratioArray = [ratio]
        # ratioUncArray = [ratioUnc]
        # etaArray = [eta]
        # etaUncArray = [etaUnc]
        # mArray = [m]
        # mUncArray = [mUnc]
        # rArray = [r]
        # rUncArray = [np.mean(rUnc)]
        # gammaArray = [gamma]
        # gammaUncArray = [np.mean(gammaUnc)]

        for item in upT:
            varray.append([])
            dVarray.append([])
            tRarray.append([])
            dtRarray.append([])
            # qArray.append([])
            # dqArray.append([])
            # ratioArray.append([])
            # ratioUncArray.append([])
            # etaArray.append([])
            # etaUncArray.append([])
            # mArray.append([])
            # mUncArray.append([])
            # rArray.append([])
            # rUncArray.append([])
            # gammaArray.append([])
            # gammaUncArray.append([])

        # print(ratioArray)
        # print(ratioUncArray)
        # print('dq:',dqArray)
        # print('dvup:',dvUpArray)
        # print('dvdown:',dvDownArray)
        # print('lencheck',len(vUp),len(dvUpArray))
        # print('eta',etaArray)
        # print('etaunc',etaUncArray)
        # print('r',rArray)
        # print('runc',rUncArray)
        # print('gamma',gammaArray)
        # print(gammaUncArray)
        # print('q',qArray)
        # print('dq',dq)
        # print('ratio',ratioArray)
        # print(ratioUncArray)
        # print('lencheck',len(qArray)==len(dqArray))
        # print('lencheck',len(ratioArray)==len(ratioUncArray))
        table1 = qLabMods.LaTeXTable(['$t_{up}$','$t_{down}$','$Voltage$','Thermistor Resistance'],
                                    [(upT,dTime,ureg.seconds),
                                    (downT,dTime,ureg.seconds),
                                    (varray,dVarray,ureg.volts),
                                    (tRarray,dtRarray,ureg.ohms)],
                                    'Drop {}: Raw Data'.format(ind+1),
                                    'DataTable',
                                    '!htp')

        # print(table1.allTogetherRaw())


    # table2 = qLabMods.LaTeXTable(['$eta$','$\gamma$','$m$','$radius$','$q_{u}$','$q/e ratio$'],
    #                              [(etaList,etaUnc,ureg.seconds),
    #                              (gammaList,gammaUncs,ureg.seconds),
    #                              (mList,mUncs,ureg.kilogram),
    #                              (rList,rUncs,ureg.meter),
    #                              (qList,dqList,ureg.meter/ureg.seconds),
    #                              (ratioList,ratioUncs,ureg.seconds)],
    #                              'All Drops: Derived Data',
    #                              'DataTable',
    #                              '!htp')
    plots(ratioList,ratioUncs)
    # print(table2.allTogetherSci())

main()
