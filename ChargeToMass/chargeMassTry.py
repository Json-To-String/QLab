import matplotlib.pyplot as plt
import numpy as np
import os
from scipy.optimize import curve_fit

def run1():
    I = .688
    dI = 0.1e-3
    V = 2.98e3
    dV = .01e3
    xBoth = np.array([9,8,7,6,5,4,3,2])
    dxPos = np.array([.01,.01,.01,.01,.01,.01,.01,.01])
    yPos = np.array([2,1.6,1.2,.9,.6,.31,.22,.04])
    dyPos = np.array([.01,.01,.05,.05,.05,.05,.05,.01])
    dxNeg = np.array([.01,.01,.01,.01,.01,.01,.01,.01])
    yNeg = np.array([-2.25,-1.8,-1.4,-1.1,-0.8,-0.6,-0.4,-.21])
    dyNeg = np.array([0.025, 0.025, 0.025, 0.025, 0.025, 0.025, 0.025, 0.025])

    return(I,dI,V,dV,xBoth,dxPos,yPos,dyPos,dxNeg,yNeg,dyNeg)

def run2():
    I = .465
    dI = 0.1e-3
    V = 2.97e3
    dV = .01e3
    xBoth = np.array([9,8,7,6,5,4,3,2])
    dxPos = np.array([.01,.01,.01,.01,.01,.01,.01,.01])
    yPos = np.array([1.3,1.0,.8,.6,.4,.22,.18,.02])
    dyPos = np.array([0.025, 0.025, 0.025, 0.025, 0.025, 0.025, 0.025, 0.025])
    dxNeg  = np.array([.01,.01,.01,.01,.01,.01,.01,.01])
    yNeg = np.array([-1.5,-1.2,-1.0,-0.8,-.58,-.4,-.25,-.2])
    dyNeg = np.array([0.025, 0.025, 0.025, 0.025, 0.025, 0.025, 0.025, 0.025])
    return(I,dI,V,dV,xBoth,dxPos,yPos,dyPos,dxNeg,yNeg,dyNeg)

def run3():
    I =.545
    dI = 0.1e-3
    V = 2.97e3
    dV = .01e3
    xBoth = np.array([9,8,7,6,5,4,3,2])
    dxPos = np.array([.01,.01,.01,.01,.01,.01,.01,.01])
    yPos = np.array([1.5,1.2,.9,.7,.5,.3,.18,.02])
    dyPos = np.array([0.025, 0.025, 0.025, 0.025, 0.025, 0.025, 0.025, 0.025])
    dxNeg = np.array([.01,.01,.01,.01,.01,.01,.01,.01])
    yNeg = np.array([-1.7,-1.4,-1.1,-.9,-.63,-.45,-.38,-.2])
    dyNeg = np.array([0.03, 0.03, 0.03, 0.03, 0.03, 0.03, 0.03, 0.03])

    return(I,dI,V,dV,xBoth,dxPos,yPos,dyPos,dxNeg,yNeg,dyNeg)

def run4():
    I = .307
    dI = 0.1e-3
    V = 2.97e3
    dV = .01e3
    xBoth = np.array([9,8,7,6,5,4,3,2])
    dxPos = np.array([.01,.01,.01,.01,.01,.01,.01,.01])
    yPos = np.array([.8,.6,.5,.38,.22,.08,.03,.001])
    dyPos = np.array([0.03, 0.03, 0.03, 0.03, 0.03, 0.03, 0.03, 0.03])
    dxNeg  = np.array([.01,.01,.01,.01,.01,.01,.01,.01])
    yNeg = np.array([-1.0,-.83,-.7,-.5,-.4,-.3,-.2,-.18])
    dyNeg = np.array([0.03, 0.03, 0.03, 0.03, 0.03, 0.03, 0.03, 0.03])

    return(I,dI,V,dV,xBoth,dxPos,yPos,dyPos,dxNeg,yNeg,dyNeg)

def run5():
    I = .455
    dI = .01e-3
    V = 2.97e3
    dV = .01e3
    xBoth = np.array([9,8,7,6,5,4,3,2])
    dxPos = np.array([.01,.01,.01,.01,.01,.01,.01,.01])
    yPos = np.array([1.2,.99,.75,.56,.4,.21,.06,.02])
    dyPos = np.array([0.025, 0.025, 0.025, 0.025, 0.025, 0.025, 0.025,
                    0.025])
    dxNeg = np.array([.01,.01,.01,.01,.01,.01,.01,.01])
    yNeg = np.array([-1.45,-1.2,-.95,-.7,-.58,-.4,-.22,-.2])
    dyNeg = np.array([0.03, 0.03, 0.03, 0.03, 0.03, 0.03, 0.03, 0.03])

    return(I,dI,V,dV,xBoth,dxPos,yPos,dyPos,dxNeg,yNeg,dyNeg)


def plot_Calc(I,dI,V,dV,xBoth,dxPos,dxNeg,yAvg,dyAvg,posB,negB,dB,runNum):
    print('\n\tRun number',runNum+1)

    # Curves to fit for both pos and neg
    def funcPos(x,r,x0,y0):
         return (y0 - np.sqrt(r**2 - (x-x0)**2))

    def funcNeg(x,r,x0,y0):
         return (y0 + np.sqrt(r**2 - (x-x0)**2))

    # Give it an interpolated interval
    xtheory = np.linspace(2,9)

    # Initial guess for curve_fit func
    p0 = [15/100,0,15/100]
    popt,pcov = curve_fit(funcPos,xBoth,yAvg,p0,sigma=dxPos)

    # Uncertainty in r comes from diagonal of cov matrix
    perr = np.sqrt(np.diag(pcov))
    # perrNeg = np.sqrt(np.diag(pcovNeg))

    r = popt[0]
    dr = perr[0]

    cM = 2*V / ((posB**2)*(r**2))
    dCMdV = 2/((posB**2)*(r**2))
    dCMdB = -4*V / ((r**2)*(posB**3))
    dCMdr = -4*V / ((posB**2)*(r**3))
    dCM = np.sqrt(
               (dCMdV**2)*dV**2
             + (dCMdB**2)*dB**2
             + (dCMdr**2)*dr**2)


    # plt.figure()
    # plt.errorbar(xBoth,yAvg,yerr=dyPos,fmt='k.')
    # plt.plot(xtheory,funcPos(xtheory,15/100,0,15/100),'b--')
    # plt.title('Guess for Theory Curve')
    # plt.xlabel('x[m]')
    # plt.ylabel('y[m]')
    # plt.show()
    #
    # plt.figure()
    # plt.errorbar(xBoth,yAvg,xerr=dxPos,yerr=dyAvg,fmt='k.',capsize =4)
    # plt.plot(xBoth,funcPos(xBoth,popt[0],popt[1],popt[2]),'k--')
    # plt.xlabel('x[m]')
    # plt.ylabel('y[m]')
    # plt.savefig('chargeMassRun{}.png'.format(runNum+1))
    # plt.show()

    return(cM,dCM,popt)


def methods(vals,cmUnc,runs):

    meanCM = np.mean(vals)
    samples = np.arange(runs)
    sdom = np.std(vals)/len(vals)

    plt.figure()
    plt.plot(samples,vals,'k.')
    plt.errorbar(samples,vals,yerr=cmUnc,fmt = 'k.',capsize = 3)
    plt.axhline(meanCM,color = 'b',linestyle = '--')
    plt.axhline(meanCM+sdom,color = 'b')
    plt.axhline(meanCM-sdom,color = 'b')
    plt.ylabel('Charge to Mass ratio $\\frac{C}{kg}$')
    plt.xlabel('Sample Number')
    plt.savefig('cmratioVsSample.png')
    # plt.show()

def main():

    # Static Variables
    muNot = 1.25663706e-6 # N/A^2
    N = 131 # turns
    dN = .1
    D = 20.8e-2 # m
    dD = 0.001 # m

    # List of runs where elements are function handles
    runList = [run1,run2,run3,run4,run5]
    runs = len(runList)

    # Preallocate empty list to later compare all c/m ratios to each other
    cmRatioList = []
    cmUncList = []

    # Loop over number of runs
    plt.figure()
    for ind in np.arange(runs):

        # Call func number at index number to call the correct number
        # Vars are redefined every pass through the loop
        I,dI,V,dV,xBoth,dxPos,yPos,dyPos,dxNeg,yNeg,dyNeg = runList[ind]()

        # Convert to m
        xBoth = xBoth/100
        dxPos = dxPos/100
        yPos = yPos/100
        dyPos = dyPos/100
        dxNeg = dxNeg/100
        yNeg = yNeg/100
        dyNeg = dyNeg/100

        yAvg = []
        for pos,neg in zip(yPos,yNeg):
            yAvg.append((pos + np.abs(neg))/2)
        dyAvg = np.std(yAvg) / np.sqrt(len(yPos)+len(yNeg))

        # B field takes in data from run
        posB = (16*muNot*N*I)/(np.sqrt(125)*D) # T
        negB = (16*muNot*N*-I)/(np.sqrt(125)*D) # T
        dB = posB*np.sqrt(((dN/N)**2) + ((dI/I)**2) + ((dD/D)**2))

        cM,dCM,popt = plot_Calc(I,dI,V,dV,xBoth,dxPos,dxNeg,yAvg,dyAvg,posB,
                                negB,dB,ind)

        cmRatioList.append(cM)
        cmUncList.append(dCM)

        def funcPos(x,r,x0,y0):
             return (y0 - np.sqrt(r**2 - (x-x0)**2))


        plt.errorbar(xBoth,yAvg,xerr=dxPos,yerr=dyAvg,fmt='k.',capsize =4)
        plt.plot(xBoth,funcPos(xBoth,popt[0],popt[1],popt[2]),'k--')
        plt.xlabel('x[m]')
        plt.ylabel('y[m]')
        # plt.savefig('chargeMassRun{}.png'.format(ind+1))
    plt.show()

    methods(cmRatioList,cmUncList,runs)

main()
