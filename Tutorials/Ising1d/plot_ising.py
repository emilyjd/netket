import numpy as np
import matplotlib.pyplot as plt
import json

plt.ion()

#N=20
exact=-1.274549484318e+00*20

#N=80
# exact=-1.273321360724e+00*80
i=0
while(i<1):
#    plt.clf()
    plt.figure()
    plt.ylabel('Energy')
    plt.xlabel('Iteration #')

    iters=[]
    energy=[]
    sigma=[]
    evar=[]
    evarsig=[]

    data=json.load(open('rmsprop_test.log'))
    for iteration in data["Output"]:
        iters.append(iteration["Iteration"])
        energy.append(iteration["Energy"]["Mean"])
        sigma.append(iteration["Energy"]["Sigma"])
        evar.append(iteration["EnergyVariance"]["Mean"])
        evarsig.append(iteration["EnergyVariance"]["Sigma"])

    nres=len(iters)
    cut=60
    if(nres>cut):

        fitx=iters[-cut:-1]
        fity=energy[-cut:-1]
        z=np.polyfit(fitx,fity,deg=0)
        p = np.poly1d(z)

        plt.xlim([0,nres])#([nres-cut,nres])
        maxval=np.max(energy)#(energy[-cut:-1])
        plt.ylim([exact-(np.abs(exact)*0.01),maxval+np.abs(maxval)*0.01])
        error=(z[0]-exact)/-exact
        plt.gca().text(0.95, 0.8, 'Relative Error : '+"{:.2e}".format(error),
        verticalalignment='bottom', horizontalalignment='right',
        color='green', fontsize=15,transform=plt.gca().transAxes)

        plt.plot(fitx,p(fitx),'b')

    plt.errorbar(iters,energy,yerr=sigma,color='red')
    plt.axhline(y=exact, xmin=0, xmax=iters[-1], linewidth=2, color = 'k',label='Exact')


    plt.legend(frameon=False)
#    plt.pause(1)
    # plt.draw()
    plt.show()
    i+=1

#plt.ioff()
#plt.show()
