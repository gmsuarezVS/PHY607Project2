import os
import matplotlib
import mpmath as mp
import numpy as np

matplotlib.use("Agg")
import matplotlib.pyplot as plt

#necessarily import the stuff from pt 1, and pt 2
#something like from "filename part1" import Particles

# the particles from part 1, and equations, rk4 from part 2

#1) define function that takes in results from part 2,
#and computes "observables" (classical ones like temperature, energy, distribution, etc.)
def potentialE(state,length,alpha)
    V=0
    z=state[:,0]+1j*state[:,1]
    for i in range(len(state)):
        for j in range(i+1,len(state)):
            #we use mp since high precision will be pretty important (rip computational time)
            V=V+alpha*mp.log(abs(mp.sin(mp.pi*(z[i]-mp.conj(z[j]))/length)/mp.sin(mp.pi*(z[i]-z[j])/length)))
    return float(V)

def totalE(state,length,alpha):
    kineticE=np.sum((state[:,2]**2+state[:,3]**2)/(2*state[:,1]**2))
    return kineticE+potentialE

def themeasure(y0,lenght,alpha,result):
    #we compute some stoof (finally)
    states=result["states"]
    E=np.empty(len(states))
    Px=np.empty(len(states))
    for i in range(len(states)):
        state=states[i]
        E[i]=totalE(state,length,alpha)
        Px[i]=np.sum(state[:,2]/state[:,1]**2)
    result["energy"]=E
    result["momentum"]=Px
    result["depth"]=np.concatenate(np.log(states[:,:,1]/y0))
    result["kinetic"]=np.concatenate(np.sum((state[:,:,2]**2+states[:,:,3]**2)/(2*state[:,:,1]**2)))
    return result

#2) plot all the heights of the particles at the end of the simulations

#3) compute collision angles, depth distribution, particle length between hits and demonstrate

#4) clean up plots/present results....

#5) compute for diffusion probability=0, the conserved energy, other deterministic features

#6) for non zero probability, since it's not deterministic (well, it is, but cause randomness is determinstic on python)
#we compute momenta over time (no longer conserved in x), average energy is temp, etc.

#^^the above are "just" def to be applied in the master script so not too worried