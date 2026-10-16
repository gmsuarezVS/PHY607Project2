#very rough outline for part 3 (no point to flesh out atm anyways)
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
#potential energy between each particle
def potentialE(state,length,alpha):
    V=0
    z=state[:,0]+1j*state[:,1]
    for i in range(len(state)):
        for j in range(i+1,len(state)):
            #we use mp since high precision will be pretty important (rip computational time)
            V=V+alpha*mp.log(abs(mp.sin(mp.pi*(z[i]-mp.conj(z[j]))/length)/mp.sin(mp.pi*(z[i]-z[j])/length)))
    return float(V)

#total energy-which needs be separate, since potential is their interaction term
def totalE(state,length,alpha):
    kineticE=np.sum((state[:,2]**2+state[:,3]**2)/(2*state[:,1]**2))
    return kineticE+potentialE

#how we "know" our code is "right" (quotes because not really in our case, but good enough)
def themeasure(y0,lenght,alpha,result):
    #we compute some stoof (finally)
    states=result["states"]
    E=np.empty(len(states))
    Px=np.empty(len(states))
    for i in range(len(states)):
        state=states[i]
        E[i]=totalE(state,length,alpha)
        #note here state[:,1] is y, state[:,2] is vx. (state[:,0] is x) ...
        Px[i]=np.sum(state[:,2]/state[:,1]**2)
    result["energy"]=E
    result["momentum"]=Px
    #here y is states[:,:,1] from part 2
    result["depth"]=np.concatenate(np.log(states[:,:,1]/y0))
    result["kinetic"]=np.concatenate(np.sum((states[:,:,2]**2+states[:,:,3]**2)/(2*states[:,:,1]**2)))
    return result

#regretting this choice of depth check but ay it is what it is
def depthmetropolis(alpha,burn,L,N,seed,sweeps,Tw,y0):
    sampler=Particles()
    np.random.seed(seed)
    state=sampler.partciles(L,N,Tw,y0)
    V=potentialE(state,L,alpha)
    depth=[]

    for sweep in range(sweeps):
        for i in range(N):
            trial=state.copy()
            trial[i,0]=np.random.uniform(0,L)
            v=np.random.uniform(0,1)
            tildeV=potentialE(trial,L,alpha)
            deltaV=tildeV-V
            if deltaV<=0 or np.random.uniform(0,1)<np.exp(-deltaV/Tw):
                state=trial
                V=tildeV
        if sweep>=burn and (sweep-burn) % sample_every ==0:
            for i in range(N):
              depth.append(np.log(state[i,1]/y0))
    return np.array(depth)

def saveresults(results,output,y0,temperature,alpha,diffuse_probability,mc_depth):
    amount=len(results)
    means=np.empty((amount,6))
    means[:]=np.nan
    sall=[]
    thetaall=[]
    lengthall=[]
    kall=[]
    finalheights=[]

    for trial in range(M):
        r=results[trial]
        s=r["depth"]
        theta=r["angles"]
        length=r["paths"]
        k=r["kinetic"]
        means[trial,0]=np.mean(s)
        if len(theta)>0:
            means[trial,1]=np.mean(np.cos(theta))
            means[trial,5]=np.mean(r["kinds"])
        if len(length)>0:
            means[trial,2]=np.mean(length)
        means[trial,3]=np.mean(k)
        means[trial,4]=np.mean(k**2)/np.mean(k)**2
        sall.extend(s)
        thetaall.extend(theta)
        lengthall.extend(length)
        kall.extend(k)
        finalheights.extend(r["states"][-1,:,1])
        E=r["energy"]
        Px=r["momentum"]
        deltaE=np.max(abs(E-E[0]))/E[0]
        deltaPx=np.max(abs(Px-Px[0]))
        print("This is run",trial+1, "The energy drift/change is:", deltaE, "x momentum:", deltaPx)

    average=np.empty(6)
    error=np.empty(6)
    #for some runs, there could be no results for certain variables so np.nan keeps it solid
    average[:]=np.nan
    error[:]=np.nan
    for i in range(6):
        values=means[:,i]
        #this removes the np.nan values
        values=values[np.isfinite(values)]
        #not taking any chances like in the last project
        if len(values)>0:
            average[j]=np.mean(values)
        if len(values)>1:
            error[j]=np.std(values,ddof=1)/np.sqrt(len(values))

#2) printing all outputs from simulation

#the 4 key statistics should go first - depth, cos(theta), path, energy, temp
    print("")

#diffusion statistics from the wall bouncing

#print out the exact comparison we wanted to see/test

#2) plot all the heights of the particles at the end of the simulations

#first depth distribution e.g. particle density vs s=log(y/y0)

#wall angles and their density e.g. theta vs cos(theta)/2

#path lengths - a histogram of lengths perhaps

#energy over time

#momenta change in x

#kinetic energy to compare to the expected equilibrium distribution at temp T

#the final particle heights over hyperbolic space (well the slice of it that we chose)

#4) clean up plots/present results....