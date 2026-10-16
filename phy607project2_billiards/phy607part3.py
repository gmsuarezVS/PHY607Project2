import matplotlib.pyplot as plt
import mpmath as mp
import numpy as np

from .phy607part1 import Particles #must update name later !!

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
    return kineticE+potentialE(state,length,alpha)

#how we "know" our code is "right" (quotes because not really in our case, but good enough)
def themeasure(y0,length,alpha,result):
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
    result["kinetic"]=np.concatenate((states[:,:,2]**2+states[:,:,3]**2)/(2*states[:,:,1]**2))
    return result

#regretting this choice of depth check but ay it is what it is
def depthmetropolis(alpha,burn,L,N,seed,sweeps,Tw,y0):
    np.random.seed(seed)
    state=Particles().particles(N,y0,L,Tw)
    V=potentialE(state,L,alpha)
    depth=[]

    for sweep in range(sweeps):
        for i in range(N):
            trial=state.copy()
            trial[i,0]=np.random.uniform(0,L)
            v=np.random.uniform(0,1)
            trial[i,1]=y0/(1-v)
            tildeV=potentialE(trial,L,alpha)
            deltaV=tildeV-V
            if deltaV<=0 or np.random.uniform(0,1)<np.exp(-deltaV/Tw):
                state=trial
                V=tildeV
        if sweep>=burn: #didn't need the sample_every idea anyways just record em all
            for i in range(N):
              depth.append(np.log(state[i,1]/y0))
    return np.array(depth)

def saveresults(results,y0,temperature,alpha,Pd,depths):
    amount=len(results)
    means=np.empty((amount,6))
    means[:]=np.nan
    sall=[]
    thetaall=[]
    lengthall=[]
    kall=[]
    finalheights=[]

    for trial in range(amount):
        r=results[trial]
        s=r["depth"]
        theta=r["angles"]
        length=r["paths"]
        k=r["kinetic"]
        means[trial,0]=np.mean(s)
        if len(theta)>0:
            means[trial,1]=np.mean(np.cos(theta))
            means[trial,5]=np.mean(r["types"])
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
            average[i]=np.mean(values)
        if len(values)>1:
            error[i]=np.std(values,ddof=1)/np.sqrt(len(values))

#2) printing all outputs from simulation, got \u00b1 online
#https://stackoverflow.com/questions/27872250/plus-minus-operator-for-python-%C2%B1
    print("Mean depth:",average[0],"\u00b1",error[0])
    print("Mean cos(angle):",average[1],"\u00b1",error[1])
    print("Mean path length:",average[2],"\u00b1",error[2])
    print("Mean kinetic energy (temperature estimate):",average[3],"\u00b1",error[3])
    print("Mean(k**2)/mean(k)**2:",average[4],"\u00b1",error[4])
    print("Diffuse fraction:",average[5],"\u00b1",error[5],"Input Pd:",Pd)
    print("Wall hits:",len(thetaall),"Complete paths:",len(lengthall))

#print out the exact comparison we wanted to see/test
    if alpha==0:
        print("Ideal mean depth, cos(angle), path: ",1,np.pi/4,np.pi)
    if alpha==0 or Pd>0:
        print("At equilibrium, the mean kinetic energy and moment ratio: ",temperature,2)
    if Pd==0:
        print("Since this is specular, go check the conservation of total energy and x momentum")
    else:
        print("The diffuse wall can/will likely change total energy and x momentum.")
    if depths is not None and len(depths)>0:
        print("Metropolis mean depth:",np.mean(depths))

#2) plot all the heights of the particles at the end of the simulations
    plt.figure() #leaving the figsize implicit for now

#first depth distribution e.g. particle density vs s=log(y/y0)
    plt.subplot(2,3,1)
    plt.title("Distribution of y depths")
    plt.hist(sall,density=True,label="computed ode")
    s=np.linspace(0,6)
    if alpha==0:
        plt.plot(s,np.exp(-s),label="The exact e^(-s)")
    if depths is not None and len(depths)>0:
        plt.hist(depths,density=True,histtype="step",label="metropolis simulation")
    plt.xlabel("log(y/y0)")
    plt.ylabel("Density of points")
    plt.legend()

#wall angles and their density e.g. theta vs cos(theta)/2
    plt.subplot(2,3,2)
    plt.title("The ingoing wall angles")
    if len(thetaall)>0:
        plt.hist(thetaall,density=True)
    theta=np.linspace(-np.pi/2,np.pi/2)
    if alpha==0 or Pd>0:
        plt.plot(theta,np.cos(theta)/2,label="Equilibrium cosine law")
        plt.legend()
    plt.xlabel("The angle (rad)")
    plt.ylabel("Density of angles")

#path lengths - a histogram of lengths perhaps
    plt.subplot(2,3,3)
    plt.title("Complete path lengths of hitting y0 to returning to y0 later")
    if len(lengthall)>0:
        plt.hist(lengthall,density=True)
    if alpha==0:
        plt.axvline(average[2],label="Measured mean")
        plt.legend()
    plt.xlabel("Path length")
    plt.ylabel("Density of paths at length L")

#energy over time
    plt.subplot(2,3,4)
    plt.title("Change in total energy")
    for trial in range(amount):
        time=results[trial]["times"]
        E=results[trial]["energy"]
        plt.plot(time,E-E[0])
    plt.xlabel("Time")
    plt.ylabel("\u0394 E")

#momenta change in x
    plt.subplot(2,3,5)
    plt.title("change in x momentum")
    for trial in range(amount):
        time=results[trial]["times"]
        Px=results[trial]["momentum"]
        plt.plot(time,Px-Px[0])
    plt.xlabel("Time")
    plt.ylabel("\u0394 Px")

#kinetic energy to compare to the expected equilibrium distribution at temp T
    plt.subplot(2,3,6)
    plt.title("Kinetic per particle")
    plt.hist(kall,density=True)
    k=np.linspace(0,6*temperature)
    if alpha==0 or Pd>0:
        plt.plot(k,np.exp(-k/temperature)/temperature,label="At Equilibrium")
        plt.legend()
    plt.xlabel("Kinetic E")
    plt.ylabel("Density of particles with energy E")

    plt.show()

#the final particle heights over hyperbolic space (well the slice of it that we chose)
    plt.figure()
    plt.title("Final particle heights")
    plt.hist(finalheights,density=True)
    if alpha==0:
        y=np.linspace(y0,max(finalheights))
        plt.plot(y,y0/y**2,label="Exact y0/y^2")
        plt.legend()
    plt.xlabel("Height y")
    plt.ylabel("Density of particles with height y")
    plt.show()
