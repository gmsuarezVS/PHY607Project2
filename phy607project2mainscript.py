import numpy as np

from phy607part1 import Particles
from phy607part2 import evolve
from phy607part3 import themeasure,depthmetropolis,saveresults

def main():
    print("Use commas to seperate answers. There are a lot of questions for variability")
    solver=("Which integration method - RK4, scipy: ")
    N=int(input("Input number of particles: "))
    Nd,dt = map(float,input("Input step quantity, step size: ").split(","))
    Nd=int(Nd)
    trials=int(input("Input number of independent trials: "))
    burn=int(input("Input starting steps to skip - this is useful to get to equilibrium; but put 0 if none)"))
    L,y0,Tw=map(float,input("Input period L, wall height y0, and temperature Tw: ").split(","))
    Pd, alpha = map(float,input("Input the diffuse probability Pd, repulsion strength alpha: ").split(","))
    seed=int(input("Input random seed for detemrinistic randomness: "))
    sweeps=int(input("Input metropolis sweep count; again put 0 if none: "))
    if sweeps>0:
        burns=int(input("Input the sweeps to skip (for same reason): "))

    results=[]
    for trial in range(trials):
        print("This is run", trial+1, "out of total", trials)
        np.random.seed(seed+trial)

        state=Particles().particles(N,y0,L,Tw)

        if burn>0:
            start=evolve(state,dt,burn,alpha,L,y0,Tw,Pd,solver)
            state=start["states"][-1]
        result=evolve(state,dt,Nd,alpha,L,y0,Tw,Pd,solver)
        result["times"]=result["times"]+burn*dt

        result=themeasure(y0,L,alpha,result)
        results.append(result)

    depths=None
    if sweeps>0:
        depthmetropolis(alpha,burns,L,N,seed+trials,sweeps,Tw,y0)
    saveresults(results,y0,Tw,alpha,Pd,depths)
    print("Check saveresults")