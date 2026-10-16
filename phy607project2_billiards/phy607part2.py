import numpy as np
from scipy.integrate import solve_ivp

from .hyperbolic_dynamics import equations, rk4_step
from .phy607part1 import Particles #need to go back and rename script 1 !!

def scipys(time,values,alpha,length):
    state=values.reshape((-1,4))
    return np.concatenate(equations(state,alpha,length))

def step_with_wall(state, dt, alpha, length, y0, Tw, Pd,solver): #order is (hopefully) right now
    # Part II calculates one proposed step
    if solver == "rk4":
        new_state = rk4_step(state, dt, equations, alpha, length)
    else:
        solution=solve_ivp(scipys,(0,dt),np.concatenate(state),args=(alpha,length))
        new_state=solution.y[:,-1].reshape(state.shape)

    hits=[]
    angles=[]
    types=[]

    # Part III checks every particle
    for i in range(len(state)):
        if new_state[i, 1] < y0:
            hits.append(i)
            angles.append(np.arctan2(new_state[i,2],-new_state[i,3]))
            new_state[i, 1] = y0

            #swapping rng random to np
            if np.random.uniform(0,1) >= Pd:
                # Specular reflection
                new_state[i, 3] = abs(new_state[i, 3])
                types.append(0) #I believe this append should be 0 for specular

            else:
                #replacing with pt 1 class now that it exists
                """# Diffuse reflection
                new_state[i, 2] = y0 * np.sqrt(Tw) * rng.normal()
                new_state[i, 3] = y0 * np.sqrt(
                    -2.0 * Tw * np.log(1.0 - rng.random())
                )"""
                new_state[i,2:4]=Particles().wallbounce_velocity(y0,Tw)
                types.append(1)

    # Periodic boundary in x
    new_state[:, 0] %= length

    return new_state,hits,angles,types

#replacing rng w np.random for now (sry)
def evolve(state, dt, steps, alpha, length, y0, Tw, Pd, solver):
    history = np.empty((steps + 1, len(state), 4))
    history[0] = state
    angles=[]
    paths=[]
    types=[]
    distance=np.zeros(len(state))
    wallhit=np.zeros(len(state))

    for step in range(steps):
        distance=distance+dt*(np.sqrt(state[:,2]**2+state[:,3]**2)/state[:,1])
        #apparently can't use type so using typE (it's fine)
        state,hits,theta,typE = step_with_wall(
            state, dt, alpha, length, y0, Tw, Pd, solver
        )
        for i in hits:
            if wallhit[i]==1:
                paths.append(distance[i])
            distance[i]=0
            wallhit[i]=1
        angles.extend(theta)
        types.extend(typE)
        history[step + 1] = state

    result={}
    result["times"]=np.arange(steps+1)*dt
    result["states"]=history
    result["angles"]=np.array(angles)
    result["paths"]=np.array(paths)
    result["types"]=np.array(types)
    return result
