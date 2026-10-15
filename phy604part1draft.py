#This is Part I of PHY607 Project 2 : Monte Carlo

import numpy as np


#This class is mostly formal, honestly, but is somewhat useful
#It makes random numbers & computes the velocity if diffuse
class Particles:

    #1) compute depth into hyperbolic space from inputs - moved into #3)

    #2) compute many uniform random numbers - the rejection sampling
    def normal(self,n):
        #accept the random numbers iff they fit under the Gaussian at that value
        random_values=np.empty(n)
        for i in range(n):
            good = False
            while not good:
              random=np.random.uniform(-5,5)
              sample=np.random.uniform(0,1)
              if sample<=np.exp(-random**2/2):
                  good=True
            value[i]=random
        return value

    #3) computes particles and their states/velocities from input

    def particles(self,n,y0,length,temperature)
        #should take in n particles, length L of space x, minimum y0, temp T at least
        if n<1 or length<=0 or y0<=0 or temperature<=0:
            print("You have to have at least one particle... and must be physical")
        state=np.empty((n,4)) #Following syntax of Di

        #using
        for i in range(n):
            state[i,0]=np.random.uniform(0,length) #x position
            state[i,1]=y0/(1-np.random.uniform(0,1)) #y position

        vx=self.normal(n)
        vy=self.normal(n)
        #the momenta formula are rescaled by y*sqrt(T)
        state[:,2]=state[:,1]*np.sqrt(temperature)*vx
        state[:,3]=state[:,1]*np.sqrt(temperature)*vy

    #4) velocity after particle hits the wall if diffuse

    def wallbounce_velocity
        #computes v_x and v_y with new formula if diffuse is chosen, then returns them
        vx=np.sqrt(temperature)*self.normal(1)[0]

    #that should be it-but this is modular enough to add more without issue later

