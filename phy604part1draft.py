#This is Part I of PHY607 Project 2 : Monte Carlo

import numpy as np


#This class is mostly formal, honestly, but is somewhat useful
#It makes random numbers!
class Particles:

    #I want this class to take in user inputs, from the main prompt,
    #and then compute various things from those input values, such as
    #periodic length L, minimum y_0, trial count, amount of steps, particle number,
    #could be integration method, temperature, diffusion probabibility, and repulsion term
    #pretty sure thats everything. could be missing some, but mostly.
    #no im definitely forgetting some stuff.... whatever we'll figure it out along the way

    #regardless, this script takes in those assuming they are done.

    #1) compute depth into hyperbolic space from inputs

    def exp

    #2) compute many uniform random numbers (the rejection sampling)

    def normal
      #input an if statement according to some distribution

    #3) computes particles and their states/velocities from input

    def particles
      #should take in n particles, length L of space x, minimum y0, temp T at least

    #4) velocity after particle hits the wall if diffuse

    def wallbounce_velocity
      #computes v_x and v_y with new formula if diffuse is chosen, then returns them

    #that should be it-but this is modular enough to add more without issue later