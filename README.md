# PHY607Project2

Di Wu and Gabe Suarez for class PHY607 in Fall 2026 under Professor Alexander Nitz

This project was to demonstrate our ability to utilize numerical differential equation methods, as well as monte carlo, as well as if we can use basic def, classes, and develop a package under the defined constraints at https://github.com/ahnitz/phy-607-2026/blob/main/class_notes/projects/project2.md. 

Please note, we did try to find a quantum mechanics project to do this with, but genuinely all the ideas fell into two categories: 1) impossibly nonlinear/prone to subtle computational errors / way too hard, 2) or just didn't evolve, since most quantum (gravity) problems are timeless. So, what better than the next best thing: scattering particles within hyperbolic space? Just like cosmological billiards, except classical. How so? Well, typically in these quantum chaos billiards systems we take our scattering space as the upper half plane H, and mod out by some discrete subgroup of PSL(2,R), Gamma, sometimes referred to as a Fuchsian group. So we're scattering in X=H/Gamma. For us Gamma is just the discrete translation subgroup which leaves us with x \ in [0,L], and we take y in [y0,infinity). 

What does this all mean??? Well, not much so far. Take N particles which interact under a Coulomb like potential V_ij for particles i and j. Give them random initial positions, and momenta, and then let them rip. 

We apply two constraints. If the particles go past 0 or L, we force them back in (well actually its modular arithmetic so it's automatic, but physically). If the particles hit the y=y0 wall, we reflect them back either specularly, like light reflecting (flip ydot -> -ydot), or with some probability Pd, we reflect back in a diffusive way, changing the velocity in x and y.

-----

Anyways if you want to run the script, you need Python 3.10 or newer. Dependencies are numpy,scipy,matplotlib, and mpmath. 

To explain the inputs a tad (though the report details them more)

rk4 vs scipy is the evolution method / ODE solver
N is particle number
Nd, and dt specify the steps and time step size
trials are the number of trials
burn is our choice of specifying how many steps from the start to not include in the analysis
L,y0,Tw are the length of the x axis, y0 is the minimum/wall height, Tw is the wall temperature
Pd, alpha are the diffusion probability, and repulsion strength
seed is just the random seed choice for deterministic randomness
sweeps is like dt but for the metropolis sweep
burns is the same as burn just for metropolis

the code inside the main phy607project2_billiards represent...

phy607main.py is the script storing the user input & calls to ...
phy607part1.py is the initial conditions of the particle
hyperbolic_dynamics.py are the forces, sets up the EOM, and general ODE solver 
phy607part2.py takes care to call to hyperbolic_dynamics.py, evolves the EOM, and takes care of wall bounces, collisions (based on wall.py original)
__init.py is just empty and required according to https://packaging.python.org/en/latest/tutorials/packaging-projects/

If you want to run this, download the whole repo, unzip it, or find some way to get to the folder which has pyproject.toml, and install the project using pip install .  (or whatever extra syntax you need)

Then you just type phy607 !

Inside of another script, you'd have to install it via
1) from phy607project2_billiards.phy607main import main
2) main()
