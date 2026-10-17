# PHY607Project2

Di Wu and Gabe Suarez for class PHY607 in Fall 2026 under Professor Alexander Nitz

This project was to demonstrate our ability to utilize numerical differential equation methods, as well as monte carlo, as well as if we can use basic def, classes, and develop a package under the defined constraints at https://github.com/ahnitz/phy-607-2026/blob/main/class_notes/projects/project2.md. 

Please note, we did try to find a quantum mechanics project to do this with, but genuinely all the ideas fell into two categories: 1) impossibly nonlinear/prone to subtle computational errors / way too hard, 2) or just didn't evolve, since most quantum (gravity) problems are timeless. So, what better than the next best thing: scattering particles within hyperbolic space? Just like cosmological billiards, except classical. How so? Well, typically in these quantum chaos billiards systems we take our scattering space as the upper half plane H, and mod out by some discrete subgroup of PSL(2,R), Gamma, sometimes referred to as a Fuchsian group. So we're scattering in X=H/Gamma. For us Gamma is just the discrete translation subgroup which leaves us with x \ in [0,L], and we take y in [y0,infinity). 

What does this all mean??? Well, not much so far. Take N particles which interact under a Coulomb like potential V_ij for particles i and j. Give them random initial positions, and momenta, and then let them rip. 

We apply two constraints. If the particles go past 0 or L, we force them back in (well actually its modular arithmetic so it's automatic, but physically). If the particles hit the y=y0 wall, we reflect them back either specularly, like light reflecting (flip ydot -> -ydot), or with some probability Pd, we reflect back in a diffusive way, changing the velocity in x and y.

-----

Anyways if you want to run the script, you need Python 3.10 or newer. Dependencies are numpy,scipy,matplotlib, and mpmath. 

Will finish later....
