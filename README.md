# SketchMM
a 2D classical molecular dynamics simulation written in python3.\
![](test.gif)
---
# Background
https://en.wikipedia.org/wiki/Molecular_dynamics\
Molecular dynamics are the use of computer simulations to numerically estimate the trajectory of a system. 

# Methodology
Molecular dynamics simulations typically have 5 main steps.

1. Initialization of particle positions and particle velocities.

2. Calculate forces and update positions and velocities according to an integration technique. 

3. Apply boundary conditions and control temperature or pressure if needed.

4. Calculate relevant statistics of the system.

5. Repeat step 2-4 until a fixed number of time steps have elapsed. 

# TODO
-Add GUI widgets for Matplotlib and graphing
-add text boxes for FPS and system values (temp, P, etc)
-switch to fixed canvas size AND fixed box size, will have to adjust pixels per unit
-add in color mapping for velocity? more red = faster particle, switch with arrows? 
-improve performance, maybe move engine to a python binding, currently at 150 particles only roughly 10 fps
-make engine benchmarking more robust, it looks ugly rn

# Reference
https://arxiv.org/pdf/2001.07089.pdf\
http://physics.weber.edu/schroeder/md/InteractiveMD.pdf\
http://www.courses.physics.helsinki.fi/fys/moldyn/lectures/L4.pdf\
https://web.northeastern.edu/afeiguin/p4840/p131spring04/node41.html\
https://arxiv.org/pdf/1111.2705.pdf\
http://fab.cba.mit.edu/classes/864.11/people/dhaval_adjodah/final_project/write_up_boltzmann.pdf