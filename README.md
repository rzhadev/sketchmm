# SketchMM
SketchMM is a 2D classical molecular dynamics simulation written in python3.
![](test.gif)
---
# Theory
-F = ma governs the dynamics of a classical molecular dynamics simulation
-Where F_a is the net force acting on a particle A, m_a is the mass of the particle, and a_a is the acceleration of the particle
-The net force F_a is the sum of all interactions provided by the negative gradient of a potential energy function. 
-LJ potential specifies the pair-wise interactions present between 2 molecules A,B of distance r_AB apart. 
-In order to speed up computation, LJ potential is shifted and truncated at a specified cutoff distance
-Periodic Boundary Conditions define the boundaries of the system. 
-The system is surrounded by periodic images (mirror images) of itself
# TODO
# Reference
https://arxiv.org/pdf/2001.07089.pdf\
http://physics.weber.edu/schroeder/md/InteractiveMD.pdf\
http://www.courses.physics.helsinki.fi/fys/moldyn/lectures/L4.pdf\
Periodic Boundary Conditions
https://web.northeastern.edu/afeiguin/p4840/p131spring04/node41.html