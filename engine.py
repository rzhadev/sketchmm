import numpy as np
from atom import Atom
# manage and handle movement of molecules
# T and P control
# link with graping system and rendering


class Engine(object):

    # elements -> np array of atom objects
    # con -> configuration object
    # iterations -> number of steps taken
    # time -> simulation time in natural units
    # running -> whether simulation is running or not
    def __init__(self, con):
        self.config = con
        self.atom_list = self.__prep()
        self.time = 0
        self.iterations = 0
        self.running = False

    # load in preset atom data
    def __prep(self):
        atoms = []
        for k, v in self.config.mapping["data"].items():
            atom = Atom(np.array(v["pos"]),
                        np.array(v["vel"]),
                        np.array(v["acc"]),
                        id=k)

            atoms.append(atom)
        return np.array(atoms)

    def debug(self):
        for atom in self.atom_list:
            print(atom)
    # run MD simulation

    def simulate(self):
        pass

    # move the simulation forward by
    # one time step using velocity verlet alg
    def step_forward(self):
        dt = self.config.mapping["parameters"]["timeStep"]
        dt_sq = self.config.mapping["parameters"]["timeStep"] * \
            self.config.mapping["parameters"]["timeStep"]
        for i in range(len(self.atom_list)):
            atom = self.atom_list[i]
            atom.pos[0] = atom.pos[0] + atom.vel[0]*dt + \
                atom.acc[0] * 0.5 * dt_sq
            atom.pos[1] = atom.pos[1] + atom.vel[1]*dt + \
                atom.acc[1] * 0.5 * dt_sq
            atom.vel[0] = atom.vel[0] + 0.5 * \
                atom.acc[0] * dt
            atom.vel[1] = atom.vel[1] + 0.5 * \
                atom.acc[1] * dt
            self.atom_list[i].pos = atom.pos
            self.atom_list[i].vel = atom.vel
        self.calculate_accelerations(self.atom_list)
        for y in range(len(self.atom_list)):
            atom = self.atom_list[y]
            atom.vel[0] = atom.vel[0] + 0.5 * \
                atom.acc[0] * dt
            atom.vel[1] = atom.vel[1] + 0.5 * \
                atom.acc[1] * dt
            self.atom_list[y].vel = atom.vel
        self.time = self.time + dt
        self.iterations = self.iterations + 1

    # use analytical derivation of L-J
    # potential to calculate interaction forces between molecules
    # walls are "hard" and follow conservation of
    # momentum (could change to springs)

    def calculate_accelerations(self, atom_list):
        epsilon = self.config.mapping["parameters"]["epsilon"]
        sigma = self.config.mapping["parameters"]["sigma"]
        mass = self.config.mapping["parameters"]["m"]
        for i in range(len(self.atom_list)-1):
            for j in range(i+1, len(self.atom_list)):
                atom1 = self.atom_list[i]
                atom2 = self.atom_list[j]
                delta_x = atom1.pos[0] - atom2.pos[0]
                delta_x_sq = delta_x * delta_x
                # add in delta_x cutoff
                delta_y = atom1.pos[1] - atom2.pos[1]
                delta_y_sq = delta_y * delta_y
                inv_r_sq = 1/(delta_x_sq + delta_y_sq)
                # used in the L-J attract/repulsive term
                pinv_r_sq = sigma/(delta_x_sq + delta_y_sq)
                term2 = pinv_r_sq * pinv_r_sq * pinv_r_sq
                term1 = term2 * term2
                forcex = delta_x * 24 * epsilon * inv_r_sq * (2*term1 - term2)
                forcey = delta_y * 24 * epsilon * inv_r_sq * (2*term1 - term2)
                self.atom_list[i].acc[0] = self.atom_list[i].acc[0] + \
                    (forcex/mass)
                self.atom_list[i].acc[1] = self.atom_list[i].acc[1] + \
                    (forcey/mass)
                self.atom_list[j].acc[0] = self.atom_list[j].acc[0] - \
                    (forcex/mass)
                self.atom_list[j].acc[1] = self.atom_list[j].acc[1] - \
                    (forcey/mass)
