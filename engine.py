import numpy as np
from config import Config
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
    def __init__(self, con_filepath="setup.yaml"):
        self.config = Config()
        self.config.load_yaml(con_filepath)
        self.atom_list = self.__prep()
        self.time = 0
        self.iterations = 0
        self.N = len(self.atom_list)
        self.running = False
        self.stats = {"potential_energy": 0.0,
                      "kinetic_energy": 0.0,
                      "energy": 0.0,
                      "pressure": 0,
                      "temperature": 0,
                      "sampling_count": 0,
                      "avg_T": 0,
                      "avg_P": 0}

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

    def compute_stats(self):
        for atom in self.atom_list:
            self.stats["kinetic_energy"] += 0.5 * \
                (atom.vel[0]*atom.vel[0] + atom.vel[1]*atom.vel[1])
        self.stats["energy"] = self.stats["kinetic_energy"] + \
            self.stats["potential_energy"]

    def debug(self):

        for atom in self.atom_list:
            print(atom)

        print(f"KE {self.stats['kinetic_energy']}")
        print(f"PE {self.stats['potential_energy']}")
        print(f"E {self.stats['energy']}")

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

            # update positions to 2nd deg accuracy
            self.atom_list[i].pos[0] += self.atom_list[i].vel[0] * \
                dt + self.atom_list[i].acc[0] * 0.5 * dt_sq
            self.atom_list[i].pos[1] += self.atom_list[i].vel[1] * \
                dt + self.atom_list[i].acc[1] * 0.5 * dt_sq

            self.atom_list[i].vel[0] += (0.5 * self.atom_list[i].acc[0] * dt)
            self.atom_list[i].vel[1] += (0.5 * self.atom_list[i].acc[1] * dt)

        self.calculate_accelerations(self.atom_list)
        for y in range(len(self.atom_list)):
            self.atom_list[y].vel[0] += (0.5 * self.atom_list[y].acc[0] * dt)
            self.atom_list[y].vel[1] += (0.5 * self.atom_list[y].acc[1] * dt)

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
        for i in range(self.N):
            for j in range(i):
                # add in force cutoff
                atom1 = self.atom_list[i]
                atom2 = self.atom_list[j]
                delta_x = atom1.pos[0] - atom2.pos[0]
                delta_x_sq = delta_x * delta_x
                delta_y = atom1.pos[1] - atom2.pos[1]
                delta_y_sq = delta_y * delta_y
                r_sq = delta_x_sq + delta_y_sq

                print(f"dist: {delta_x_sq + delta_y_sq}")

                inv_r_sq = 1/(r_sq)
                pinv_r_sq = sigma/(r_sq)
                # repulsive term
                term2 = pinv_r_sq * pinv_r_sq * pinv_r_sq
                # attractive term
                term1 = term2 * term2

                self.stats["potential_energy"] += 4 * \
                    epsilon * (term1 - term2)

                forcex = delta_x * (24 * epsilon * ((2*term1) - term2) * inv_r_sq)   # noqa
                forcey = delta_y * (24 * epsilon * ((2*term1) - term2) * inv_r_sq)  # noqa

                self.atom_list[i].acc[0] += (forcex/mass)
                self.atom_list[i].acc[1] += (forcey/mass)
                self.atom_list[j].acc[0] -= (forcex/mass)
                self.atom_list[j].acc[1] -= (forcey/mass)
