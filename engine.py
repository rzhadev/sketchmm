import numpy as np
from config import Config
from atom import Atom
from data import Collector
# manage and handle movement of molecules
# T and P control
# link with graping system and rendering


class Engine(object):

    # elements -> np array of atom objects
    # N -> number of molecules
    # con -> configuration object
    # iterations -> number of steps taken
    # time -> simulation time in natural units
    # running -> whether simulation is running or not
    # stats -> store system quantities
    def __init__(self, con_filepath="setup.yaml"):
        self.config = Config()
        self.config.load_yaml(con_filepath)
        self.atom_list = self.__prep()
        self.time = 0.0
        self.iterations = 0
        self.N = len(self.atom_list)
        self.running = False
        self.stats = {"potential_energy": 0.0,
                      "kinetic_energy": 0.0,
                      "energy": 0.0,
                      "pressure": 0.0,
                      "temperature": 0.0,
                      "sampling_count": 0,
                      "avg_T": 0.0,
                      "avg_P": 0.0}
        self.dist = []
        self.times = []

    # load in preset atom data from yaml file
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
        self.stats["kinetic_energy"] = 0.0
        mass = self.config.mapping["parameters"]["m"]
        for atom in self.atom_list:
            self.stats["kinetic_energy"] += 0.5 * \
                ((mass*atom.vel[0]*atom.vel[0]) +
                 (mass*atom.vel[1]*atom.vel[1]))
        self.stats["energy"] = self.stats["kinetic_energy"] + \
            self.stats["potential_energy"]

    def debug(self):

        for atom in self.atom_list:
            print(atom)

        # print(f"KE {self.stats['kinetic_energy']}")
        # print(f"PE {self.stats['potential_energy']}")
        print(f"E {self.stats['energy']}\n")

    # run MD simulation
    def simulate(self):
        self.running = True
        while self.running:
            self.step_forward()
            self.compute_stats()
            self.debug()

    # run a fixed number of steps
    def fixed_steps_simulation(self, steps):
        for _ in range(steps):
            self.step_forward()
            self.compute_stats()
            self.debug()

    # move the simulation forward by one time step
    # using velocity verlet alg

    def step_forward(self):
        self.iterations += 1
        print(f"time: {self.time}")
        print(f"iterations: {self.iterations}")
        dt = self.config.mapping["parameters"]["timeStep"]
        dt_sq = dt * dt
        for i in range(self.N):

            # update positions to 2nd deg accuracy
            self.atom_list[i].pos[0] += self.atom_list[i].vel[0] * \
                dt + self.atom_list[i].acc[0] * 0.5 * dt_sq
            self.atom_list[i].pos[1] += self.atom_list[i].vel[1] * \
                dt + self.atom_list[i].acc[1] * 0.5 * dt_sq

            self.atom_list[i].vel[0] += (0.5 * self.atom_list[i].acc[0] * dt)
            self.atom_list[i].vel[1] += (0.5 * self.atom_list[i].acc[1] * dt)

        self.update_accelerations()
        # update velocities by full time step
        for y in range(self.N):
            self.atom_list[y].vel[0] += (0.5 * self.atom_list[y].acc[0] * dt)
            self.atom_list[y].vel[1] += (0.5 * self.atom_list[y].acc[1] * dt)

        self.time += dt

    # use analytical derivation of L-J
    # potential to calculate interaction forces between molecules
    # use F=ma to find accelerations

    def update_accelerations(self):
        self.stats["potential_energy"] = 0.0
        epsilon = self.config.mapping["parameters"]["epsilon"]
        sigma = self.config.mapping["parameters"]["sigma"]
        mass = self.config.mapping["parameters"]["m"]

        # set F=0 at cutoff distance
        cutoff = 3.5
        cutoff_sq = cutoff * cutoff

        # subtract cutoff potential energy to avoid discontinuity
        cut_pE = 4 * epsilon * \
            (np.power((sigma/cutoff), 12) - np.power((sigma/cutoff), 6))

        for i in range(self.N):
            for j in range(i):
                atom1 = self.atom_list[i]
                atom2 = self.atom_list[j]
                delta_x = atom1.pos[0] - atom2.pos[0]
                delta_x_sq = delta_x * delta_x
                # do not calculate if distance is too big
                if delta_x_sq < cutoff_sq:
                    delta_y = atom1.pos[1] - atom2.pos[1]
                    delta_y_sq = delta_y * delta_y
                    if delta_y_sq < cutoff_sq:
                        r_sq = delta_x_sq + delta_y_sq
                        if r_sq < cutoff_sq:
                            # print(f"delta x: {delta_x}")
                            # print(f"delta y: {delta_y}")
                            print(f"dist: {np.power(r_sq, 0.5)}")
                            self.dist.append(np.power(r_sq, 0.5))
                            
                            self.times.append(self.time)
                            inv_r_sq = 1.0 / (r_sq)
                            pinv_r_sq = sigma/(r_sq)
                            # attractive term
                            term2 = pinv_r_sq * pinv_r_sq * pinv_r_sq
                            # repulsive term
                            term1 = term2 * term2

                            self.stats["potential_energy"] += 4.0 * \
                                epsilon * (term1 - term2) - cut_pE

                            forcex = delta_x * 24.0 * epsilon * ((2.0 * term1) - term2) * inv_r_sq  # noqa
                            forcey = delta_y * 24.0 * epsilon * ((2.0 * term1) - term2) * inv_r_sq  # noqa
                            # print(f"forcex : {forcex}")
                            # print(f"forcey : {forcey}")
                            # print(f"force mag: {np.power(forcex * forcex + forcey * forcey, 0.5)}")  # noqa

                            self.atom_list[i].acc[0] = (forcex/mass)
                            self.atom_list[i].acc[1] = (forcey/mass)
                            self.atom_list[j].acc[0] = -(forcex/mass)
                            self.atom_list[j].acc[1] = -(forcey/mass)
