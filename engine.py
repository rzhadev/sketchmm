#!/usr/bin/env pypy

import numpy as np
import time
import yaml
import os
import os.path


class Atom():
    # np arrays with length 2 to store coords
    def __init__(self, p, v, a, id, c=0):
        self.name = id
        self.pos = p
        self.vel = v
        self.acc = a
        # adjust color based on velocity
        self.color = c

    def __str__(self):
        return f" id:{self.name}, pos:{self.pos}, vel:{self.vel}"


class Engine(object):

    def __init__(self):
        self.__load_yaml()
        self.N = len(self.atom_list)
        self.running = False
        self.stats = {"pE": 0.0,
                      "kE": 0.0,
                      "E": 0.0,
                      "instantP": 0.0,
                      "instantT": 0.0,
                      "sampleCount": 0,
                      "avgT": 0.0,
                      "avgP": 0.0,
                      "totalT": 0.0,
                      "totalP": 0.0,
                      "last_sample": 0.0,
                      "time": 0.0,
                      "iterations": 0}

    def __load_yaml(self, setup="./setup.yaml", data="./data.yaml"):
        dir_path = os.path.dirname(os.path.realpath(__file__))
        setup_path = os.path.join(dir_path, setup)
        data_path = os.path.join(dir_path, data)
        try:
            with open(setup_path, 'r') as f:
                self.config = yaml.load(f, Loader=yaml.SafeLoader)
                print("Loaded configuration from %s" % setup_path)
            with open(data_path, 'r') as f:
                data = yaml.load(f, Loader=yaml.SafeLoader)
                atoms = []
                for k, v in data.items():
                    atom = Atom(np.array(v["pos"]),
                                np.array(v["vel"]),
                                np.array(v["acc"]),
                                id=k)
                    atoms.append(atom)
                print("Loaded configuration from %s" % data_path)
                self.atom_list = atoms

        except FileNotFoundError and ImportError:
            print("File could not be found.")

    def compute_stats(self):
        sampleInt = self.config["sampleInterval"]
        self.stats["kE"] = 0.0
        mass = self.config["parameters"]["m"]
        kB = self.config["parameters"]["kB"]
        for atom in self.atom_list:
            self.stats["kE"] += 0.5 * \
                ((mass*atom.vel[0]*atom.vel[0]) +
                 (mass*atom.vel[1]*atom.vel[1]))
        self.stats["E"] = self.stats["kE"] + \
            self.stats["pE"]
        self.stats["instantT"] = self.stats["kE"] / (self.N * kB)

        if (self.stats["time"] - self.stats["last_sample"]) >= sampleInt:
            self.stats["totalT"] += self.stats["instantT"]
            self.stats["totalP"] += self.stats["instantP"]
            self.stats["sampleCount"] += 1
            self.stats["avgT"] = self.stats["totalT"] / \
                self.stats["sampleCount"]
            self.stats["avgP"] = self.stats["totalP"] / \
                self.stats["sampleCount"]
            self.stats["last_sample"] = self.stats["time"]

    def debug(self):

        for atom in self.atom_list:
            print(atom)

        print(f"E: {self.stats['E']}")
        print(f"avgT: {self.stats['avgT']}")
        print(f"avgP: {self.stats['avgP']}\n")

    def simulate(self):
        self.running = True
        start_time = time.time()
        while self.running:
            self.step_forward()
            self.debug()
        end_time = time.time()
        print(f"{end_time-start_time} seconds elapsed")

    # run a fixed number of steps
    def fixed_steps_simulation(self, steps):
        start_time = time.time()
        for _ in range(steps):
            self.step_forward()
            self.debug()
        end_time = time.time()
        print(f"{end_time-start_time} seconds elapsed")

    # move the simulation forward by one time step
    # using velocity verlet alg

    def step_forward(self):
        self.stats["iterations"] += 1
        print(f"time: {self.stats['time']:.4f}")
        print(f"iterations: {self.stats['iterations']}")
        dt = self.config["parameters"]["timeStep"]
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
        self.stats["time"] += dt
        self.apply_boundary()
        self.compute_stats()

    # check for wall collisions
    # apply "hard" wall with conservation of momentum
    def apply_boundary(self, buffer=0.5):
        boundary = self.config["boxSize"]
        wallForce = 0.0
        for i in range(self.N):
            atom = self.atom_list[i]
            # check if atom has hit a wall
            leftx = atom.pos[0] < buffer
            rightx = atom.pos[0] > (boundary - buffer)
            lefty = atom.pos[1] < buffer
            righty = atom.pos[1] > (boundary - buffer)
            # apply conservation of momentum,
            # forces applied on different walls will cancel out
            if leftx or rightx:
                self.atom_list[i].vel[0] = -(self.atom_list[i].vel[0])
                if leftx:
                    wallForce -= atom.acc[0]
                elif rightx:
                    wallForce += atom.acc[0]
            if lefty or righty:
                self.atom_list[i].vel[1] = -(self.atom_list[i].vel[1])
                if lefty:
                    wallForce -= atom.acc[1]
                elif righty:
                    wallForce += atom.acc[1]
        # pressure = force per unit area
        self.stats["instantP"] = wallForce / (4 * boundary)

    # use analytical derivation of L-J
    # potential to calculate interaction forces between molecules
    # use F=ma to find accelerations
    def update_accelerations(self):
        self.stats["pE"] = 0.0
        epsilon = self.config["parameters"]["epsilon"]
        sigma = self.config["parameters"]["sigma"]
        mass = self.config["parameters"]["m"]

        # set F=0 at cutoff distance
        cutoff = 3.5
        cutoff_sq = cutoff * cutoff

        # subtract cutoff potential energy to avoid discontinuity
        cut_pE = 4 * epsilon * \
            (np.power((sigma/cutoff), 12) - np.power((sigma/cutoff), 6))

        # N^2 double loop for force calculation
        # could switch to N alg like cell lists
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
                            inv_r_sq = 1.0 / (r_sq)
                            pinv_r_sq = sigma/(r_sq)
                            # attractive term
                            term2 = pinv_r_sq * pinv_r_sq * pinv_r_sq
                            # repulsive term
                            term1 = term2 * term2

                            self.stats["pE"] += 4.0 * \
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
