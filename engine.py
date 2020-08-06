#!/usr/bin/env pypy

import numpy as np
import time
import yaml
import os
import os.path


class Engine(object):

    def __init__(self, setup="./setup.yaml", data="./data.xyz"):
        self.N = 0
        self.DIM = 2
        self.config = None
        self.pos = None
        self.vel = None
        self.acc = None
        self.id = []
        self.sim_name = ""
        self.running = False
        self.pE = 0.0
        self.kE = 0.0
        self.E = 0.0
        self.instantP = 0.0
        self.instantT = 0.0
        self.sampleCount = 0
        self.avgT = 0.0
        self.avgP = 0.0
        self.totalT = 0.0
        self.totalP = 0.0
        self.last_sample_time = 0.0
        self.time = 0.0
        self.iterations = 0
        self.__rev(setup, data)

    def __rev(self, setup, data):
        dir_path = os.path.dirname(os.path.realpath(__file__))
        setup_path = os.path.join(dir_path, setup)
        data_path = os.path.join(dir_path, data)
        try:
            with open(setup_path, 'r') as f:
                self.config = yaml.load(f, Loader=yaml.SafeLoader)
                print("Loaded configuration from %s" % setup_path)
                f.close()
        except FileNotFoundError:
            print("File not found at %s" % setup_path)
        try:
            with open(data_path, 'r') as f:
                self.N = int(f.readline())
                self.pos = np.zeros((self.N, self.DIM))
                self.vel = np.zeros((self.N, self.DIM))
                self.acc = np.zeros((self.N, self.DIM))
                self.sim_name = f.readline()
                data = [i.split() for i in f.readlines()]
                if len(data) != self.N:
                    raise SizeMissmatch
                for i in range(self.N):
                    self.id.append(data[i][0])
                    self.pos[i] = np.array(
                        (float(data[i][1]), float(data[i][2])))
                    #self.vel[i] = np.random.randn(self.DIM)
                    #self.acc[i] = np.random.randn(self.DIM)
                    self.vel[i] = np.array([0, 0])
                    self.acc[i] = np.array([0, 0])

                print("Loaded configuration from %s" % data_path)
                f.close()

        except FileNotFoundError:
            print("File not found at %s" % data_path)

    def compute_stats(self):
        sampleInt = self.config["sampleInterval"]
        self.kE = 0.0
        mass = self.config["parameters"]["m"]
        kB = self.config["parameters"]["kB"]
        for i in range(self.N):

            self.kE += 0.5 * \
                ((mass*self.vel[i][0]*self.vel[i][0]) +
                 (mass*self.vel[i][1]*self.vel[i][1]))

        self.E = self.kE + self.pE
        self.instantT = self.kE / (self.N * kB)
        if (self.time - self.last_sample_time) >= sampleInt:
            self.totalT += self.instantT
            self.totalP += self.instantP
            self.sampleCount += 1
            self.avgT = self.totalT / \
                self.sampleCount
            self.avgP = self.totalP / \
                self.sampleCount
            self.last_sample_time = self.time

    def debug(self):
        print(f"E: {self.E}")
        print(f"avgT: {self.avgT}")
        print(f"avgP: {self.avgP}\n")

    def simulate(self):
        self.running = True
        start_time = time.time()
        while self.running:
            self.step_forward()
            self.debug()
        end_time = time.time()
        print(f"{end_time-start_time} seconds elapsed")

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
        self.iterations += 1
        print(f"time: {self.time:.4f}")
        print(f"iterations: {self.iterations}")
        dt = self.config["parameters"]["timeStep"]
        dt_sq = dt * dt
        for i in range(self.N):
            # update positions to 2nd deg accuracy
            self.pos[i][0] += self.vel[i][0] * \
                dt + self.acc[i][0] * 0.5 * dt_sq
            self.pos[i][1] += self.vel[i][1] * \
                dt + self.acc[i][1] * 0.5 * dt_sq
            self.vel[i][0] += (0.5 * self.acc[i][0] * dt)
            self.vel[i][1] += (0.5 * self.acc[i][1] * dt)
        self.update_accelerations()
        # update velocities by full time step
        for y in range(self.N):
            self.vel[y][0] += (0.5 * self.acc[y][0] * dt)
            self.vel[y][1] += (0.5 * self.acc[y][1] * dt)
        self.time += dt
        # self.apply_boundary()
        self.compute_stats()

    # change to periodic boundary
    # check for wall collisions
    # apply "hard" wall with conservation of momentum
    """ def apply_boundary(self, buffer=0.5):
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
                self.vel[i][0] = -(self.vel[i][0])
                if leftx:
                    wallForce -= atom.acc[0]
                elif rightx:
                    wallForce += atom.acc[0]
            if lefty or righty:
                self.vel[i][1] = -(self.vel[i][1])
                if lefty:
                    wallForce -= atom.acc[1]
                elif righty:
                    wallForce += atom.acc[1]
        # pressure = force per unit area
        self.stats["instantP"] = wallForce / (4 * boundary)
    """
# use analytical derivation of L-J
# potential to calculate interaction forces between molecules
# use F=ma to find accelerations

    def update_accelerations(self):
        self.pE = 0.0
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
        for i in range(self.N-1):
            for j in range(i+1, self.N):
                delta_x = self.pos[i][0] - self.pos[j][0]
                delta_x_sq = delta_x * delta_x
                # do not calculate if distance is too big
                if delta_x_sq < cutoff_sq:
                    delta_y = self.pos[i][1] - self.pos[j][1]
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
                            self.pE += 4.0 * \
                                epsilon * (term1 - term2) - cut_pE
                            forcex = delta_x * 24.0 * epsilon * ((2.0 * term1) - term2) * inv_r_sq  # noqa
                            forcey = delta_y * 24.0 * epsilon * ((2.0 * term1) - term2) * inv_r_sq  # noqa
                            # print(f"forcex : {forcex}")
                            # print(f"forcey : {forcey}")
                            # print(f"force mag: {np.power(forcex * forcex + forcey * forcey, 0.5)}")  # noqa
                            self.acc[i][0] = (forcex/mass)
                            self.acc[i][1] = (forcey/mass)
                            self.acc[j][0] = -(forcex/mass)
                            self.acc[j][1] = -(forcey/mass)


class SizeMissmatch(Exception):
    pass


class DimensionMissmatch(Exception):
    pass


if __name__ == "__main__":
    e = Engine()
    e.fixed_steps_simulation(1000)
