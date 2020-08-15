#!/usr/bin/env pypy

import numpy as np
import time
import yaml
import sys
import os
from scipy.stats import maxwell
import os.path


class Engine(object):

    # load configuration settings, or default settings (parse from command line)
    # initialize particles using lattice/predefined positions
    # sample gaussian splitribution
    # according to maxwell boltzmann splitribution velocities
    def __init__(self):
        self.config = {}
        self.N = 0
        self.DIM = 2
        self.pos = np.zeros((0, 0), float)
        self.vel = np.zeros((0, 0), float)
        self.acc = np.zeros((0, 0), float)
        self.pE = 0.0
        self.kE = 0.0
        self.E = 0.0
        self.instantT = 0.0
        self.avgT = 0.0
        self.totalT = 0.0
        self.last_sample = 0.0
        self.sample_count = 0.0
        self.time = 0.0
        self.iterations = 0
    # load from config from yaml file

    def load_settings(self, settings_file='./settings.yaml'):
        dir_path = os.path.dirname(os.path.realpath(__file__))
        settings_path = os.path.join(dir_path, settings_file)
        try:
            f = open(settings_path, 'r')
            self.config = yaml.load(f, Loader=yaml.SafeLoader)
            f.close()
        except FileNotFoundError:
            print("Settings file not found.")

    # optionally load data from .xyz file
    def load_positions(self, data_file='./data.xyz'):
        dir_path = os.path.dirname(os.path.realpath(__file__))
        data_path = os.path.join(dir_path, data_file)
        try:
            f = open(data_path, 'r')
            if self.N != int(f.readline()):
                raise NSizeMissmatch
            else:
                self.N = int(f.readline())
            print(f"{f.readline()}")
            data = np.asarray([line.split() for line in f.readlines()])
            self.pos = np.zeros((self.N, self.DIM))
            self.pos[:, 0] = data[:, 1].astype(float)
            self.pos[:, 1] = data[:, 2].astype(float)
            if self.pos.shape != (self.N, self.DIM):
                raise NSizeMissmatch
            f.close()
        except FileNotFoundError:
            print("Data file not found.")

    # randomly initialize particle velocities according to distribution
    # use lattice structure for positions
    def init_simulation(self):
        self.N = self.config["N"]
        self.pos = np.zeros((self.N, self.DIM))
        self.vel = np.zeros((self.N, self.DIM))
        self.acc = np.zeros((self.N, self.DIM))

        # sample maxwell boltzmann splitribution for target temp
        # scale param a = sqrt((kB * T)/m) for maxwell boltzmann
        kB = self.config["parameters"]["kB"]
        m = self.config["parameters"]["m"]
        T = self.config["targetTemp"]
        scale = np.sqrt((kB*T/m))
        mags = maxwell.rvs(size=self.N, scale=scale, loc=0.0)
        theta = np.random.random(self.N) * 2 * np.pi
        self.vel[:, 0] = mags * np.cos(theta)
        self.vel[:, 1] = mags * np.sin(theta)

        # lattice generation
        # -for a box of L x L dimensions,
        # there are sqrt(N) boxes that cover an L length side
        # with an L/sqrt(N) distance between boxes
        boxes = int(round(np.power(self.N, (1.0 / self.DIM))))
        split = self.config["boxSize"] / float(boxes)
        index = 0
        for i in range(boxes):
            x = (0.5 + i) * split
            for j in range(boxes):
                if index < self.N:
                    y = (0.5 + j) * split
                    self.pos[index] = np.asarray((x, y))
                    index += 1

    # # calculate energy, potential energy
    # kinetic energy, instant temperature, average temperature
    # of the system
    # (only calculate T if within the sample interval)

    def stats(self):
        pass

    # move simulation one time step further according to
    # velocity verlet
    # apply periodic boundary conditions
    # and rescale velocities if within the interval
    # rescale by factor of sqrt(targetT/avgT)
    def step_forward(self):
        pass

    # calculate pair-wise interactions between molecules
    # or within periodic images of the system (if greater than cutoff)
    # apply forces to all particles using newton's second law
    def apply_force_field(self):
        self.pE = 0.0
        self.acc = np.zeros((self.N, self.DIM))
        eps = self.config["parameters"]["epsilon"]
        sig = self.config["parameters"]["sigma"]
        mass = self.config["parameters"]["mass"]
        # distance for LJ cutoff
        cutoff = 2.5 * sig
        box_length = self.config["boxSize"]
        # truncate potential energy at cut off distance
        cut_pE = 4 * eps * \
            (np.power((sig/cutoff), 12) - np.power((sig/cutoff), 6))
        for i in range(self.N - 1):
            for j in range(i+1, self.N):
                # apply minimum image convention
                # maximum separation between particles is L/2
                # use location of periodic image of particle j if
                # distance is > L/2.
                delta_x = self.pos[i][0] - self.pos[j][0]
                if delta_x > box_length/2:
                    delta_x -= box_length
                elif delta_x <= -box_length/2:
                    delta_x += box_length
                delta_y = self.pos[i][1] - self.pos[j][1]
                if delta_y > box_length/2:
                    delta_y -= box_length
                elif delta_y <= -box_length/2:
                    delta_y += box_length

                delta_x_sq = delta_x ** 2
                delta_y_sq = delta_y ** 2
                dist_sq = delta_x_sq + delta_y_sq
                if dist_sq < (cutoff * cutoff):
                    inv_dist_sq = 1.0 / dist_sq
                    para_inv = sig / dist_sq
                    attra_term = para_inv ** 3
                    repul_term = attractive_term ** 2
                    self.pE += 4.0 * eps * \
                        (repul_term - attra_term) - cut_pE

                    # negative gradient of LJ potential
                    forcex = delta_x * 24.0 * eps * \
                        ((2.0 * repul_term) - attra_term) * inv_dist_sq
                    forcey = delta_y * 24.0 * eps * \
                        ((2.0 * repul_term) - attra_term) * inv_dist_sq
                    self.acc[i][0] += forcex/mass
                    self.acc[i][1] += forcey/mass
                    self.acc[j][0] -= forcex/mass
                    self.acc[j][1] -= forcey/mass


class NSizeMissmatch(Exception):
    pass


if __name__ == "__main__":
    e = Engine()
    e.load_settings()
    e.init_simulation()
