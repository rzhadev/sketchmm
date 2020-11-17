import time
import matplotlib.pyplot as plt
import numpy as np
from engine import Engine


# benchmark the performance of the force field functions
def perft():
    pass


# maxwell boltzmann probability density function in 2 dimensions

def maxwellboltzmann(v, m, k, T):
    return ((m * v) / (k * T)) * np.exp((-m * v ** 2)/(2 * k * T))


if __name__ == "__main__":
    perft()
