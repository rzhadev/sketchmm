import time
import matplotlib.pyplot as plt
import numpy as np
from engine import Engine


# benchmark the performance of the force field functions
def perft():
    e = Engine()
    e.config["boxSize"] = 50
    N = np.arange(2, 5000, step=10)
    loop_time = np.zeros((len(N), 2))
    cell_time = np.zeros((len(N), 2))
    for i in range(len(N)):
        e.config["N"] = N[i]
        e.reset()
        start = time.time()
        e.apply_force_field1()
        end = time.time()
        e.apply_force_field2()
        end1 = time.time()
        cell_time[i] = np.asarray((N[i], end-start))
        loop_time[i] = np.asarray((N[i], end1-end))
        print(i)
    return loop_time, cell_time


# maxwell boltzmann probability density function in 2 dimensions


def maxwellboltzmann(v, m, k, T):
    return ((m * v) / (k * T)) * np.exp((-m * v ** 2)/(2 * k * T))


if __name__ == "__main__":
    perft()
