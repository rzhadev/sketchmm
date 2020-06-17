import numpy as np


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
        return f" id {self.name}, pos {self.pos}, vel {self.vel}, acc {self.acc}"  # noqa
