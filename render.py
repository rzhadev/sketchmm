import os
import tkinter as tk
from engine import Engine


class Renderer():
    def __init__(self, root, size, steps_per_frame):
        self.engine = Engine()
        self.size = size
        self.steps_per_frame = steps_per_frame
        self.canvas = tk.Canvas(root, width=size, height=size)
        self.canvas.pack()
        # get the size of 1 natural unit in pixels
        self.unit_pixel = (size / self.engine.config["boxSize"])
        self.radius = self.unit_pixel / 2
        self.atoms = []
        self.add_atoms()
        self.draw()

    def add_atoms(self):
        for i in range(self.engine.N):
            posx = self.engine.pos[i][0] * self.unit_pixel
            posy = (self.engine.pos[i][1] * self.unit_pixel)
            # center the oval at posx,posy
            leftx = posx - self.radius
            lefty = posy - self.radius
            rightx = posx + self.radius
            righty = posy + self.radius
            # add fill color to match velocity
            self.atoms.append(
                self.canvas.create_oval(
                    leftx, lefty, rightx, righty, fill="red")
            )

    def draw(self):
        for _ in range(self.steps_per_frame):
            self.engine.step_forward()
            self.engine.debug()
        for i in range(self.engine.N):
            posx = self.engine.pos[i][0] * self.unit_pixel
            posy = (self.engine.pos[i][1] * self.unit_pixel)
            self.radius = self.unit_pixel/2
            leftx = posx - self.radius
            lefty = posy - self.radius
            rightx = posx + self.radius
            righty = posy + self.radius
            self.canvas.coords(self.atoms[i], leftx, lefty, rightx, righty)
        self.canvas.after(10, self.draw)
