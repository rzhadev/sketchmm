import tkinter as tk
import os
import threading
from engine import Engine
from render import Renderer

CANVAS_DIM = 800
root = tk.Tk()
root.title("Sketch Molecular Dynamics")
root.resizable(False, False)
renderer = Renderer(root, CANVAS_DIM, 1)
root.mainloop()
