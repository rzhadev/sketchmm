from tkinter import *
from engine import Engine

canvas_width = 1200
canvas_height = 1200

top = Tk()
w = Canvas(top, width=canvas_width, height=canvas_height)

w.pack()

y = int(canvas_height / 2)
w.create_line(0, y, canvas_width, y, fill="#476042")
mainloop()
