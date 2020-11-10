import pygame
import numpy as np
pygame.init()

from engine import Engine 

engine = Engine()
engine.load_settings("settings.yaml")
pxPerU = 20
CANVASLENGTH = pxPerU * engine.config["boxSize"]
screen = pygame.display.set_mode([CANVASLENGTH,CANVASLENGTH])

running = True
while running:

    # Did the user click the window close button?
    for event in pygame.event.get():
        if event.type == pygame.QUIT:
            running = False

    # Fill the background with white
    screen.fill((255, 255, 255))

        
    for i in range(engine.config["N"]):
        px = engine.pos[i][0] * pxPerU
        py = CANVASLENGTH - engine.pos[i][1] * pxPerU
        unit_vel = (engine.vel / np.linalg.norm(engine.vel)) * 20
        pygame.draw.circle(screen, (0, 0, 0), (int(round(px)), int(round(py))), int(round(pxPerU/2)))
        pygame.draw.line(screen, (255,0,0), (int(round(px)), int(round(py))), (int(round(px+unit_vel[i][0]*5)), int(round(py+unit_vel[i][1]*5))))
        
    # Flip the display
    for i in range(5):
        engine.step_forward()

    pygame.display.flip()

# Done! Time to quit.
pygame.quit()