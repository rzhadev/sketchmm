import pygame
import pygame_gui as pygg
import numpy as np
pygame.init()

from engine import Engine 

engine = Engine()
engine.load_settings("settings.yaml")
#engine.debugging = True
pxPerU = 20
RAD = int(round(pxPerU/2))
CANVASLENGTH = pxPerU * engine.config["boxSize"]
screen = pygame.display.set_mode((CANVASLENGTH,CANVASLENGTH))
manager = pygg.UIManager((CANVASLENGTH, CANVASLENGTH))
clock = pygame.time.Clock()

running = True
while running:
    time_delta = clock.tick(60)/1000.0
    for event in pygame.event.get():
        if event.type == pygame.QUIT:
            running = False
        manager.process_events(event)
    
    screen.fill((255, 255, 255))
    unit_vel = (engine.vel / np.linalg.norm(engine.vel)) * pxPerU
    unit_vel[:, 1] *= -1
    for i in range(engine.config["N"]):
        px = int(round(engine.pos[i][0] * pxPerU))
        py = int(round(CANVASLENGTH - engine.pos[i][1] * pxPerU))
        # reverse y-direction of the arrow (because Y axis is reversed)
        vx = int(round(px+unit_vel[i][0]*10))
        vy = int(round(py+unit_vel[i][1]*10))
        pygame.draw.circle(screen, (0, 0, 0), (px, py), RAD)
        pygame.draw.line(screen, (255, 0, 0), (px, py), (vx, vy), 3)
        
    for i in range(5):
        engine.step_forward()

    print(f"fps: {clock.get_fps()}", end='\r')
    manager.update(time_delta)
    manager.draw_ui(screen)
    pygame.display.flip()

pygame.quit()