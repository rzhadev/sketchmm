import pygame
import pygame_gui as pygg
import numpy as np
from engine import Engine 

pygame.init()
engine = Engine()
engine.load_settings("settings.yaml")
#engine.debugging = True
CANVASLENGTH = 800
pxPerU = int(round(CANVASLENGTH / engine.config["boxSize"]))
RAD = int(round(pxPerU/2))
MAGFACTOR = 2

screen = pygame.display.set_mode((1200,CANVASLENGTH))
pygame.display.set_caption("SketchMM")
manager = pygg.UIManager((1200, CANVASLENGTH), theme_path="theme.json")
clock = pygame.time.Clock()
start_button = pygg.elements.UIButton(relative_rect=pygame.Rect((840, 10), (170, 50)), 
                                        text='Start', 
                                        manager=manager)
step_button = pygg.elements.UIButton(relative_rect=pygame.Rect((1020, 10), (170, 50)),
                                        text="Step",
                                        manager=manager)
info_block = pygg.elements.ui_label.UILabel(relative_rect=pygame.Rect((840, 70), (350, 50)), 
                                        text="FPS: 0.0    Time: 0.0    Iterations: 0.0", 
                                        manager=manager)

panel = pygame.Surface((400, 800))


simulating = False
running = True
while running:
    time_delta = clock.tick(60)/1000.0
    for event in pygame.event.get():
        if event.type == pygame.QUIT:
            running = False 

        if event.type == pygame.USEREVENT:
            if event.user_type == pygg.UI_BUTTON_PRESSED:
                if event.ui_element == start_button:
                    if (simulating):
                        start_button.set_text("Resume")
                    else:
                        start_button.set_text("Pause")
                    simulating = not simulating

                elif event.ui_element == step_button:
                    engine.step_forward()

        manager.process_events(event)

    screen.fill((255, 255, 255))
    panel.fill((169, 169, 169))
    fact = np.linalg.norm(engine.vel)
    norm_vel = (engine.vel / fact)
    for i in range(engine.config["N"]):
        px = int(round(engine.pos[i][0] * pxPerU))
        py = int(round(CANVASLENGTH - engine.pos[i][1] * pxPerU))
        # reverse y-direction of the arrow 
        vx = int(round(px+norm_vel[i][0]*MAGFACTOR*pxPerU))
        vy = int(round(py+norm_vel[i][1]*MAGFACTOR*pxPerU*-1))
        pygame.draw.circle(screen, (0, 0, 0), (px, py), RAD)
        pygame.draw.line(screen, (255, 0, 0), (px, py), (vx, vy), 3)
    
    if simulating:
        for i in range(5):
            engine.step_forward()
    info_block.set_text(f"FPS: {clock.get_fps():.1f}   Time: {engine.time:.2f}  Iterations: {engine.iterations}")
    screen.blit(panel, (830, 0))
    manager.update(time_delta)
    manager.draw_ui(screen)
    pygame.display.flip()

pygame.quit()