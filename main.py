from engine import Engine
from config import Config
if __name__ == "__main__":
    c = Config()
    c.load_yaml("setup.yaml")

    e = Engine(c)

    e.debug()
    e.step_forward()
    e.debug()
    
