from engine import Engine
from config import Config
if __name__ == "__main__":
    c = Config()
    c.load_yaml("setup.yaml")

    e = Engine(c)

    for i in range(1, 6):
        print(f"step {i}")
        e.step_forward()
        e.debug()
