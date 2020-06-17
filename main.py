from engine import Engine
from config import Config
if __name__ == "__main__":
    e = Engine("setup.yaml")
    for i in range(100):
        print(f"step {i}")
        e.step_forward()
        e.compute_stats()
        e.debug()
        print('\n')
