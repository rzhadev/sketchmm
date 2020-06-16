import yaml
import os
# natural units of
# epsilon=sigma=kB=m=1, dt=1microsecond, using natural units


class Config():

    def __init__(self, s=1, e=1, k=1, mass=1, dt=1e-4):
        self.mapping = {"parameters": {"sigma": s, "epsilon": e,
                                       "kB": k, "m": mass,
                                       "timeStep": dt},
                        "data": {
                            "atom1": {"pos": [0.0, 0.0], "vel": [0.0, 0.0], "acc": [0.0, 0.0]},  # noqa
                            "atom2": {"pos": [10.0, 10.0], "vel": [0.0, 0.0], "acc": [0.0, 0.0]},  # noqa
                        }
        }

    def to_yaml(self, filename):
        with open(filename, 'w') as f:
            try:
                data = yaml.dump(self.mapping, f)
                print(f"Saved config to {filename}")
            except ImportError:
                print("Failed saving config file")

    def load_yaml(self, filename):
        with open(filename, 'r') as f:
            try:
                data = yaml.load(f, Loader=yaml.SafeLoader)
                self.mapping = data
                print(f"Loaded data from {filename}")
            except ImportError:
                print("Failed loading config file")


if __name__ == "__main__":
    con = Config()
    con.to_yaml('setup.yaml')