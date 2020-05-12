import numpy as np
import os,json

r = os.path.realpath(os.path.dirname(__file__))
with open(os.path.join(r,"boundaryConditions/temperature.json"), "r") as f:
	t = json.load(f)

t1 = t["SOUTH"]["value"]
t2 = t["NORTH"]["value"]

# ΔT sinh(pi y) sinh(pi x) + T1
def analyticalSolution_XY(x,y):
    return (t2-t1)*np.sinh(np.pi*y)*np.sin(np.pi*x)/np.sinh(np.pi) + t1
