import os,json

r = os.path.realpath(os.path.dirname(__file__))
with open(os.path.join(r,"boundaryConditions","temperature.json"), "r") as f:
	t = json.load(f)
with open(os.path.join(r, "properties.json"), "r") as f:
	p = json.load(f)

t1 = t["WEST"]["value"]
t2 = t["EAST"]["value"]
k  = p["BODY"]["Conductivity"]
q  = p["BODY"]["HeatGeneration"]

def analyticalSolution_X(x):
    return -0.5*x*x*q/k + (t2-t1 + 0.5*q/k)*x + t1
