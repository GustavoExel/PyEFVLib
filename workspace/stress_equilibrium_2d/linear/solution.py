from matplotlib import pyplot as plt, colors, cm
import os, json, pandas as pd, numpy as np

with open(os.path.join(os.path.dirname(__file__), "properties.json"), "r") as f:
	properties = json.load(f)
with open(os.path.join(os.path.dirname(__file__), "boundaryConditions/v.json"), "r") as f:
	v_condition = json.load(f)

resultsPath = "/home/gustavoe/Documents/Sinmec/HTRelated/PyEFVLib/results/stress_equilibrium/linear/Results.csv"
resultsData = pd.read_csv(resultsPath)

top_stress = v_condition["North"]["value"]
shearModulus = properties["Body"]["ShearModulus"]
poissonsRatio = properties["Body"]["PoissonsRatio"]
lameParameter = 2*shearModulus*poissonsRatio/(1-2*poissonsRatio)
density = properties["Body"]["Density"]
gravity = properties["Body"]["Gravity"]
height = 1.0

Y, X = resultsData["Y"], resultsData["X"]
fieldValues = resultsData[ resultsData.columns[-1] ]
Y, fieldValues = zip(*sorted([(y,val) for x,y,val in zip(X,Y,fieldValues) if abs(x-0.5)<0.1]))
Y, fieldValues = np.array(Y), np.array(fieldValues)

a_vals=Y*(top_stress+density*gravity*(height-Y/2.0))/(2.0*shearModulus+lameParameter)
plt.figure()
plt.scatter(Y,1000*fieldValues, marker='X', color='r', label="Resultados Numéricos")
plt.plot(Y,1000*a_vals, color='k', label="Solução Analítica")
plt.xlabel("Y (m)")
plt.ylabel("v (mm)")
plt.legend()	
plt.title("Deslocamento em y")

plt.show()
