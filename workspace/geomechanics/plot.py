import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import json
import sys, os

dirname = os.path.realpath( os.path.dirname(__file__) )
sys.path.append(dirname + "/../../../../../solgeom")
from solgeom.Terzaghi import Solution

filePath = dirname+"/../../../../results/geomechanics/Results.csv"

df = pd.read_csv(filePath)
X = df["X"]
Y = df["Y"]

numberOfTimeSteps = int(df.keys()[-1].split(" - ")[0].replace("TimeStep",""))
numberOfSamples = 20

# stepInTimeSamples = np.round(np.linspace(1, numberOfTimeSteps-1, numberOfSamples)).astype(int)
# stepInTimeSamples = np.concatenate(( np.round(np.linspace(1, int(0.1*numberOfTimeSteps), numberOfSamples//3)).astype(int), np.round(np.linspace(1, numberOfTimeSteps-1, 2*numberOfSamples//3)).astype(int)[1:] ))
stepInTimeSamples = [1, 10, 50, 125, 250, 499]

mm = 1000.
kPa = 1/1000.

height = 6.0
load = 1.0e+5
gravity = 0.0
timeStep = 20

rock = json.load(open(dirname+"/../../../solgeom/examples/solid.json", "r"))
fluid = json.load(open(dirname+"/../../../solgeom/examples/fluid.json", "r"))

terza = Solution(height, load, rock, fluid , gravity)
z = terza.getPositionValues()

fig, ax = plt.subplots(2, 2, figsize=(8,7))
fig.subplots_adjust(left=0.070, right=0.975, top=0.970, bottom=0.065, hspace=0.235, wspace=0.300)

# Plot pressure profiles ---------------------------------------------------
times_1 = [N*timeStep for N in stepInTimeSamples]
for idx,t in zip(stepInTimeSamples,times_1):
	p = terza.getPressureValuesConstTime(t)
	p_n = df[f"TimeStep{idx} - p"]
	ax[0][0].plot(p*kPa, z, label=f"t = {timeStep*idx}s")
	ax[0][0].scatter(p_n*kPa, Y, color="k", marker=".",linewidth=0.5)
ax[0][0].set_xlabel("Pressure (kPa)", size=12)
ax[0][0].set_ylabel("Height (m)", size=12)
ax[0][0].grid(True)
# ax[0][0].legend()
# --------------------------------------------------------------------------

# Plot displacement profiles -----------------------------------------------
for idx,t in zip(stepInTimeSamples,times_1):
	w = terza.getDisplacementValuesConstTime(t)
	w_n = df[f"TimeStep{idx} - u_y"]
	ax[0][1].plot(w*mm, z)
	ax[0][1].scatter(w_n*mm, Y, color="k", marker=".",linewidth=0.5)
ax[0][1].set_xlabel("Displacement (mm)", size=12)
ax[0][1].set_ylabel("Height (m)", size=12)
ax[0][1].grid(True)
# --------------------------------------------------------------------------

# Plot bottom pressure over time -------------------------------------------
# times_2 = np.linspace(0, 2000., 100)
times_2 = [timeStep*s for s in range(1,1+numberOfTimeSteps)]
p_a = terza.getPressureValuesAtPosition(0.0, times_2) # Bottom pressure (z=0.0)
p_n = np.array([ np.average([p for p,y in zip(df[f"TimeStep{s} - p"],Y) if y<0.05]) for s in range(1,1+numberOfTimeSteps) ])
ax[1][0].plot(times_2, p_a*kPa)
ax[1][0].scatter(times_2, p_n*kPa, color="k", marker=".", linewidth=0.5)
ax[1][0].set_ylim((0, 55))
ax[1][0].set_xlabel("Time (s)", size=12)
ax[1][0].set_ylabel("Bottom Pressure (kPa)", size=12)
ax[1][0].grid(True)
# --------------------------------------------------------------------------

# Plot top displacement over time ------------------------------------------
w_a = terza.getDisplacementValuesAtPosition(height, times_2) # Top displacement (z=height)
w_n = np.array([ np.average([w for w,y in zip(df[f"TimeStep{s} - u_y"],Y) if y>height-0.05]) for s in range(1,1+numberOfTimeSteps) ])
ax[1][1].plot(times_2, w_a*mm)
ax[1][1].scatter(times_2, w_n*mm, color="k", marker=".", linewidth=0.5)
ax[1][1].set_xlabel("Time (s)", size=12)
ax[1][1].set_ylabel("Top Displacement (mm)", size=12)
ax[1][1].grid(True)
# --------------------------------------------------------------------------

plt.show()