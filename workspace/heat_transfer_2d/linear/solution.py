import os, sys, json, pandas as pd, matplotlib.pyplot as plt

with open(os.path.join(os.path.dirname(__file__), "properties.json"), "r") as f:
	properties = json.load(f)
a = 0.5 * properties["Body"]["HeatGeneration"] / properties["Body"]["Conductivity"]

def analyticalTemperature(t1, t2, a, x):
	return -a*x**2 + (t2-t1+a)*x + t1

results = os.path.join(os.path.dirname(__file__), *3*[os.path.pardir], "results", "heat_transfer_2d", "linear", "Results.csv")
print(results)
data = pd.read_csv(results)

X, Y= data["X"], data["Y"]
nT = data[data.columns[-1]]
aT = analyticalTemperature(300, 350, a, X)

X,nT,aT = zip(*sorted([(x,nt,at) for x,y,nt,at in zip(X,Y,nT,aT) if abs(y-0.5)<0.1]))

plt.scatter(X, nT, marker='X', color='r', label="Resultados Numéricos")
plt.plot(X, aT, color='k', label="Solução Analítica")

plt.ylabel("Temperature (K)")
plt.xlabel("X (m)")
plt.legend()
plt.show()