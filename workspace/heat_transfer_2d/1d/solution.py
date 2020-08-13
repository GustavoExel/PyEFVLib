import os, sys, json, pandas as pd, matplotlib.pyplot as plt

results = os.path.join(os.path.dirname(__file__), *3*[os.path.pardir], "results", "heat_transfer_2d", "1d", "Results.csv")
data = pd.read_csv(results)

X, Y= data["X"], data["Y"]
nT = data[data.columns[-1]]

X,nT=zip(*sorted([(x,nt) for x,nt in zip(X,nT)]))
# X,nT = zip(*sorted([(x,nt) for x,y,nt in zip(X,Y,nT) if abs(y-0.5)<0.1]))	

plt.plot(X, nT, marker='X', color='r', label="Resultados Numéricos")

plt.ylabel("Temperature (K)")
plt.xlabel("X (m)")
plt.legend()
plt.show()