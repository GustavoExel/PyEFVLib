import numpy as np
import scipy.interpolate
import matplotlib.pyplot as plt
from matplotlib import cm
from matplotlib.colors import ListedColormap as CM, Normalize

def colorplot(X, Y, D):
	Xi, Yi = np.meshgrid( np.linspace(min(X), max(X), len(X)), np.linspace(min(Y), max(Y), len(Y)) )
	Di = scipy.interpolate.griddata((X,Y), D, (Xi,Yi), method="linear")
	plt.pcolormesh(Xi,Yi,Di, shading="auto", cmap=CM( cm.get_cmap("RdBu",64)(np.linspace(1,0,64)) )) # Makes BuRd instead of RdBu
	plt.colorbar()
	plt.show()

if __name__ == "__main__":
	X = list( np.random.uniform(0,1,25**2) ) + [0,1,1,0]
	Y = list( np.random.uniform(0,1,25**2) ) + [0,0,1,1]
	D = np.array([x**2+(y-0.5)**2 for (x,y) in zip (X,Y)])
	colorplot(X,Y,D)
