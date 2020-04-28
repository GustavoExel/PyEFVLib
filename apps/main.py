import sys,os
sys.path.append(os.path.join(os.path.dirname(__file__), os.path.pardir))
from libs.simulation.HeatTransfer2D import HeatTransfer2D
from matplotlib import pyplot as plt, cm
from scipy.interpolate import griddata
import numpy as np

s = HeatTransfer2D()


X,Y = zip(*[v.getCoordinates()[:-1] for v in s.grid.vertices])
nT = s.numericalTemperature

Xi, Yi = np.meshgrid( np.linspace(min(X), max(X), len(X)), np.linspace(min(Y), max(Y), len(Y)) )
nTi = griddata((X,Y), nT, (Xi,Yi), method='linear')

plt.pcolor(Xi,Yi,nTi, cmap=CM( cm.get_cmap("RdBu",256)(np.linspace(1,0,256)) ))
plt.title("Numerical Temperature")
plt.colorbar()
plt.show()
