from libs.geometry.MSHReader import MSHReader
from libs.geometry.Grid import Grid
import time
from libs.simulation.HeatTransfer2D import HeatTransfer2D

# r = MSHReader("/home/gustavoe/Documents/Sinmec/GTRelated/GridReader/meshes/test.msh")#QuadPlate.msh")
# g = Grid(r.getData())

s = HeatTransfer2D("/home/gustavoe/Documents/Sinmec/GTRelated/GridReader/meshes/Square.msh")

# import numpy as np
# from matplotlib import pyplot as plt
# from mpl_toolkits import mplot3d
# from scipy.interpolate import griddata

# X,Y = zip(*[v.getCoordinates()[:-1] for v in s.grid.vertices])
# nT = s.numerical

# Xi, Yi = np.meshgrid( np.linspace(min(X), max(X), len(X)), np.linspace(min(Y), max(Y), len(Y)) )
# nTi = griddata((X,Y), nT, (Xi,Yi), method='linear')

# # ax1 = plt.axes(projection='3d')					# plt.pcolor(Xi,Yi,aTi, cmap='RdBu')
# # ax1.plot_surface(Xi,Yi,nTi, cmap='RdBu')		# plt.colorbar()
# plt.pcolor(Xi,Yi,nTi, cmap="RdBu")
# plt.title("Numerical")
# plt.show()

