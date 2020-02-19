from libs.simulation.HeatTransfer2D import HeatTransfer2D
from matplotlib import pyplot as plt
import h5py

s = HeatTransfer2D()
# Read benchmark (including solidProperties and boundaryConditions(for printing))

# Define the Analytical Solution
def analyticalSolution(x, y):
	# Considering NORTH and SOUTH isolated, 20K in the WEST, and 50K in the EAST
	return 30*x + 20

# Read Results.cgns
file = h5py.File(s.cgnsSaver.outputPath, "r+")
zone = file["BASE/ZONE"]
lastTS = max([int(key.strip("TimeStep")) for key in zone.keys() if "TimeStep" in key])

X = zone["GridCoordinates/CoordinateX/ data"][()]
Y = zone["GridCoordinates/CoordinateY/ data"][()]
nT = zone[f"TimeStep{str(lastTS)}/numerical temperature/ data"][()]


# Compute error
aT = [analyticalSolution(x,y) for x,y in zip(X,Y)]
volumes = [vertex.volume for vertex in s.grid.vertices]
euclideanError = sum([v * (at-nt)**2 for v, at, nt in zip(volumes, aT, nT)]) ** 0.5
maximumError = max([abs(at-nt) for at, nt in zip(aT, nT)])
print("\t{}e_{}v".format(s.grid.elements.size, s.grid.vertices.size))


# Plot
print("\tEuclidean Error : {:.5e}".format( euclideanError ) )
print("\tMaximum Error   : {:.5e}".format( maximumError ) )

data = list(zip(aT, nT, X))
data.sort(key=lambda triple:triple[2])
aT, nT, X = zip(*data)
del data

plt.plot(X, nT, label="Numerical Solution", color='k')
plt.scatter(X, aT, label="Analytical Solution", color='r', marker='X')
plt.xlabel("X (m)")
plt.ylabel("Temperature (K)")
plt.legend()
plt.show()
