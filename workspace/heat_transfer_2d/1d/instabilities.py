import sys,os
sys.path.append(os.path.join(os.path.dirname(__file__), *3*[os.path.pardir]))
from apps.heat_transfer import heatTransfer
from PyEFVLib import ProblemData, MSHReader, Grid
from matplotlib import pyplot as plt
import numpy as np

# Setting model
model = "workspace/heat_transfer_2d/1d"
if len(sys.argv)>1 and not "-" in sys.argv[1]: model=sys.argv[1]

problemData = ProblemData(model)

reader = MSHReader(problemData.meshFilePath)
grid = Grid(reader.getData())
problemData.setGrid(grid)
problemData.read()

# Calculating parameters
k, rho, cp = problemData.propertyData.get(0, "Conductivity"), problemData.propertyData.get(0, "Density"), problemData.propertyData.get(0, "HeatCapacity")
h2 = sum([element.volume for element in grid.elements])/grid.elements.size
alpha = k / (rho * cp)

print(f"h^2 = {h2}")
print(f"h^2/alpha = {h2/alpha}")


df=1
di=0.05
# df=0.01
# di=0.001
for i in np.arange(0.001,df,di):
	X,Y = zip(*[v.getCoordinates()[:-1] for v in grid.vertices])
	finalTemperatureField=heatTransfer(
		model 	  = model,
		libraryPath = problemData.libraryPath,
		outputPath = problemData.outputFilePath,
		fileName="Results",
		extension = "csv",

		grid 	  = grid,
		propertyData = problemData.propertyData,

		initialValues = problemData.initialValues,
		neumannBoundaries = problemData.neumannBoundaries,
		dirichletBoundaries = problemData.dirichletBoundaries,

		timeStep  = (h2/alpha) * i,
		finalTime = problemData.finalTime,
		maxNumberOfIterations = problemData.maxNumberOfIterations,
		tolerance = problemData.tolerance,
		verbosity = False
	)

	X, finalTemperatureField = zip(*sorted([(x,t) for x,t in zip(X, finalTemperatureField)]))
	plt.plot(X, finalTemperatureField, label=f"dt = (h^2/alpha) * {i:.04f}", color=(i/df,i/df,1-i/df))

plt.xlabel("X [m]")
plt.ylabel("Temperature [K]")
plt.legend()
plt.show()