import sys,os
sys.path.append(os.path.join(os.path.dirname(__file__), *3*[os.path.pardir]))
from apps.heat_transfer import heatTransfer
from PyEFVLib import ProblemData, MSHReader, Grid

model = "workspace/heat_transfer_2d/1d"
if len(sys.argv)>1 and not "-" in sys.argv[1]: model=sys.argv[1]

problemData = ProblemData(model)

reader = MSHReader(problemData.paths["Grid"])
grid = Grid(reader.getData())
problemData.setGrid(grid)
problemData.read()

heatTransfer(
	model 	  = model,
	extension = "csv",
	grid 	  = grid,
	timeStep  = problemData.timeStep,
	outputPath = problemData.paths["Output"],
	libraryPath = problemData.libraryPath,
	initialValues = problemData.initialValues,
	maxNumberOfIterations = problemData.maxNumberOfIterations,
	propertyData = problemData.propertyData,
	neumannBoundaries = problemData.neumannBoundaries,
	dirichletBoundaries = problemData.dirichletBoundaries,
	finalTime = problemData.finalTime,
	tolerance = problemData.tolerance,
	fileName="Results"
)
