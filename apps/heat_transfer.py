import sys,os
sys.path.append(os.path.join(os.path.dirname(__file__), os.path.pardir))
from PyEFVLib import MSHReader, Grid, ProblemData, CgnsSaver, CsvSaver
import numpy as np
from scipy import sparse
import scipy.sparse.linalg
import time

def heatTransfer(
			libraryPath,			# PyEFVLib path
			outputPath,				# Results directory path (Ex.: "results/heat_transfer_2d/...")
			extension,				# Extension type. Either "csv" or "cgns"

			grid,					# Object of class Grid
			propertyData,			# List of dictionaries containing the properties

			initialValues,			# Dictionary whose keys are the field names, and values are the field values
			neumannBoundaries,		# Dictionary whose keys are the field names, and values are objects of the class NeumannBoundaryCondition
			dirichletBoundaries,	# Dictionary whose keys are the field names, and values are objects of the class DirichletBoundaryCondition

			timeStep,				# Floating point number indicating the timeStep used in the simulation (constant)
			finalTime,				# The time at which, if reached, the simulation stops. If None, then it is not used.
			maxNumberOfIterations,	# Number of iterations at which, if reached, the simulation stops. If None, then it is not used.
			tolerance,				# The value at which, if the maximum difference between field values reach, the simulation stops. If None, then it is not used.
			
			fileName="Results",		# File name
			transient=True,			# If False, the transient term is not added to the equation, and it's solved in one iteration
			verbosity=True 			# If False does not print iteration info
			):

	#-------------------------SETTINGS----------------------------------------------
	print(f"timeStep = {timeStep}")
	print(f"finalTime = {finalTime}")
	print(f"maxNumberOfIterations = {maxNumberOfIterations}")
	print(f"tolerance = {tolerance}")
	print(f"transient = {transient}")

	initialTime = time.time()

	dimension = grid.dimension
	currentTime = 0.0

	savers = {"cgns": CgnsSaver, "csv": CsvSaver}
	saver = savers[extension](grid, outputPath, libraryPath, fileName=fileName)

	temperatureField = np.repeat(initialValues["temperature"], grid.vertices.size)
	prevTemperatureField = np.repeat(initialValues["temperature"], grid.vertices.size)

	coords,matrixVals = [], []
	difference = 0.0
	iteration = 0
	converged = False

	def add(i, j, val):
		coords.append((i,j))
		matrixVals.append(val)

	#-------------------------------------------------------------------------------
	#-------------------------SIMULATION MAIN LOOP----------------------------------
	#-------------------------------------------------------------------------------
	while not converged:
		if maxNumberOfIterations != None and iteration > maxNumberOfIterations:
			break
		#-------------------------ADD TO LINEAR SYSTEM------------------------------
		independent = np.zeros(grid.vertices.size)

		# Generation Term
		for region in grid.regions:
			heatGeneration = propertyData[region.handle]["HeatGeneration"]
			for element in region.elements:
				local = 0
				for vertex in element.vertices:
					independent[vertex.handle] += element.subelementVolumes[local] * heatGeneration
					local += 1

		# Diffusion Term
		if iteration == 0:
			for region in grid.regions:
				conductivity = propertyData[region.handle]["Conductivity"]
				for element in region.elements:
					for innerFace in element.innerFaces:
						diffusiveFlux = conductivity * np.matmul( np.transpose(innerFace.globalDerivatives) , innerFace.area.getCoordinates()[:dimension] )
						backwardVertexHandle = element.vertices[element.shape.innerFaceNeighborVertices[innerFace.local][0]].handle
						forwardVertexHandle = element.vertices[element.shape.innerFaceNeighborVertices[innerFace.local][1]].handle
						
						i=0
						for vertex in element.vertices:
							coefficient = -1.0 * diffusiveFlux[i]
							add(backwardVertexHandle, vertex.handle, coefficient)
							add(forwardVertexHandle, vertex.handle, -coefficient)
							i+=1

		# Transient Term
		if transient:	# If user knows that the accumulation term is irrelevant to the problem
			for region in grid.regions:
				density = propertyData[region.handle]["Density"]
				heatCapacity = propertyData[region.handle]["HeatCapacity"]
				accumulation = density * heatCapacity / timeStep

				for element in region.elements:
					local = 0
					for vertex in element.vertices:
						independent[vertex.handle] += element.subelementVolumes[local] * accumulation * prevTemperatureField[vertex.handle]
						if iteration == 0:
							add(vertex.handle, vertex.handle, element.subelementVolumes[local] * accumulation)						
						local += 1

		# Neumann Boundary Condition
		for bCondition in neumannBoundaries["temperature"]:
			for facet in bCondition.boundary.facets:
				for outerFace in facet.outerFaces:
					independent[outerFace.vertex.handle] -= bCondition.getValue(outerFace.handle) * np.linalg.norm(outerFace.area.getCoordinates())

		# Dirichlet Boundary Condition
		for bCondition in dirichletBoundaries["temperature"]:
			for vertex in bCondition.boundary.vertices:
				independent[vertex.handle] = bCondition.getValue(vertex.handle)
		if iteration == 0:
			for bCondition in dirichletBoundaries["temperature"]:
				for vertex in bCondition.boundary.vertices:
					matrixVals = [val for coord, val in zip(coords, matrixVals) if coord[0] != vertex.handle]
					coords 	   = [coord for coord in coords if coord[0] != vertex.handle]
					add(vertex.handle, vertex.handle, 1.0)

		#-------------------------SOLVE LINEAR SYSTEM-------------------------------
		if iteration == 0:
			matrix = sparse.coo_matrix( (matrixVals, zip(*coords)) )
			matrix = sparse.csc_matrix( matrix )
			inverseMatrix = sparse.linalg.inv( matrix )
		temperatureField = inverseMatrix * independent

		#-------------------------PRINT ITERATION DATA------------------------------
		if iteration > 0 and verbosity:
			print("{:>9}\t{:>14.2e}\t{:>14.2e}\t{:>14.2e}".format(iteration, currentTime, timeStep, difference))

		#-------------------------INCREMENT TIME------------------------------------
		currentTime += timeStep

		#-------------------------SAVE RESULTS--------------------------------------
		# saver.timeSteps	= np.append(saver.timeSteps,  currentTime)
		# saver.fields  = np.vstack([saver.fields, temperatureField])
		saver.save("temperature field", temperatureField, currentTime)

		#-------------------------CHECK CONVERGENCE---------------------------------
		converged = False
		difference = max([abs(temp-oldTemp) for temp, oldTemp in zip(temperatureField, prevTemperatureField)])
		prevTemperatureField = temperatureField
		if finalTime != None and currentTime > finalTime:
			converged = True
		elif iteration > 0 and tolerance != None:
			converged = difference < tolerance
		#-------------------------INCREMENT ITERATION-------------------------------
		iteration += 1   


	#-------------------------------------------------------------------------------
	#-------------------------AFTER END OF MAIN LOOP ITERATION------------------------
	#-------------------------------------------------------------------------------
	finalSimulationTime = time.time()
	if verbosity:
		print("Ended Simultaion, elapsed {:.2f}s".format(finalSimulationTime-initialTime))

	saver.finalize()
	if verbosity:
		print("Saved file: elapsed {:.2f}s".format(time.time()-finalSimulationTime))

		print("\n\t\033[1;35mresult:\033[0m", saver.outputPath, "\n")

	return temperatureField

if __name__ == "__main__":
	if "--help" in sys.argv:
		print("\npython apps/heat_transfer_2d.py workspace_file for opening a described model in workspace\n")
		print("-p\t for permanent regime (without the accumulation term)")
		print("-g\t for show results graphicaly")
		print("-s\t for verbosity 0")
		print("-1d\t compare 1d analytical with numerical solution along a graph")
		print("-2d\t show 2d analytical solution colorplot. Useful for really discrepant differences")
		print("--extension=cgns for output file in cgns extension\n")
		exit(0)
	
	model = "workspace/heat_transfer_2d/linear"
	if len(sys.argv)>1 and not "-" in sys.argv[1]: model=sys.argv[1]

	problemData = ProblemData(model)

	reader = MSHReader(problemData.paths["Grid"])
	grid = Grid(reader.getData())
	problemData.setGrid(grid)
	problemData.read()

	if not "-s" in sys.argv:
		for key,path in zip( ["input", "output", "grids"] , [os.path.join(problemData.libraryPath,"workspace",model) , problemData.paths["Output"], problemData.paths["Grid"]] ):
			print("\t\033[1;35m{}\033[0m\n\t\t{}\n".format(key, path))
		print("\t\033[1;35msolid\033[0m")
		for region in grid.regions:
			print("\t\t\033[36m{}\033[0m".format(region.name))
			for _property in problemData.propertyData[region.handle].keys():
				print("\t\t\t{}   : {}".format(_property, problemData.propertyData[region.handle][_property]))
			print("")
		print("\n{:>9}\t{:>14}\t{:>14}\t{:>14}".format("Iteration", "CurrentTime", "TimeStep", "Difference"))

	finalTemperatureField = heatTransfer(
		libraryPath = problemData.libraryPath,
		outputPath = problemData.paths["Output"],
		extension = "csv" if not "--extension=cgns" in sys.argv else "cgns",
		
		grid 	  = grid,
		propertyData = problemData.propertyData,
		
		initialValues = problemData.initialValues,
		neumannBoundaries = problemData.neumannBoundaries,
		dirichletBoundaries = problemData.dirichletBoundaries,

		timeStep  = problemData.timeStep,
		finalTime = problemData.finalTime,
		maxNumberOfIterations = problemData.maxNumberOfIterations,
		tolerance = problemData.tolerance,
		
		transient = not "-p" in sys.argv,
		verbosity = not "-s" in sys.argv
	)

	#-------------------------------------------------------------------------------
	#-------------------------SHOW RESULTS GRAPHICALY-------------------------------
	#-------------------------------------------------------------------------------
	if "-g" in sys.argv:
		import matplotlib
		from matplotlib import pyplot as plt
		from matplotlib import cm
		from matplotlib.colors import ListedColormap as CM, Normalize
		from scipy.interpolate import griddata

		X,Y = zip(*[v.getCoordinates()[:-1] for v in grid.vertices])

		Xi, Yi = np.meshgrid( np.linspace(min(X), max(X), len(X)), np.linspace(min(Y), max(Y), len(Y)) )
		nTi = griddata((X,Y), finalTemperatureField, (Xi,Yi), method="linear")

		plt.pcolor(Xi,Yi,nTi, cmap=CM( cm.get_cmap("RdBu",64)(np.linspace(1,0,64)) )) # Makes BuRd instead of RdBu
		plt.title("Numerical Temperature")
		plt.colorbar()
		plt.show()

	if "--paraview" in sys.argv:
		os.system(f"/usr/bin/paraview {saver.outputPath}")