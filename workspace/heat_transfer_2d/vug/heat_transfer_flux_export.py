import sys,os
sys.path.append(os.path.join(os.path.dirname(__file__), *[os.path.pardir]*3))
from PyEFVLib import MSHReader, Grid, ProblemData, CgnsSaver, CsvSaver
import numpy as np
import pandas as pd
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
		verbosity=True, 			# If False does not print iteration info
		color=True
	):
	#-------------------------SETTINGS----------------------------------------------
	initialTime = time.time()

	dimension = grid.dimension
	currentTime = 0.0

	savers = {"cgns": CgnsSaver, "csv": CsvSaver}
	saver = savers[extension](grid, outputPath, libraryPath, fileName=fileName)
	fluxDF = pd.DataFrame( {
		"Element": [elemIdx for elemIdx in range(grid.elements.size) for facetIdx in range(4)],
		"Facet": grid.elements.size*[0,1,2,3],
		"X1": [element.vertices[i].x for element in grid.elements for i in [0,1,2,3] ],
		"Y1": [element.vertices[i].y for element in grid.elements for i in [0,1,2,3] ],
		"X2": [element.vertices[i].x for element in grid.elements for i in [1,2,3,0] ],
		"Y2": [element.vertices[i].y for element in grid.elements for i in [1,2,3,0] ]

	} )
	fluxDF.to_csv("fluxes.csv", index=False)
	del fluxDF

	temperatureField = np.repeat(0.0, grid.vertices.size)
	prevTemperatureField = initialValues["temperature"].copy()

	coords,matrixVals = [], []
	difference = 0.0
	iteration = 0
	converged = False

	def print_purple(text, end="\n"):
		print(f"\n\t{text}", end=end)
	def add(i, j, val):
		coords.append((i,j))
		matrixVals.append(val)

	#-------------------------------------------------------------------------------
	#-------------------------SIMULATION MAIN LOOP----------------------------------
	#-------------------------------------------------------------------------------
	while not converged:
		if maxNumberOfIterations != None and iteration > maxNumberOfIterations:
			break
		if iteration > 1 and not transient:
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

		facetFlux = np.zeros((grid.elements.size,4,2))
		for region in grid.regions:
			conductivity = propertyData[region.handle]["Conductivity"]
			for elementIdx, element in enumerate(region.elements):
				for innerFaceIdx, innerFace in enumerate(element.innerFaces):
					temperatureVector = np.array([temperatureField[vertex.handle] for vertex in element.vertices])

					heatFlux = -conductivity * np.matmul( innerFace.globalDerivatives, temperatureVector )

					facetFlux[elementIdx][innerFaceIdx] = heatFlux
		fluxDF = pd.read_csv("fluxes.csv")
		fluxDF[f"time_step_{iteration} - q''x"] = np.array([facetData[0] for elementData in facetFlux for facetData in elementData])
		fluxDF[f"time_step_{iteration} - q''y"] = np.array([facetData[1] for elementData in facetFlux for facetData in elementData])
		fluxDF[f"time_step_{iteration} - q''"] = np.array([np.linalg.norm(facetData) for elementData in facetFlux for facetData in elementData])
		fluxDF.to_csv("fluxes.csv", index=False)
		del fluxDF

		#-------------------------PRINT ITERATION DATA------------------------------
		if iteration > 0 and verbosity:
			print("{:>9}\t{:>14.2e}\t{:>14.2e}\t{:>14.2e}".format(iteration, currentTime, timeStep, difference))

		#-------------------------INCREMENT TIME------------------------------------
		currentTime += timeStep

		#-------------------------SAVE RESULTS--------------------------------------
		saver.save("temperature", temperatureField, currentTime)

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

		[print, print_purple][color]("\n\tresult: ", end="")
		print(os.path.realpath(saver.outputPath), "\n")

	return temperatureField

if __name__ == "__main__":
	model = "workspace/heat_transfer_2d/vug"
	problemData = ProblemData(model)
	reader = MSHReader(problemData.paths["Grid"])
	grid = Grid(reader.getData())
	problemData.setGrid(grid)
	problemData.read()

	K = [1e3, 1e12, 1e22]
	for k in K:
	# if True:
	# 	k=1
		problemData.propertyData[0]["Conductivity"] = 1
		problemData.propertyData[1]["Conductivity"] = k
		
		# open("fluxes.csv", "w").close()

		finalTemperatureField = heatTransfer(
			libraryPath = problemData.libraryPath,
			outputPath = problemData.paths["Output"],
			extension = "csv" if not "--extension=cgns" in sys.argv else "cgns",
			
			grid 	  = grid,
			propertyData = problemData.propertyData,
			
			# initialValues = problemData.initialValues,
			initialValues = { "temperature": np.repeat( problemData.initialValues["temperature"], grid.vertices.size ) },
			neumannBoundaries = problemData.neumannBoundaries,
			dirichletBoundaries = problemData.dirichletBoundaries,

			timeStep  = problemData.timeStep,
			finalTime = problemData.finalTime,
			maxNumberOfIterations = problemData.maxNumberOfIterations,
			tolerance = problemData.tolerance,
			
			transient = not "-p" in sys.argv,
			verbosity = not "-s" in sys.argv,
		)

		os.rename("fluxes.csv", f"fluxos - vec\\results\\k1=1, k2={k:.0e} - 30x30.csv")