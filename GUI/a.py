import sys,os
sys.path.append(os.path.join(os.path.dirname(__file__), os.path.pardir))
from PyEFVLib import MSHReader, Grid, ProblemData, CgnsSaver, CsvSaver
import numpy as np, pandas as pd
from scipy import sparse
import scipy.sparse.linalg
from matplotlib import pyplot as plt
import time

model = "workspace/heat_transfer_2d/linear"

problemData = ProblemData(model)

reader = MSHReader(problemData.paths["Grid"])
grid = Grid(reader.getData())
problemData.setGrid(grid)
problemData.read()

libraryPath = problemData.libraryPath
outputPath = problemData.paths["Output"]

propertyData = problemData.propertyData

# initialValues = problemData.initialValues,
initialValues = { "temperature": np.repeat( problemData.initialValues["temperature"], grid.vertices.size ) }
neumannBoundaries = problemData.neumannBoundaries
dirichletBoundaries = problemData.dirichletBoundaries

timeStep  = problemData.timeStep
finalTime = problemData.finalTime
maxNumberOfIterations = problemData.maxNumberOfIterations
tolerance = problemData.tolerance

transient = False

#-------------------------SETTINGS----------------------------------------------
initialTime = time.time()

dimension = grid.dimension
currentTime = 0.0

saver = CsvSaver(grid, outputPath, libraryPath, fileName="Results")

temperatureField = np.repeat(0.0, grid.vertices.size)
prevTemperatureField = initialValues["temperature"].copy()

coords,matrixVals = [], []
difference = 0.0
iteration = 0
converged = False

def print_purple(text, end="\n"):
	print(f"\n\t\033[1;35m{text}\033[0m", end=end)
def add(i, j, val):
	# VAL = round(val,2) + 1000
	VAL = val
	coords.append((i,j))
	matrixVals.append(VAL)

#-------------------------------------------------------------------------------
#-------------------------SIMULATION MAIN LOOP----------------------------------
#-------------------------------------------------------------------------------
it = 0
while not converged and it == 0:
	# print("it = ", it)
	it += 1
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
		M = matrix.toarray()
		for i,line in enumerate(M):
			if not [x for x in line if x != 0.0 and abs(x) > 0.1]:
				print(i,line)
		inverseMatrix = sparse.linalg.inv( matrix )
	temperatureField = inverseMatrix * independent

	#-------------------------INCREMENT TIME------------------------------------
	currentTime += timeStep

	#-------------------------SAVE RESULTS--------------------------------------
	# saver.timeSteps	= np.append(saver.timeSteps,  currentTime)
	# saver.fields  = np.vstack([saver.fields, temperatureField])
	saver.save("temperature", temperatureField, currentTime)

	#-------------------------CHECK CONVERGENCE---------------------------------
	converged = False
	difference = max([abs(temp-oldTemp) for temp, oldTemp in zip(temperatureField, prevTemperatureField)])
	print(difference)
	prevTemperatureField = temperatureField
	if finalTime != None and currentTime > finalTime:
		converged = True
	elif iteration > 0 and tolerance != None:
		converged = difference < tolerance
	#-------------------------INCREMENT ITERATION-------------------------------
	iteration += 1   


d = pd.DataFrame(matrix.toarray())
d["independent"] = independent
d.to_csv("temp.csv")

#-------------------------------------------------------------------------------
#-------------------------AFTER END OF MAIN LOOP ITERATION------------------------
#-------------------------------------------------------------------------------
saver.finalize()
