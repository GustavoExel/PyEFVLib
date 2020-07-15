import sys,os
sys.path.append(os.path.join(os.path.dirname(__file__), os.path.pardir))

from PyEFVLib import MSHReader, Grid, ProblemData, CsvSaver
import numpy as np
from scipy import sparse
import scipy.sparse.linalg

#-------------------------SETTINGS----------------------------------------------
problemData = ProblemData('heat_transfer_2d')

reader = MSHReader(problemData.paths["Grid"])
grid = Grid(reader.getData())
problemData.grid = grid
problemData.read()

numberOfVertices = grid.vertices.size
timeStep = problemData.timeStep
currentTime = 0.0

saver = CsvSaver(grid, problemData.paths["Output"], problemData.libraryPath)

temperatureField = np.repeat(problemData.initialValue, numberOfVertices)
prevTemperatureField = np.repeat(problemData.initialValue["temperature"], numberOfVertices)

coords,matrixVals = [], []

difference = 0.0
iteration = 0
converged = False

def add(i, j, val):
	global coords, matrixVals
	coords.append((i,j))
	matrixVals.append(val)

print("{:>9}\t{:>14}\t{:>14}\t{:>14}".format("Iteration", "CurrentTime", "TimeStep", "Difference"))
#-------------------------------------------------------------------------------
#-------------------------SIMULATION MAIN LOOP----------------------------------
#-------------------------------------------------------------------------------
while not converged and iteration < problemData.maxNumberOfIterations:
	#-------------------------ADD TO LINEAR SYSTEM------------------------------
	independent = np.zeros(numberOfVertices)

	# Generation Term
	for region in grid.regions:
		heatGeneration = problemData.propertyData[region.handle]["HeatGeneration"]
		for element in region.elements:
			local = 0
			for vertex in element.vertices:
				independent[vertex.handle] += element.subelementVolumes[local] * heatGeneration
				local += 1

	# Diffusion Term
	if iteration == 0:
		for region in grid.regions:
			conductivity = problemData.propertyData[region.handle]["Conductivity"]
			for element in region.elements:
				for innerFace in element.innerFaces:
					diffusiveFlux = conductivity * np.matmul( np.transpose(innerFace.globalDerivatives) , innerFace.area.getCoordinates()[:-1] )
					backwardVertexHandle = element.vertices[element.shape.innerFaceNeighborVertices[innerFace.local][0]].handle
					forwardVertexHandle = element.vertices[element.shape.innerFaceNeighborVertices[innerFace.local][1]].handle
					
					i=0
					for vertex in element.vertices:
						coefficient = -1.0 * diffusiveFlux[i]
						add(backwardVertexHandle, vertex.handle, coefficient)
						add(forwardVertexHandle, vertex.handle, -coefficient)
						i+=1

	# Transient Term
	for region in grid.regions:
		density = problemData.propertyData[region.handle]["Density"]
		heatCapacity = problemData.propertyData[region.handle]["HeatCapacity"]
		accumulation = density * heatCapacity / timeStep

		for element in region.elements:
			local = 0
			for vertex in element.vertices:
				independent[vertex.handle] += element.subelementVolumes[local] * accumulation * prevTemperatureField[vertex.handle]
				if iteration == 0:
					add(vertex.handle, vertex.handle, element.subelementVolumes[local] * accumulation)						
				local += 1

	# Neumann Boundary Condition
	for bCondition in problemData.neumannBoundaries["temperature"]:
		for facet in bCondition.boundary.facets:
			for outerFace in facet.outerFaces:
				independent[outerFace.vertex.handle] -= bCondition.getValue(outerFace.handle) * np.linalg.norm(outerFace.area.getCoordinates())

	# Dirichlet Boundary Condition
	for bCondition in problemData.dirichletBoundaries["temperature"]:
		for vertex in bCondition.boundary.vertices:
			independent[vertex.handle] = bCondition.getValue(vertex.handle)
	if iteration == 0:
		for bCondition in problemData.dirichletBoundaries["temperature"]:
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
	if iteration > 0:
		print("{:>9}\t{:>14e}\t{:>14e}\t{:>14e}".format(iteration, currentTime, timeStep, difference))

	#-------------------------INCREMENT TIME------------------------------------
	currentTime += timeStep

	#-------------------------SAVE RESULTS--------------------------------------
	saver.save('temperature field', temperatureField, currentTime)

	#-------------------------CHECK CONVERGENCE---------------------------------
	converged = False
	difference = max([abs(temp-oldTemp) for temp, oldTemp in zip(temperatureField, prevTemperatureField)])
	prevTemperatureField = temperatureField
	if currentTime > problemData.finalTime:
		converged = True
	elif iteration > 0:
		converged = difference < problemData.tolerance

	#-------------------------INCREMENT ITERATION-------------------------------
	iteration += 1   

#-------------------------------------------------------------------------------
#-------------------------AFTER END OF MAIN LOOP ITERATION------------------------
#-------------------------------------------------------------------------------
saver.finalize()

print("\n\t\033[1;35mresult:\033[0m", saver.outputPath, '\n')
#-------------------------------------------------------------------------------
#-------------------------SHOW RESULTS GRAPHICALY-------------------------------
#-------------------------------------------------------------------------------
from matplotlib import pyplot as plt, colors, cm
from scipy.interpolate import griddata

X,Y = zip(*[v.getCoordinates()[:-1] for v in grid.vertices])

Xi, Yi = np.meshgrid( np.linspace(min(X), max(X), len(X)), np.linspace(min(Y), max(Y), len(Y)) )
nTi = griddata((X,Y), temperatureField, (Xi,Yi), method='linear')

plt.pcolor(Xi,Yi,nTi, cmap=colors.ListedColormap( cm.get_cmap("RdBu",256)(np.linspace(1,0,256)) ))
plt.title("Numerical Temperature")
plt.colorbar()	
plt.show()