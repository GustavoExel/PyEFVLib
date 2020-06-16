import sys,os
sys.path.append(os.path.join(os.path.dirname(__file__), os.path.pardir))

from PyEFVLib import MSHReader, Grid, ProblemData, CgnsSaver
import numpy as np

#-------------------------SETTINGS----------------------------------------------
problemData = ProblemData('heat_transfer_2d')

reader = MSHReader(problemData.paths["Grid"])
grid = Grid(reader.getData())
problemData.grid = grid
problemData.read()

timeStep = problemData.timeStep
currentTime = 0.0

cgnsSaver = CgnsSaver(grid, problemData.paths["Output"], problemData.libraryPath)

temperatureField = np.repeat(problemData.initialValue, grid.vertices.size)
prevTemperatureField = np.repeat(problemData.initialValue, grid.vertices.size)

matrix = np.zeros([grid.vertices.size, grid.vertices.size])
difference = 0.0
iteration = 0
converged = False

print("{:>9}\t{:>14}\t{:>14}\t{:>14}".format("Iteration", "CurrentTime", "TimeStep", "Difference"))
#-------------------------------------------------------------------------------
#-------------------------SIMULATION MAIN LOOP----------------------------------
#-------------------------------------------------------------------------------
while not converged and iteration < problemData.maxNumberOfIterations:
	#-------------------------ADD TO LINEAR SYSTEM------------------------------
	independent = np.zeros(grid.vertices.size)

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
						matrix[backwardVertexHandle][vertex.handle] += coefficient
						matrix[forwardVertexHandle][vertex.handle] -= coefficient
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
						matrix[vertex.handle][vertex.handle] += element.subelementVolumes[local] * accumulation
				local += 1

	# Neumann Boundary Condition
	for bCondition in problemData.neumannBoundaries:
		for facet in bCondition.boundary.facets:
			for outerFace in facet.outerFaces:
				independent[outerFace.vertex.handle] -= bCondition.getValue(outerFace.handle) * np.linalg.norm(outerFace.area.getCoordinates())

	# Dirichlet Boundary Condition
	for bCondition in problemData.dirichletBoundaries:
		for vertex in bCondition.boundary.vertices:
			independent[vertex.handle] = bCondition.getValue(vertex.handle)
	if iteration == 0:
		for bCondition in problemData.dirichletBoundaries:
			for vertex in bCondition.boundary.vertices:
				matrix[vertex.handle] = np.zeros(grid.vertices.size)
				matrix[vertex.handle][vertex.handle] = 1.0

	#-------------------------SOLVE LINEAR SYSTEM-------------------------------
	temperatureField = np.linalg.solve(matrix, independent)

	#-------------------------PRINT ITERATION DATA------------------------------
	if iteration > 0:
		print("{:>9}\t{:>14e}\t{:>14e}\t{:>14e}".format(iteration, currentTime, timeStep, difference))

	#-------------------------INCREMENT TIME------------------------------------
	currentTime += timeStep

	#-------------------------SAVE RESULTS--------------------------------------
	cgnsSaver.save('temperature field', temperatureField, currentTime)

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
cgnsSaver.finalize()

print("\n\t\033[1;35mresult:\033[0m", problemData.paths["Output"]+"Results.cgns", '\n')
#-------------------------------------------------------------------------------
#-------------------------SHOW RESULTS GRAPHICALY-------------------------------
#-------------------------------------------------------------------------------
from matplotlib import pyplot as plt, colors, cm
from scipy.interpolate import griddata

X,Y = zip(*[v.getCoordinates()[:-1] for v in grid.vertices])

Xi, Yi = np.meshgrid( np.linspace(min(X), max(X), len(X)), np.linspace(min(Y), max(Y), len(Y)) )
nTi = griddata((X,Y), temperatureField, (Xi,Yi), method='linear')

plt.pcolor(Xi,Yi,nTi, cmap=colors.ListedColormap( cm.get_cmap("RdBu",256)(np.linspace(1,0,256)) ))
# plt.pcolor(Xi,Yi,nTi, cmap="RdBu")
plt.title("Numerical Temperature")
plt.colorbar()	
plt.show()