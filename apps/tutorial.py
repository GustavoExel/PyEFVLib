import sys,os
sys.path.append(os.path.join(os.path.dirname(__file__), os.path.pardir))
from libs.geometry.MSHReader import MSHReader
from libs.geometry.Grid import Grid
from libs.simulation.ProblemData import ProblemData
from libs.simulation.CgnsSaver import CgnsSaver
import numpy as np
import time

#-------------------------SETTINGS----------------------------------------------
initialTime = time.time()
problemData = ProblemData('heat_transfer_1d')

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


for key,path in zip( ["input", "output", "grids"] , [problemData.libraryPath+"/benchmark/heat_transfer_2d/" , problemData.paths["Output"], problemData.paths["Grid"]] ):
	print(f"\t\033[1;35m{key}\033[0m\n\t\t{path}\n")
print(f"\t\033[1;35msolid\033[0m")
for region in grid.regions:
	print(f"\t\t\033[36m{region.name}\033[0m")
	for _property in problemData.propertyData[region.handle].keys():
		print(f"\t\t\t{_property}   : {problemData.propertyData[region.handle][_property]}")
	print("")
print("\n{:>9}\t{:>14}\t{:>14}\t{:>14}".format("Iteration", "CurrentTime", "TimeStep", "Difference"))

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
			for vertex in element.vertices:
				independent[vertex.handle] = vertex.volume * heatGeneration

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

		for vertex in grid.vertices:
			independent[vertex.handle] += vertex.volume * accumulation * prevTemperatureField[vertex.handle]
			if iteration == 0:
					matrix[vertex.handle][vertex.handle] += vertex.volume * accumulation

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
	cgnsSaver.timeSteps	= np.append(cgnsSaver.timeSteps,  currentTime)
	cgnsSaver.timeFields  = np.vstack([cgnsSaver.timeFields, temperatureField])

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
finalSimulationTime = time.time()
print("Ended Simultaion, elapsed {:.2f}s".format(finalSimulationTime-initialTime))

cgnsSaver.finalize()
print("Saved file: elapsed {:.2f}s".format(time.time()-finalSimulationTime))

print("\n\t\033[1;35mresult:\033[0m", problemData.paths["Output"]+"Results.cgns", '\n')
#-------------------------------------------------------------------------------
#-------------------------SHOW RESULTS GRAPHICALY-------------------------------
#-------------------------------------------------------------------------------
import matplotlib
from matplotlib import pyplot as plt
from matplotlib import cm
from matplotlib.colors import ListedColormap as CM, Normalize
from scipy.interpolate import griddata
from workspace.heat_transfer_1d.analyticalSolution import analyticalSolution_X

X,Y = zip(*[v.getCoordinates()[:-1] for v in grid.vertices])

Xi, Yi = np.meshgrid( np.linspace(min(X), max(X), len(X)), np.linspace(min(Y), max(Y), len(Y)) )
nTi = griddata((X,Y), temperatureField, (Xi,Yi), method='linear')

plt.pcolor(Xi,Yi,nTi, cmap=CM( cm.get_cmap("RdBu",256)(np.linspace(1,0,256)) )) # Makes BuRd instead of RdBu
plt.title("Numerical Temperature")
plt.colorbar()	

plt.figure()

x=np.linspace(0,1,15)
X,T = zip(*sorted(zip(X,temperatureField), key=lambda t:t[0]))
plt.plot(X,T, color='k', label='Numerical Results')
plt.scatter(x,[analyticalSolution_X(xx) for xx in x], marker='X', color='r', label='Analytical Results')
plt.xlabel("X")
plt.ylabel("Temperature")
plt.legend()

plt.show()

err = [abs(T - analyticalSolution_X(v.getCoordinates()[0])) for v,T in zip(grid.vertices,temperatureField)]
print("Max Error: {:.3e}".format(max(err)))
print("Mean Error: {:.3e}".format(sum(err)/len(err)))
print("Euclidean Error: {:.3e}".format(np.sqrt(sum([v.volume*e**2 for v,e in zip(grid.vertices, err)]))))