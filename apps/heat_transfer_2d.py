import sys,os
sys.path.append(os.path.join(os.path.dirname(__file__), os.path.pardir))
from PyEFVLib import MSHReader, Grid, ProblemData, CgnsSaver
import numpy as np
import time

if '--help' in sys.argv:
	print('\npython main.py workspace_file for opening a described model in workspace\n')
	print('-p\t for permanent regime (without the accumulation term)')
	print('-g\t for show results graphicaly')
	print('-s\t for verbosity 0')
	print('-1d\t compare 1d analytical with numerical solution along a graph')
	print('-2d\t show 2d analytical solution colorplot. Useful for really discrepant differences\n')
	exit(0)


model = 'sine_distribution'
if len(sys.argv)>1 and not '-' in sys.argv[1]: model=sys.argv[1]
#-------------------------SETTINGS----------------------------------------------
initialTime = time.time()
problemData = ProblemData(model)

reader = MSHReader(problemData.paths["Grid"])
grid = Grid(reader.getData())
problemData.grid = grid
problemData.read()

timeStep = problemData.timeStep
currentTime = 0.0

cgnsSaver = CgnsSaver(grid, problemData.paths["Output"], problemData.libraryPath, ['temperature field'])

temperatureField = np.repeat(problemData.initialValue["temperature"], grid.vertices.size)
prevTemperatureField = np.repeat(problemData.initialValue["temperature"], grid.vertices.size)

matrix = np.zeros([grid.vertices.size, grid.vertices.size])
difference = 0.0
iteration = 0
converged = False

if not '-s' in sys.argv:
	for key,path in zip( ["input", "output", "grids"] , [problemData.libraryPath+"/workspace/"+model , problemData.paths["Output"], problemData.paths["Grid"]] ):
		print("\t\033[1;35m{}\033[0m\n\t\t{}\n".format(key, path))
	print("\t\033[1;35msolid\033[0m")
	for region in grid.regions:
		print("\t\t\033[36m{}\033[0m".format(region.name))
		for _property in problemData.propertyData[region.handle].keys():
			print("\t\t\t{}   : {}".format(_property, problemData.propertyData[region.handle][_property]))
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
	if not '-p' in sys.argv:	# If user knows that the accumulation term is irrelevant to the problem
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
				matrix[vertex.handle] = np.zeros(grid.vertices.size)
				matrix[vertex.handle][vertex.handle] = 1.0

	#-------------------------SOLVE LINEAR SYSTEM-------------------------------
	temperatureField = np.linalg.solve(matrix, independent)

	#-------------------------PRINT ITERATION DATA------------------------------
	if iteration > 0 and not '-s' in sys.argv:
		print("{:>9}\t{:>14e}\t{:>14e}\t{:>14e}".format(iteration, currentTime, timeStep, difference))

	#-------------------------INCREMENT TIME------------------------------------
	currentTime += timeStep

	#-------------------------SAVE RESULTS--------------------------------------
	# cgnsSaver.timeSteps	= np.append(cgnsSaver.timeSteps,  currentTime)
	# cgnsSaver.fields  = np.vstack([cgnsSaver.fields, temperatureField])
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
finalSimulationTime = time.time()
if not '-s' in sys.argv:
	print("Ended Simultaion, elapsed {:.2f}s".format(finalSimulationTime-initialTime))

cgnsSaver.finalize()
if not '-s' in sys.argv:
	print("Saved file: elapsed {:.2f}s".format(time.time()-finalSimulationTime))

	print("\n\t\033[1;35mresult:\033[0m", problemData.paths["Output"]+"Results.cgns", '\n')
#-------------------------------------------------------------------------------
#-------------------------SHOW RESULTS GRAPHICALY-------------------------------
#-------------------------------------------------------------------------------
if '-g' in sys.argv:
	import matplotlib
	from matplotlib import pyplot as plt
	from matplotlib import cm
	from matplotlib.colors import ListedColormap as CM, Normalize
	from scipy.interpolate import griddata

	X,Y = zip(*[v.getCoordinates()[:-1] for v in grid.vertices])

	Xi, Yi = np.meshgrid( np.linspace(min(X), max(X), len(X)), np.linspace(min(Y), max(Y), len(Y)) )
	nTi = griddata((X,Y), temperatureField, (Xi,Yi), method='linear')

	plt.pcolor(Xi,Yi,nTi, cmap=CM( cm.get_cmap("RdBu",64)(np.linspace(1,0,64)) )) # Makes BuRd instead of RdBu
	plt.title("Numerical Temperature")
	plt.colorbar()	

	if '-1d' in sys.argv:
		from workspace.heat_transfer_1d.analyticalSolution import analyticalSolution_X
		plt.figure()

		x=np.linspace(0,1,200)
		X,T = zip(*[(x,t) for x,y,t in zip(X,Y,temperatureField) if y==1])
		plt.plot(x,[analyticalSolution_X(xx) for xx in x], color='k', label='Analytical Results')
		plt.scatter(X,T, marker='X', color='r', label='Numerical Results')
		plt.xlabel("X (m)")
		plt.ylabel("Temperature (K)")
		plt.legend()

		err = [abs(T - analyticalSolution_X(v.getCoordinates()[0])) for v,T in zip(grid.vertices,temperatureField)]
		print("Max Error: {:.3e}".format(max(err)))
		print("Mean Error: {:.3e}".format(sum(err)/len(err)))
		print("Euclidean Error: {:.3e}".format(np.sqrt(sum([v.volume*e**2 for v,e in zip(grid.vertices, err)]))))

	elif '-2d' in sys.argv:
		from workspace.sine_distribution.analyticalSolution import analyticalSolution_XY
		plt.figure()
		aT  = [analyticalSolution_XY(x,y) for x,y in zip(X,Y)]
		aTi = griddata((X,Y), aT, (Xi,Yi), method='linear')
		plt.pcolor(Xi,Yi,aTi, cmap=CM( cm.get_cmap("RdBu",64)(np.linspace(1,0,64)) ))
		plt.title("Analytical Temperature")
		plt.colorbar()

		err = [abs(T - analyticalSolution_XY(*v.getCoordinates()[:-1])) for v,T in zip(grid.vertices,temperatureField)]
		print("Max Error: {:.3e}".format(max(err)))
		print("Mean Error: {:.3e}".format(sum(err)/len(err)))
		print("Euclidean Error: {:.3e}".format(np.sqrt(sum([v.volume*e**2 for v,e in zip(grid.vertices, err)]))))


	plt.show()

if '-1d' in sys.argv and not '-g' in sys.argv:
	from workspace.heat_transfer_1d.analyticalSolution import analyticalSolution_X
	err = [abs(T - analyticalSolution_X(v.getCoordinates()[0])) for v,T in zip(grid.vertices,temperatureField)]
	print("Max Error: {:.3e}".format(max(err)))
	print("Mean Error: {:.3e}".format(sum(err)/len(err)))
	print("Euclidean Error: {:.3e}".format(np.sqrt(sum([v.volume*e**2 for v,e in zip(grid.vertices, err)]))))


if '-2d' in sys.argv and not '-g' in sys.argv:
	from workspace.sine_distribution.analyticalSolution import analyticalSolution_XY
	err = [abs(T - analyticalSolution_XY(*v.getCoordinates()[:-1])) for v,T in zip(grid.vertices,temperatureField)]
	print("Max Error: {:.3e}".format(max(err)))
	print("Mean Error: {:.3e}".format(sum(err)/len(err)))
	print("Euclidean Error: {:.3e}".format(np.sqrt(sum([v.volume*e**2 for v,e in zip(grid.vertices, err)]))))
