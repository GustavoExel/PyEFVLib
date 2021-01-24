import sys,os,shutil
sys.path.append(os.path.join(os.path.dirname(__file__), *[os.path.pardir]*3))
from PyEFVLib import MSHReader, Grid, ProblemData, CgnsSaver, CsvSaver, VtuSaver, VtmSaver
import numpy as np
import pandas as pd
import scipy as sp
from scipy import sparse
import scipy.sparse.linalg
import time

def heatTransfer(libraryPath, outputPath, extension, grid, propertyData, initialValues, neumannBoundaries, dirichletBoundaries, timeStep, finalTime, maxNumberOfIterations, tolerance, fileName="Results", transient=True, verbosity=True, color=True, exportFluxes=True, fluxesOutputPath="FLUXES-EXP.csv"):
	#-------------------------SETTINGS----------------------------------------------
	initialTime = time.time()

	dimension = grid.dimension
	currentTime = 0.0

	savers = {"cgns": CgnsSaver, "csv": CsvSaver, "vtu": VtuSaver, "vtm": VtmSaver}
	saver = savers[extension](grid, outputPath, libraryPath, fileName=fileName)
	if exportFluxes:
		fluxDF = pd.DataFrame( {
			"Element": [elemIdx for elemIdx in range(grid.elements.size) for facetIdx in range(4)],
			"Facet": grid.elements.size*[0,1,2,3],
			"X1": [element.vertices[i].x for element in grid.elements for i in [0,1,2,3] ],
			"Y1": [element.vertices[i].y for element in grid.elements for i in [0,1,2,3] ],
			"X2": [element.vertices[i].x for element in grid.elements for i in [1,2,3,0] ],
			"Y2": [element.vertices[i].y for element in grid.elements for i in [1,2,3,0] ]

		} )
		fluxDF.to_csv(fluxesOutputPath, index=False)
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

	#-------------------------SIMULATION MAIN LOOP----------------------------------
	while not converged:
		print(f"timeStep = {iteration}")
		if maxNumberOfIterations != None and iteration > maxNumberOfIterations:
			break
		if iteration > 1 and not transient:
			break
		#-------------------------ADD TO LINEAR SYSTEM------------------------------
		independent = np.zeros(grid.vertices.size)

		# Generation Term
		for region in grid.regions:
			heatGeneration = propertyData.get(region.handle, "HeatGeneration")
			for element in region.elements:
				local = 0
				for vertex in element.vertices:
					independent[vertex.handle] += element.subelementVolumes[local] * heatGeneration
					local += 1

		# Diffusion Term
		if iteration == 0:
			for region in grid.regions:
				conductivity = propertyData.get(region.handle, "Conductivity")
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
				density = propertyData.get(region.handle, "Density")
				heatCapacity = propertyData.get(region.handle, "HeatCapacity")
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
					independent[outerFace.vertex.handle] += bCondition.getValue(outerFace.handle) * np.linalg.norm(outerFace.area.getCoordinates())

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

		#----------------------------EXPORT FLUXES----------------------------------
		if exportFluxes:
			facetFlux = np.zeros((grid.elements.size,4,2))
			for region in grid.regions:
				conductivity = propertyData.get(region.handle, "Conductivity")
				for elementIdx, element in enumerate(region.elements):
					temperatureVector = np.array([temperatureField[vertex.handle] for vertex in element.vertices])
					for innerFaceIdx, innerFace in enumerate(element.innerFaces):
						heatFlux = -conductivity * np.matmul( innerFace.globalDerivatives, temperatureVector )
						# Como queremos o fluxo na direção da face, precisamos achar o produto interno do fluxo de calor com a área da face, dividir pelo quadrado da norma da área da face, e multiplicar pelo vetor área da face
						normFace = innerFace.area.getCoordinates()[:-1] / np.linalg.norm( innerFace.area.getCoordinates()[:-1] )
						heatFlux = normFace * np.dot(heatFlux, normFace)

						facetFlux[element.handle][innerFaceIdx] = heatFlux
			fluxDF = pd.read_csv(fluxesOutputPath)
			fluxDF[f"time_step_{iteration} - q''x"] = np.array([facetData[0] for elementData in facetFlux for facetData in elementData])
			fluxDF[f"time_step_{iteration} - q''y"] = np.array([facetData[1] for elementData in facetFlux for facetData in elementData])
			fluxDF[f"time_step_{iteration} - q''"] = np.array([np.linalg.norm(facetData) for elementData in facetFlux for facetData in elementData])
			fluxDF.to_csv(fluxesOutputPath, index=False)
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

	#-------------------------AFTER END OF MAIN LOOP ITERATION------------------------
	saver.finalize()

	return

if __name__ == "__main__":
	def init(size):
		global model,problemData,reader,grid
		model = "workspace/heat_transfer_2d/vug"
		problemData = ProblemData(model)
		problemData.meshFilePath = "/".join( problemData.meshFilePath.split(os.path.sep)[:-1] + [ f"vug {size:.0f}x{size:.0f}.msh" ] )
		reader = MSHReader(problemData.meshFilePath)
		grid = Grid(reader.getData())
		problemData.setGrid(grid)
		problemData.read()

	def run(exportFluxes,twoRegion,meshSize,k,transient,verbosity,outputPath=None):
		print( f"SETTINGS: exportFluxes = {exportFluxes}, twoRegion = {twoRegion}, meshSize = {meshSize}, k = {k:.0e}, transient = {transient}, verbosity = {verbosity}" )

		init(meshSize)
		problemData.propertyData.set(0, "Conductivity", 1 if twoRegion else k)
		problemData.propertyData.set(1, "Conductivity", k)
		if outputPath!=None: problemData.outputFilePath = outputPath

		# open("fluxes.csv", "w").close()

		finalTemperatureField = heatTransfer(
			libraryPath = problemData.libraryPath,
			outputPath = problemData.outputFilePath,
			extension = "vtm",#"csv" if not "--extension=cgns" in sys.argv else "cgns",
			
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
			
			transient = transient,
			verbosity = verbosity,
			exportFluxes = exportFluxes,
			fluxesOutputPath = "FLUXES-EXP.csv"
		)

		if exportFluxes:
			twoRegionName = f"fluxos - vec\\results\\fluxos-k1=1,k2={k:.0e}-{meshSize}x{meshSize}.csv"
			oneRegionName = f"fluxos - vec\\results\\fluxos-k1={k:.0e},k2={k:.0e}-{meshSize}x{meshSize}.csv"
			fileName = twoRegionName if twoRegion else oneRegionName

			try: os.remove(fileName)
			except: pass
			os.rename("FLUXES-EXP.csv", fileName)

	#---------------------SETTINGS----------------------------------
	outputPath = os.path.realpath( os.path.dirname( __file__ ) )
	exportFluxes = True
	twoRegion = True
	meshSize = 30
	# K = [1]
	K = [1,1e3,1e12,1e14,1e15,1e16,1e22]
	transient = True
	verbosity = False

	# run(exportFluxes,twoRegion,meshSize,K,transient,verbosity)

	if os.path.isdir("fluxos - vec"): shutil.rmtree("fluxos - vec")
	os.makedirs("fluxos - vec"); os.makedirs(os.path.join("fluxos - vec","results")); os.makedirs(os.path.join("fluxos - vec","images"))

	for meshSize in [15,30]:
		for k in K:
			pass
			run(exportFluxes,twoRegion,meshSize,k,transient,verbosity)
			# os.rename("Results.csv", os.path.join("temperatureFields", f"temperatureField-k1=1,k2={k:.0e}-{meshSize}x{meshSize}.csv"))
