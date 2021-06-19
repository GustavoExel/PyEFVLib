print("""
	PROBLEMA DE MANDEL:
	(EXEMPLO DO RELATÓRIO)
		2D
		COLUNA DE 1m x 1m
		Malha de Triângulos
		610 Nós
		1134 Elementos

		1MPa aplicado em cima
		
		timeStep = 10s
		finalTime = 1000s
		100 time steps
""")

import sys,os
dirname = os.path.realpath( os.path.dirname(__file__) )
sys.path.append(os.path.join(dirname, *[os.path.pardir]*3))
sys.path.append(os.path.join(dirname, *[os.path.pardir]*4, "solgeom"))
import PyEFVLib
from solgeom.Mandel import Solution
import numpy as np
import json

load = 1.0e+6
gravity = 0.0 #+

timeStep = 1
finalTime = 100

length = 1.0
height = 1.0
meshFilePath = "{MESHES}/msh/2D/Fine.msh"

run = True
plot = True

times = np.logspace(-2,2,50)
timeSteps = times[1:] - times[:-1]


def geomechanics(problemData):
	propertyData 	 = problemData.propertyData
	timeStep 		 = problemData.timeStep
	grid 			 = problemData.grid
	numberOfVertices = grid.numberOfVertices
	dimension 		 = grid.dimension

	csvSaver = PyEFVLib.CsvSaver(grid, problemData.outputFilePath, problemData.libraryPath)
	meshioSaver = PyEFVLib.MeshioSaver(grid, problemData.outputFilePath, problemData.libraryPath, extension="xdmf")

	uField    = np.repeat(0.0, dimension*numberOfVertices)
	oldUField = np.concatenate((problemData.initialValues["u_x"], problemData.initialValues["u_y"], problemData.initialValues["u_z"])) if dimension==3 else np.concatenate((problemData.initialValues["u_x"], problemData.initialValues["u_y"]))

	pField    = np.repeat(0.0, numberOfVertices)
	oldPField = problemData.initialValues["p"].copy()

	matrix 		= np.zeros(((1+dimension)*numberOfVertices, (1+dimension)*numberOfVertices))

	def getTransposedVoigtArea(innerFace):
		Sx, Sy, Sz = innerFace.area.getCoordinates()
		return np.array([[Sx,0,Sy],[0,Sy,Sx]]) if dimension==2 else np.array([[Sx,0,0,Sy,0,Sz],[0,Sy,0,Sx,Sz,0],[0,0,Sz,0,Sy,Sx]])

	def getVoigtGradientOperator(globalDerivatives):
		if len(globalDerivatives) == 2:
			Nx,Ny = globalDerivatives
			zero=np.zeros(Nx.size)
			return np.array([[Nx,zero],[zero,Ny],[Ny,Nx]])

		if len(globalDerivatives) == 3:
			Nx,Ny,Nz = globalDerivatives
			zero=np.zeros(Nx.size)
			return np.array([[Nx,zero,zero],[zero,Ny,zero],[zero,zero,Nz],[Ny,Nx,zero],[zero,Nz,Ny],[Nz,zero,Nx]])

	def assembleMatrix():
		matrix 		= np.zeros(((1+dimension)*numberOfVertices, (1+dimension)*numberOfVertices))
		for region in grid.regions:
			nu         = propertyData.get(region.handle, "nu")
			G          = propertyData.get(region.handle, "G")
			cs         = propertyData.get(region.handle, "cs")
			phi        = propertyData.get(region.handle, "phi")
			k          = propertyData.get(region.handle, "k")
			cf         = propertyData.get(region.handle, "cf")
			mu         = propertyData.get(region.handle, "mu")
			rhos       = propertyData.get(region.handle, "rhos")
			rhof       = propertyData.get(region.handle, "rhof")

			g = np.array([0.0, 0.0, 0.0])[:dimension]
			lame = 2*G*nu/(1-2*nu)
			Ce = np.array([[2*G+lame, lame, 0], [lame, 2*G+lame, 0], [0, 0, G]]) if dimension==2 else np.array([[2*G+lame, lame, lame, 0, 0, 0], [lame, 2*G+lame, lame, 0, 0, 0], [lame, lame, 2*G+lame, 0, 0, 0], [0, 0, 0, G, 0, 0], [0, 0, 0, 0, G, 0], [0, 0, 0, 0, 0, G]])
			rho = phi * rhof + (1-phi) * rhos
			K = 2*G*(1 + nu) / 3*(1-2*nu)
			cb = 1 / K
			alpha = 1 - cs / cb
			S = (phi * cf + (alpha-phi) * cs)

			for element in region.elements:
				# S * (1/timeStep) * p
				for vertex in element.vertices:
					matrix[vertex.handle+(dimension)*numberOfVertices][vertex.handle+(dimension)*numberOfVertices] += vertex.getSubElementVolume(element) * S * (1/timeStep) 

				# (-1) * alpha * grad(p)
				for face in element.innerFaces:
					area = face.area.getCoordinates()[:dimension]
					m = len(element.vertices)
					transposedVoigtArea = getTransposedVoigtArea(face)
					shapeFunctions = face.getShapeFunctions()
					identityShapeFunctionMatrix = np.array([shapeFunctions, shapeFunctions, shapeFunctions, np.zeros(m), np.zeros(m), np.zeros(m)]) if dimension == 3 else np.array([shapeFunctions, shapeFunctions, np.zeros(m)])
					matrixCoefficients = (-1) * alpha * np.matmul(transposedVoigtArea, identityShapeFunctionMatrix)
					backwardsHandle, forwardHandle = face.getNeighborVerticesHandles()
					for coord in range(dimension):
						for local, vertex in enumerate(element.vertices):
							matrix[backwardsHandle+(coord+0)*numberOfVertices][vertex.handle+(dimension)*numberOfVertices] += matrixCoefficients[coord][local]
							matrix[forwardHandle+(coord+0)*numberOfVertices][vertex.handle+(dimension)*numberOfVertices] += -matrixCoefficients[coord][local]

				# (-1) * k * mu * grad(p)
				for face in element.innerFaces:
					area = face.area.getCoordinates()[:dimension]
					matrixCoefficients = (-1) * k * mu * np.matmul( area.T, face.globalDerivatives )
					backwardsHandle, forwardHandle = face.getNeighborVerticesHandles()
					for local, vertex in enumerate(element.vertices):
						matrix[backwardsHandle+(dimension)*numberOfVertices][vertex.handle+(dimension)*numberOfVertices] += matrixCoefficients[local]
						matrix[forwardHandle+(dimension)*numberOfVertices][vertex.handle+(dimension)*numberOfVertices] += -matrixCoefficients[local]

				# Ce * grad_s(u)
				for face in element.innerFaces:
					area = face.area.getCoordinates()[:dimension]
					transposedVoigtArea = getTransposedVoigtArea(face)
					voigtGradientOperator = getVoigtGradientOperator(face.globalDerivatives)
					coeff = Ce
					matrixCoefficients = np.einsum("ij,jk,kmn->imn", transposedVoigtArea, coeff, voigtGradientOperator)
					backwardsHandle, forwardHandle = face.getNeighborVerticesHandles()
					for i in range(dimension):
						for j in range(dimension):
							for local, vertex in enumerate(element.vertices):
								matrix[backwardsHandle+(i+0)*numberOfVertices][vertex.handle+(j+0)*numberOfVertices] += matrixCoefficients[i][j][local]
								matrix[forwardHandle+(i+0)*numberOfVertices][vertex.handle+(j+0)*numberOfVertices] += -matrixCoefficients[i][j][local]

				# alpha * (1/timeStep) * u
				for face in element.faces:
					area = face.area.getCoordinates()[:dimension]
					shapeFunctions = face.getShapeFunctions()
					for coord in range(dimension):
						for local, vertex in enumerate(element.vertices):
							if type(face) == PyEFVLib.InnerFace:
								backwardsHandle, forwardHandle = face.getNeighborVerticesHandles()

								matrix[backwardsHandle+(dimension)*numberOfVertices][vertex.handle+(coord+0)*numberOfVertices] += alpha * (1/timeStep) * shapeFunctions[local] * area[coord]
								matrix[forwardHandle+(dimension)*numberOfVertices][vertex.handle+(coord+0)*numberOfVertices] += -alpha * (1/timeStep) * shapeFunctions[local] * area[coord]

							elif type(face) == PyEFVLib.OuterFace:
								matrix[face.vertex.handle+(dimension)*numberOfVertices][vertex.handle+(coord+0)*numberOfVertices] += alpha * (1/timeStep) * shapeFunctions[local] * area[coord]

		# Dirichlet Boundary Conditions
		for bCondition in problemData.dirichletBoundaries["u_x"]:
			for vertex in bCondition.boundary.vertices:
				matrix[vertex.handle+(0)*numberOfVertices] = np.zeros((1+dimension)*numberOfVertices)
				matrix[vertex.handle+(0)*numberOfVertices][vertex.handle+(0)*numberOfVertices] = 1.0

		for bCondition in problemData.dirichletBoundaries["u_y"]:
			for vertex in bCondition.boundary.vertices:
				matrix[vertex.handle+(1)*numberOfVertices] = np.zeros((1+dimension)*numberOfVertices)
				matrix[vertex.handle+(1)*numberOfVertices][vertex.handle+(1)*numberOfVertices] = 1.0

		if dimension == 3:
			for bCondition in problemData.dirichletBoundaries["u_z"]:
				for vertex in bCondition.boundary.vertices:
					matrix[vertex.handle+(2)*numberOfVertices] = np.zeros((1+dimension)*numberOfVertices)
					matrix[vertex.handle+(2)*numberOfVertices][vertex.handle+(2)*numberOfVertices] = 1.0

		for bCondition in problemData.dirichletBoundaries["p"]:
			for vertex in bCondition.boundary.vertices:
				matrix[vertex.handle+(dimension)*numberOfVertices] = np.zeros((1+dimension)*numberOfVertices)
				matrix[vertex.handle+(dimension)*numberOfVertices][vertex.handle+(dimension)*numberOfVertices] = 1.0

		# Inverse Matrix
		# inverseMatrix = np.linalg.inv(matrix)
		return matrix

	def assembleIndependent(currentTime):
		independent = np.zeros((1+dimension)*numberOfVertices)

		for region in grid.regions:
			nu         = propertyData.get(region.handle, "nu")
			G          = propertyData.get(region.handle, "G")
			cs         = propertyData.get(region.handle, "cs")
			phi        = propertyData.get(region.handle, "phi")
			k          = propertyData.get(region.handle, "k")
			cf         = propertyData.get(region.handle, "cf")
			mu         = propertyData.get(region.handle, "mu")
			rhos       = propertyData.get(region.handle, "rhos")
			rhof       = propertyData.get(region.handle, "rhof")

			g = np.array([0.0, 0.0, 0.0])[:dimension]
			lame = 2*G*nu/(1-2*nu)
			Ce = np.array([[2*G+lame, lame, 0], [lame, 2*G+lame, 0], [0, 0, G]]) if dimension==2 else np.array([[2*G+lame, lame, lame, 0, 0, 0], [lame, 2*G+lame, lame, 0, 0, 0], [lame, lame, 2*G+lame, 0, 0, 0], [0, 0, 0, G, 0, 0], [0, 0, 0, 0, G, 0], [0, 0, 0, 0, 0, G]])
			rho = phi * rhof + (1-phi) * rhos
			K = 2*G*(1 + nu) / 3*(1-2*nu)
			cb = 1 / K
			alpha = 1 - cs / cb
			S = (phi * cf + (alpha-phi) * cs)

			for element in region.elements:
				# (-1) * rho * g
				for vertex in element.vertices:
					for coord in range(dimension):
						independent[vertex.handle+coord*numberOfVertices] += vertex.getSubElementVolume(element) * (-1) * rho * g[coord]

				# S * (1/timeStep) * p_old
				for vertex in element.vertices:
					independent[vertex.handle+(dimension)*numberOfVertices] += vertex.getSubElementVolume(element) * S * (1/timeStep) * oldPField[vertex.handle]

				# (-1) * k * mu * rho * g
				for innerFace in element.innerFaces:
					area = innerFace.area.getCoordinates()[:dimension]
					coefficient = (-1) * k * mu * rho * g
					coefficient = np.matmul(area.T, coefficient)

					backwardsHandle, forwardHandle = innerFace.getNeighborVerticesHandles()
					independent[backwardsHandle+(dimension)*numberOfVertices] += coefficient
					independent[forwardHandle+(dimension)*numberOfVertices]   -= coefficient

				# alpha * (1/timeStep) * u_old
				for face in element.faces:
					area = face.area.getCoordinates()[:dimension]
					shapeFunctions = face.getShapeFunctions()

					for coord in range(dimension):
						for local, vertex in enumerate(element.vertices):
							if type(face) == PyEFVLib.InnerFace:
								backwardsHandle, forwardHandle = face.getNeighborVerticesHandles()
								independent[backwardsHandle+(dimension)*numberOfVertices] += alpha * (1/timeStep) * shapeFunctions[local] * area[coord] * oldUField[vertex.handle + (coord)*numberOfVertices]
								independent[forwardHandle+(dimension)*numberOfVertices] -= alpha * (1/timeStep) * shapeFunctions[local] * area[coord] * oldUField[vertex.handle + (coord)*numberOfVertices]

							elif type(face) == PyEFVLib.OuterFace:
								independent[face.vertex.handle+(dimension)*numberOfVertices] += alpha * (1/timeStep) * shapeFunctions[local] * area[coord] * oldUField[vertex.handle + (coord)*numberOfVertices]

		# Neumann Boundary Condition
		for bCondition in problemData.neumannBoundaries["u_x"]:
			for facet in bCondition.boundary.facets:
				for outerFace in facet.outerFaces:
					independent[outerFace.vertex.handle+(0)*numberOfVertices] += bCondition.getValue(outerFace.handle) * np.linalg.norm(outerFace.area.getCoordinates())

		for bCondition in problemData.neumannBoundaries["u_y"]:
			for facet in bCondition.boundary.facets:
				for outerFace in facet.outerFaces:
					independent[outerFace.vertex.handle+(1)*numberOfVertices] += bCondition.getValue(outerFace.handle) * np.linalg.norm(outerFace.area.getCoordinates())

		if dimension == 3:
			for bCondition in problemData.neumannBoundaries["u_z"]:
				for facet in bCondition.boundary.facets:
					for outerFace in facet.outerFaces:
						independent[outerFace.vertex.handle+(2)*numberOfVertices] += bCondition.getValue(outerFace.handle) * np.linalg.norm(outerFace.area.getCoordinates())

		for bCondition in problemData.neumannBoundaries["p"]:
			for facet in bCondition.boundary.facets:
				for outerFace in facet.outerFaces:
					independent[outerFace.vertex.handle+(dimension)*numberOfVertices] += bCondition.getValue(outerFace.handle) * np.linalg.norm(outerFace.area.getCoordinates())

		# Dirichlet Boundary Condition
		for bCondition in problemData.dirichletBoundaries["u_x"]:
			for vertex in bCondition.boundary.vertices:
				independent[vertex.handle+(0)*numberOfVertices] = bCondition.getValue(vertex.handle)

		for bCondition in problemData.dirichletBoundaries["u_y"]:
			for vertex in bCondition.boundary.vertices:
				if bCondition.boundary.name == "North":
					independent[vertex.handle+(1)*numberOfVertices] = mandel.getVertDisplacementValue(height, currentTime)
					# print(currentTime, mandel.getVertDisplacementValue(height, currentTime))
				else:
					independent[vertex.handle+(1)*numberOfVertices] = bCondition.getValue(vertex.handle)

		if dimension == 3:
			for bCondition in problemData.dirichletBoundaries["u_z"]:
				for vertex in bCondition.boundary.vertices:
					independent[vertex.handle+(2)*numberOfVertices] = bCondition.getValue(vertex.handle)

		for bCondition in problemData.dirichletBoundaries["p"]:
			for vertex in bCondition.boundary.vertices:
				independent[vertex.handle+(dimension)*numberOfVertices] = bCondition.getValue(vertex.handle)

		return independent

	tolerance = problemData.tolerance
	difference = 2*tolerance
	iteration = 0
	currentTime = 0.0
	converged = False

	rock = json.load(open(dirname+"/../../../../solgeom/examples/solid.json", "r"))
	fluid = json.load(open(dirname+"/../../../../solgeom/examples/fluid.json", "r"))
	mandel = Solution(length, height, load, rock, fluid)


	try:
		while not converged:
			timeStep = timeSteps[iteration]

			currentTime += timeStep
			matrix = assembleMatrix()
			independent = assembleIndependent( currentTime )

			results = np.linalg.solve(matrix, independent)
			uField = results[(0)*numberOfVertices:(0+dimension)*numberOfVertices]
			pField = results[(dimension)*numberOfVertices:(dimension+1)*numberOfVertices]

			difference = max( max(abs(pField - oldPField)), max(abs(uField - oldUField)) )

			oldPField = pField.copy()
			oldUField = uField.copy()

			iteration += 1

			csvSaver.save("u_x", uField[0*numberOfVertices:1*numberOfVertices], currentTime)
			csvSaver.save("u_y", uField[1*numberOfVertices:2*numberOfVertices], currentTime)
			if dimension == 3:
				csvSaver.save("u_z", uField[2*numberOfVertices:3*numberOfVertices], currentTime)
			csvSaver.save("p", pField, currentTime)

			meshioSaver.save("u_x", uField[0*numberOfVertices:1*numberOfVertices], currentTime)
			meshioSaver.save("u_y", uField[1*numberOfVertices:2*numberOfVertices], currentTime)
			if dimension == 3:
				meshioSaver.save("u_z", uField[2*numberOfVertices:3*numberOfVertices], currentTime)
			meshioSaver.save("p", pField, currentTime)

			print("{:>9}	{:>14.2e}	{:>14.2e}	{:>14.2e}".format(iteration, currentTime, timeStep, difference))
			converged = ( difference <= tolerance ) or ( currentTime >= problemData.finalTime ) or ( iteration >= problemData.maxNumberOfIterations )

			if iteration >= problemData.maxNumberOfIterations:
				break
	except: # Allows for Ctrl+C to interrupt the simulation and save the current results
		pass

	csvSaver.finalize()
	meshioSaver.finalize()



def geomechanics(problemData):
	propertyData 	 = problemData.propertyData
	timeStep 		 = problemData.timeStep
	grid 			 = problemData.grid
	numberOfVertices = grid.numberOfVertices
	dimension 		 = grid.dimension

	csvSaver = PyEFVLib.CsvSaver(grid, problemData.outputFilePath, problemData.libraryPath)
	meshioSaver = PyEFVLib.MeshioSaver(grid, problemData.outputFilePath, problemData.libraryPath, extension="xdmf")

	nextUField = np.repeat(0.0, dimension*numberOfVertices)
	uField    = np.repeat(0.0, dimension*numberOfVertices)
	prevUField = np.concatenate((problemData.initialValues["u_x"], problemData.initialValues["u_y"], problemData.initialValues["u_z"])) if dimension==3 else np.concatenate((problemData.initialValues["u_x"], problemData.initialValues["u_y"]))

	nextPField    = np.repeat(0.0, numberOfVertices)
	pField    = np.repeat(0.0, numberOfVertices)
	prevPField = problemData.initialValues["p"].copy()

	def getTransposedVoigtArea(innerFace):
		Sx, Sy, Sz = innerFace.area.getCoordinates()
		return np.array([[Sx,0,Sy],[0,Sy,Sx]]) if dimension==2 else np.array([[Sx,0,0,Sy,0,Sz],[0,Sy,0,Sx,Sz,0],[0,0,Sz,0,Sy,Sx]])

	def getVoigtGradientOperator(globalDerivatives):
		if len(globalDerivatives) == 2:
			Nx,Ny = globalDerivatives
			zero=np.zeros(Nx.size)
			return np.array([[Nx,zero],[zero,Ny],[Ny,Nx]])

		if len(globalDerivatives) == 3:
			Nx,Ny,Nz = globalDerivatives
			zero=np.zeros(Nx.size)
			return np.array([[Nx,zero,zero],[zero,Ny,zero],[zero,zero,Nz],[Ny,Nx,zero],[zero,Nz,Ny],[Nz,zero,Nx]])

	def assembleMatrix(timeStep):
		massMatrix 		= np.zeros((numberOfVertices, numberOfVertices))
		geoMatrix 		= np.zeros((dimension*numberOfVertices, dimension*numberOfVertices))

		for region in grid.regions:
			nu         = propertyData.get(region.handle, "nu")
			G          = propertyData.get(region.handle, "G")
			cs         = propertyData.get(region.handle, "cs")
			phi        = propertyData.get(region.handle, "phi")
			k          = propertyData.get(region.handle, "k")
			cf         = propertyData.get(region.handle, "cf")
			mu         = propertyData.get(region.handle, "mu")
			rhos       = propertyData.get(region.handle, "rhos")
			rhof       = propertyData.get(region.handle, "rhof")

			g = np.array([0.0, 0.0, 0.0])[:dimension]
			lame = 2*G*nu/(1-2*nu)
			Ce = np.array([[2*G+lame, lame, 0], [lame, 2*G+lame, 0], [0, 0, G]]) if dimension==2 else np.array([[2*G+lame, lame, lame, 0, 0, 0], [lame, 2*G+lame, lame, 0, 0, 0], [lame, lame, 2*G+lame, 0, 0, 0], [0, 0, 0, G, 0, 0], [0, 0, 0, 0, G, 0], [0, 0, 0, 0, 0, G]])
			rho = phi * rhof + (1-phi) * rhos
			K = 2*G*(1 + nu) / 3*(1-2*nu)
			cb = 1 / K
			alpha = 1 - cs / cb
			S = (phi * cf + (alpha-phi) * cs)

			for element in region.elements:
				# S * (1/timeStep) * p
				for vertex in element.vertices:
					massMatrix[vertex.handle][vertex.handle] += vertex.getSubElementVolume(element) * S * (1/timeStep) 
					massMatrix[vertex.handle][vertex.handle] += vertex.getSubElementVolume(element) * cb* (1/timeStep) * alpha**2

				for face in element.innerFaces:
					area = face.area.getCoordinates()[:dimension]

					transposedVoigtArea = getTransposedVoigtArea(face)
					voigtGradientOperator = getVoigtGradientOperator(face.globalDerivatives)

					matrixCoefficientsH = (-1) * k * mu * np.matmul( area.T, face.globalDerivatives )
					matrixCoefficientsK = np.einsum("ij,jk,kmn->imn", transposedVoigtArea, Ce, voigtGradientOperator)

					backwardsHandle, forwardHandle = face.getNeighborVerticesHandles()
					for local, vertex in enumerate(element.vertices):
						# (-1) * k * mu * grad(p)
						massMatrix[backwardsHandle][vertex.handle] += matrixCoefficientsH[local]
						massMatrix[forwardHandle][vertex.handle] += -matrixCoefficientsH[local]

						# Ce * grad_s(u)
						for i in range(dimension):
							for j in range(dimension):
								geoMatrix[backwardsHandle+(i+0)*numberOfVertices][vertex.handle+(j+0)*numberOfVertices] += matrixCoefficientsK[i][j][local]
								geoMatrix[forwardHandle+(i+0)*numberOfVertices][vertex.handle+(j+0)*numberOfVertices] += -matrixCoefficientsK[i][j][local]

		# Dirichlet Boundary Conditions
		for bCondition in problemData.dirichletBoundaries["u_x"]:
			for vertex in bCondition.boundary.vertices:
				geoMatrix[vertex.handle+(0)*numberOfVertices] = np.zeros(dimension*numberOfVertices)
				geoMatrix[vertex.handle+(0)*numberOfVertices][vertex.handle+(0)*numberOfVertices] = 1.0

		for bCondition in problemData.dirichletBoundaries["u_y"]:
			for vertex in bCondition.boundary.vertices:
				geoMatrix[vertex.handle+(1)*numberOfVertices] = np.zeros(dimension*numberOfVertices)
				geoMatrix[vertex.handle+(1)*numberOfVertices][vertex.handle+(1)*numberOfVertices] = 1.0

		if dimension == 3:
			for bCondition in problemData.dirichletBoundaries["u_z"]:
				for vertex in bCondition.boundary.vertices:
					geoMatrix[vertex.handle+(2)*numberOfVertices] = np.zeros(dimension*numberOfVertices)
					geoMatrix[vertex.handle+(2)*numberOfVertices][vertex.handle+(2)*numberOfVertices] = 1.0

		for bCondition in problemData.dirichletBoundaries["p"]:
			for vertex in bCondition.boundary.vertices:
				massMatrix[vertex.handle] = np.zeros(numberOfVertices)
				massMatrix[vertex.handle][vertex.handle] = 1.0

		# Inverse Matrix
		inverseMassMatrix = np.linalg.inv(massMatrix)
		inverseGeoMatrix = np.linalg.inv(geoMatrix)
		return inverseMassMatrix, inverseGeoMatrix

	def assembleMassIndependent(prevPField, prevUField, pField, uField, timeStep):
		massIndependent = np.zeros(numberOfVertices)

		for region in grid.regions:
			nu         = propertyData.get(region.handle, "nu")
			G          = propertyData.get(region.handle, "G")
			cs         = propertyData.get(region.handle, "cs")
			phi        = propertyData.get(region.handle, "phi")
			k          = propertyData.get(region.handle, "k")
			cf         = propertyData.get(region.handle, "cf")
			mu         = propertyData.get(region.handle, "mu")
			rhos       = propertyData.get(region.handle, "rhos")
			rhof       = propertyData.get(region.handle, "rhof")

			g = np.array([0.0, 0.0, 0.0])[:dimension]
			lame = 2*G*nu/(1-2*nu)
			Ce = np.array([[2*G+lame, lame, 0], [lame, 2*G+lame, 0], [0, 0, G]]) if dimension==2 else np.array([[2*G+lame, lame, lame, 0, 0, 0], [lame, 2*G+lame, lame, 0, 0, 0], [lame, lame, 2*G+lame, 0, 0, 0], [0, 0, 0, G, 0, 0], [0, 0, 0, 0, G, 0], [0, 0, 0, 0, 0, G]])
			rho = phi * rhof + (1-phi) * rhos
			K = 2*G*(1 + nu) / 3*(1-2*nu)
			cb = 1 / K
			alpha = 1 - cs / cb
			S = (phi * cf + (alpha-phi) * cs)

			for element in region.elements:
				# S * (1/timeStep) * p_prev
				for vertex in element.vertices:
					massIndependent[vertex.handle] += vertex.getSubElementVolume(element) * S * (1/timeStep) * prevPField[vertex.handle]
					massIndependent[vertex.handle] += vertex.getSubElementVolume(element) * cb* (1/timeStep) * pField[vertex.handle] * alpha**2


				for face in element.innerFaces:
					area = face.area.getCoordinates()[:dimension]

					shapeFunctions = face.getShapeFunctions()

					# (-1) * k * mu * rho * g
					coefficientH = (-1) * k * mu * rho * g
					coefficientH = np.matmul(area.T, coefficientH)

					backwardsHandle, forwardHandle = face.getNeighborVerticesHandles()

					if type(face) == PyEFVLib.InnerFace:
						# (-1) * k * mu * rho * g
						massIndependent[backwardsHandle] += coefficientH
						massIndependent[forwardHandle]   -= coefficientH

					for coord in range(dimension):
						for local, vertex in enumerate(element.vertices):
							if type(face) == PyEFVLib.InnerFace:
								# alpha * (1/timeStep) * u_prev
								massIndependent[backwardsHandle] += alpha * (1/timeStep) * shapeFunctions[local] * area[coord] * prevUField[vertex.handle + (coord)*numberOfVertices]
								massIndependent[forwardHandle] -= alpha * (1/timeStep) * shapeFunctions[local] * area[coord] * prevUField[vertex.handle + (coord)*numberOfVertices]
								# alpha * (1/timeStep) * u
								massIndependent[backwardsHandle] -= alpha * (1/timeStep) * shapeFunctions[local] * area[coord] * uField[vertex.handle+(coord+0)*numberOfVertices]
								massIndependent[forwardHandle] -= -alpha * (1/timeStep) * shapeFunctions[local] * area[coord] * uField[vertex.handle+(coord+0)*numberOfVertices]

							elif type(face) == PyEFVLib.OuterFace:
								# alpha * (1/timeStep) * u_prev
								massIndependent[face.vertex.handle] += alpha * (1/timeStep) * shapeFunctions[local] * area[coord] * prevUField[vertex.handle + (coord)*numberOfVertices]
								# alpha * (1/timeStep) * u
								massIndependent[face.vertex.handle] -= alpha * (1/timeStep) * shapeFunctions[local] * area[coord] * uField[vertex.handle+(coord+0)*numberOfVertices]


		# Neumann Boundary Condition
		for bCondition in problemData.neumannBoundaries["p"]:
			for facet in bCondition.boundary.facets:
				for outerFace in facet.outerFaces:
					massIndependent[outerFace.vertex.handle] += bCondition.getValue(outerFace.handle) * np.linalg.norm(outerFace.area.getCoordinates())

		# Dirichlet Boundary Condition
		for bCondition in problemData.dirichletBoundaries["p"]:
			for vertex in bCondition.boundary.vertices:
				massIndependent[vertex.handle] = bCondition.getValue(vertex.handle)

		return massIndependent

	def assembleGeoIndependent(nextPField, timeStep):
		geoIndependent = np.zeros(dimension*numberOfVertices)

		for region in grid.regions:
			nu         = propertyData.get(region.handle, "nu")
			G          = propertyData.get(region.handle, "G")
			cs         = propertyData.get(region.handle, "cs")
			phi        = propertyData.get(region.handle, "phi")
			k          = propertyData.get(region.handle, "k")
			cf         = propertyData.get(region.handle, "cf")
			mu         = propertyData.get(region.handle, "mu")
			rhos       = propertyData.get(region.handle, "rhos")
			rhof       = propertyData.get(region.handle, "rhof")

			g = np.array([0.0, 0.0, 0.0])[:dimension]
			lame = 2*G*nu/(1-2*nu)
			Ce = np.array([[2*G+lame, lame, 0], [lame, 2*G+lame, 0], [0, 0, G]]) if dimension==2 else np.array([[2*G+lame, lame, lame, 0, 0, 0], [lame, 2*G+lame, lame, 0, 0, 0], [lame, lame, 2*G+lame, 0, 0, 0], [0, 0, 0, G, 0, 0], [0, 0, 0, 0, G, 0], [0, 0, 0, 0, 0, G]])
			rho = phi * rhof + (1-phi) * rhos
			K = 2*G*(1 + nu) / 3*(1-2*nu)
			cb = 1 / K
			alpha = 1 - cs / cb
			S = (phi * cf + (alpha-phi) * cs)

			for element in region.elements:
				# (-1) * rho * g
				for vertex in element.vertices:
					for coord in range(dimension):
						geoIndependent[vertex.handle+coord*numberOfVertices] += vertex.getSubElementVolume(element) * (-1) * rho * g[coord]

				# (-1) * alpha * grad(p)
				for face in element.innerFaces:
					area = face.area.getCoordinates()[:dimension]
					m = len(element.vertices)
					transposedVoigtArea = getTransposedVoigtArea(face)
					shapeFunctions = face.getShapeFunctions()
					identityShapeFunctionMatrix = np.array([shapeFunctions, shapeFunctions, shapeFunctions, np.zeros(m), np.zeros(m), np.zeros(m)]) if dimension == 3 else np.array([shapeFunctions, shapeFunctions, np.zeros(m)])
					matrixCoefficients = (-1) * alpha * np.matmul(transposedVoigtArea, identityShapeFunctionMatrix)
					backwardsHandle, forwardHandle = face.getNeighborVerticesHandles()
					for coord in range(dimension):
						for local, vertex in enumerate(element.vertices):
							geoIndependent[backwardsHandle+(coord+0)*numberOfVertices] -= matrixCoefficients[coord][local] * nextPField[vertex.handle]
							geoIndependent[forwardHandle+(coord+0)*numberOfVertices] -= -matrixCoefficients[coord][local] * nextPField[vertex.handle]

		# Neumann Boundary Condition
		for bCondition in problemData.neumannBoundaries["u_x"]:
			for facet in bCondition.boundary.facets:
				for outerFace in facet.outerFaces:
					geoIndependent[outerFace.vertex.handle+(0)*numberOfVertices] += bCondition.getValue(outerFace.handle) * np.linalg.norm(outerFace.area.getCoordinates())

		for bCondition in problemData.neumannBoundaries["u_y"]:
			for facet in bCondition.boundary.facets:
				for outerFace in facet.outerFaces:
					geoIndependent[outerFace.vertex.handle+(1)*numberOfVertices] += bCondition.getValue(outerFace.handle) * np.linalg.norm(outerFace.area.getCoordinates())

		if dimension == 3:
			for bCondition in problemData.neumannBoundaries["u_z"]:
				for facet in bCondition.boundary.facets:
					for outerFace in facet.outerFaces:
						geoIndependent[outerFace.vertex.handle+(2)*numberOfVertices] += bCondition.getValue(outerFace.handle) * np.linalg.norm(outerFace.area.getCoordinates())

		# Dirichlet Boundary Condition
		for bCondition in problemData.dirichletBoundaries["u_x"]:
			for vertex in bCondition.boundary.vertices:
				geoIndependent[vertex.handle+(0)*numberOfVertices] = bCondition.getValue(vertex.handle)

		for bCondition in problemData.dirichletBoundaries["u_y"]:
			for vertex in bCondition.boundary.vertices:
				if bCondition.boundary.name == "North":
					geoIndependent[vertex.handle+(1)*numberOfVertices] = mandel.getVertDisplacementValue(height, currentTime)
				else:
					geoIndependent[vertex.handle+(1)*numberOfVertices] = bCondition.getValue(vertex.handle)

		if dimension == 3:
			for bCondition in problemData.dirichletBoundaries["u_z"]:
				for vertex in bCondition.boundary.vertices:
					geoIndependent[vertex.handle+(2)*numberOfVertices] = bCondition.getValue(vertex.handle)

		return geoIndependent

	tolerance = problemData.tolerance
	difference = 2*tolerance
	currentTime = 0.0
	converged = False
	timeIteration = 0

	rock = json.load(open(dirname+"/../../../../solgeom/examples/solid.json", "r"))
	fluid = json.load(open(dirname+"/../../../../solgeom/examples/fluid.json", "r"))
	mandel = Solution(length, height, load, rock, fluid)

	while not converged:
		timeStep = timeSteps[timeIteration]
		difference = 2*tolerance

		pField = prevPField.copy()
		uField = prevUField.copy()

		innerIteration = 0

		inverseMassMatrix, inverseGeoMatrix = assembleMatrix(timeStep)

		while difference >= tolerance or innerIteration<=1:
			massIndependent = assembleMassIndependent(prevPField, prevUField, pField, uField, timeStep)
			nextPField = np.matmul(inverseMassMatrix, massIndependent)
			
			geoIndependent = assembleGeoIndependent(nextPField, timeStep)
			nextUField = np.matmul(inverseGeoMatrix, geoIndependent)
			
			difference = max( max(abs(nextPField - pField)), max(abs(nextUField - uField)) )
			print(difference)

			pField = nextPField.copy()
			uField = nextUField.copy()

			innerIteration += 1

		timeIteration += 1

		difference = max( max(abs(pField - prevPField)), max(abs(uField - prevUField)) )
		
		prevPField = pField.copy()
		prevUField = uField.copy()

		currentTime += timeStep

		csvSaver.save("u_x", uField[0*numberOfVertices:1*numberOfVertices], currentTime)
		csvSaver.save("u_y", uField[1*numberOfVertices:2*numberOfVertices], currentTime)
		if dimension == 3:
			csvSaver.save("u_z", uField[2*numberOfVertices:3*numberOfVertices], currentTime)
		csvSaver.save("p", pField, currentTime)

		meshioSaver.save("u_x", uField[0*numberOfVertices:1*numberOfVertices], currentTime)
		meshioSaver.save("u_y", uField[1*numberOfVertices:2*numberOfVertices], currentTime)
		if dimension == 3:
			meshioSaver.save("u_z", uField[2*numberOfVertices:3*numberOfVertices], currentTime)
		meshioSaver.save("p", pField, currentTime)

		print("{:>9}	{:>14.2e}	{:>14.2e}	{:>14.2e}".format(innerIteration, currentTime, timeStep, difference))
		converged = ( difference <= tolerance ) or ( currentTime >= problemData.finalTime ) or ( timeIteration >= problemData.maxNumberOfIterations )

		if timeIteration >= problemData.maxNumberOfIterations:
			break

	csvSaver.finalize()
	meshioSaver.finalize()




def main_run():
	problemData = PyEFVLib.ProblemData(
		meshFilePath = meshFilePath,
		outputFilePath = "{RESULTS}/geomechanics",
		numericalSettings = PyEFVLib.NumericalSettings( timeStep = timeStep, tolerance=0.25, maxNumberOfIterations = 50-1 ),
		propertyData = PyEFVLib.PropertyData({
			'Body':
			{
				'nu': 0.2,
				'G': 6000000000.0,
				'cs': 0.0,
				'phi': 0.19,
				'k': 1.9e-15,
				'cf': 3.0303e-10,
				'mu': 1000.0,
				'rhos': 2700.0,
				'rhof': 1000.0,
			},
		}),
		boundaryConditions = PyEFVLib.BoundaryConditions({
			'u_x': {
				'InitialValue': 0.0,
				'West': { 'condition' : PyEFVLib.Dirichlet, 'type' : PyEFVLib.Constant, 'value' : 0.0 },
				'East': { 'condition' : PyEFVLib.Neumann, 'type' : PyEFVLib.Constant, 'value' : 0.0 },
				'South': { 'condition' : PyEFVLib.Neumann, 'type' : PyEFVLib.Constant, 'value' : 0.0 },
				'North': { 'condition' : PyEFVLib.Neumann, 'type' : PyEFVLib.Constant, 'value' : 0.0 },
			},
			'u_y': {
				'InitialValue': 0.0,
				'West': { 'condition' : PyEFVLib.Neumann, 'type' : PyEFVLib.Constant, 'value' : 0.0 },
				'East': { 'condition' : PyEFVLib.Neumann, 'type' : PyEFVLib.Constant, 'value' : 0.0 },
				'South': { 'condition' : PyEFVLib.Dirichlet, 'type' : PyEFVLib.Constant, 'value' : 0.0 },
				'North': { 'condition' : PyEFVLib.Dirichlet, 'type' : PyEFVLib.Constant, 'value' : 0.0 },
			},
			'p': {
				'InitialValue': 0.0,
				'West': { 'condition' : PyEFVLib.Neumann, 'type' : PyEFVLib.Constant, 'value' : 0.0 },
				'East': { 'condition' : PyEFVLib.Dirichlet, 'type' : PyEFVLib.Constant, 'value' : 0.0 },
				'South': { 'condition' : PyEFVLib.Neumann, 'type' : PyEFVLib.Constant, 'value' : 0.0 },
				'North': { 'condition' : PyEFVLib.Neumann, 'type' : PyEFVLib.Constant, 'value' : 0.0 },
			},
		}),
	)

	geomechanics( problemData )

if run and __name__ == '__main__':
	main_run()


################################################################################################################################################################
################################################################################################################################################################
################################################################################################################################################################
################################################################################################################################################################
################################################################################################################################################################
################################################################################################################################################################
################################################################################################################################################################
################################################################################################################################################################




import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import json
import sys, os

dirname = os.path.realpath( os.path.dirname(__file__) )
sys.path.append(dirname + "/../../../../solgeom")
from solgeom.Mandel import Solution

filePath = dirname+"/../../../results/geomechanics/Results.csv"

df = pd.read_csv(filePath)
X = df["X"]
Y = df["Y"]

numberOfTimeSteps = int(df.keys()[-1].split(" - ")[0].replace("TimeStep",""))

stepInTimeSamples = [5, int(numberOfTimeSteps*.8), int(numberOfTimeSteps*.9), numberOfTimeSteps-1]

mm = 1000.
kPa = 1/1000.

def main_plot():
	solid = json.load(open(dirname+"/../../../../solgeom/examples/solid.json", "r"))
	fluid = json.load(open(dirname+"/../../../../solgeom/examples/fluid.json", "r"))

	mandel = Solution(length, height, load, solid, fluid)
	x = mandel.getXPositionValues()
	y = mandel.getYPositionValues()

	fig, ax = plt.subplots(2, 3, figsize=(16,7))
	fig.subplots_adjust(left=0.060, right=0.940, top=0.970, bottom=0.085, hspace=0.235, wspace=0.300)

	# Plot pressure profiles ---------------------------------------------------
	# times_1 = [0.01, 1., 10., 30., 80.]
	times_1 = [times[N] for N in stepInTimeSamples]
	for idx,t in zip(stepInTimeSamples,times_1):
		p_a = mandel.getPressureValuesConstTime(t)
		p_n = np.array([p for p,x,y in zip(df[f"TimeStep{idx} - p"],X,Y) if y==0.0])
		x_n = np.array([x for p,x,y in zip(df[f"TimeStep{idx} - p"],X,Y) if y==0.0])

		ax[0][0].plot(x, p_a*kPa)
		ax[0][0].scatter(x_n, p_n*kPa, color="r", marker=".",linewidth=0.5)
	ax[0][0].set_xlabel("Length (m)", size=12)
	ax[0][0].set_ylabel("Pressure (kPa)", size=12)
	ax[0][0].grid(True)
	# --------------------------------------------------------------------------

	# Plot pressure over time --------------------------------------------------
	times_2 = np.logspace(-2, 3., 100)
	# times_n = timeStep * np.arange(1,numberOfTimeSteps+1)
	times_n = times[:numberOfTimeSteps]
	p_a = mandel.getPressureValuesAtPosition(0.0, times_2)
	vtx = [idx for idx,(x,y) in enumerate(zip(X,Y)) if x==0.0 and y==0.0][0]
	p_n = np.array([df[f"TimeStep{s} - p"][vtx] for s in range(1,numberOfTimeSteps+1)])

	# ax[1][0].plot(times_2, p_a*kPa)
	ax[1][0].semilogx(times_2, p_a*kPa)
	ax[1][0].scatter(times_n, p_n*kPa, marker='.', color='r')
	ax[1][0].set_xlabel("Time (s)", size=12)
	ax[1][0].set_ylabel("Pressure at %s"%(r"$x=0$"), size=12)
	ax[1][0].grid(True)
	# --------------------------------------------------------------------------

	# Plot horizontal displacement profiles ------------------------------------
	for idx,t in zip(stepInTimeSamples,times_1):
		u_a = mandel.getHorDisplacementValuesConstTime(t)
		u_n = df[f"TimeStep{idx} - u_x"]
		ax[0][1].plot(x, u_a*mm, label=f"t = {t:.02f}s")
		ax[0][1].scatter(X, u_n*mm, marker='.', color='r')
	ax[0][1].set_xlabel("x (m)", size=12)
	ax[0][1].set_ylabel("Hor. displacement (mm)", size=12)
	ax[0][1].grid(True)
	ax[0][1].legend()
	# --------------------------------------------------------------------------

	# Plot horizontal displacement at x=L --------------------------------------
	u_a = mandel.getHorDisplacementValuesAtPosition(length, times_2)
	vtx = [idx for idx,(x,y) in enumerate(zip(X,Y)) if x==length and y==0.0][0]
	u_n = np.array([df[f"TimeStep{s} - u_x"][vtx] for s in range(1,numberOfTimeSteps+1)])
	# ax[1][1].plot(times_2, u_a*mm)
	ax[1][1].semilogx(times_2, u_a*mm)
	ax[1][1].scatter(times_n, u_n*mm, marker='.', color='r')
	ax[1][1].set_xlabel("Time (s)", size=12)
	ax[1][1].set_ylabel("Hor. displacement at %s"%(r"$x=L$"), size=12)
	ax[1][1].grid(True)
	# --------------------------------------------------------------------------

	# Plot vertical displacement profiles --------------------------------------
	for idx,t in zip(stepInTimeSamples,times_1):
		v_a = mandel.getVertDisplacementValuesConstTime(t)
		v_n = df[f"TimeStep{idx} - u_y"]
		ax[0][2].plot(v_a*mm, y)
		ax[0][2].scatter(v_n*mm, Y, marker='.', color='r')

	ax[0][2].set_xlabel("Vert. displacement (mm)", size=12)
	ax[0][2].set_ylabel("y (m)", size=12)
	ax[0][2].grid(True)
	# --------------------------------------------------------------------------

	# Plot vertical displacement at y=H ----------------------------------------
	v_a = mandel.getVertDisplacementValuesAtPosition(height, times_2)
	vtx = [idx for idx,(x,y) in enumerate(zip(X,Y)) if x==0.0 and y==height][0]
	v_n = np.array([df[f"TimeStep{s} - u_y"][vtx] for s in range(1,numberOfTimeSteps+1)])

	ax[1][2].semilogx(times_2, v_a*mm)
	# ax[1][2].plot(times_2, v_a*mm)
	ax[1][2].scatter(times_n, v_n*mm, marker='.', color='r')
	ax[1][2].set_xlabel("Time (s)", size=12)
	ax[1][2].set_ylabel("Vert. displacement at %s"%(r"$y=H$"), size=12)
	ax[1][2].grid(True)
	# --------------------------------------------------------------------------





	plt.legend()
	plt.show()



if plot and __name__ == '__main__':
	main_plot()