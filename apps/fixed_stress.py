import sys,os
sys.path.append(os.path.join(os.path.dirname(__file__), os.pardir))
import PyEFVLib
import numpy as np

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

	def assembleMatrix():
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

				# (-1) * k * mu * grad(p)
				for face in element.innerFaces:
					area = face.area.getCoordinates()[:dimension]
					matrixCoefficients = (-1) * k * mu * np.matmul( area.T, face.globalDerivatives )
					backwardsHandle, forwardHandle = face.getNeighborVerticesHandles()
					for local, vertex in enumerate(element.vertices):
						massMatrix[backwardsHandle][vertex.handle] += matrixCoefficients[local]
						massMatrix[forwardHandle][vertex.handle] += -matrixCoefficients[local]

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
								geoMatrix[backwardsHandle+(i+0)*numberOfVertices][vertex.handle+(j+0)*numberOfVertices] += matrixCoefficients[i][j][local]
								geoMatrix[forwardHandle+(i+0)*numberOfVertices][vertex.handle+(j+0)*numberOfVertices] += -matrixCoefficients[i][j][local]

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

	def assembleMassIndependent(prevPField, prevUField, pField, uField):
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

				# (-1) * k * mu * rho * g
				for innerFace in element.innerFaces:
					area = innerFace.area.getCoordinates()[:dimension]
					coefficient = (-1) * k * mu * rho * g
					coefficient = np.matmul(area.T, coefficient)

					backwardsHandle, forwardHandle = innerFace.getNeighborVerticesHandles()
					massIndependent[backwardsHandle] += coefficient
					massIndependent[forwardHandle]   -= coefficient

				# alpha * (1/timeStep) * u_prev
				for face in element.faces:
					area = face.area.getCoordinates()[:dimension]
					shapeFunctions = face.getShapeFunctions()

					for coord in range(dimension):
						for local, vertex in enumerate(element.vertices):
							if type(face) == PyEFVLib.InnerFace:
								backwardsHandle, forwardHandle = face.getNeighborVerticesHandles()
								massIndependent[backwardsHandle] += alpha * (1/timeStep) * shapeFunctions[local] * area[coord] * prevUField[vertex.handle + (coord)*numberOfVertices]
								massIndependent[forwardHandle] -= alpha * (1/timeStep) * shapeFunctions[local] * area[coord] * prevUField[vertex.handle + (coord)*numberOfVertices]

							elif type(face) == PyEFVLib.OuterFace:
								massIndependent[face.vertex.handle] += alpha * (1/timeStep) * shapeFunctions[local] * area[coord] * prevUField[vertex.handle + (coord)*numberOfVertices]

				# alpha * (1/timeStep) * u
				for face in element.faces:
					area = face.area.getCoordinates()[:dimension]
					shapeFunctions = face.getShapeFunctions()
					for coord in range(dimension):
						for local, vertex in enumerate(element.vertices):
							if type(face) == PyEFVLib.InnerFace:
								backwardsHandle, forwardHandle = face.getNeighborVerticesHandles()

								massIndependent[backwardsHandle] -= alpha * (1/timeStep) * shapeFunctions[local] * area[coord] * uField[vertex.handle+(coord+0)*numberOfVertices]
								massIndependent[forwardHandle] -= -alpha * (1/timeStep) * shapeFunctions[local] * area[coord] * uField[vertex.handle+(coord+0)*numberOfVertices]

							elif type(face) == PyEFVLib.OuterFace:
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

	def assembleGeoIndependent(nextPField):
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
					m = element.numberOfVertices
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

	inverseMassMatrix, inverseGeoMatrix = assembleMatrix()

	while not converged:
		difference = 2*tolerance

		pField = prevPField.copy()
		uField = prevUField.copy()

		iteration = 0

		while difference >= tolerance or iteration<=1:
			massIndependent = assembleMassIndependent(prevPField, prevUField, pField, uField)
			nextPField = np.matmul(inverseMassMatrix, massIndependent)
			
			geoIndependent = assembleGeoIndependent(nextPField)
			nextUField = np.matmul(inverseGeoMatrix, geoIndependent)
			
			difference = max( max(abs(nextPField - pField)), max(abs(nextUField - uField)) )

			pField = nextPField.copy()
			uField = nextUField.copy()
			
			iteration += 1

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

		print("{:>9}	{:>14.2e}	{:>14.2e}	{:>14.2e}".format(iteration, currentTime, timeStep, difference))
		converged = ( difference <= tolerance ) or ( currentTime >= problemData.finalTime ) or ( iteration >= problemData.maxNumberOfIterations )

		if iteration >= problemData.maxNumberOfIterations:
				break

	csvSaver.finalize()
	meshioSaver.finalize()

def main():
	problemData = PyEFVLib.ProblemData(
		meshFilePath = "{MESHES}/msh/2D/Square.msh",
		outputFilePath = "{RESULTS}/geomechanics",
		numericalSettings = PyEFVLib.NumericalSettings( timeStep = 20, finalTime = 100, maxNumberOfIterations = 600 ),
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
				'East': { 'condition' : PyEFVLib.Dirichlet, 'type' : PyEFVLib.Constant, 'value' : 0.0 },
				'South': { 'condition' : PyEFVLib.Neumann, 'type' : PyEFVLib.Constant, 'value' : 0.0 },
				'North': { 'condition' : PyEFVLib.Neumann, 'type' : PyEFVLib.Constant, 'value' : 0.0 },
			},
			'u_y': {
				'InitialValue': 0.0,
				'West': { 'condition' : PyEFVLib.Neumann, 'type' : PyEFVLib.Constant, 'value' : 0.0 },
				'East': { 'condition' : PyEFVLib.Neumann, 'type' : PyEFVLib.Constant, 'value' : 0.0 },
				'South': { 'condition' : PyEFVLib.Dirichlet, 'type' : PyEFVLib.Constant, 'value' : 0.0 },
				'North': { 'condition' : PyEFVLib.Neumann, 'type' : PyEFVLib.Constant, 'value' : 1e5 },
			},
			'p': {
				'InitialValue': 0.0,
				'West': { 'condition' : PyEFVLib.Neumann, 'type' : PyEFVLib.Constant, 'value' : 0.0 },
				'East': { 'condition' : PyEFVLib.Neumann, 'type' : PyEFVLib.Constant, 'value' : 0.0 },
				'South': { 'condition' : PyEFVLib.Neumann, 'type' : PyEFVLib.Constant, 'value' : 0.0 },
				'North': { 'condition' : PyEFVLib.Dirichlet, 'type' : PyEFVLib.Constant, 'value' : 0.0 },
			},
		}),
	)

	geomechanics( problemData )

if __name__ == '__main__':
	main()
