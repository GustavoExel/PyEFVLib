# print("""
# 	STRIP FOOTING
# """)

load = 1.0e+5
gravity = 0.0 #+

width = 16.0
height = 8.0
load_width = 1.0

t_load = 440750.0

meshName = "strip_footing_eq_spacing"
# meshName = "strip_footing_uneq_spacing_x"
# meshName = "strip_footing_uneq_spacing_xy"
meshFilePath = f"{{MESHES}}/msh/2D/{meshName}.msh"

run = True
plot = False

import sys,os
sys.path.append(os.path.join(os.path.dirname(__file__), *[os.path.pardir]*3))
import PyEFVLib
import numpy as np

# times = np.logspace(-1,2,100) - 1e-1
# timeSteps = times[1:]-times[:-1]

dts = [4407.5,8815.1,22038,44075,88151,220380,440750,881510,2203800,4407500,8815100,22038000,44075000]
Ns = [10,5,6,5,5,6,5,5,6,5,5,6,5]

timeSteps = [t for dt,N in zip(dts,Ns) for t in N*[dt]]
times = [sum(timeSteps[:i]) for i in range(len(timeSteps))]

# A carga deve ser aplicada até 440750s, mas temos apenas 440753.5.
# Contudo devemos aplicar até t<440750s pois o tempo atual leva até o próximo +dt

def geomechanics(problemData):
	print(meshName)
	propertyData 	 = problemData.propertyData
	timeStep 		 = problemData.timeStep
	grid 			 = problemData.grid
	numberOfVertices = grid.numberOfVertices
	dimension 		 = grid.dimension

	csvSaver = PyEFVLib.CsvSaver(grid, problemData.outputFilePath, problemData.libraryPath, fileName=meshName)
	meshioSaver = PyEFVLib.MeshioSaver(grid, problemData.outputFilePath, problemData.libraryPath, extension="xdmf", fileName=meshName)

	uField    = np.repeat(0.0, dimension*numberOfVertices)
	oldUField = np.concatenate((problemData.initialValues["u_x"], problemData.initialValues["u_y"], problemData.initialValues["u_z"])) if dimension==3 else np.concatenate((problemData.initialValues["u_x"], problemData.initialValues["u_y"]))

	pField    = np.repeat(0.0, numberOfVertices)
	oldPField = problemData.initialValues["p"].copy()

	matrix 		= np.zeros(((1+dimension)*numberOfVertices, (1+dimension)*numberOfVertices))

	loadBoundary = [b for b in grid.boundaries if b.name=="Load"][0]
	lastFacet = max(loadBoundary.facets,key=lambda f:f.centroid.x)
	of0,of1 = sorted(lastFacet.outerFaces,key=lambda o:o.centroid.x)
	v0,v1 = of0.vertex,of1.vertex
	takeFromBoth = (lastFacet.centroid.x > load_width)
	if takeFromBoth:
		r0 = (lastFacet.centroid.x-load_width)/(lastFacet.centroid.x-v0.x)
		r1 = 1.0
	else:
		r0 = 0.0
		r1 = (v1.x-load_width)/(v1.x-lastFacet.centroid.x)

	global lastLoadForce
	lastLoadForce = 0.0

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
					m = element.vertices.size
					transposedVoigtArea = getTransposedVoigtArea(face)
					shapeFunctions = face.getShapeFunctions()
					identityShapeFunctionMatrix = np.array([shapeFunctions, shapeFunctions, shapeFunctions, np.zeros(m), np.zeros(m), np.zeros(m)]) if dimension == 3 else np.array([shapeFunctions, shapeFunctions, np.zeros(m)])
					matrixCoefficients = (-1) * alpha * np.matmul(transposedVoigtArea, identityShapeFunctionMatrix)
					backwardsHandle, forwardHandle = face.getNeighborVerticesHandles()
					for coord in range(dimension):
						for local, vertex in enumerate(element.vertices):
							matrix[backwardsHandle+(coord+0)*numberOfVertices][vertex.handle+(dimension)*numberOfVertices] += matrixCoefficients[coord][local]
							matrix[forwardHandle+(coord+0)*numberOfVertices][vertex.handle+(dimension)*numberOfVertices] += -matrixCoefficients[coord][local]

				# (-1) * (k / mu) * grad(p)
				for face in element.innerFaces:
					area = face.area.getCoordinates()[:dimension]
					matrixCoefficients = (-1) * (k / mu) * np.matmul( area.T, face.globalDerivatives )
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
		# return inverseMatrix
		return matrix

	def assembleIndependent():
		global lastLoadForce

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

				# (-1) * (k / mu) * rho * g
				for innerFace in element.innerFaces:
					area = innerFace.area.getCoordinates()[:dimension]
					coefficient = (-1) * (k / mu) * rho * g
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
					independent[outerFace.vertex.handle+(1)*numberOfVertices] += bCondition.getValue(outerFace.handle, time=currentTime) * np.linalg.norm(outerFace.area.getCoordinates())
		
		bCondition = [bc for bc in problemData.neumannBoundaries["u_y"] if bc.boundary.name=="Load"][0]
		lastLoadForce = bCondition.getValue(of0.handle, time=currentTime)
		
		if __custom:
			independent[of0.vertex.handle+(1)*numberOfVertices] -= lastLoadForce * r0 * np.linalg.norm(of0.area.getCoordinates())
			independent[of1.vertex.handle+(1)*numberOfVertices] -= lastLoadForce * r1 * np.linalg.norm(of1.area.getCoordinates())

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

	print("{:>9}	{:>14}	{:>14}	{:>14}	{:>14}".format("Iteration", "Current Time", "Time Step", "Difference", "Load"))

	while not converged:
		timeStep = timeSteps[iteration]
		# currentTime = times[iteration]

		matrix = assembleMatrix()
		independent = assembleIndependent()

		results = np.linalg.solve(matrix, independent)
		uField = results[(0)*numberOfVertices:(0+dimension)*numberOfVertices]
		pField = results[(dimension)*numberOfVertices:(dimension+1)*numberOfVertices]

		difference = max( max(abs(pField - oldPField)), max(abs(uField - oldUField)) )

		oldPField = pField.copy()
		oldUField = uField.copy()

		currentTime += timeStep
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

		print("{:>9}	{:>14.2e}	{:>14.2e}	{:>14.2e}	{:>14.2e}".format(iteration, currentTime, timeStep, difference, lastLoadForce))
		converged = ( difference <= tolerance ) or ( currentTime >= problemData.finalTime ) or ( iteration >= problemData.maxNumberOfIterations )
		converged = (iteration > 1) and converged

		if iteration >= problemData.maxNumberOfIterations:
			break

	csvSaver.finalize()
	meshioSaver.finalize()

def main_run():
	problemData = PyEFVLib.ProblemData(
		meshFilePath = meshFilePath,
		outputFilePath = f"{{RESULTS}}/geomechanics/2strip_footing - {['std','custom'][int(__custom)]}/{meshName}",
		numericalSettings = PyEFVLib.NumericalSettings( timeStep = None, maxNumberOfIterations = len(times)-1 ),
		propertyData = PyEFVLib.PropertyData({
			'Body':
			{
				'nu': 0.3,
				'G': 7.6923e+05,
				'cs': 0.0,
				'phi': 1.0,
				'k': 1.1798e-17,
				'cf': 0.0,
				'mu': 1e-3,
				'rhos': 1702.0,
				'rhof': 1000.0,
			},
		}),
		boundaryConditions = PyEFVLib.BoundaryConditions({
			'u_x': {
				'InitialValue': 0.0,
				'West': { 'condition' : PyEFVLib.Dirichlet, 'type' : PyEFVLib.Constant, 'value' : 0.0 },
				'East': { 'condition' : PyEFVLib.Dirichlet, 'type' : PyEFVLib.Constant, 'value' : 0.0 },
				'South': { 'condition' : PyEFVLib.Dirichlet, 'type' : PyEFVLib.Constant, 'value' : 0.0 }, # ????
				'North': { 'condition' : PyEFVLib.Neumann, 'type' : PyEFVLib.Constant, 'value' : 0.0 },
				'Load': { 'condition' : PyEFVLib.Neumann, 'type' : PyEFVLib.Constant, 'value' : 0.0 },
			},
			'u_y': {
				'InitialValue': 0.0,
				'West': { 'condition' : PyEFVLib.Neumann, 'type' : PyEFVLib.Constant, 'value' : 0.0 },
				'East': { 'condition' : PyEFVLib.Neumann, 'type' : PyEFVLib.Constant, 'value' : 0.0 },
				'South': { 'condition' : PyEFVLib.Dirichlet, 'type' : PyEFVLib.Constant, 'value' : 0.0 },
				'North': { 'condition' : PyEFVLib.Neumann, 'type' : PyEFVLib.Constant, 'value' : 0.0 },
				# 'Load': { 'condition' : PyEFVLib.Neumann, 'type' : PyEFVLib.Variable, 'value' : f"{load}" },
				'Load': { 'condition' : PyEFVLib.Neumann, 'type' : PyEFVLib.Variable, 'value' : f"{load} * (t/{t_load}) if t < {t_load} else {load}" },
			},
			'p': {
				'InitialValue': 0.0,
				'West': { 'condition' : PyEFVLib.Neumann, 'type' : PyEFVLib.Constant, 'value' : 0.0 },
				'East': { 'condition' : PyEFVLib.Neumann, 'type' : PyEFVLib.Constant, 'value' : 0.0 },
				'South': { 'condition' : PyEFVLib.Neumann, 'type' : PyEFVLib.Constant, 'value' : 0.0 },
				'North': { 'condition' : PyEFVLib.Dirichlet, 'type' : PyEFVLib.Constant, 'value' : 0.0 },
				'Load': { 'condition' : PyEFVLib.Neumann, 'type' : PyEFVLib.Constant, 'value' : 0.0 },
			},
		}),
	)

	geomechanics( problemData )

if run and __name__ == '__main__':
	__custom = True
	meshName = "strip_footing_super_fine"
	meshFilePath = f"{{MESHES}}/msh/2D/{meshName}.msh"
	main_run()
	# for __custom in [True,False]:
	# 	for meshName in ["2strip_footing_eq_spacing","2strip_footing_uneq_spacing_x","2strip_footing_uneq_spacing_xy"]:
	# 		meshFilePath = f"{{MESHES}}/msh/2D/{meshName}.msh"
	# 		main_run()


import matplotlib.pyplot as plt
import pandas as pd

kPa = 1000
mm = 1e-3

dirname = os.path.realpath( os.path.dirname(__file__) )
def main_plot():
	filePath = dirname+f"/../../../results/geomechanics/strip_footing/{meshName}/{meshName}.csv"
	df = pd.read_csv(filePath)
	X = df["X"]
	Y = df["Y"]
	numberOfTimeSteps = int(df.keys()[-1].split(" - ")[0].replace("TimeStep",""))

	# Pressures

	fig1, ax1 = plt.subplots(2,2)
	fig1.subplots_adjust(hspace=0.4)

	# time = timeStep*np.arange(1,numberOfTimeSteps+1)
	time = times[1:1+numberOfTimeSteps]
	dst2 = lambda x,y,x0,y0:(x-x0)**2+(y-y0)**2
	idxA = min([(idx,dst2(x,y,0.0,height)) for idx,(x,y) in enumerate(zip(X,Y))],key=lambda t:t[1])[0]
	idxB = min([(idx,dst2(x,y,0.0,height-load_width)) for idx,(x,y) in enumerate(zip(X,Y))],key=lambda t:t[1])[0]
	idxC = min([(idx,dst2(x,y,load_width,height-load_width)) for idx,(x,y) in enumerate(zip(X,Y))],key=lambda t:t[1])[0]
	idxD = min([(idx,dst2(x,y,load_width,height)) for idx,(x,y) in enumerate(zip(X,Y))],key=lambda t:t[1])[0]

	pA = np.array([df[f"TimeStep{s} - p"][idxA] for s in range(1,numberOfTimeSteps+1)])
	ax1[0][0].scatter(time, pA/kPa, marker='.', color='k')
	ax1[0][0].set_title("A")
	ax1[0][0].set_ylabel("Pressure (kPa)")
	ax1[0][0].set_xlabel("Time (s)")
	ax1[0][0].semilogx()

	pB = np.array([df[f"TimeStep{s} - p"][idxB] for s in range(1,numberOfTimeSteps+1)])
	ax1[1][0].scatter(time, pB/kPa, marker='.', color='k')
	ax1[1][0].set_title("B")
	ax1[1][0].set_ylabel("Pressure (kPa)")
	ax1[1][0].set_xlabel("Time (s)")
	ax1[1][0].semilogx()

	pC = np.array([df[f"TimeStep{s} - p"][idxC] for s in range(1,numberOfTimeSteps+1)])
	ax1[1][1].scatter(time, pC/kPa, marker='.', color='k')
	ax1[1][1].set_title("C")
	ax1[1][1].set_ylabel("Pressure (kPa)")
	ax1[1][1].set_xlabel("Time (s)")
	ax1[1][1].semilogx()

	pD = np.array([df[f"TimeStep{s} - p"][idxD] for s in range(1,numberOfTimeSteps+1)])
	ax1[0][1].scatter(time, pD/kPa, marker='.', color='k')
	ax1[0][1].set_title("D")
	ax1[0][1].set_ylabel("Pressure (kPa)")
	ax1[0][1].set_xlabel("Time (s)")
	ax1[0][1].semilogx()


	# X Displacements

	fig2, ax2 = plt.subplots(2,2)
	fig2.subplots_adjust(hspace=0.4)

	uA = np.array([df[f"TimeStep{s} - u_x"][idxA] for s in range(1,numberOfTimeSteps+1)])
	ax2[0][0].scatter(time, uA/mm, marker='.', color='k')
	ax2[0][0].set_title("A")
	ax2[0][0].set_ylabel("X Displacement (mm)")
	ax2[0][0].set_xlabel("Time (s)")
	ax2[0][0].set_ylim((-14e-4,14e-4))	
	ax2[0][0].semilogx()

	uB = np.array([df[f"TimeStep{s} - u_x"][idxB] for s in range(1,numberOfTimeSteps+1)])
	ax2[1][0].scatter(time, uB/mm, marker='.', color='k')
	ax2[1][0].set_title("B")
	ax2[1][0].set_ylabel("X Displacement (mm)")
	ax2[1][0].set_xlabel("Time (s)")
	ax2[1][0].set_ylim((-14e-4,14e-4))	
	ax2[1][0].semilogx()

	uC = np.array([df[f"TimeStep{s} - u_x"][idxC] for s in range(1,numberOfTimeSteps+1)])
	ax2[1][1].scatter(time, uC/mm, marker='.', color='k')
	ax2[1][1].set_title("C")
	ax2[1][1].set_ylabel("X Displacement (mm)")
	ax2[1][1].set_xlabel("Time (s)")
	ax2[1][1].semilogx()

	uD = np.array([df[f"TimeStep{s} - u_x"][idxD] for s in range(1,numberOfTimeSteps+1)])
	ax2[0][1].scatter(time, uD/mm, marker='.', color='k')
	ax2[0][1].set_title("D")
	ax2[0][1].set_ylabel("X Displacement (mm)")
	ax2[0][1].set_xlabel("Time (s)")
	ax2[0][1].semilogx()




	# Y Displacements

	fig3, ax3 = plt.subplots(2,2)
	fig3.subplots_adjust(hspace=0.4)

	vA = np.array([df[f"TimeStep{s} - u_y"][idxA] for s in range(1,numberOfTimeSteps+1)])
	ax3[0][0].scatter(time, vA/mm, marker='.', color='k')
	ax3[0][0].set_title("A")
	ax3[0][0].set_ylabel("Y Displacement (mm)")
	ax3[0][0].set_xlabel("Time (s)")
	ax3[0][0].semilogx()

	vB = np.array([df[f"TimeStep{s} - u_y"][idxB] for s in range(1,numberOfTimeSteps+1)])
	ax3[1][0].scatter(time, vB/mm, marker='.', color='k')
	ax3[1][0].set_title("B")
	ax3[1][0].set_ylabel("Y Displacement (mm)")
	ax3[1][0].set_xlabel("Time (s)")
	ax3[1][0].semilogx()

	vC = np.array([df[f"TimeStep{s} - u_y"][idxC] for s in range(1,numberOfTimeSteps+1)])
	ax3[1][1].scatter(time, vC/mm, marker='.', color='k')
	ax3[1][1].set_title("C")
	ax3[1][1].set_ylabel("Y Displacement (mm)")
	ax3[1][1].set_xlabel("Time (s)")
	ax3[1][1].semilogx()

	vD = np.array([df[f"TimeStep{s} - u_y"][idxD] for s in range(1,numberOfTimeSteps+1)])
	ax3[0][1].scatter(time, vD/mm, marker='.', color='k')
	ax3[0][1].set_title("D")
	ax3[0][1].set_ylabel("Y Displacement (mm)")
	ax3[0][1].set_xlabel("Time (s)")
	ax3[0][1].semilogx()



	plt.show()

if plot and __name__ == "__main__":
	# main_plot()

	sys.path.append(dirname+f"/../../../results/geomechanics/strip_footing")
	import plot