print("""
	TESTE
""")

load = 1.0e+6
gravity = -9.81 #+

timeStep = 10
finalTime = 1000

height = 6.0
meshFilePath = "{MESHES}/msh/2D/column6.msh"

run = True
plot = True

stepInTimeSamples = [1, 2, 10, 25, 50, 99]

import sys,os
sys.path.append(os.path.join(os.path.dirname(__file__), *[os.path.pardir]*3))
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

	oldUField = np.concatenate((problemData.initialValues["u_x"], problemData.initialValues["u_y"], problemData.initialValues["u_z"])) if dimension==3 else np.concatenate((problemData.initialValues["u_x"], problemData.initialValues["u_y"]))
	uField    = oldUField.copy()

	oldPField = problemData.initialValues["p"].copy()
	pField    = oldPField.copy()

	matrix 		= np.zeros(((1+dimension)*numberOfVertices, (1+dimension)*numberOfVertices))

	meshioSaver.save("u_x", oldUField[0*numberOfVertices:1*numberOfVertices], 0.0)
	meshioSaver.save("u_y", oldUField[1*numberOfVertices:2*numberOfVertices], 0.0)
	if dimension == 3:
		meshioSaver.save("u_z", oldUField[2*numberOfVertices:3*numberOfVertices], 0.0)
	meshioSaver.save("p", oldPField, 0.0)

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

			g = np.array([0.0, -gravity, 0.0])[:dimension]
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
		inverseMatrix = np.linalg.inv(matrix)
		return inverseMatrix

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

			g = np.array([0.0, -gravity, 0.0])[:dimension]
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

	inverseMatrix = assembleMatrix()

	while not converged:
		independent = assembleIndependent(currentTime)

		results = np.matmul(inverseMatrix, independent)
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

		print("{:>9}	{:>14.2e}	{:>14.2e}	{:>14.2e}".format(iteration, currentTime, timeStep, difference))
		converged = ( difference <= tolerance ) or ( currentTime >= problemData.finalTime ) or ( iteration >= problemData.maxNumberOfIterations )

		if iteration >= problemData.maxNumberOfIterations:
			break

	csvSaver.finalize()
	meshioSaver.finalize()

	return pField, uField[0*numberOfVertices:1*numberOfVertices], uField[1*numberOfVertices:2*numberOfVertices]

def main_run():
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
	})

	boundaryConditionsDict = {
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
			'North': { 'condition' : PyEFVLib.Neumann, 'type' : PyEFVLib.Constant, 'value' : load },
		},
		'p': {
			'InitialValue': 0.0,
			'West': { 'condition' : PyEFVLib.Neumann, 'type' : PyEFVLib.Constant, 'value' : 0.0 },
			'East': { 'condition' : PyEFVLib.Neumann, 'type' : PyEFVLib.Constant, 'value' : 0.0 },
			'South': { 'condition' : PyEFVLib.Neumann, 'type' : PyEFVLib.Constant, 'value' : 0.0 },
			'North': { 'condition' : PyEFVLib.Neumann, 'type' : PyEFVLib.Constant, 'value' : 0.0 },
		},
	}

	problemData = PyEFVLib.ProblemData(
		meshFilePath = meshFilePath,
		outputFilePath = "{RESULTS}/geomechanics",
		numericalSettings = PyEFVLib.NumericalSettings( timeStep = 5*timeStep, tolerance=1e-2, maxNumberOfIterations = 200 ),
		propertyData = propertyData,
		boundaryConditions = PyEFVLib.BoundaryConditions(boundaryConditionsDict),
	)

	initialPField, initialUXField, initialUYField = geomechanics( problemData )

	boundaryConditionsDict['p']['InitialValue'] = initialPField
	boundaryConditionsDict['u_x']['InitialValue'] = initialUXField
	boundaryConditionsDict['u_y']['InitialValue'] = initialUYField
	boundaryConditionsDict['p']['North']['condition'] = PyEFVLib.Dirichlet

	problemData = PyEFVLib.ProblemData(
		meshFilePath = meshFilePath,
		outputFilePath = "{RESULTS}/geomechanics",
		numericalSettings = PyEFVLib.NumericalSettings( timeStep = timeStep, tolerance=1e-2, finalTime = finalTime, maxNumberOfIterations = 600 ),
		propertyData = propertyData,
		boundaryConditions = PyEFVLib.BoundaryConditions(boundaryConditionsDict),
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
from solgeom.Terzaghi import Solution

filePath = dirname+"/../../../results/geomechanics/Results.csv"

df = pd.read_csv(filePath)
X = df["X"]
Y = df["Y"]

numberOfTimeSteps = int(df.keys()[-1].split(" - ")[0].replace("TimeStep",""))

mm = 1000.
kPa = 1/1000.

def main_plot():
	rock = json.load(open(dirname+"/../../../../solgeom/examples/solid.json", "r"))
	fluid = json.load(open(dirname+"/../../../../solgeom/examples/fluid.json", "r"))

	terza = Solution(height, load, rock, fluid , gravity)
	z = terza.getPositionValues()

	fig, ax = plt.subplots(2, 2, figsize=(8,7))
	fig.subplots_adjust(left=0.070, right=0.975, top=0.970, bottom=0.065, hspace=0.235, wspace=0.300)

	# Plot pressure profiles ---------------------------------------------------
	times_1 = [N*timeStep for N in stepInTimeSamples]
	for idx,t in zip(stepInTimeSamples,times_1):
		p = terza.getPressureValuesConstTime(t)
		p_n = df[f"TimeStep{idx} - p"]
		ax[0][0].plot(p*kPa, z, label=f"t = {timeStep*idx}s")
		ax[0][0].scatter(p_n*kPa, Y, color="k", marker=".",linewidth=0.5)
	ax[0][0].set_xlabel("Pressure (kPa)", size=12)
	ax[0][0].set_ylabel("Height (m)", size=12)
	ax[0][0].grid(True)
	# ax[0][0].legend()
	# --------------------------------------------------------------------------

	# Plot displacement profiles -----------------------------------------------
	for idx,t in zip(stepInTimeSamples,times_1):
		w = terza.getDisplacementValuesConstTime(t)
		w_n = df[f"TimeStep{idx} - u_y"]
		ax[0][1].plot(w*mm, z)
		ax[0][1].scatter(w_n*mm, Y, color="k", marker=".",linewidth=0.5)
	ax[0][1].set_xlabel("Displacement (mm)", size=12)
	ax[0][1].set_ylabel("Height (m)", size=12)
	ax[0][1].grid(True)
	# --------------------------------------------------------------------------

	# Plot bottom pressure over time -------------------------------------------
	# times_2 = np.linspace(0, 2000., 100)
	times_2 = [timeStep*s for s in range(1,1+numberOfTimeSteps)]
	p_a = terza.getPressureValuesAtPosition(0.0, times_2) # Bottom pressure (z=0.0)
	p_n = np.array([ np.average([p for p,y in zip(df[f"TimeStep{s} - p"],Y) if y<0.05]) for s in range(1,1+numberOfTimeSteps) ])
	ax[1][0].plot(times_2, p_a*kPa)
	ax[1][0].scatter(times_2, p_n*kPa, color="k", marker=".", linewidth=0.5)
	ax[1][0].set_xlabel("Time (s)", size=12)
	ax[1][0].set_ylabel("Bottom Pressure (kPa)", size=12)
	ax[1][0].grid(True)
	# --------------------------------------------------------------------------

	# Plot top displacement over time ------------------------------------------
	w_a = terza.getDisplacementValuesAtPosition(height, times_2) # Top displacement (z=height)
	w_n = np.array([ np.average([w for w,y in zip(df[f"TimeStep{s} - u_y"],Y) if y>height-0.05]) for s in range(1,1+numberOfTimeSteps) ])
	ax[1][1].plot(times_2, w_a*mm)
	ax[1][1].scatter(times_2, w_n*mm, color="k", marker=".", linewidth=0.5)
	ax[1][1].set_xlabel("Time (s)", size=12)
	ax[1][1].set_ylabel("Top Displacement (mm)", size=12)
	ax[1][1].grid(True)
	# --------------------------------------------------------------------------

	plt.show()

if plot and __name__ == '__main__':
	main_plot()