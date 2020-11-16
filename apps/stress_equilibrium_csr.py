import sys,os
sys.path.append(os.path.join(os.path.dirname(__file__), os.path.pardir))

from PyEFVLib import MSHReader, Grid, ProblemData, CsvSaver, CgnsSaver
from PyEFVLib.simulation import LinearSystem as ls
import numpy as np
from scipy import sparse
from scipy.sparse.linalg import spsolve
import scipy.sparse.linalg

#-------------------------SETTINGS----------------------------------------------
def stressEquilibrium(
		libraryPath,			# PyEFVLib path
		outputPath,				# Results directory path (Ex.: "results/heat_transfer_2d/...")
		extension,				# Extension type. Either "csv" or "cgns"

		grid,					# Object of class Grid
		propertyData,			# List of dictionaries containing the properties

		neumannBoundaries,		# Dictionary whose keys are the field names, and values are objects of the class NeumannBoundaryCondition
		dirichletBoundaries,	# Dictionary whose keys are the field names, and values are objects of the class DirichletBoundaryCondition
		boundaryConditions,		# List of dictionaries whose keys are field names and values are BoundaryCondition objects

		fileName="Results",		# File name
		verbosity=True 			# If False does not print iteration info
	):

	# from PyEFVLib.boundaryConditionPrinter import stressEquilibriumBoundaryConditionsPrinter
	# stressEquilibriumBoundaryConditionsPrinter(problemData.boundaryConditions)

	savers = {"csv": CsvSaver, "cgns": CgnsSaver}
	saver = savers[extension](grid, outputPath, libraryPath, fileName=fileName)

	currentTime = 0.0
	numberOfVertices = grid.vertices.size
	displacements = np.repeat(0.0, 2*numberOfVertices)

	#---------------------------HELPER FUNCTIONS------------------------------------


	def getConstitutiveMatrix(region):
		shearModulus = propertyData[region.handle]["ShearModulus"]
		poissonsRatio = propertyData[region.handle]["PoissonsRatio"]

		lameParameter=2*shearModulus*poissonsRatio/(1-2*poissonsRatio)
		constitutiveMatrix = np.array([[lameParameter*(1-poissonsRatio)/poissonsRatio,lameParameter 							   ,0			],
									   [lameParameter								 ,lameParameter*(1-poissonsRatio)/poissonsRatio,0			],
									   [0			 								 ,0											   ,shearModulus]])
		return constitutiveMatrix

	def getTransposedVoigtArea(face):
		Sx, Sy, _ = face.area.getCoordinates()
		return np.array([[Sx,0,Sy],[0,Sy,Sx]])

	def getVoigtGradientOperator(globalDerivatives):
		Nx,Ny = globalDerivatives
		zero=np.zeros(Nx.size)
		return np.array([[Nx,zero],[zero,Ny],[Ny,Nx]])

	def getOuterFaceGlobalDerivatives(outerFace):
		localDerivatives = outerFace.facet.element.shape.vertexShapeFunctionDerivatives[ outerFace.vertexLocalIndex ]
		return outerFace.facet.element.getGlobalDerivatives(localDerivatives)

	#-------------------------ADD TO LINEAR SYSTEM------------------------------
	ls_csr = ls.LinearSystemCSR(grid.stencil, 2)
	ls_csr.initialize()

	U = lambda handle: handle + numberOfVertices * 0
	V = lambda handle: handle + numberOfVertices * 1

	# Gravity Term
	for region in grid.regions:
		density = propertyData[region.handle]["Density"]
		gravity = propertyData[region.handle]["Gravity"]
		for element in region.elements:
			for local, vertex in enumerate(element.vertices):
				ls_csr.addValueToRHS(V(vertex.handle), - density * gravity * element.subelementVolumes[local])

	# Stress Term
	for region in grid.regions:
		constitutiveMatrix = getConstitutiveMatrix(region)
		for element in region.elements:
			for innerFace in element.innerFaces:
				transposedVoigtArea = getTransposedVoigtArea(innerFace)
				voigtGradientOperator = getVoigtGradientOperator(innerFace.globalDerivatives)

				matrixCoefficient = np.einsum("ij,jk,kmn->imn", transposedVoigtArea, constitutiveMatrix, voigtGradientOperator)
				backwardVertexHandle = element.vertices[element.shape.innerFaceNeighborVertices[innerFace.local][0]].handle
				forwardVertexHandle = element.vertices[element.shape.innerFaceNeighborVertices[innerFace.local][1]].handle

				for local, vertex in enumerate(element.vertices):
					for neighborVertex in [backwardVertexHandle, forwardVertexHandle]:
						ls_csr.addValueToMatrix(U(neighborVertex), U(vertex.handle), matrixCoefficient[0][0][local])
						ls_csr.addValueToMatrix(U(neighborVertex), V(vertex.handle), matrixCoefficient[0][1][local])
						ls_csr.addValueToMatrix(V(neighborVertex), U(vertex.handle), matrixCoefficient[1][0][local])
						ls_csr.addValueToMatrix(V(neighborVertex), V(vertex.handle), matrixCoefficient[1][1][local])

						matrixCoefficient *= -1

	# Boundary Conditions
	for bc in boundaryConditions:
		boundary=bc["u"].boundary
		uBoundaryType = bc["u"].__type__
		vBoundaryType = bc["v"].__type__


		if uBoundaryType == "NEUMANN":
			for facet in boundary.facets:
				for outerFace in facet.outerFaces:
					ls_csr.addValueToRHS(U(outerFace.vertex.handle), - bc["u"].getValue(outerFace.handle) * np.linalg.norm(outerFace.area.getCoordinates()))

		if vBoundaryType == "NEUMANN":
			for facet in boundary.facets:
				for outerFace in facet.outerFaces:
					ls_csr.addValueToRHS(V(outerFace.vertex.handle), - bc["v"].getValue(outerFace.handle) * np.linalg.norm(outerFace.area.getCoordinates()))

		if uBoundaryType == "DIRICHLET":
			for vertex in boundary.vertices:
				ls_csr.setValueToRHS(U(vertex.handle), bc["u"].getValue(vertex.handle))
				ls_csr.matZeroRow(U(vertex.handle), 1.0)


		if vBoundaryType == "DIRICHLET":
			for vertex in boundary.vertices:
				ls_csr.setValueToRHS(V(vertex.handle), bc["v"].getValue(vertex.handle))
				ls_csr.matZeroRow(V(vertex.handle), 1.0)



	#-------------------------SOLVE LINEAR SYSTEM-------------------------------
	displacements = spsolve(ls_csr.matrix, ls_csr.rhs)

	#-------------------------SAVE RESULTS--------------------------------------
	saver.save('u', displacements[:numberOfVertices], currentTime)
	saver.save('v', displacements[numberOfVertices:], currentTime)
	saver.finalize()

	print("\n\t\033[1;35mresult:\033[0m", saver.outputPath, '\n')

	return displacements


if __name__ == "__main__":
	model = "workspace/stress_equilibrium_2d/linear"

	problemData = ProblemData(model)

	reader = MSHReader(problemData.paths["Grid"])
	grid = Grid(reader.getData())
	grid.buildStencil()
	problemData.setGrid(grid)
	problemData.read()

	displacements = stressEquilibrium(
		libraryPath = problemData.libraryPath,
		outputPath = problemData.paths["Output"],
		extension = "csv" if not "--extension=cgns" in sys.argv else "cgns",

		grid 	  = grid,
		propertyData = problemData.propertyData,

		neumannBoundaries = problemData.neumannBoundaries,
		dirichletBoundaries = problemData.dirichletBoundaries,
		boundaryConditions = problemData.boundaryConditions,

		verbosity = not "-s" in sys.argv
	)

	#-------------------------------------------------------------------------------
	#-------------------------SHOW RESULTS GRAPHICALY-------------------------------
	#-------------------------------------------------------------------------------
	from matplotlib import pyplot as plt, colors, cm
	def show_1d(fieldValues, name):
		top_stress = problemData.boundaryConditionData["v"]["North"]["value"]
		shearModulus = problemData.propertyData[0]["ShearModulus"]
		poissonsRatio = problemData.propertyData[0]["PoissonsRatio"]
		lameParameter = 2*shearModulus*poissonsRatio/(1-2*poissonsRatio)
		density = problemData.propertyData[0]["Density"]
		gravity = problemData.propertyData[0]["Gravity"]
		height = 1.0

		y, vals = zip(*[ (vertex.getCoordinates()[1], val) for vertex, val in zip(grid.vertices, fieldValues) if 0.1 > np.abs(vertex.getCoordinates()[0]-0.5)])
		y, vals = zip(*( sorted( zip(y, vals), key=lambda p:p[0] ) ))
		y, vals = np.array(y), np.array(vals)

		a_vals=y*(top_stress+density*gravity*(height-y/2))/(2*shearModulus+lameParameter)
		plt.figure()
		plt.scatter(1000*vals, y, marker='.', color='k', label="Resultados Numéricos")
		plt.plot(1000*a_vals, y, label="Solução Analítica")
		plt.ylabel("X (m)")
		plt.xlabel("v (mm)")
		plt.grid(True)
		plt.legend()
		plt.title(name)

	show_1d(displacements[grid.vertices.size:], "Deslocamento em y")
	plt.show()
