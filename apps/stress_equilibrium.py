import sys,os
sys.path.append(os.path.join(os.path.dirname(__file__), os.path.pardir))

from PyEFVLib import MSHReader, Grid, ProblemData, CsvSaver, CgnsSaver
import numpy as np
from scipy import sparse
import scipy.sparse.linalg

#-------------------------SETTINGS----------------------------------------------
problemData = ProblemData( "stress_equilibrium_2d/linear" )

reader = MSHReader(problemData.paths["Grid"])
grid = Grid(reader.getData())
problemData.grid = grid
problemData.read()

# from PyEFVLib.boundaryConditionPrinter import stressEquilibriumBoundaryConditionsPrinter
# stressEquilibriumBoundaryConditionsPrinter(problemData.boundaryConditions)

saver = CsvSaver(grid, problemData.paths["Output"], problemData.libraryPath)

currentTime = 0.0
numberOfVertices = grid.vertices.size
displacements = np.repeat(problemData.initialValue["u"], 2*numberOfVertices)
prevDisplacements = np.repeat(problemData.initialValue["u"], numberOfVertices)

coords,matrixVals = [], []
#---------------------------HELPER FUNCTIONS------------------------------------
def add(i, j, val):
	global coords, matrixVals
	coords.append((i,j))
	matrixVals.append(val)

def getConstitutiveMatrix(region):
	shearModulus = problemData.propertyData[region.handle]["ShearModulus"]
	poissonsRatio = problemData.propertyData[region.handle]["PoissonsRatio"]

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
independent = np.zeros(2*numberOfVertices)

U = lambda handle: handle + numberOfVertices * 0
V = lambda handle: handle + numberOfVertices * 1

# Gravity Term
for region in grid.regions:
	density = problemData.propertyData[region.handle]["Density"]
	gravity = problemData.propertyData[region.handle]["Gravity"]
	for element in region.elements:
		local = 0
		for vertex in element.vertices:
			independent[V(vertex.handle)] += - density * gravity * element.subelementVolumes[local]
			local += 1

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
			
			local=0
			for vertex in element.vertices:
				for neighborVertex in [backwardVertexHandle, forwardVertexHandle]:
					add( U(neighborVertex), U(vertex.handle), matrixCoefficient[0][0][local] )
					add( U(neighborVertex), V(vertex.handle), matrixCoefficient[0][1][local] )
					add( V(neighborVertex), U(vertex.handle), matrixCoefficient[1][0][local] )
					add( V(neighborVertex), V(vertex.handle), matrixCoefficient[1][1][local] )

					matrixCoefficient *= -1
				local+=1

# Boundary Conditions
for bc in problemData.boundaryConditions:
	boundary=bc["u"].boundary
	uBoundaryType = bc["u"].__type__
	vBoundaryType = bc["v"].__type__


	if uBoundaryType == "NEUMANN":
		for facet in boundary.facets:
			for outerFace in facet.outerFaces:
				independent[U(outerFace.vertex.handle)] -= bc["u"].getValue(outerFace.handle) * np.linalg.norm(outerFace.area.getCoordinates())

	if vBoundaryType == "NEUMANN":
		for facet in boundary.facets:
			for outerFace in facet.outerFaces:
				independent[V(outerFace.vertex.handle)] -= bc["v"].getValue(outerFace.handle) * np.linalg.norm(outerFace.area.getCoordinates())

	if uBoundaryType == "DIRICHLET":
		for vertex in boundary.vertices:
			independent[vertex.handle] = bc["u"].getValue(vertex.handle)
			matrixVals = [val for coord, val in zip(coords, matrixVals) if coord[0] != vertex.handle]
			coords 	   = [coord for coord in coords if coord[0] != vertex.handle]
			add(vertex.handle, vertex.handle, 1.0)


	if vBoundaryType == "DIRICHLET":
		for vertex in boundary.vertices:
			independent[vertex.handle+numberOfVertices] = bc["v"].getValue(vertex.handle)
			matrixVals = [val for coord, val in zip(coords, matrixVals) if coord[0] != vertex.handle+numberOfVertices]
			coords 	   = [coord for coord in coords if coord[0] != vertex.handle+numberOfVertices]
			add(vertex.handle+numberOfVertices, vertex.handle+numberOfVertices, 1.0)



#-------------------------SOLVE LINEAR SYSTEM-------------------------------
matrix = sparse.coo_matrix( (matrixVals, zip(*coords)) )
matrix = sparse.csc_matrix( matrix )
displacements = scipy.sparse.linalg.spsolve(matrix, independent)

#-------------------------SAVE RESULTS--------------------------------------
saver.save('u', displacements[:numberOfVertices], currentTime)
saver.save('v', displacements[numberOfVertices:], currentTime)
saver.finalize()

print("\n\t\033[1;35mresult:\033[0m", saver.outputPath, '\n')
# os.system("/usr/bin/paraview %sResults.cgns" % problemData.paths["Output"])

from matplotlib import pyplot as plt, colors, cm
def show_1d(fieldValues, name):
	top_stress = problemData.boundaryConditionData["v"]["NORTH"]["value"]
	shearModulus = problemData.propertyData[region.handle]["ShearModulus"]
	poissonsRatio = problemData.propertyData[region.handle]["PoissonsRatio"]
	lameParameter = 2*shearModulus*poissonsRatio/(1-2*poissonsRatio)
	density = problemData.propertyData[region.handle]["Density"]
	gravity = problemData.propertyData[region.handle]["Gravity"]
	height = 1.0

	y, vals = zip(*[ (vertex.getCoordinates()[1], val) for vertex, val in zip(grid.vertices, fieldValues) if 0.1 > np.abs(vertex.getCoordinates()[0]-0.5)])
	y, vals = zip(*( sorted( zip(y, vals), key=lambda p:p[0] ) ))
	y, vals = np.array(y), np.array(vals)
	
	a_vals=y*(top_stress+density*gravity*(height-y/2))/(2*shearModulus+lameParameter)
	plt.figure()
	plt.scatter(y,1000*vals, marker='.', color='k', label="Resultados Numéricos")
	plt.plot(y,1000*a_vals, label="Solução Analítica")
	plt.xlabel("X (m)")
	plt.ylabel("v (mm)")
	plt.legend()	
	plt.title(name)

# show_1d(displacements[numberOfVertices:], "Deslocamento em y")
# plt.show()
