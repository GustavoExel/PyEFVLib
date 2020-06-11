import sys,os
sys.path.append(os.path.join(os.path.dirname(__file__), os.path.pardir))

from PyEFVLib import MSHReader, Grid, ProblemData, CgnsSaver
import numpy as np

#-------------------------SETTINGS----------------------------------------------
problemData = ProblemData( "stress_equilibrium" )

reader = MSHReader(problemData.paths["Grid"])
grid = Grid(reader.getData())
problemData.grid = grid
problemData.read()

from PyEFVLib.boundaryConditionPrinter import stressEquilibriumBoundaryConditionsPrinter
stressEquilibriumBoundaryConditionsPrinter(problemData.boundaryConditions)

cgnsSaver = CgnsSaver(grid, problemData.paths["Output"], problemData.libraryPath)

currentTime = 0.0
numberOfVertices = grid.vertices.size
displacements = np.repeat(problemData.initialValue["u"], 2*numberOfVertices)
prevDisplacements = np.repeat(problemData.initialValue["u"], numberOfVertices)

matrix = np.zeros([2*numberOfVertices, 2*numberOfVertices])
#---------------------------HELPER FUNCTIONS------------------------------------
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
				matrix[U(backwardVertexHandle)][U(vertex.handle)] += matrixCoefficient[0][0][local]
				matrix[U(backwardVertexHandle)][V(vertex.handle)] += matrixCoefficient[0][1][local]
				matrix[V(backwardVertexHandle)][U(vertex.handle)] += matrixCoefficient[1][0][local]
				matrix[V(backwardVertexHandle)][V(vertex.handle)] += matrixCoefficient[1][1][local]
				matrix[U(forwardVertexHandle)][U(vertex.handle)]  -= matrixCoefficient[0][0][local]
				matrix[U(forwardVertexHandle)][V(vertex.handle)]  -= matrixCoefficient[0][1][local]
				matrix[V(forwardVertexHandle)][U(vertex.handle)]  -= matrixCoefficient[1][0][local]
				matrix[V(forwardVertexHandle)][V(vertex.handle)]  -= matrixCoefficient[1][1][local]

				local+=1

# Boundary Conditions
for bc in problemData.boundaryConditions:
	boundary=bc["u"].boundary
	uBoundaryType = bc["u"].__type__
	vBoundaryType = bc["v"].__type__

	if uBoundaryType != vBoundaryType:
		dirichletLabel = "u" if uBoundaryType=="DIRICHLET" else "v"
		neumannLabel = "u" if uBoundaryType=="NEUMANN" else "v"

		dirichletIndex = 0 if uBoundaryType=="DIRICHLET" else 1
		neumannIndex = 0 if uBoundaryType=="NEUMANN" else 1

		D=lambda vertex: int( vertex.handle + numberOfVertices * dirichletIndex )
		N=lambda vertex: int( vertex.handle + numberOfVertices * neumannIndex )

		for vertex in boundary.vertices:
			matrix[D(vertex)] = np.zeros(2*numberOfVertices)
			matrix[D(vertex)][D(vertex)] = 1.0
			independent[D(vertex)] = bc[dirichletLabel].getValue(vertex.handle)

			matrix[N(vertex)] = np.zeros(2*numberOfVertices)
			independent[N(vertex)] = 0.0

		for facet in boundary.facets:
			shearModulus = problemData.propertyData[region.handle]["ShearModulus"]
			constitutiveMatrix = getConstitutiveMatrix(facet.element.region)
			for outerFace in facet.outerFaces:
				globalDerivatives = getOuterFaceGlobalDerivatives(outerFace)
				Nx, Ny = globalDerivatives

				Sx,Sy,_ = outerFace.area.getCoordinates()

				normalArea = [Sx,Sy][neumannIndex]
				shearArea = [Sx,Sy][dirichletIndex]

				transposedVoigtArea = getTransposedVoigtArea(outerFace)
				voigtGradientOperator = getVoigtGradientOperator(globalDerivatives)
				coefficient = np.einsum("ij,jk,kmn->imn", transposedVoigtArea, constitutiveMatrix, voigtGradientOperator)

				local=0
				for vertex in facet.element.vertices:
					matrix[N(outerFace.vertex)][U(vertex.handle)] += coefficient[neumannIndex][0][local]
					matrix[N(outerFace.vertex)][V(vertex.handle)] += coefficient[neumannIndex][1][local]
					matrix[N(outerFace.vertex)][U(vertex.handle)] -= shearArea * shearModulus * Ny[local]
					matrix[N(outerFace.vertex)][V(vertex.handle)] -= shearArea * shearModulus * Nx[local]
					independent[N(outerFace.vertex)] += normalArea * bc[neumannLabel].getValue( outerFace.vertex.handle )
					local+=1


	if uBoundaryType == "NEUMANN" and vBoundaryType == "NEUMANN":
		for facet in boundary.facets:
			for outerFace in facet.outerFaces:
				independent[U(outerFace.vertex.handle)] = bc["u"].getValue(outerFace.handle) * np.linalg.norm(outerFace.area.getCoordinates())
				independent[V(outerFace.vertex.handle)] = bc["v"].getValue(outerFace.handle) * np.linalg.norm(outerFace.area.getCoordinates())

	if uBoundaryType == "DIRICHLET" and vBoundaryType == "DIRICHLET":
		for vertex in boundary.vertices:
			matrix[U(vertex.handle)] = np.zeros(2*numberOfVertices)
			matrix[V(vertex.handle)] = np.zeros(2*numberOfVertices)
			matrix[U(vertex.handle)][U(vertex.handle)] = 1.0
			matrix[V(vertex.handle)][V(vertex.handle)] = 1.0

			independent[U(vertex.handle)] = bc["u"].getValue(vertex.handle)
			independent[V(vertex.handle)] = bc["v"].getValue(vertex.handle)



	# if uBoundaryType == "NEUMANN":
	# 	for facet in boundary.facets:
	# 		for outerFace in facet.outerFaces:
	# 			independent[U(outerFace.vertex.handle)] = bc["u"].getValue(outerFace.handle) * np.linalg.norm(outerFace.area.getCoordinates())

	# if vBoundaryType == "NEUMANN":
	# 	for facet in boundary.facets:
	# 		for outerFace in facet.outerFaces:
	# 			independent[V(outerFace.vertex.handle)] = bc["v"].getValue(outerFace.handle) * np.linalg.norm(outerFace.area.getCoordinates())

	# if uBoundaryType == "DIRICHLET":
	# 	for vertex in boundary.vertices:
	# 		independent[vertex.handle] = bc["u"].getValue(vertex.handle)
	# 		matrix[vertex.handle] = np.zeros(2*numberOfVertices)
	# 		matrix[vertex.handle][vertex.handle] = 1.0

	# if vBoundaryType == "DIRICHLET":
	# 	for vertex in boundary.vertices:
	# 		independent[vertex.handle+numberOfVertices] = bc["v"].getValue(vertex.handle)
	# 		matrix[vertex.handle+numberOfVertices] = np.zeros(2*numberOfVertices)
	# 		matrix[vertex.handle+numberOfVertices][vertex.handle+numberOfVertices] = 1.0


#-------------------------SOLVE LINEAR SYSTEM-------------------------------
displacements = np.linalg.solve(matrix, independent)

#-------------------------SAVE RESULTS--------------------------------------
cgnsSaver.save('u', displacements[:numberOfVertices], currentTime)
cgnsSaver.save('v', displacements[numberOfVertices:], currentTime)
cgnsSaver.finalize()

print("\n\t\033[1;35mresult:\033[0m", problemData.paths["Output"]+"Results.cgns", '\n')
# os.system("/usr/bin/paraview %sResults.cgns" % problemData.paths["Output"])

from matplotlib import pyplot as plt, colors, cm
def show_1d(fieldValues, name):
	top_stress = problemData.boundaryConditionData["v"]["NORTH"]["value"]
	shearModulus = problemData.propertyData[region.handle]["ShearModulus"]
	poissonsRatio = problemData.propertyData[region.handle]["PoissonsRatio"]
	lameParameter = 2*shearModulus*poissonsRatio/(1-2*poissonsRatio)

	y, vals = zip(*[ (vertex.getCoordinates()[1], val) for vertex, val in zip(grid.vertices, fieldValues) if 0.1 > np.abs(vertex.getCoordinates()[0]-0.5)])
	y, vals = zip(*( sorted( zip(y, vals), key=lambda p:p[0] ) ))
	y, vals = np.array(y), np.array(vals)
	
	a_vals=y*top_stress/(2*shearModulus+lameParameter)
	print("sum(vals): ", sum(vals)) 
	print("max(vals): ", max(vals)) 
	print("min(vals): ", min(vals)) 
	print("avg(vals): ", sum(vals)/len(vals))
	plt.figure()
	plt.scatter(y,vals, marker='.', color='k')
	plt.plot(y, a_vals)
	plt.legend()	
	plt.title(name)

show_1d(displacements[numberOfVertices:], "Y displacement (v) along beam")
plt.show()
