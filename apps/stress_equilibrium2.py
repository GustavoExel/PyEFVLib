import sys,os
sys.path.append(os.path.join(os.path.dirname(__file__), os.path.pardir))

from PyEFVLib import MSHReader, Grid, ProblemData, CgnsSaver
import numpy as np

models = ["stress_equilibrium", "stress_1d", "test_1d"]
#-------------------------SETTINGS----------------------------------------------
problemData = ProblemData( models[1] )

reader = MSHReader(problemData.paths["Grid"])
grid = Grid(reader.getData())
problemData.grid = grid
problemData.read()

cgnsSaver = CgnsSaver(grid, problemData.paths["Output"], problemData.libraryPath)

currentTime = 0.0
numberOfVertices = grid.vertices.size
displacements = np.repeat(problemData.initialValue["u"], 2*numberOfVertices)
prevDisplacements = np.repeat(problemData.initialValue["u"], numberOfVertices)

matrix = np.zeros([2*numberOfVertices, 2*numberOfVertices])
difference = 0.0
iteration = 0
converged = False

#-------------------------LINEAR SYSTEM FUNCTIONS---------------------------
def getVoigtArea(face):
	Sx, Sy, Sz = face.area.getCoordinates()
	return np.array([[Sx,0.0],[0.0,Sy],[Sy,Sx]])

def getConstitutiveMatrix(region):
	G = problemData.propertyData[region.handle]["ShearModulus"]
	poissonsRatio = problemData.propertyData[region.handle]["PoissonsRatio"]

	l=2*G*poissonsRatio/(1-2*poissonsRatio)
	b=l*(1-poissonsRatio)/poissonsRatio

	return np.array([[b,l,0],[l,b,0],[0,0,G]])

def getVoigtGradientOperator(globalDerivatives):
	Nx, Ny = globalDerivatives
	zero = np.zeros( Nx.size )
	return np.array([[Nx,zero],[zero,Ny],[Ny,Nx]])

#-------------------------ADD TO LINEAR SYSTEM------------------------------
U=lambda vertex: vertex + numberOfVertices * 0
V=lambda vertex: vertex + numberOfVertices * 1
independent = np.zeros(2*numberOfVertices)

# Gravity Term
for region in grid.regions:
	volumetricForce = problemData.propertyData[region.handle]["Density"] * problemData.propertyData[region.handle]["Gravity"]
	for element in region.elements:
		local=0
		for vertex in element.vertices:
			independent[V(vertex.handle)] += -volumetricForce * element.subelementVolumes[local]
			local+=1

# Stress Term
for region in grid.regions:
	constitutiveMatrix = getConstitutiveMatrix(region)
	for element in region.elements:
		for innerFace in element.innerFaces:
			voigtArea = getVoigtArea(innerFace)
			voigtGradientOperator = getVoigtGradientOperator(innerFace.globalDerivatives)
			
			coefficient = np.einsum("ji,jk,kmn->imn", voigtArea, constitutiveMatrix, voigtGradientOperator)
			i=0
			for vertex in element.vertices:
				neighborIndexes = [	element.vertices[element.shape.innerFaceNeighborVertices[innerFace.local][0]].handle,
									element.vertices[element.shape.innerFaceNeighborVertices[innerFace.local][1]].handle ]

				for neighborIndex in neighborIndexes:
					matrix[U(neighborIndex)][U(vertex.handle)] += coefficient[0][0][i]
					matrix[U(neighborIndex)][V(vertex.handle)] += coefficient[0][1][i]
					matrix[V(neighborIndex)][U(vertex.handle)] += coefficient[1][0][i]
					matrix[V(neighborIndex)][V(vertex.handle)] += coefficient[1][1][i]
					coefficient *= -1.0
				i+=1

# Boundary Conditions
# for bc in problemData.boundaryConditions:
# 	boundary=bc["u"].boundary
# 	uBoundaryType = bc["u"].__type__
# 	vBoundaryType = bc["v"].__type__

# 	if uBoundaryType == vBoundaryType == "DIRICHLET":
# 		for vertex in boundary.vertices:
# 			matrix[U(vertex.handle)] = np.zeros(2*numberOfVertices)
# 			matrix[V(vertex.handle)] = np.zeros(2*numberOfVertices)
# 			matrix[U(vertex.handle)][U(vertex.handle)] = 1.0
# 			matrix[V(vertex.handle)][V(vertex.handle)] = 1.0
# 			independent[U(vertex.handle)] = u_boundaryCondition
# 			independent[V(vertex.handle)] = v_boundaryCondition

# 	if uBoundaryType == vBoundaryType == "NEUMANN":
# 		pass
# 		for vertex in boundary.vertices:
# 			matrix[U(vertex.handle)] = np.zeros(2*numberOfVertices)
# 			matrix[V(vertex.handle)] = np.zeros(2*numberOfVertices)
# 			independent[U(vertex.handle)] = 0.0
# 			independent[V(vertex.handle)] = 0.0

# 		for facet in boundary.facets: # Same as for element in boundary.elements (but belongs to boundary)
# 			element=facet.element
# 			for outerFace in facet.outerFaces: # Same as for vertex in element.vertices (but belongs to boundary)
# 				vertex=outerFace.vertex

# 				sx,sy,_ = outerFace.area.getCoordinates()

# 				Tx = bc["u"].getValue(vertex.handle)
# 				Ty = bc["v"].getValue(vertex.handle)
# 				xNormalStress=sx*(Tx*sx+Ty*sy)/(sx**2+sy**2)
# 				yNormalStress=sy*(Tx*sx+Ty*sy)/(sx**2+sy**2)
# 				shearStress=((Tx-xNormalStress)**2+(Ty-yNormalStress)**2)**(1/2)

# 				transposedVoigtArea = getVoigtArea(outerFace).T
# 				constitutiveMatrix = getConstitutiveMatrix(element.region)
# 				globalDerivatives = element.getGlobalDerivatives(element.shape.vertexShapeFunctionDerivatives[ list(element.vertices).index(vertex) ])
# 				voigtGradientOperator = getVoigtGradientOperator(globalDerivatives)

# 				matrixCoefficient=np.einsum( "ij,jk,kmn->imn", transposedVoigtArea, constitutiveMatrix, voigtGradientOperator )
# 				independentCoefficient=transposedVoigtArea @ [xNormalStress, yNormalStress, shearStress]

# 				local=0
# 				for vertex in element.vertices:
# 					matrix[U(vertex.handle)][U(vertex.handle)] += matrixCoefficient[0][0][local]
# 					matrix[U(vertex.handle)][V(vertex.handle)] += matrixCoefficient[0][1][local]
# 					matrix[V(vertex.handle)][U(vertex.handle)] += matrixCoefficient[1][0][local]
# 					matrix[V(vertex.handle)][V(vertex.handle)] += matrixCoefficient[1][1][local]

# 					independent[U(vertex.handle)] += independentCoefficient[0]
# 					independent[V(vertex.handle)] += independentCoefficient[1]

# 					local += 1

# 	else:
# 		dirichletLabel = "u" if uBoundaryType=="DIRICHLET" else "v"
# 		neumannLabel = "u" if uBoundaryType=="NEUMANN" else "v"

# 		dirichletIndex = 0 if uBoundaryType=="DIRICHLET" else 1
# 		neumannIndex = 0 if uBoundaryType=="NEUMANN" else 1

# 		D=lambda vertex: int( vertex.handle + numberOfVertices * dirichletIndex )
# 		N=lambda vertex: int( vertex.handle + numberOfVertices * neumannIndex )

# 		for vertex in boundary.vertices:
# 			matrix[D(vertex)] = np.zeros(2*numberOfVertices)
# 			matrix[D(vertex)][D(vertex)] = 1.0
# 			independent[D(vertex)] = bc[dirichletLabel].value

# 			matrix[N(vertex)] = np.zeros(2*numberOfVertices)
# 			independent[N(vertex)] = 0.0

# 		for facet in boundary.facets:
# 			element=facet.element
# 			shearModulus = problemData.propertyData[element.region.handle]["ShearModulus"]
# 			for outerFace in facet.outerFaces:
# 				vertex=outerFace.vertex

# 				sx,sy,_ = outerFace.area.getCoordinates()
# 				Tx = bc["u"].getValue(vertex.handle) if neumannLabel=="u" else 0.0
# 				Ty = bc["v"].getValue(vertex.handle) if neumannLabel=="v" else 0.0

# 				# Note that the way we are calculating the normal stress, the only case
# 				# the result is going to be accurate is when the boundary is orthogonal
# 				# with the tension applied. Otherwise prefer to use the same boundary
# 				# condition formulation (Dirichlet/Neumann) in order to get an accurate Results
# 				normalArea = [sx, sy][ neumannIndex ]
# 				normalStress=normalArea*(Tx*sx+Ty*sy)/(sx**2+sy**2)

# 				transposedVoigtArea = getVoigtArea(outerFace).T
# 				constitutiveMatrix = getConstitutiveMatrix(facet.element.region)
# 				globalDerivatives = element.getGlobalDerivatives(element.shape.vertexShapeFunctionDerivatives[ list(element.vertices).index(vertex) ])
# 				voigtGradientOperator = getVoigtGradientOperator(globalDerivatives)

# 				matrixCoefficient = np.einsum("ij,jk,kmn->imn", transposedVoigtArea, constitutiveMatrix, voigtGradientOperator)

# 				local=0
# 				for vertex in facet.element.vertices:
# 					matrix[N(vertex)][U(vertex.handle)] += matrixCoefficient[neumannIndex][0][local]
# 					matrix[N(vertex)][V(vertex.handle)] += matrixCoefficient[neumannIndex][1][local]

# 					matrix[N(vertex)][U(vertex.handle)] -= sy * shearModulus
# 					matrix[N(vertex)][V(vertex.handle)] -= sx * shearModulus

# 					independent[N(vertex)] += normalArea * normalStress

# 					local += 1


#-------------------------SOLVE LINEAR SYSTEM-------------------------------
displacements = np.linalg.solve(matrix, independent)

#-------------------------SAVE RESULTS--------------------------------------
# ACTUALLY DISPLACEMENTS WILL BE ONE MATRIX
cgnsSaver.save('u', displacements[:numberOfVertices], 0.0)
cgnsSaver.save('v', displacements[numberOfVertices:], 0.0)
cgnsSaver.finalize()

print("\n\033[1;35mresult:\033[0m", problemData.paths["Output"]+"Results.cgns", '\n')
#-------------------------------------------------------------------------------
#-------------------------SHOW RESULTS GRAPHICALY-------------------------------
#-------------------------------------------------------------------------------
if "-g" in sys.argv:
	from matplotlib import pyplot as plt, colors, cm
	from scipy.interpolate import griddata

	X,Y = zip(*[v.getCoordinates()[:-1] for v in grid.vertices])
	Xi, Yi = np.meshgrid( np.linspace(min(X), max(X), len(X)), np.linspace(min(Y), max(Y), len(Y)) )
	def show(fieldValues, name):
		plt.figure()
		gridValues = griddata((X,Y), fieldValues, (Xi,Yi), method='linear')
		plt.pcolor(Xi,Yi,gridValues, cmap="RdBu")#, cmap=colors.ListedColormap( cm.get_cmap("RdBu",256)(np.linspace(1,0,256)) ))
		plt.title(name)
		plt.colorbar()

	def show_1d(fieldValues, name):
		top_stress = problemData.boundaryConditionData["v"]["NORTH"]["value"]
		shear_modulus = problemData.propertyData[region.handle]["ShearModulus"]
		poissonsRatio = problemData.propertyData[region.handle]["PoissonsRatio"]
		lameParameter = 2*shearModulus*poissonsRatio/(1-2*poissonsRatio)

		y, vals = zip(*[ (vertex.getCoordinates()[1], val) for vertex, val in zip(grid.vertices, fieldValues) if 0.1 > np.abs(vertex.getCoordinates()[0]-0.5)])
		y, vals = zip(*( sorted( zip(y, vals), key=lambda p:p[0] ) ))
		y, vals = np.array(y), np.array(vals)
		
		a_vals=	-y*top_stress/(2*shear_modulus+lameParameter)
		print("sum(vals): ", sum(vals)) 
		print("max(vals): ", max(vals)) 
		print("min(vals): ", min(vals)) 
		print("avg(vals): ", sum(vals)/len(vals))
		plt.figure()
		plt.scatter(y,vals, marker='.', color='k')
		plt.plot(y, a_vals)
		plt.legend()	
		plt.title(name)

	# show(displacements[:numberOfVertices], "X displacement (u)")
	# show(displacements[numberOfVertices:], "Y displacement (v)")
	# show(normal_stain_x, "X normal strain (exx)")
	# show(normal_stain_y, "Y normal strain (eyy)")
	# show(shear_strain, "Shear strain (γxy)")
	show_1d(displacements[numberOfVertices:], "Y displacement (v) along beam")

	plt.show()
else:
	try:
		print("Opening Paraview...")
		os.system("/usr/bin/paraview %sResults.cgns" % problemData.paths["Output"])
	except:
		print("Could not launch /usr/bin/paraview")