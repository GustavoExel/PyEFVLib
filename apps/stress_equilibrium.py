
import sys,os
sys.path.append(os.path.join(os.path.dirname(__file__), os.path.pardir))

from PyEFVLib import MSHReader, Grid, ProblemData, CgnsSaver
import numpy as np

models = ["stress_equilibrium", "stress_1d"]
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

def show3dArray(arr):
	t=""
	ll=[]
	i=0
	for l in arr:
		lll=[]
		for e in l:
			ttt=', '.join([str(round(ee, 3)) for ee in e])
			ttt="[%s]" % ttt
			lll.append(ttt)
		tt=", ".join(lll)
		tt="[%s]" % tt
		if i != arr.shape[1]:
			tt += ','
		if i != 0:
			tt = "  "+tt

		ll.append(tt)
		i+=1
	t='\n'.join(ll)
	return "[%s]" % t




def computeLocalMatrixCoefficient(innerFace, shearModulus, poissonsRatio):
	Sx, Sy, Sz = innerFace.area.getCoordinates()

	if '--debug' in sys.argv: Sx, Sy = 1 , 2
	if '--debug' in sys.argv: innerFace.globalDerivatives = np.array([ [1,2,3],[4,5,6] ])
	if '--debug' in sys.argv: shearModulus, poissonsRatio = 10,1/3
	# v=l/2(G+l)
	voigtAreaMatrix = np.array([[Sx, 0],[0, Sy],[Sy, Sx]])

	lameParameter=l=2*shearModulus*poissonsRatio/(1-2*poissonsRatio)
	b=l*(1-poissonsRatio)/poissonsRatio
	G=shearModulus

	Nx, Ny = innerFace.globalDerivatives
	zero=np.zeros(Nx.size)
	voigtGradientOperator = np.array([[Nx, zero],[zero, Ny],[Ny, Nx]])

	constitutiveMatrix = np.array([[lameParameter*(1-poissonsRatio)/poissonsRatio,lameParameter 							   ,0			],
								   [lameParameter								 ,lameParameter*(1-poissonsRatio)/poissonsRatio,0			],
								   [0			 								 ,0											   ,shearModulus]])
	
	if '--debug' in sys.argv:
		print('-'*100)
		print("𝜆 = ", round(lameParameter))
		print("𝛽 = ", b)
		print("G = ", shearModulus)
		print("Nx = ", Nx)
		print("Ny = ", Ny, end='\n\n')
		print("voigtAreaMatrix =\n", voigtAreaMatrix, end='\n\n')
		print("voigtGradientOperator =\n", show3dArray(voigtGradientOperator), end='\n\n')
		print("constitutiveMatrix =\n", constitutiveMatrix, end='\n\n')
		print('-'*100)

	coefficient = np.zeros([2,2,innerFace.element.vertices.size])
	coefficient[0][0] = np.matmul( np.array([ Sx*b, Sy*G ]), innerFace.globalDerivatives )
	coefficient[0][1] = np.matmul( np.array([ Sy*G, Sx*l ]), innerFace.globalDerivatives )
	coefficient[1][0] = np.matmul( np.array([ Sy*l, Sx*G ]), innerFace.globalDerivatives )
	coefficient[1][1] = np.matmul( np.array([ Sx*G, Sy*b ]), innerFace.globalDerivatives )

	c1=np.einsum("ki,kj,jdn->jd", voigtAreaMatrix, constitutiveMatrix, voigtGradientOperator)
	print(c1.shape)
	c1=np.zeros([2,2,innerFace.element.vertices.size])
	for i in range(voigtAreaMatrix.shape[1]): 					# dimension
		for j in range(constitutiveMatrix.shape[1]): 			# N_var
			for k in range(3):									# N_var
				for d in range(voigtGradientOperator.shape[1]): # dimension
					for n in range(voigtGradientOperator.shape[2]):
						c1[i][d] += voigtAreaMatrix[k][i] * constitutiveMatrix[k][j] * voigtGradientOperator[j][d][n]
	print(c1)
	if '--debug' in sys.argv:
		print("Not Voigt")
		print(show3dArray(c1))
		print("\nVoigt:")
		print(show3dArray(coefficient))
		print("\nAre them equal?:")
		try:
			print(c1==coefficient)
		except:
			print("Different sizes")

	return c1


def computeLocalMatrixCoefficient2(innerFace, G, poissonsRatio):
	Sx, Sy, Sz = innerFace.area.getCoordinates()
	l=2*shearModulus*poissonsRatio/(1-2*poissonsRatio)
	b=l*(1-poissonsRatio)/poissonsRatio

	coefficient = np.zeros([2,2,innerFace.element.vertices.size])
	coefficient[0][0] = np.matmul( np.array([ Sx*b, Sy*G ]), innerFace.globalDerivatives )
	coefficient[0][1] = np.matmul( np.array([ Sy*G, Sx*l ]), innerFace.globalDerivatives )
	coefficient[1][0] = np.matmul( np.array([ Sy*l, Sx*G ]), innerFace.globalDerivatives )
	coefficient[1][1] = np.matmul( np.array([ Sx*G, Sy*b ]), innerFace.globalDerivatives )

	#--------------------
	Sx, Sy, Sz = innerFace.area.getCoordinates()
	Nx, Ny = innerFace.globalDerivatives
	zero=np.zeros(Nx.size)
	lameParameter=2*shearModulus*poissonsRatio/(1-2*poissonsRatio)

	voigtAreaMatrix = np.array([[Sx, 0],[0, Sy],[Sy, Sx]])
	voigtGradientOperator = np.array([[Nx, zero],[zero, Ny],[Ny, Nx]])
	constitutiveMatrix = np.array([[lameParameter*(1-poissonsRatio)/poissonsRatio,lameParameter 							   ,0			],
								   [lameParameter								 ,lameParameter*(1-poissonsRatio)/poissonsRatio,0			],
								   [0			 								 ,0											   ,shearModulus]])

	# coefficient = np.einsum('ik,kjn->ijn',np.matmul(voigtAreaMatrix.T, constitutiveMatrix),voigtGradientOperator)
	
	
	c1=np.einsum("ki,kj,jdn->jd", voigtAreaMatrix, constitutiveMatrix, voigtGradientOperator)
	c1=np.zeros([2,2,innerFace.element.vertices.size])
	for i in range(voigtAreaMatrix.shape[1]): # dimension
		for j in range(constitutiveMatrix.shape[1]): # N_var
			for k in range(3):	# N_var
				for d in range(voigtGradientOperator.shape[1]): # dimension
					for n in range(voigtGradientOperator.shape[2]):
						c1[i][d] += voigtAreaMatrix[k][i] * constitutiveMatrix[k][j] * voigtGradientOperator[j][d][n]

	if '--debug' in sys.argv:
		print("Not Voigt")
		print(c1)
		print("\nVoigt:")
		print(coefficient)
		print("\nAre them equal?:")
		try:
			print(c1/coefficient)
		except:
			print("Different sizes")

	return coefficient

#-------------------------ADD TO LINEAR SYSTEM------------------------------
independent = np.zeros(2*numberOfVertices)

# Gravity Term
for region in grid.regions:
	volumetricForce = problemData.propertyData[region.handle]["Density"] * problemData.propertyData[region.handle]["Gravity"]
	for element in region.elements:
		local = 0
		for vertex in element.vertices:
			independent[vertex.handle + numberOfVertices] += -element.subelementVolumes[local] * volumetricForce
			local += 1

# Stress Term
for region in grid.regions:
	shearModulus = problemData.propertyData[region.handle]["ShearModulus"]
	poissonsRatio = problemData.propertyData[region.handle]["PoissonsRatio"]
	
	for element in region.elements:
		for innerFace in element.innerFaces:
			backwardVertexHandle = element.vertices[element.shape.innerFaceNeighborVertices[innerFace.local][0]].handle
			forwardVertexHandle = element.vertices[element.shape.innerFaceNeighborVertices[innerFace.local][1]].handle
			coefficient = computeLocalMatrixCoefficient(innerFace, shearModulus, poissonsRatio)

			i=0
			for vertex in element.vertices:
				matrix[backwardVertexHandle][vertex.handle] += coefficient[0][0][i]											# eq: u, var: u
				matrix[forwardVertexHandle][vertex.handle] -= coefficient[0][0][i]											# eq: u, var: u
				matrix[backwardVertexHandle][vertex.handle + numberOfVertices] += coefficient[0][1][i]						# eq: u, var: v
				matrix[forwardVertexHandle][vertex.handle + numberOfVertices] -= coefficient[0][1][i]						# eq: u, var: v
				matrix[backwardVertexHandle + numberOfVertices][vertex.handle] += coefficient[1][0][i]						# eq: v, var: u
				matrix[forwardVertexHandle + numberOfVertices][vertex.handle] -= coefficient[1][0][i]						# eq: v, var: u
				matrix[backwardVertexHandle + numberOfVertices][vertex.handle + numberOfVertices] += coefficient[1][1][i]	# eq: v, var: v
				matrix[forwardVertexHandle + numberOfVertices][vertex.handle + numberOfVertices] -= coefficient[1][1][i]	# eq: v, var: v
				i+=1
				if '--debug' in sys.argv: exit(0)



# Neumann Boundary Condition
for bCondition in problemData.neumannBoundaries["u"]:
	for facet in bCondition.boundary.facets:
		for outerFace in facet.outerFaces:
			independent[outerFace.vertex.handle] = bCondition.getValue(outerFace.handle) * np.linalg.norm(outerFace.area.getCoordinates())

for bCondition in problemData.neumannBoundaries["v"]:
	for facet in bCondition.boundary.facets:
		for outerFace in facet.outerFaces:
			independent[outerFace.vertex.handle+numberOfVertices] = bCondition.getValue(outerFace.handle) * np.linalg.norm(outerFace.area.getCoordinates())

# Dirichlet Boundary Condition
for bCondition in problemData.dirichletBoundaries["u"]:
	for vertex in bCondition.boundary.vertices:
		independent[vertex.handle] = bCondition.getValue(vertex.handle)
for bCondition in problemData.dirichletBoundaries["u"]:
	for vertex in bCondition.boundary.vertices:
		matrix[vertex.handle] = np.zeros(2*numberOfVertices)
		matrix[vertex.handle][vertex.handle] = 1.0

for bCondition in problemData.dirichletBoundaries["v"]:
	for vertex in bCondition.boundary.vertices:
		independent[vertex.handle+numberOfVertices] = bCondition.getValue(vertex.handle)
for bCondition in problemData.dirichletBoundaries["v"]:
	for vertex in bCondition.boundary.vertices:
		matrix[vertex.handle+numberOfVertices] = np.zeros(2*numberOfVertices)
		matrix[vertex.handle+numberOfVertices][vertex.handle+numberOfVertices] = 1.0

#-------------------------SOLVE LINEAR SYSTEM-------------------------------
displacements = np.linalg.solve(matrix, independent)

normal_stain_x = np.zeros(numberOfVertices)
normal_stain_y = np.zeros(numberOfVertices)
shear_strain = np.zeros(numberOfVertices)
for element in grid.elements:
	element_u = np.array([ displacements[vertex.handle] for vertex in element.vertices ])
	element_v = np.array([ displacements[vertex.handle+numberOfVertices] for vertex in element.vertices ])

	for local in range(element.vertices.size):
		localDerivatives = element.shape.vertexShapeFunctionDerivatives[local]
		globalDerivatives = np.matmul(np.linalg.inv(element.getTransposedJacobian(localDerivatives)) , np.transpose(localDerivatives))

		u_grad = np.matmul( globalDerivatives, element_u )
		v_grad = np.matmul( globalDerivatives, element_v )

		normal_stain_x[element.vertices[local].handle] += u_grad[0]
		normal_stain_y[element.vertices[local].handle] += v_grad[1]
		shear_strain[element.vertices[local].handle] += u_grad[1] + v_grad[0]

#-------------------------SAVE RESULTS--------------------------------------
# ACTUALLY DISPLACEMENTS WILL BE ONE MATRIX
cgnsSaver.save('u', displacements[:numberOfVertices], 0.0)
cgnsSaver.save('v', displacements[numberOfVertices:], 0.0)
cgnsSaver.save('exx', normal_stain_x, 0.0)
cgnsSaver.save('eyy', normal_stain_y, 0.0)
cgnsSaver.save('exy', shear_strain, 0.0)

#-------------------------CHECK CONVERGENCE---------------------------------
converged = False
difference = max([abs(temp-oldTemp) for temp, oldTemp in zip(displacements, prevDisplacements)])
prevDisplacements = displacements

#-------------------------------------------------------------------------------
#-------------------------AFTER END OF MAIN LOOP ITERATION------------------------
#-------------------------------------------------------------------------------
cgnsSaver.finalize()

print("\n\t\033[1;35mresult:\033[0m", problemData.paths["Output"]+"Results.cgns", '\n')
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