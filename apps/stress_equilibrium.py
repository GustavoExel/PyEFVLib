import sys,os
sys.path.append(os.path.join(os.path.dirname(__file__), os.path.pardir))

from PyEFVLib import MSHReader, Grid, ProblemData, CgnsSaver
import numpy as np

#-------------------------SETTINGS----------------------------------------------
problemData = ProblemData('stress_equilibrium')

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

print("{:>9}\t{:>14}\t{:>14}\t{:>14}".format("Iteration", "CurrentTime", "TimeStep", "Difference"))

def computeLocalMatrixCoefficient(innerFace, G, poissonsRatio):
	Sx, Sy, Sz = innerFace.area.getCoordinates()
	l=2*shearModulus*poissonsRatio/(1-poissonsRatio)
	b=l*(-1+1/poissonsRatio)

	coefficient = np.zeros([2,2,innerFace.element.vertices.size])
	coefficient[0][0] = np.matmul( np.array([ Sx*b, Sy*G ]), innerFace.globalDerivatives )
	coefficient[0][1] = np.matmul( np.array([ Sy*G, Sx*l ]), innerFace.globalDerivatives )
	coefficient[1][0] = np.matmul( np.array([ Sy*l, Sx*G ]), innerFace.globalDerivatives )
	coefficient[1][1] = np.matmul( np.array([ Sx*G, Sy*b ]), innerFace.globalDerivatives )

	return coefficient

#-------------------------------------------------------------------------------
#-------------------------SIMULATION MAIN LOOP----------------------------------
#-------------------------------------------------------------------------------
while not converged and iteration < problemData.maxNumberOfIterations:
	#-------------------------ADD TO LINEAR SYSTEM------------------------------
	independent = np.zeros(2*numberOfVertices)

	# Gravity Term
	for region in grid.regions:
		# heatGeneration = problemData.propertyData[region.handle]["HeatGeneration"]
		for element in region.elements:
			local = 0
			for vertex in element.vertices:
				# independent[vertex.handle] += element.subelementVolumes[local] * heatGeneration
				local += 1

	# Stress Term
	if iteration == 0:
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


	# Neumann Boundary Condition
	for bCondition in problemData.neumannBoundaries["u"]:
		for facet in bCondition.boundary.facets:
			for outerFace in facet.outerFaces:
				independent[outerFace.vertex.handle] = bCondition.getValue(outerFace.handle) * np.linalg.norm(outerFace.area.getCoordinates())

	for bCondition in problemData.neumannBoundaries["v"]:
		for facet in bCondition.boundary.facets:
			for outerFace in facet.outerFaces:
				independent[outerFace.vertex.handle] = bCondition.getValue(outerFace.handle) * np.linalg.norm(outerFace.area.getCoordinates())

	# Dirichlet Boundary Condition
	for bCondition in problemData.dirichletBoundaries["u"]:
		for vertex in bCondition.boundary.vertices:
			independent[vertex.handle] = bCondition.getValue(vertex.handle)
	if iteration == 0:
		for bCondition in problemData.dirichletBoundaries["u"]:
			for vertex in bCondition.boundary.vertices:
				matrix[vertex.handle] = np.zeros(2*numberOfVertices)
				matrix[vertex.handle][vertex.handle] = 1.0

	for bCondition in problemData.dirichletBoundaries["v"]:
		for vertex in bCondition.boundary.vertices:
			independent[vertex.handle+numberOfVertices] = bCondition.getValue(vertex.handle)
	if iteration == 0:
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
			localDerivatives = element.shape.subelementShapeFunctionDerivatives[local]
			globalDerivatives = np.matmul(np.linalg.inv(element.getTransposedJacobian(localDerivatives)) , np.transpose(localDerivatives))

			u_grad = np.matmul( globalDerivatives, element_u )
			v_grad = np.matmul( globalDerivatives, element_v )

			normal_stain_x[element.vertices[local].handle] += u_grad[0]
			normal_stain_y[element.vertices[local].handle] += v_grad[1]
			shear_strain[element.vertices[local].handle] += u_grad[1] + v_grad[0]

	#-------------------------PRINT ITERATION DATA------------------------------
	if iteration > 0:
		print("{:>9}\t{:>14e}\t{:>14e}\t{:>14e}".format(iteration, currentTime, 0.0, difference))

	#-------------------------SAVE RESULTS--------------------------------------
	# ACTUALLY DISPLACEMENTS WILL BE ONE MATRIX
	cgnsSaver.save('u', displacements[:numberOfVertices], 0.0)
	cgnsSaver.save('v', displacements[numberOfVertices:], 0.0)
	cgnsSaver.save('exx', normal_stain_x, 0.0)
	cgnsSaver.save('eyy', normal_stain_y, 0.0)

	#-------------------------CHECK CONVERGENCE---------------------------------
	converged = False
	difference = max([abs(temp-oldTemp) for temp, oldTemp in zip(displacements, prevDisplacements)])
	prevDisplacements = displacements
	if iteration > 0:
		converged = difference < problemData.tolerance

	#-------------------------INCREMENT ITERATION-------------------------------
	iteration += 1   

#-------------------------------------------------------------------------------
#-------------------------AFTER END OF MAIN LOOP ITERATION------------------------
#-------------------------------------------------------------------------------
cgnsSaver.finalize()

print("\n\t\033[1;35mresult:\033[0m", problemData.paths["Output"]+"Results.cgns", '\n')
#-------------------------------------------------------------------------------
#-------------------------SHOW RESULTS GRAPHICALY-------------------------------
#-------------------------------------------------------------------------------
from matplotlib import pyplot as plt, colors, cm
from scipy.interpolate import griddata

X,Y = zip(*[v.getCoordinates()[:-1] for v in grid.vertices])
Xi, Yi = np.meshgrid( np.linspace(min(X), max(X), 2*len(X)), np.linspace(min(Y), max(Y), 2*len(Y)) )
def show(fieldValues, name):
	plt.figure()
	gridValues = griddata((X,Y), fieldValues, (Xi,Yi), method='linear')
	plt.pcolor(Xi,Yi,gridValues, cmap=colors.ListedColormap( cm.get_cmap("RdBu",256)(np.linspace(1,0,256)) ))
	plt.title(name)
	plt.colorbar()

show(displacements[:numberOfVertices], "X displacement (u)")
show(displacements[numberOfVertices:], "Y displacement (v)")
show(normal_stain_x, "X normal strain (exx)")
show(normal_stain_y, "Y normal strain (eyy)")

plt.show()