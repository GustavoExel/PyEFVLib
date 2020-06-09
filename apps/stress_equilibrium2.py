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

#-------------------------ADD TO LINEAR SYSTEM------------------------------
independent = np.zeros(2*numberOfVertices)

# Gravity Term

# Stress Term

# Neumann Boundary Condition
# Dirichlet Boundary Condition
for bc in problemData.boundaryConditions:
	boundary=bc["u"].boundary

	if bc["u"].__type__ == "DIRICHLET" and bc["v"].__type__ == "DIRICHLET":
		for vertex in bc.boundary.vertices:
			# Assigns the displacement to both u and v
			matrix[vertex.handle] = 				  np.zeros(2*numberOfVertices)
			matrix[vertex.handle][vertex.handle] = 1.0
			matrix[vertex.handle+numberOfVertices] = np.zeros(2*numberOfVertices)
			matrix[vertex.handle+numberOfVertices][vertex.handle+numberOfVertices] = 1.0
			independent[vertex.handle] = 				  bc["u"].getValue(vertex.handle)
			independent[vertex.handle+numberOfVertices] = bc["v"].getValue(vertex.handle)

	elif bc["u"].__type__ == "NEUMANN" and bc["v"].__type__ == "NEUMANN":
		for vertex in boundary.vertices:
			# Clears the gravity term in both u and v in order to assign the tension
			independent[vertex.handle]=0.0
			independent[vertex.handle+numberOfVertices]=0.0
		for facet in boundary.facets:
			for outerFace in facet.outerFaces:
				# Tx isn't necessarily xNormalStress because the boundary may not be
				# parallel with the x or y axis, and that's why we're calculating
				# normal and shear stresses. In theory at least it shouldn't affect when 
				# the boundaries are parallel with the axis
				sx, sy, sz = outerFace.area.getCoordinates()
				Tx = bc["u"].getValue(outerFace.vertex.handle)
				Ty = bc["v"].getValue(outerFace.vertex.handle)
				xNormalStress=sx*(Tx*sx+Ty*sy)/(sx**2+sy**2)
				yNormalStress=sy*(Tx*sx+Ty*sy)/(sx**2+sy**2)
				shearStress=np.linalg.norm(np.array([Tx,Ty])-np.array([xNormalStress, yNormalStress]))

				independent[outerFace.vertex.handle] += sx*xNormalStress+sy*shearStress
				independent[outerFace.vertex.handle+numberOfVertices] += sx*shearStress+sy*yNormalStress

	elif bc["u"].__type__ == "DIRICHLET" and bc["v"].__type__ == "NEUMANN":
		for vertex in boundary.vertices:
			# Assigns the displacement in that vertex
			matrix[vertex.handle] = np.zeros(2*numberOfVertices)
			matrix[vertex.handle][vertex.handle] = 1.0
			independent[vertex.handle] = bc["u"].getValue(vertex.handle)

			# Clears the gravity term in order to keep only the b.c. summatory
			independent[vertex.handle+numberOfVertices] = 0.0

		for facet in boundary.facets:
			for outerFcae in facet.outerFaces:
				sx, sy, sz = outerFace.area.getCoordinates()
				Tx = bc["u"].getValue(outerFace.vertex.handle)
				Ty = bc["v"].getValue(outerFace.vertex.handle)
				yNormalStress=sy*(Tx*sx+Ty*sy)/(sx**2+sy**2)

				independent["v"] += Sy*yNormalStress

		for innerFace in outerFace.element.innerFaces:
			Nx,Ny=innerFace.globalShapeFunctionDerivatives.T[vertex]
			matrix["u"][outerFace][vertex] -= Sy*shearModulus*Ny
			matrix["v"][outerFace][vertex] -= Sy*shearModulus*Nx

	elif bc["u"].__type__ == "NEUMANN" and bc["v"].__type__ == "DIRICHLET":
		for vertex in boundary.vertices:
			matrix[vertex.handle+numberOfVertices] = np.zeros(2*numberOfVertices)
			matrix[vertex.handle+numberOfVertices][vertex.handle+numberOfVertices] = 1.0
			independent[vertex.handle+numberOfVertices] = bc["v"].getValue(vertex.handle)
			independent[vertex.handle] = 0.0

		for facet in boundary.facets:
			for outerFace in facet.outerFaces:
				Sx, Sy, Sz = outerFace.area.getCoordinates()
				Tx, Ty = bc["u"].getValue(outerFace.vertex.handle), bc["v"].getValue(outerFace.vertex.handle)
				xNormalStress=Sx*(Tx*Sx+Ty*Sy)/(Sx**2+Sy**2)

				independent[outerFace.vertex.handle] += Sx*xNormalStress

			for vertex in facet.element.vertices:
				Nx,Ny=outerFace.globalShapeFunctionDerivatives[vertex]
				matrix["u"][outerFace][vertex] -= Sy*shearModulus*Ny
				matrix["v"][outerFace][vertex] -= Sy*shearModulus*Nx


-------------------------SOLVE LINEAR SYSTEM-------------------------------
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
		# os.system("/usr/bin/paraview %sResults.cgns" % problemData.paths["Output"])
	except:
		print("Could not launch /usr/bin/paraview")