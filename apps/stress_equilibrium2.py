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
for boundary in boundaries:
    if boundary["u"].type == "Dirichlet" && boundary["v"].type == "Dirichlet":
		for vertex in boundary.vertices:
			matrix["u"][vertex] = [0] * (2*numberOfVertices)
			matrix["v"][vertex] = [0] * (2*numberOfVertices)
			independent["u"][vertex] = u_boundaryCondition
			independent["v"][vertex] = v_boundaryCondition

    Sx, Sy = outerFace.area
    elif boundary["u"].type == "Neumann" && boundary["v"].type == "Neumann":
    	for vertex in boundary.vertices:
    		independent["u"][vertex]=0.0
    		independent["v"][vertex]=0.0
		for outerFace in boundary.outerFaces:
			xNormalStress=Sx*(Tx*Sx+Ty*Sy)/(Sx**2+Sy**2)
			yNormalStress=Sy*(Tx*Sx+Ty*Sy)/(Sx**2+Sy**2)
			shearStress=sqrt((Tx-xNormalStress)**2+(Ty-yNormalStress)**2)

			independent["u"] += Sx*xNormalStress+Sy*shearStress
			independent["v"] += Sx*shearStress+Sy*yNormalStress

    elif boundary["u"].type == "Dirichlet" && boundary["v"].type == "Neumann":
    	for vertex in boundary.vertices:
			matrix["u"][vertex] = [0] * (2*numberOfVertices)
			independent["u"][vertex] = u_boundaryCondition
			independent["v"][vertex] = 0.0

		for outerFace in boundary.outerFaces:
			yNormalStress=Sy*(Tx*Sx+Ty*Sy)/(Sx**2+Sy**2)
			independent["v"] += Sy*yNormalStress

			for vertex in outerFace.element.vertices:
				Nx,Ny=outerFace.globalShapeFunctionDerivatives[vertex]
				matrix["u"][outerFace][vertex] -= Sy*shearModulus*Ny
				matrix["v"][outerFace][vertex] -= Sy*shearModulus*Nx

    elif boundary["u"].type == "Neumann" && boundary["v"].type == "Dirichlet":
    	for vertex in boundary.vertices:
			matrix["v"][vertex] = [0] * (2*numberOfVertices)
			independent["v"][vertex] = v_boundaryCondition
			independent["u"][vertex] = 0.0

		for outerFace in boundary.outerFaces:
			xNormalStress=Sx*(Tx*Sx+Ty*Sy)/(Sx**2+Sy**2)
			independent["u"] += Sx*xNormalStress
			for vertex in outerFace.element.vertices:
				Nx,Ny=outerFace.globalShapeFunctionDerivatives[vertex]
				matrix["u"][outerFace][vertex] -= Sy*shearModulus*Ny
				matrix["v"][outerFace][vertex] -= Sy*shearModulus*Nx


#-------------------------SOLVE LINEAR SYSTEM-------------------------------
# displacements = np.linalg.solve(matrix, independent)

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