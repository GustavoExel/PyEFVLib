import sys,os
pyEFVLibPath = os.path.join(os.path.dirname(__file__), os.path.pardir, os.path.pardir)
workspacePath = os.path.join(os.path.dirname(__file__), os.path.pardir)
sys.path += [pyEFVLibPath, workspacePath]

import PyEFVLib
from colorplot import colorplot
from matplotlib import pyplot as plt
import numpy as np

def reconstruct2D(grid,Fx,Fy,Fxx,Fyy):
	aDivs = [Fxx(vertex.x, vertex.y) + Fyy(vertex.x, vertex.y) for vertex in grid.vertices]
	nDivs = grid.numberOfVertices * [0,]

	for element in grid.elements:
		elementXFieldVector = np.array( [Fx(vertex.x, vertex.y) for vertex in element.vertices] )
		elementYFieldVector = np.array( [Fy(vertex.x, vertex.y) for vertex in element.vertices] )

		for innerFace in element.innerFaces:
			backwardVertex = element.vertices[ element.shape.innerFaceNeighborVertices[innerFace.local][0] ]
			forwardVertex = element.vertices[ element.shape.innerFaceNeighborVertices[innerFace.local][1] ]

			shapeFunctionValues = innerFace.getShapeFunctions()

			xValAtIP = np.dot(elementXFieldVector, shapeFunctionValues)
			yValAtIP = np.dot(elementYFieldVector, shapeFunctionValues)


			area = innerFace.area.getCoordinates()[:-1]

			nDivs[backwardVertex.handle] += np.dot( (xValAtIP, yValAtIP), area )
			nDivs[forwardVertex.handle] -= np.dot( (xValAtIP, yValAtIP), area )

		for outerFace in element.outerFaces:
			shapeFunctionValues = element.shape.outerFaceShapeFunctionValues[outerFace.facet.elementLocalIndex][outerFace.local]

			xValAtIP = np.dot(elementXFieldVector, shapeFunctionValues)
			yValAtIP = np.dot(elementYFieldVector, shapeFunctionValues)

			area = outerFace.area.getCoordinates()[:-1]

			nDivs[outerFace.vertex.handle] += np.dot( (xValAtIP, yValAtIP), area )


	nDivs = [ nDiv/vertex.volume for vertex, nDiv in zip(grid.vertices, nDivs) ]
	
	# boundaryVertices = [ vertex.handle for boundary in grid.boundaries for vertex in boundary.vertices ]
	# aDivs, nDivs = zip(*[ (div, r1div) for div,r1div,vertex in zip(aDivs,nDivs,grid.vertices) if vertex.handle not in boundaryVertices ])

	aDivs, nDivs = np.array(aDivs), np.array(nDivs)

	print("\nReconstructed values at integration points")
	print(f"Max difference = {max(abs(aDivs-nDivs)) :.4f}, Field range = [{min(aDivs) :.4f}, {max(aDivs) :.4f}]\t| {100*(max(abs(aDivs-nDivs)))/(max(abs(aDivs))):.2f}%")



def reconstruct3D(grid,Fx,Fy,Fz,Fxx,Fyy,Fzz):
	aDivs = [Fxx(vertex.x, vertex.y, vertex.z) + Fyy(vertex.x, vertex.y, vertex.z) + Fzz(vertex.x, vertex.y, vertex.z) for vertex in grid.vertices]
	nDivs = grid.numberOfVertices * [0,]

	for element in grid.elements:
		elementXFieldVector = np.array( [Fx(vertex.x, vertex.y, vertex.z) for vertex in element.vertices] )
		elementYFieldVector = np.array( [Fy(vertex.x, vertex.y, vertex.z) for vertex in element.vertices] )
		elementZFieldVector = np.array( [Fz(vertex.x, vertex.y, vertex.z) for vertex in element.vertices] )

		for innerFace in element.innerFaces:
			backwardVertex = element.vertices[ element.shape.innerFaceNeighborVertices[innerFace.local][0] ]
			forwardVertex = element.vertices[ element.shape.innerFaceNeighborVertices[innerFace.local][1] ]

			shapeFunctionValues = innerFace.getShapeFunctions()

			xValAtIP = np.dot(elementXFieldVector, shapeFunctionValues)
			yValAtIP = np.dot(elementYFieldVector, shapeFunctionValues)
			zValAtIP = np.dot(elementZFieldVector, shapeFunctionValues)

			area = innerFace.area.getCoordinates()

			nDivs[backwardVertex.handle] += np.dot( (xValAtIP, yValAtIP, zValAtIP), area )
			nDivs[forwardVertex.handle] -= np.dot( (xValAtIP, yValAtIP, zValAtIP), area )

		for outerFace in element.outerFaces:
			shapeFunctionValues = element.shape.outerFaceShapeFunctionValues[outerFace.facet.elementLocalIndex][outerFace.local]

			xValAtIP = np.dot(elementXFieldVector, shapeFunctionValues)
			yValAtIP = np.dot(elementYFieldVector, shapeFunctionValues)
			zValAtIP = np.dot(elementZFieldVector, shapeFunctionValues)

			area = outerFace.area.getCoordinates()

			nDivs[outerFace.vertex.handle] += np.dot( (xValAtIP, yValAtIP, zValAtIP), area )


	nDivs = [ nDiv/vertex.volume for vertex, nDiv in zip(grid.vertices, nDivs) ]
	
	boundaryVertices = [ vertex.handle for boundary in grid.boundaries for vertex in boundary.vertices ]
	aDivs, nDivs = zip(*[ (div, r1div) for div,r1div,vertex in zip(aDivs,nDivs,grid.vertices) if vertex.handle not in boundaryVertices ])

	aDivs, nDivs = np.array(aDivs), np.array(nDivs)

	print("\nReconstructed values at integration points")
	print(f"Max difference = {max(abs(aDivs-nDivs)) :.4f}, Field range = [{min(aDivs) :.4f}, {max(aDivs) :.4f}]\t| {100*(max(abs(aDivs-nDivs)))/(max(abs(aDivs))):.2f}%")



if __name__ == "__main__":
	for meshName in ["Fine.msh", "10x10.msh"]:
		grid = PyEFVLib.read( os.path.join(pyEFVLibPath, "meshes", "msh", "2D", meshName) )
		print("------------------------------------\n", meshName)

		Fx = lambda x,y: x+y
		Fy = lambda x,y: x*y
		Fxx = lambda x,y: 1
		Fyy = lambda x,y: x
		print("\n__________________")
		print("(x+y, x*y)")
		reconstruct2D(grid,Fx,Fy,Fxx,Fyy)

		Fx = lambda x,y: x**2 + y**2
		Fy = lambda x,y: np.log(y+1) + 1/(x+1)
		Fxx = lambda x,y: 2*x
		Fyy = lambda x,y: 1/(y+1)
		print("\n__________________")
		print("(x^2 + y^2, ln(y) + 1/x)")
		reconstruct2D(grid,Fx,Fy,Fxx,Fyy)

		Fx = lambda x,y: -np.exp(-x) * np.sin(y) + x
		Fy = lambda x,y: np.exp(-x) * np.cos(y)
		Fxx = lambda x,y: np.exp(-x) * np.sin(y) + 1
		Fyy = lambda x,y: -np.exp(-x) * np.sin(y)
		print("\n__________________")
		print("(-exp(-x) * sin(y) + x, exp(-x) * cos(y))")
		reconstruct2D(grid,Fx,Fy,Fxx,Fyy)

	for meshName in ["Hexas.msh", "Prisms.msh", "Tetras.msh", "Pyrams.msh"]:
		grid = PyEFVLib.read( os.path.join(pyEFVLibPath, "meshes", "msh", "3D", meshName) )
		print("------------------------------------\n", meshName)
		
		Fx = lambda x,y,z: 2*x
		Fy = lambda x,y,z: 2*y
		Fz = lambda x,y,z: 2*z
		Fxx = lambda x,y,z: 2
		Fyy = lambda x,y,z: 2
		Fzz = lambda x,y,z: 2
		print("\n__________________")
		print("(2x, 2y, 2z)")
		reconstruct3D(grid,Fx,Fy,Fz,Fxx,Fyy,Fzz)

		Fx = lambda x,y,z: y/(z+1) + x
		Fy = lambda x,y,z: x/(z+1) + y
		Fz = lambda x,y,z: -0.1*x*y/(z+1)**2 + 0.1*z
		Fxx = lambda x,y,z: 1
		Fyy = lambda x,y,z: 1
		Fzz = lambda x,y,z: 0.1*(2*x*y)/(z+1)**3 + 0.1
		print("\n__________________")
		print("(y/(z+1), x/(z+1), -xy/(z+1)^2)")
		reconstruct3D(grid,Fx,Fy,Fz,Fxx,Fyy,Fzz)