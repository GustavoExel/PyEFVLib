import sys,os
pyEFVLibPath = os.path.join(os.path.dirname(__file__), os.path.pardir, os.path.pardir)
workspacePath = os.path.join(os.path.dirname(__file__), os.path.pardir)
sys.path += [pyEFVLibPath, workspacePath]

import PyEFVLib
from colorplot import colorplot
from matplotlib import pyplot as plt
import numpy as np

def reconstruct2D(grid,F,Fx,Fy,cplot=False,diff=False,quiver=False):
	X,Y = zip(*[v.getCoordinates()[:-1] for v in grid.vertices])
	Xip,Yip = zip( *[innerFace.centroid.getCoordinates()[:-1] for element in grid.elements for innerFace in element.innerFaces] )

	X,Y=np.array(X),np.array(Y)
	Xip,Yip=np.array(Xip),np.array(Yip)

	fieldAtVertices = F(X,Y)
	fieldAtIP = F(Xip,Yip)
	
	gradFieldAtVertices = np.array(list(zip( Fx(X,Y), Fy(X,Y) )))
	rGradFieldAtVertices = grid.numberOfVertices * [(0,0),]

	for element in grid.elements:
		elementFieldVector = np.array( [fieldAtVertices[vertex.handle] for vertex in element.vertices] )
		for innerFaceIndex, innerFace in enumerate(element.innerFaces):
			backwardVertex = element.vertices[ element.shape.innerFaceNeighborVertices[innerFaceIndex][0] ]
			forwardVertex = element.vertices[ element.shape.innerFaceNeighborVertices[innerFaceIndex][1] ]

			shapeFunctionValues = element.shape.innerFaceShapeFunctionValues[innerFaceIndex]
			pressureAtIP = np.dot(elementFieldVector, shapeFunctionValues)

			area = innerFace.area.getCoordinates()[:-1]

			rGradFieldAtVertices[backwardVertex.handle] += pressureAtIP * area	
			rGradFieldAtVertices[forwardVertex.handle] -= pressureAtIP * area	


	for facet in grid.facets:
		elementFieldVector = np.array( [fieldAtVertices[vertex.handle] for vertex in facet.element.vertices] )
		for outerFace in facet.outerFaces:
			shapeFunctionValues = facet.element.shape.outerFaceShapeFunctionValues[facet.elementLocalIndex][outerFace.local]

			pressureAtIP = np.dot(elementFieldVector, shapeFunctionValues)
			area = outerFace.area.getCoordinates()[:-1]

			rGradFieldAtVertices[outerFace.vertex.handle] += pressureAtIP * area	



	rGradFieldAtVertices = [ grad/vertex.volume for vertex, grad in zip( grid.vertices, rGradFieldAtVertices ) ]

	# boundaryVertices = [ vertex.handle for boundary in grid.boundaries for vertex in boundary.vertices ]
	# gradFieldAtVertices, rGradFieldAtVertices = zip(*[ (grad,rgrad) for grad,rgrad,vertex in zip(gradFieldAtVertices,rGradFieldAtVertices,grid.vertices) if vertex.handle not in boundaryVertices ])

	# print("\nSum of values in the control volume contour")
	gNFaV = np.array([np.linalg.norm(grad) for grad in gradFieldAtVertices])
	r2gNFaV = np.array([np.linalg.norm(grad) for grad in rGradFieldAtVertices])
	print(f"Max difference = {max(abs(gNFaV-r2gNFaV)) :>8.4f}, Field range = [{min(gNFaV) :.4f}, {max(gNFaV) :.4f}]\t| {100*max(abs(gNFaV-r2gNFaV))/max(gNFaV):.6f}%")

def reconstruct3D(grid,F,Fx,Fy,Fz):
	X,Y,Z = zip(*[v.getCoordinates() for v in grid.vertices])
	Xip,Yip,Zip = zip( *[innerFace.centroid.getCoordinates() for element in grid.elements for innerFace in element.innerFaces] )

	X,Y,Z=np.array(X),np.array(Y),np.array(Z)
	Xip,Yip,Zip=np.array(Xip),np.array(Yip),np.array(Zip)

	fieldAtVertices = F(X,Y,Z)
	gradFieldAtVertices = np.array(list(zip( Fx(X,Y,Z), Fy(X,Y,Z), Fz(X,Y,Z) )))
	fieldAtIP = F(Xip,Yip,Zip)
	
	rGradFieldAtVertices = grid.numberOfVertices * [(0,0,0),]


	for element in grid.elements:
		elementFieldVector = np.array( [fieldAtVertices[vertex.handle] for vertex in element.vertices] )
		for innerFaceIndex, innerFace in enumerate(element.innerFaces):
			backwardVertex = element.vertices[ element.shape.innerFaceNeighborVertices[innerFaceIndex][0] ]
			forwardVertex = element.vertices[ element.shape.innerFaceNeighborVertices[innerFaceIndex][1] ]

			shapeFunctionValues = element.shape.innerFaceShapeFunctionValues[innerFaceIndex]
			pressureAtIP = np.dot(elementFieldVector, shapeFunctionValues)

			area = innerFace.area.getCoordinates()

			rGradFieldAtVertices[backwardVertex.handle] += pressureAtIP * area	
			rGradFieldAtVertices[forwardVertex.handle] -= pressureAtIP * area	

	for facet in grid.facets:
		elementFieldVector = np.array( [fieldAtVertices[vertex.handle] for vertex in facet.element.vertices] )
		for outerFace in facet.outerFaces:
			shapeFunctionValues = facet.element.shape.outerFaceShapeFunctionValues[facet.elementLocalIndex][outerFace.local]

			pressureAtIP = np.dot(elementFieldVector, shapeFunctionValues)
			area = outerFace.area.getCoordinates()

			rGradFieldAtVertices[outerFace.vertex.handle] += pressureAtIP * area	


	rGradFieldAtVertices = [ grad/vertex.volume for vertex, grad in zip( grid.vertices, rGradFieldAtVertices ) ]

	# boundaryVertices = [ vertex.handle for boundary in grid.boundaries for vertex in boundary.vertices ]
	# gradFieldAtVertices, rGradFieldAtVertices = zip(*[ (grad,rgrad) for grad,rgrad,vertex in zip(gradFieldAtVertices,rGradFieldAtVertices,grid.vertices) if vertex.handle not in boundaryVertices ])

	gNFaV = np.array([np.linalg.norm(grad) for grad in gradFieldAtVertices])
	rgNFaV = np.array([np.linalg.norm(grad) for grad in rGradFieldAtVertices])
	print(f"Max difference = {max(abs(gNFaV-rgNFaV)) :>8.4f}, Field range = [{min(gNFaV) :.4f}, {max(gNFaV) :.4f}]\t| {100*max(abs(gNFaV-rgNFaV))/max(gNFaV):.2f}%")

if __name__ == "__main__":
	for meshName in ["Fine.msh", "10x10.msh"]:
		grid = PyEFVLib.read( os.path.join(pyEFVLibPath, "meshes", "msh", "2D", meshName) )
		print("------------------------------------\n", meshName)

		F = lambda x,y: 3*np.power(x,2) + np.sin(2*y)
		Fx = lambda x,y: 6*x
		Fy = lambda x,y: 2*np.cos(2*y)
		print("\n3(x^2) + sin(2y)")
		reconstruct2D(grid,F,Fx,Fy,cplot=False,diff=False,quiver=False)

		F = lambda x,y: 3*np.power(x,2) + 4*np.power(y,2)
		Fx = lambda x,y: 6*x
		Fy = lambda x,y: 8*y
		print("\n3(x^2) + 4(y^2)")
		reconstruct2D(grid,F,Fx,Fy)

		F = lambda x,y: np.exp(-x) * np.sin(y)
		Fx = lambda x,y: -np.exp(-x) * np.sin(y)
		Fy = lambda x,y: np.exp(-x) * np.cos(y)
		print("\nexp(-x) * sin(y)")
		reconstruct2D(grid,F,Fx,Fy)

	for meshName in ["Tetras.msh", "Pyrams.msh", "Hexas.msh", "Prisms.msh"]:
		grid = PyEFVLib.read( os.path.join(pyEFVLibPath, "meshes", "msh", "3D", meshName) )
		print("------------------------------------\n", meshName)
		
		F = lambda x,y,z: np.power(x,2) + np.power(y,2) + np.power(z,2)
		Fx = lambda x,y,z: 2*x
		Fy = lambda x,y,z: 2*y
		Fz = lambda x,y,z: 2*z
		print("\n(x^2) + (y^2) + (z^2)")
		reconstruct3D(grid,F,Fx,Fy,Fz)

		F = lambda x,y,z: (x*y)/(z+1)
		Fx = lambda x,y,z: y/(z+1)
		Fy = lambda x,y,z: x/(z+1)
		Fz = lambda x,y,z: -x*y/(z+1)**2
		print("\n(xy)/z")
		reconstruct3D(grid,F,Fx,Fy,Fz)


# def f(x,y,show=False):
# 	grid = PyEFVLib.read( "../../meshes/msh/2D/Square.msh" )

# 	VERTEX = min(grid.vertices, key=lambda v:abs(v.x-x)+abs(v.y-y))
# 	elements = VERTEX.elements
# 	IP = [ innerFace for element in elements for iFidx, innerFace in enumerate(element.innerFaces) if list(element.vertices).index(VERTEX) in element.shape.innerFaceNeighborVertices[iFidx] ]
# 	IPDir = [ 1 if ip.element.shape.innerFaceNeighborVertices[list(ip.element.innerFaces).index(ip)][0]==list(ip.element.vertices).index(VERTEX) else -1 for ip in IP ]
# 	IPCoords = [ ip.centroid.getCoordinates()[:-1] for ip in IP ]
# 	IPAreas = [ ipDir * ip.area.getCoordinates()[:-1] for ip,ipDir in zip(IP, IPDir) ]

# 	F = lambda x,y: 3*np.power(x,2) + np.sin(2*y)
# 	Fx = lambda x,y: 6*x
# 	Fy = lambda x,y: 2*np.cos(2*y)
# 	FIP = [ F(x,y) for x,y in IPCoords ]

# 	grad = [ Fx( VERTEX.x, VERTEX.y ), Fy( VERTEX.x, VERTEX.y ) ]
# 	rgrad = sum([ fip*area for fip,area in zip(FIP, IPAreas) ])/VERTEX.volume

# 	print(f"error = {100*np.linalg.norm(np.array(grad)-np.array(rgrad))/max(np.linalg.norm(grad), np.linalg.norm(rgrad)):.02f}%")

# 	if show:
# 		for element in VERTEX.elements:
# 			Xe,Ye = zip(*[vertex.getCoordinates()[:-1] for vertex in element.vertices])
# 			plt.plot(Xe+(Xe[0],), Ye+(Ye[0],), color='k')
# 			plt.scatter(Xe, Ye, color='r')

# 		Xip, Yip = zip(*IPCoords)
# 		plt.scatter(Xip, Yip, color='g')
# 		sx, sy = zip(*IPAreas)
# 		plt.quiver(Xip, Yip, sx, sy)

# 		plt.scatter(VERTEX.x, VERTEX.y, color='b')
# 		plt.show()
# if __name__ == "y__main__":
# 	f(0.7,0.7,False)