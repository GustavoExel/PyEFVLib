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
	gradFieldAtIP = np.array(list(zip( Fx(Xip,Yip), Fy(Xip,Yip) )))

	rFieldAtIP = []
	rGradFieldAtIP = []

	for element in grid.elements:
		fieldVector = np.array( [fieldAtVertices[vertex.handle] for vertex in element.vertices] )
		for innerFace in element.innerFaces:
			grad = np.matmul( innerFace.globalDerivatives, fieldVector )
			rGradFieldAtIP.append(grad)

			val = np.dot( fieldVector, innerFace.getShapeFunctions() )
			rFieldAtIP.append(val)


	FxIP, FyIP = zip(*gradFieldAtIP); FxIP, FyIP = np.array(FxIP), np.array(FyIP)
	rFxIP, rFyIP = zip(*rGradFieldAtIP); rFxIP, rFyIP = np.array(rFxIP), np.array(rFyIP)
	rFieldAtIP = np.array(rFieldAtIP)

	# colorplot(Xip,Yip,FyIP-rFyIP)

	print(f"Max difference between Fx field = {max(abs(FxIP-rFxIP)) :.4f}, Field range = [{min(FxIP) :.4f}, {max(FxIP) :.4f}]\t| {100*(max(abs(FxIP-rFxIP)))/(max(abs(FxIP))):.2f}%")
	print(f"Max difference between Fy field = {max(abs(FyIP-rFyIP)) :.4f}, Field range = [{min(FyIP) :.4f}, {max(FyIP) :.4f}]\t| {100*(max(abs(FyIP-rFyIP)))/(max(abs(FyIP))):.2f}%")

	print(f"Max difference between F field = {max(abs(fieldAtIP-rFieldAtIP)) :.4f}, Field range = [{min(fieldAtIP) :.4f}, {max(fieldAtIP) :.4f}]\t\t| {100*(max(abs(fieldAtIP-rFieldAtIP)))/(max(abs(fieldAtIP))):.2f}%")

	if cplot:
		colorplot(Xip,Yip,rFieldAtIP)

	if diff:
		colorplot(Xip,Yip,fieldAtIP-rFieldAtIP)

	if quiver:
		plt.quiver(Xip,Yip,rFxIP,rFyIP)
		plt.scatter(X,Y,marker='.',linewidths=0.5,color='k')
		plt.show()

def reconstruct3D(grid,F,Fx,Fy,Fz):
	X,Y,Z = zip(*[v.getCoordinates() for v in grid.vertices])
	Xip,Yip,Zip = zip( *[innerFace.centroid.getCoordinates() for element in grid.elements for innerFace in element.innerFaces] )

	X,Y,Z=np.array(X),np.array(Y),np.array(Z)
	Xip,Yip,Zip=np.array(Xip),np.array(Yip),np.array(Zip)

	fieldAtVertices = F(X,Y,Z)

	fieldAtIP = F(Xip,Yip,Zip)
	gradFieldAtIP = np.array(list(zip( Fx(Xip,Yip,Zip), Fy(Xip,Yip,Zip), Fz(Xip,Yip,Zip) )))

	rFieldAtIP = []
	rGradFieldAtIP = []

	for element in grid.elements:
		fieldVector = np.array( [fieldAtVertices[vertex.handle] for vertex in element.vertices] )
		for innerFace in element.innerFaces:
			grad = np.matmul( innerFace.globalDerivatives, fieldVector )
			rGradFieldAtIP.append(grad)

			val = np.dot( fieldVector, innerFace.getShapeFunctions() )
			rFieldAtIP.append(val)


	FxIP, FyIP, FzIP = zip(*gradFieldAtIP); FxIP, FyIP, FzIP = np.array(FxIP), np.array(FyIP), np.array(FzIP)
	rFxIP, rFyIP, rFzIP = zip(*rGradFieldAtIP); rFxIP, rFyIP, rFzIP = np.array(rFxIP), np.array(rFyIP), np.array(rFzIP)
	rFieldAtIP = np.array(rFieldAtIP)

	print(f"Max difference between Fx field = {max(abs(FxIP-rFxIP)) :.4f}, Field range = [{min(FxIP) :.4f}, {max(FxIP) :.4f}]\t| {100*(max(abs(FxIP-rFxIP)))/(max(abs(FxIP))):.2f}%")
	print(f"Max difference between Fy field = {max(abs(FyIP-rFyIP)) :.4f}, Field range = [{min(FyIP) :.4f}, {max(FyIP) :.4f}]\t| {100*(max(abs(FyIP-rFyIP)))/(max(abs(FyIP))):.2f}%")
	print(f"Max difference between Fz field = {max(abs(FzIP-rFzIP)) :.4f}, Field range = [{min(FzIP) :.4f}, {max(FzIP) :.4f}]\t| {100*(max(abs(FzIP-rFzIP)))/(max(abs(FzIP))):.2f}%")

	print(f"Max difference between F field = {max(abs(fieldAtIP-rFieldAtIP)) :.4f}, Field range = [{min(fieldAtIP) :.4f}, {max(fieldAtIP) :.4f}]\t\t| {100*(max(abs(fieldAtIP-rFieldAtIP)))/(max(abs(fieldAtIP))):.2f}%")

if __name__ == "__main__":
	for meshName in ["Fine.msh", "10x10.msh"]:
		grid = PyEFVLib.read( os.path.join(pyEFVLibPath, "meshes", meshName) )
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

	for meshName in ["Hexas.msh", "Pyrams.msh"]:
		grid = PyEFVLib.read( os.path.join(pyEFVLibPath, "meshes", "3DGeometries", meshName) )
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
