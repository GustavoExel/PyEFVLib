import numpy as np
import subprocess, os, sys
from PyEFVLib.geometry.Shape import Triangle, Quadrilateral, Tetrahedron, Hexahedron, Prism, Pyramid
from PyEFVLib.simulation.Saver import Saver

class VtuSaver(Saver):
	def __init__(self, grid, outputPath, basePath, fileName="Results", **kwargs): 
		Saver.__init__(self, grid, outputPath, basePath, 'vtu', fileName)

	def finalize(self):
		shapesDict	 = { Triangle:5,Quadrilateral:9,Tetrahedron:10,Hexahedron:12,Prism:13,Pyramid:14 }
		connectivity = [ [ vertex.handle for vertex in element.vertices ] for element in self.grid.elements ]
		shapes		 = [ shapesDict[element.shape] for element in self.grid.elements ]
		offsets		 = [ len(conn) for conn in connectivity ]
		offsets		 = [ sum(offsets[:i+1]) for i in range( len(offsets) ) ]

		if not os.path.isdir(os.path.dirname(self.outputPath)):
			os.makedirs(os.path.dirname(self.outputPath))

		if os.path.isfile(self.outputPath):
			os.remove(self.outputPath)

		# Write file
		with open(self.outputPath, "w") as file:
			file.write("<?xml version=\"1.0\"?>\n<VTKFile type=\"UnstructuredGrid\" version=\"0.1\" byte_order=\"LittleEndian\">\n\t<UnstructuredGrid>\n")
			file.write(f"\t\t<Piece NumberOfPoints=\"{ self.grid.numberOfVertices }\" NumberOfCells=\"{ self.grid.elements.size }\">\n\t\t\t<Points>\n")
			file.write("\t\t\t\t<DataArray type=\"Float64\" Name=\"Points\" NumberOfComponents=\"3\" format=\"ascii\">\n")
			file.write( "".join( [ f"\t\t\t\t\t{c:.8f}\n" for vertex in self.grid.vertices for c in vertex.getCoordinates() ] ) )
			file.write("\t\t\t\t</DataArray>\n\t\t\t</Points>\n\t\t\t<Cells>\n")
			file.write("\t\t\t\t<DataArray type=\"Int32\" Name=\"connectivity\" format=\"ascii\">\n")
			file.write( "\t\t\t\t\t" + "".join( [ f"{v:.0f} " for c in connectivity for v in c ] ) + "\n")
			file.write("\t\t\t\t</DataArray>\n")
			file.write("\t\t\t\t<DataArray type=\"Int32\" Name=\"offsets\" format=\"ascii\">\n")
			file.write( "\t\t\t\t\t" + "".join( [ f"{o:.0f} " for o in offsets ] ) + "\n")
			file.write("\t\t\t\t</DataArray>\n")
			file.write("\t\t\t\t<DataArray type=\"Int32\" Name=\"types\" format=\"ascii\">\n")
			file.write( "\t\t\t\t\t" + "".join( [ f"{s:.0f} " for s in shapes ] ) + "\n" )
			file.write("\t\t\t\t</DataArray>\n")
			file.write("\t\t\t</Cells>\n")
			file.write("\t\t\t<PointData>\n")
			for fieldName in self.fields.keys():
				file.write(f"\t\t\t\t<DataArray type=\"Float64\" Name=\"{fieldName}\" format=\"ascii\">\n")
				file.write("\t\t\t\t\t" + "".join([ f"{d:.15f} " for d in self.fields[fieldName][-1] ]) + "\n")
				file.write("\t\t\t\t</DataArray>\n")
			file.write("\t\t\t</PointData>\n\t\t</Piece>\n\t</UnstructuredGrid>\n</VTKFile>\n")

		self.finalized = True