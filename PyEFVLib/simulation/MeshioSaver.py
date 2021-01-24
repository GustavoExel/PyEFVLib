import numpy as np
import subprocess, os, sys
import meshio
from PyEFVLib.simulation.Saver import Saver
from PyEFVLib.geometry.Shape import Triangle, Quadrilateral, Tetrahedron, Hexahedron, Prism, Pyramid

class MeshioSaver(Saver):
	# line triangle quad tetra pyramid wedge hexahedron
	# Formats:
	# msh mdpa ply stl vtk vtu xdmf xmf h5m med inp mesh meshb bdf fem nas obj off post post.gz dato dato.gz su2 svg dat tec ugrid wkt 
	def __init__(self, grid, outputPath, basePath, extension, fileName="Results", **kwargs): 
		Saver.__init__(self, grid, outputPath, basePath, extension, fileName)
		if not os.path.exists(outputPath):
			os.makedirs( outputPath )

	def finalize(self):
		self.points = np.array( [v.getCoordinates() for v in self.grid.vertices] )
		
		meshioShapes   = ["triangle", "quad", "tetra", "pyramid", "wedge", "hexahedron"]
		pyEFVLibShapes = [Triangle, Quadrilateral, Tetrahedron, Pyramid, Prism, Hexahedron]

		self.cells = [ ( shape , np.array([[vertex.handle for vertex in element.vertices] for element in self.grid.elements if element.shape == shapeClass], dtype=np.uint64) ) for shape, shapeClass in zip(meshioShapes, pyEFVLibShapes) ]
		self.cells = [ cell for cell in self.cells if cell[1].size ]

		# Two separate writers because the only that supports time series data is xdmf
		if self.extension == "xdmf":
			self.xdmfWrite()
		else:
			self.meshioWrite()
	
		self.finalized = True
	
	def meshioWrite(self):
		data  = { fieldName : self.fields[fieldName][-1] for fieldName in self.fields }
		
		meshioMesh = meshio.Mesh( self.points, self.cells, point_data=data )
		meshioMesh.write( self.outputPath )

	def xdmfWrite(self):
		with meshio.xdmf.TimeSeriesWriter(self.outputPath) as writer:
			writer.write_points_cells(self.points, self.cells)
			for idx, timeStep, in enumerate(self.timeSteps):
				data  = { fieldName : self.fields[fieldName][idx] for fieldName in self.fields }
				writer.write_data(timeStep, point_data=data)