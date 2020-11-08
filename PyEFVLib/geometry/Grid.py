import numpy as np
from PyEFVLib.geometry.Shape import Triangle, Quadrilateral, Tetrahedron, Hexahedron, Prism, Pyramid
from PyEFVLib.geometry.Vertex import Vertex
from PyEFVLib.geometry.Element import Element
from PyEFVLib.geometry.Region import Region
from PyEFVLib.geometry.Boundary import Boundary, BoundaryBuilder
from PyEFVLib.geometry.GridData import GridData

class Grid:
	def __init__(self, gridData):
		if gridData.__class__ != GridData:
			raise Exception("Grid argument must be of class GridData")
		self.gridData = gridData
		self.dimension = gridData.dimension
		self.build()

	def build(self):
		self.buildVertices()
		self.buildElements()
		self.buildRegions()
		self.buildBoundaries()

	def buildVertices(self):
		handle = 0
		self.vertices = np.array([])
		for coord in self.gridData.vertices:
			vertex = Vertex(coord, handle)
			self.vertices = np.append(self.vertices, vertex)
			handle += 1

	def buildElements(self):
		self.innerFaceCounter = 0
		handle = 0
		self.elements = np.array([])
		for iElem in self.gridData.elementsConnectivities:
			elem = Element([self.vertices[iVertex] for iVertex in iElem], self, handle)
			self.elements = np.append(self.elements, elem)
			handle += 1

	def buildRegions(self):
		handle = 0
		self.regions = np.array([])
		for iRegion, regionName in zip(self.gridData.regionsElementsIndexes, self.gridData.regionsNames):
			region = Region([self.elements[iElem] for iElem in iRegion], regionName, self, handle)
			self.regions = np.append(self.regions, region)
			handle += 1

	def buildBoundaries(self):
		self.boundaries = np.array([])
		BoundaryBuilder(self)

	def getShapes(self):
		return [Triangle, Quadrilateral, Tetrahedron, Hexahedron, Prism, Pyramid]

	def buildStencil(self):
		nVertices = len(self.vertices)
		self.stencil = [[] for i in range(nVertices)]
		for element in self.elements:
			localHandle = 0
			for v1 in element.vertices:
				for v2 in element.vertices[localHandle:]:
					if not v2.handle in self.stencil[v1.handle]:		self.stencil[v1.handle].append(v2.handle)
					if not v1.handle in self.stencil[v2.handle]:		self.stencil[v2.handle].append(v1.handle)
				localHandle += 1
