import numpy as np
from PyEFVLib.Shape import Triangle, Quadrilateral, Tetrahedron, Hexahedron, Prism, Pyramid
from PyEFVLib.Vertex import Vertex
from PyEFVLib.Element import Element
from PyEFVLib.Region import Region
from PyEFVLib.Facet import Facet
from PyEFVLib.Boundary import Boundary
from PyEFVLib.GridData import GridData

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
		self.vertices = []
		for coord in self.gridData.verticesCoordinates:
			vertex = Vertex(coord, handle)
			self.vertices.append(vertex)
			handle += 1
		self.numberOfVertices = len(self.vertices)

	def buildElements(self):
		self.innerFaceCounter = 0
		handle = 0
		self.elements = []
		for elementConnectivity in self.gridData.elementsConnectivities:
			element = Element(self, elementConnectivity, handle)
			self.elements.append(element)
			handle += 1
	shapes = [Triangle, Quadrilateral, Tetrahedron, Hexahedron, Prism, Pyramid]

	def buildRegions(self):
		handle = 0
		self.regions = []
		for regionElementsIndexes, regionName in zip(self.gridData.regionsElementsIndexes, self.gridData.regionsNames):
			region = Region(self, regionElementsIndexes, regionName, handle)
			self.regions.append(region)
			handle += 1

	def buildBoundaries(self):
		self.outerFaceCounter = 0
		handle = 0
		self.facets = []
		for facetConnectivity in self.gridData.boundariesConnectivities:
			facet = Facet(self, facetConnectivity, handle)
			self.facets.append(facet)
			handle += 1

		handle = 0
		self.boundaries = []
		for boundaryFacetsIdexes, boundaryName in zip(self.gridData.boundariesIndexes, self.gridData.boundariesNames):
			boundary = Boundary(self, boundaryFacetsIdexes, boundaryName, handle)
			self.boundaries.append(boundary)
			handle += 1

	# def buildStencil(self):
	# 	nVertices = len(self.vertices)
	# 	self.stencil = [[] for i in range(nVertices)]
	# 	for element in self.elements:
	# 		localHandle = 0
	# 		for v1 in element.vertices:
	# 			for v2 in element.vertices[localHandle:]:
	# 				if not v2.handle in self.stencil[v1.handle]:		self.stencil[v1.handle].append(v2.handle)
	# 				if not v1.handle in self.stencil[v2.handle]:		self.stencil[v2.handle].append(v1.handle)
	# 			localHandle += 1
