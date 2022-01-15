import numpy as np
from PyEFVLib.Shape import Triangle, Quadrilateral, Tetrahedron, Hexahedron, Prism, Pyramid
from PyEFVLib.Vertex import Vertex
from PyEFVLib.Element import Element
from PyEFVLib.Region import Region
from PyEFVLib.Facet import Facet
from PyEFVLib.Boundary import Boundary
from PyEFVLib.GridData import GridData

class Grid:
	"""
	Grid is the main geometric entity in PyEFVLib, and it contains all the other entities 
	necessary to solve a problem.

	Example:
		>>> import PyEFVLib
		>>> grid = PyEFVLib.read("mesh.msh")

	Attributes
	----------
	vertices	: list[PyEFVLib.Vertex]
	 	A list of the grid's vertices
	elements	: list[PyEFVLib.Element]
	 	A list of the grid's elements
	regions		: list[PyEFVLib.Region]
	 	A list of the grid's regions
	boundaries	: list[PyEFVLib.Boundary]
	 	A list of the grid's boundaries
	dimension	: int
	 	the dimension of the mesh (2 or 3)
	gridData	: PyEFVLib.GridData
	 	A data structure that holds the information provided by mesh file
	"""
	def __init__(self, gridData):
		if gridData.__class__ != GridData:
			raise Exception("Grid argument must be of class GridData")
		self.gridData = gridData
		self.dimension = gridData.dimension
		self.correctedForNegativeVolume = False
		self.__build()

	def __build(self):
		self.__buildVertices()
		self.__buildElements()
		self.__buildRegions()
		self.__buildBoundaries()

	def __buildVertices(self):
		handle = 0
		self.vertices = []
		for coord in self.gridData.verticesCoordinates:
			vertex = Vertex(coord, handle)
			self.vertices.append(vertex)
			handle += 1
		self.numberOfVertices = len(self.vertices)

	def __buildElements(self):
		self.correctForNegativeVolume = False
		self.innerFaceCounter = 0
		handle = 0
		self.elements = []
		for elementConnectivity in self.gridData.elementsConnectivities:
			element = Element(self, elementConnectivity, handle)
			self.elements.append(element)
			handle += 1
		if self.correctForNegativeVolume and not self.correctedForNegativeVolume:
			print("Negative volumes detected. Correcting for them.\nBe sure to fix the mesh in order to avoid inefficiencies during the mesh reading.")
			self.correctedForNegativeVolume = True
			self.__buildElements()
		elif self.correctedForNegativeVolume and self.correctForNegativeVolume:
			raise Exception("Negative volumes were found and couldn't be fixed.")
	shapes = [Triangle, Quadrilateral, Tetrahedron, Hexahedron, Prism, Pyramid]

	def __buildRegions(self):
		handle = 0
		self.regions = []
		for regionElementsIndexes, regionName in zip(self.gridData.regionsElementsIndexes, self.gridData.regionsNames):
			region = Region(self, regionElementsIndexes, regionName, handle)
			self.regions.append(region)
			handle += 1

	def __buildBoundaries(self):
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