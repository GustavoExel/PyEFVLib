import numpy as np
from libs.geometry.Shape import Triangle, Quadrilateral
from libs.geometry.Vertex import Vertex
from libs.geometry.Element import Element
from libs.geometry.Region import Region
from libs.geometry.Boundary import Boundary, BoundaryBuilder

class Grid:
	def __init__(self, gridData):
		self.gridData = gridData
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
		for iElem in self.gridData.elemConnectivity:
			elem = Element([self.vertices[iVertex] for iVertex in iElem], self, handle)
			self.elements = np.append(self.elements, elem)
			handle += 1

	def buildRegions(self):
		handle = 0
		self.regions = np.array([])
		for iRegion, regionName in zip(self.gridData.regionElements, self.gridData.regionNames):
			region = Region([self.elements[iElem] for iElem in iRegion], regionName, handle)
			region.setGrid(self)
			self.regions = np.append(self.regions, region)
			handle += 1

	def buildBoundaries(self):
		self.boundaries = np.array([])
		BoundaryBuilder(self)

	def getVertices(self):
		return self.vertices

	def getElements(self):
		return self.elements

	def getRegions(self):
		return self.regions

	def getShapes(self):
		return [Triangle, Quadrilateral]