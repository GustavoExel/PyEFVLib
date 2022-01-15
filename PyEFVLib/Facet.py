from PyEFVLib.Point import Point
from PyEFVLib.OuterFace import OuterFace
import numpy as np

class Facet:
	"""
	Facet mediates the distribution of attributes between outerFace, boundary, and element.

	Attributes
	----------
	handle				: int
	 	The facet index on the grid
	boundary			: PyEFVLib.Boundary
	 	The boundary to which facet belongs
	element				: PyEFVLib.Element
	 	The element to which facet belongs
	vertices			: list[PyEFVLib.Vertex]
	 	The vertices that lie on facet
	outerFaces			: list[PyEFVLib.OuterFace]
	 	The outerFaces that make facet
	area				: PyEFVLib.Point
	 	Facet's area
	centroid			: PyEFVLib.Point
	 	Facet's centroid
	elementLocalIndex	: int
	 	Its local index on element
	boundaryLocalIndex	: int
	 	Its local index on boundary
	"""

	def __init__(self, grid, verticesIndexes, handle):
		self.vertices = [grid.vertices[vertexIndex] for vertexIndex in verticesIndexes]
		self.handle = handle

		self.__computeAreaVector()
		self.__computeCentroid()
		self.__buildOuterFaces()

	def _setBoundary(self, boundary, boundaryLocalIndex):
		self.boundary = boundary
		self.boundaryLocalIndex = boundaryLocalIndex

	def _setElement(self, element, elementLocalIndex):
		self.element = element
		self.elementLocalIndex = elementLocalIndex

	def __computeAreaVector(self):
		if len(self.vertices) == 2:
			self.area = Point( self.vertices[0].y-self.vertices[1].y , self.vertices[1].x-self.vertices[0].x, 0.0 ) 

		if len(self.vertices) == 3:
			v0,v1,v2 = [vertex.getCoordinates() for vertex in self.vertices]
			self.area = Point( *np.cross(v1-v0,v2-v0)/2 )

		if len(self.vertices) == 4:
			CM = self.vertices[2] - self.vertices[0]
			LR = self.vertices[3] - self.vertices[1]
			x = 0.5 * (CM.y*LR.z - CM.z*LR.y)
			y = 0.5 * (CM.z*LR.x - CM.x*LR.z)
			z = 0.5 * (CM.x*LR.y - CM.y*LR.x)
			self.area = Point(x, y, z)

	def __computeCentroid(self):
		self.centroid = Point(*sum([vertex.getCoordinates() for vertex in self.vertices])/len(self.vertices))

	def __buildOuterFaces(self):
		self.outerFaces = []
		for vertex in self.vertices:
			outerFace = OuterFace(vertex, self)
			self.outerFaces.append(outerFace)