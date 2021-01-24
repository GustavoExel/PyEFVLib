from PyEFVLib.geometry.Point import Point
from PyEFVLib.geometry.OuterFace import OuterFace
import numpy as np

class Facet:
	def __init__(self, grid, verticesIndexes, handle):
		self.vertices = np.array([grid.vertices[vertexIndex] for vertexIndex in verticesIndexes])
		self.handle = handle

		self.computeAreaVector()
		self.computeCentroid()
		self.buildOuterFaces()

	def setBoundary(self, boundary, boundaryLocalIndex):
		self.boundary = boundary
		self.boundaryLocalIndex = boundaryLocalIndex

	def setElement(self, element, elementLocalIndex):
		self.element = element
		self.elementLocalIndex = elementLocalIndex

	def computeAreaVector(self):
		if self.vertices.size == 2:
			self.area = Point( self.vertices[0].y-self.vertices[1].y , self.vertices[1].x-self.vertices[0].x, 0.0 ) 

		if self.vertices.size == 3:
			v0,v1,v2 = [vertex.getCoordinates() for vertex in self.vertices]
			self.area = Point( *np.cross(v1-v0,v2-v0)/2 )

		if self.vertices.size == 4:
			CM = self.vertices[2] - self.vertices[0]
			LR = self.vertices[3] - self.vertices[1]
			x = 0.5 * (CM.y*LR.z - CM.z*LR.y)
			y = 0.5 * (CM.z*LR.x - CM.x*LR.z)
			z = 0.5 * (CM.x*LR.y - CM.y*LR.x)
			self.area = Point(x, y, z)

	def computeCentroid(self):
		self.centroid = Point(*sum([vertex.getCoordinates() for vertex in self.vertices])/self.vertices.size)

	def buildOuterFaces(self):
		self.outerFaces = np.array([])
		for vertex in self.vertices:
			outerFace = OuterFace(vertex, self)
			self.outerFaces = np.append(self.outerFaces, outerFace)