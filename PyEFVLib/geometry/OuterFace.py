from PyEFVLib.geometry.Point import Point
import numpy as np

class OuterFace:
	def __init__(self, vertex, facet):
		self.vertex = vertex
		self.facet = facet

	def setLocal(self, local):
		self.local = local

	def computeCentroid(self):
		shapeFunctionValues = self.facet.element.shape.outerFaceShapeFunctionValues[self.facet.elementLocalIndex][self.local]
		elementVerticesCoords = [vertex.getCoordinates()[:self.facet.element.shape.dimension] for vertex in self.facet.element.vertices]

		self.centroid = Point(*np.dot(shapeFunctionValues, elementVerticesCoords))

	def computeAreaVector(self):
		self.area = self.facet.area / self.facet.vertices.size