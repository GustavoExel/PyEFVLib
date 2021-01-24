import numpy as np
from PyEFVLib.geometry.Point import Point

class InnerFace:
	def __init__(self, element, handle, local):
		self.handle = handle
		self.local = local
		self.element = element

		self.evaluateCentroid()
		self.calculateGlobalDerivatives()

	def evaluateCentroid(self):
		shapeFunctionValues = self.element.shape.innerFaceShapeFunctionValues[self.local]
		coords = np.array([0.0,0.0,0.0])
		for vertex, weight in zip(self.element.vertices, shapeFunctionValues):
			coords += weight * vertex.getCoordinates()
		self.centroid = Point(*coords)

	def setArea(self, area):
		self.area = area

	def calculateGlobalDerivatives(self):
		derivatives = self.element.shape.innerFaceShapeFunctionDerivatives[self.local]
		self.globalDerivatives = np.matmul(np.linalg.inv(self.element.getTransposedJacobian(derivatives)) , np.transpose(derivatives))

	def getNeighborVertices(self):
		backwardVertex = self.element.vertices[ self.element.shape.innerFaceNeighborVertices[self.local][0] ]
		forwardVertex = self.element.vertices[ self.element.shape.innerFaceNeighborVertices[self.local][1] ]
		return backwardVertex, forwardVertex
