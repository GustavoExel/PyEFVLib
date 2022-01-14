import numpy as np
from PyEFVLib.Point import Point

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

	def getShapeFunctions(self):
		return self.element.shape.innerFaceShapeFunctionValues[self.local]

	def getNeighborVerticesLocals(self):
		return self.element.shape.innerFaceNeighborVertices[self.local][:2]

	def getNeighborVerticesHandles(self):
		backwardLocal, forwardLocal = self.getNeighborVerticesLocals()
		return self.element.vertices[backwardLocal].handle, self.element.vertices[forwardLocal].handle

"""
		m = self.element.numberOfVertices
		d = self.element.grid.dimension
		transposedJacobian = [[0.0 for _ in range(d)] for _ in range(d)]
		invTransposedJacobian = [[0.0 for _ in range(d)] for _ in range(d)]
		self.globalDerivatives = [[0.0 for _ in range(m)] for _ in range(d)]

		innerFaceShapeFunctionDerivatives = self.element.shape.innerFaceShapeFunctionDerivatives[self.local]

		for i in range(d):
			for j in range(d):
				# transposedJacobian[i][j] = 0.0
				for k in range(m):
					transposedJacobian[i][j] += innerFaceShapeFunctionDerivatives[k][i] * self.element.vertices[k].getCoordinates()[j]

		jacobianDeterminant = 0.0
		if d==2:
			jacobianDeterminant = transposedJacobian[0][0]*transposedJacobian[1][1]-transposedJacobian[0][1]*transposedJacobian[1][0]

			invTransposedJacobian[0][0] = transposedJacobian[1][1] / jacobianDeterminant
			invTransposedJacobian[0][1] = -transposedJacobian[0][1] / jacobianDeterminant
			invTransposedJacobian[1][0] = -transposedJacobian[1][0] / jacobianDeterminant
			invTransposedJacobian[1][1] = transposedJacobian[0][0] / jacobianDeterminant
		elif d==3:
			jacobianDeterminant = transposedJacobian[0][0]*transposedJacobian[1][1]*transposedJacobian[2][2]+transposedJacobian[0][1]*transposedJacobian[1][2]*transposedJacobian[2][0]+transposedJacobian[1][0]*transposedJacobian[2][1]*transposedJacobian[0][2]-transposedJacobian[0][2]*transposedJacobian[1][1]*transposedJacobian[2][0]-transposedJacobian[0][1]*transposedJacobian[1][0]*transposedJacobian[2][2]-transposedJacobian[1][2]*transposedJacobian[2][1]*transposedJacobian[0][0]

			invTransposedJacobian[0][0] = ( transposedJacobian[1][1]*transposedJacobian[2][2]-transposedJacobian[1][2]*transposedJacobian[2][1] ) / jacobianDeterminant
			invTransposedJacobian[0][1] = ( transposedJacobian[0][2]*transposedJacobian[2][1]-transposedJacobian[0][1]*transposedJacobian[2][2] ) / jacobianDeterminant
			invTransposedJacobian[0][2] = ( transposedJacobian[0][1]*transposedJacobian[1][2]-transposedJacobian[0][2]*transposedJacobian[1][1] ) / jacobianDeterminant
			invTransposedJacobian[1][0] = ( transposedJacobian[1][2]*transposedJacobian[2][0]-transposedJacobian[1][0]*transposedJacobian[2][2] ) / jacobianDeterminant
			invTransposedJacobian[1][1] = ( transposedJacobian[0][0]*transposedJacobian[2][2]-transposedJacobian[0][2]*transposedJacobian[2][0] ) / jacobianDeterminant
			invTransposedJacobian[1][2] = ( transposedJacobian[0][2]*transposedJacobian[1][0]-transposedJacobian[0][0]*transposedJacobian[1][2] ) / jacobianDeterminant
			invTransposedJacobian[2][0] = ( transposedJacobian[1][0]*transposedJacobian[2][1]-transposedJacobian[1][1]*transposedJacobian[2][0] ) / jacobianDeterminant
			invTransposedJacobian[2][1] = ( transposedJacobian[0][1]*transposedJacobian[2][0]-transposedJacobian[0][0]*transposedJacobian[2][1] ) / jacobianDeterminant
			invTransposedJacobian[2][2] = ( transposedJacobian[0][0]*transposedJacobian[1][1]-transposedJacobian[0][1]*transposedJacobian[1][0] ) / jacobianDeterminant

		for i in range(d):
			for j in range(m):
				for k in range(d):
					self.globalDerivatives[i][j] += invTransposedJacobian[i][k] * innerFaceShapeFunctionDerivatives[j][k]

"""