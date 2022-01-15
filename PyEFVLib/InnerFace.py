import numpy as np
from PyEFVLib.Point import Point

class InnerFace:
	"""
	The InnerFace is a very important class, because it composes the control surface of
	each control volume, and all surface intgrals in the conservation equations will be
	computed in the innerFaces.

	Attributes
	----------
	handle				: int
		the innerFace index in the grid
	local				: int
		the innerFace local index in the element
	element				: int
		the innerFace element
	area				: PyEFVLib.Point
		the innerFace area
	centroid			: PyEFVLib.Point
		the innerFace centroid
	globalDerivatives	: np.darray[d, m]
		The gradient opeerator evaluated at the centroid of the innerFace.
	"""
	def __init__(self, element, handle, local):
		self.handle = handle
		self.local = local
		self.element = element

		self.__evaluateCentroid()
		self.__calculateGlobalDerivatives()

	def __evaluateCentroid(self):
		shapeFunctionValues = self.element.shape.innerFaceShapeFunctionValues[self.local]
		coords = np.array([0.0,0.0,0.0])
		for vertex, weight in zip(self.element.vertices, shapeFunctionValues):
			coords += weight * vertex.getCoordinates()
		self.centroid = Point(*coords)

	def __calculateGlobalDerivatives(self):
		derivatives = self.element.shape.innerFaceShapeFunctionDerivatives[self.local]
		self.globalDerivatives = np.linalg.inv(self.element.getTransposedJacobian(derivatives)) @ derivatives.T

	def getNeighborVertices(self):
		"""
		Returns the vertices that are right next to the innerFace. That means that the
		innerFace belongs to their control surface. The difference between the backward
		and forward vertex, is that the area vector points away from the backward vertex
		and towards the forward vertex. When integrating in the control surface the
		default is to integrate with the area vector pointing away from the center, that's
		why the sign is flipped for the forward vertex.
		"""
		backwardVertex = self.element.vertices[ self.element.shape.innerFaceNeighborVertices[self.local][0] ]
		forwardVertex = self.element.vertices[ self.element.shape.innerFaceNeighborVertices[self.local][1] ]
		return backwardVertex, forwardVertex

	def getShapeFunctions(self):
		"""
		Returns the shape functions evaluated at the innerFace centroid. The order in which
		they apper refers to the order of the vertices in the element.vertices list.
		"""
		return self.element.shape.innerFaceShapeFunctionValues[self.local]

	def getNeighborVerticesLocals(self):
		"""
		Returns the local indices in the element.vertices list of the vertices that are right
		next to the innerFace. That means that the innerFace belongs to their control surface.
		The difference between the backward and forward vertex, is that the area vector points
		away from the backward vertex and towards the forward vertex. When integrating in the
		control surface the default is to integrate with the area vector pointing away from the
		center, that's why the sign is flipped for the forward vertex.
		"""
		return self.element.shape.innerFaceNeighborVertices[self.local][:2]

	def getNeighborVerticesHandles(self):
		"""
		Returns the global indices in the grid.vertices list of the vertices that are right
		next to the innerFace. That means that the innerFace belongs to their control surface.
		The difference between the backward and forward vertex, is that the area vector points
		away from the backward vertex and towards the forward vertex. When integrating in the
		control surface the default is to integrate with the area vector pointing away from the
		center, that's why the sign is flipped for the forward vertex.
		"""
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