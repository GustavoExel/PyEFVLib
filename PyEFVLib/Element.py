import numpy as np
from PyEFVLib.Point import Point
from PyEFVLib.InnerFace import InnerFace

class Element:
	def __init__(self, grid, verticesIndexes, handle):
		self.handle = handle
		self.grid = grid
		self.vertices = [grid.vertices[vertexIndex] for vertexIndex in verticesIndexes]

		for vertex in self.vertices:
			vertex.addElement(self)
		
		self.innerFaces = []
		self.outerFaces = []
		self.faces      = []

		self.tellShape()
		self.buildInnerFaces()
		self.buildSubelement()

	def tellShape(self):
		for shape in self.grid.shapes:
			if shape._is(self):
				self.shape = shape
				return
		raise Exception("This element has not been registered in Grid yet")

	def buildInnerFaces(self):
		for i in range(self.shape.numberOfInnerFaces):
			innerFace = InnerFace(self, self.grid.innerFaceCounter, i)
			innerFace.setArea( self.shape.getInnerFaceAreaVector(i, self.centroid, self.vertices) )

			self.innerFaces.append(innerFace)
			self.faces.append(innerFace)

		self.grid.innerFaceCounter += self.shape.numberOfInnerFaces

	def buildSubelement(self):
		self.subelementVolumes = []
		self.volume = 0.0
		for local in range(self.numberOfVertices):
			shapeFunctionDerivatives = self.shape.subelementShapeFunctionDerivatives[local]
			volume = self.shape.subelementTransformedVolumes[local] * np.linalg.det(self.getTransposedJacobian(shapeFunctionDerivatives))

			self.volume += volume
			self.vertices[local].volume += volume
			self.subelementVolumes.append(volume)

	def setRegion(self, region):
		self.region = region

	def setOuterFace(self, outerFace):
		self.outerFaces.append(outerFace)
		self.faces.append(outerFace)

	@property
	def centroid(self):
		if hasattr(self, "_centroid"):
			return self._centroid
		else:
			self._centroid = sum( self.vertices ) / self.numberOfVertices
			return self._centroid

	@property
	def numberOfVertices(self):
		return len(self.vertices)

	def getTransposedJacobian(self, shapeFunctionDerivatives):	# shapeFunctionDerivatives must be already a numpy array
		vertices = np.array([vertex.getCoordinates()[:self.shape.dimension] for vertex in self.vertices])
		return np.matmul(np.transpose(shapeFunctionDerivatives), vertices)

	def getGlobalDerivatives(self, derivatives):
		return np.matmul(np.linalg.inv(self.getTransposedJacobian(derivatives)) , np.transpose(derivatives))
