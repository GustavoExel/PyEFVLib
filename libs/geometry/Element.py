import numpy as np
from libs.geometry.Point import Point
from libs.geometry.InnerFace import InnerFace

class Element:
	def __init__(self, vertices, grid, handle):
		self.handle = handle
		self.grid = grid
		self.vertices = np.array(vertices)

		for vertex in vertices:
			vertex.addElement(self)

		self.tellShape()
		self.buildInnerFaces()
		self.buildSubelement()
		
	def tellShape(self):
		for shape in self.grid.getShapes():
			if shape._is(self):
				self.shape = shape(self)
				return
		raise Exception("Either this element isn\'t 2-dimensional or isn\'t registered yet")

	def buildInnerFaces(self):
		self.innerFaces = np.array([])
		centroid = Point(*sum([v.getCoordinates() for v in self.vertices])/self.vertices.size)
		for i in range(self.shape.numberOfInnerFaces):
			innerFace = InnerFace(self, self.grid.innerFaceCounter, i)
			innerFace.area = self.shape.getInnerFaceAreaVector(i,centroid, self.vertices)

			self.innerFaces = np.append(self.innerFaces, innerFace)

		self.grid.innerFaceCounter += self.shape.numberOfInnerFaces

	def buildSubelement(self):
		self.subelementVolumes = []
		self.volume = 0.0
		for local in range(self.vertices.size):
			shapeFunctionDerivatives = self.shape.subelementShapeFunctionDerivatives[local]
			volume = self.shape.subelementTransformedVolumes[local] * np.linalg.det(self.getTransposeJacobian(shapeFunctionDerivatives))

			self.volume += volume 
			self.vertices[local].volume += volume
			self.subelementVolumes.append(volume)

	def setRegion(self, region):
		self.region = region

	def getTransposeJacobian(self, shapeFunctionDerivatives):	# shapeFunctionDerivatives must be already a numpy array
		vertices = np.array([[vertex.getCoordinates()[k] for k in range(self.shape.dimension)] for vertex in self.vertices])
		return np.matmul(np.transpose(shapeFunctionDerivatives), vertices)