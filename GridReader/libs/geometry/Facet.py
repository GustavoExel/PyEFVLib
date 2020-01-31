import numpy as np

class Facet:
	def __init__(self, element, elementLocalIndex, handle):
		self.element = element
		self.elementLocalIndex = elementLocalIndex
		self.handle = handle
		self.vertices = np.array([])
		self.outerFaces = np.array([])

	def addVertex(self, vertex):
		self.vertices = np.append(self.vertices, vertex)