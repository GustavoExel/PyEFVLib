import numpy as np
from libs.geometry.Point import Point

class Vertex(Point):
	def __init__(self, coordinates, handle):
		Point.__init__(self, *coordinates)
		self.handle = handle
		self.elements = np.array([])
		self.volume = 0.0

	def addElement(self, element):
		self.elements = np.append(self.elements, element)

	def getElements(self):
		return self.elements

	def getIndedx(self):
		return self.handle