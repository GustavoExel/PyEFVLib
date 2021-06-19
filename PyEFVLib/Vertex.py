import numpy as np
from PyEFVLib.Point import Point

class Vertex(Point):
	def __init__(self, coordinates, handle):
		Point.__init__(self, *coordinates)
		self.handle = handle
		self.elements = []
		self.volume = 0.0

	def addElement(self, element):
		self.elements.append(element)

	def getLocal(self, element):
		return list(element.vertices).index(self)

	def getInnerFaces(self):
		return [ innerFace for element in self.elements for innerFace in element.innerFaces if self in innerFace.getNeighborVertices() ]

	def getSubElementVolume(self, element):
		return element.subelementVolumes[ self.getLocal(element) ]