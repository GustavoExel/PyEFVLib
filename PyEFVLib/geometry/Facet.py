from PyEFVLib.geometry.Point import Point
import numpy as np

class Facet:
	def __init__(self, element, elementLocalIndex, handle):
		self.element = element
		self.elementLocalIndex = elementLocalIndex
		self.handle = handle
		self.area = Point(0.0, 0.0, 0.0)
		self.vertices = np.array([])
		self.outerFaces = np.array([])

	def addVertex(self, vertex):
		self.vertices = np.append(self.vertices, vertex)

	def addOuterFace(self, outerFace):
		self.outerFaces = np.append(self.outerFaces, outerFace)
