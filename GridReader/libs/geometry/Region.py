import numpy as np

class Region:
	def __init__(self, elements, name, handle):
		self.handle = handle
		self.name = name
		self.elements = np.array(elements)
		for element in elements:
			element.setRegion(self)

	def setGrid(self, grid):
		self.grid = grid

	def getGrid(self):
		return self.grid

	def getIndex(self):
		return self.handle