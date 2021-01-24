import numpy as np

class Region:
	def __init__(self, grid, elementsIndexes, name, handle):
		self.name = name
		self.grid = grid
		self.handle = handle
		self.elements = np.array([grid.elements[elementIndex] for elementIndex in elementsIndexes])
		for element in self.elements:
			element.setRegion(self)