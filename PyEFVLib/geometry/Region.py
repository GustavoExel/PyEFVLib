import numpy as np

class Region:
	def __init__(self, elements, name, grid, handle):
		self.name = name
		self.grid = grid
		self.handle = handle
		self.elements = np.array(elements)
		for element in elements:
			element.setRegion(self)