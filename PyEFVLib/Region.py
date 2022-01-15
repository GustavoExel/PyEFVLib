import numpy as np

class Region:
	"""
	The region class serves to store physical properties, such that we can have
	different regions with different properties across the domain. One important
	observation is that the equations are applied to each vertex, but the properties
	are stored in the elements (which belong to region).

	Attributes
	----------
	name	 : str
		Identifies the region in order to differentiate different regions
	handle	 : int
		Indexes all regions belonging to grid.
	vertices : list[PyEFVLib.Vertex]
		A list with the region vertices
	elements	 : list[PyEFVLib.Element]
		A list with the region elements
	"""
	def __init__(self, grid, elementsIndexes, name, handle):
		self.name = name
		self.grid = grid
		self.handle = handle
		self.elements = [grid.elements[elementIndex] for elementIndex in elementsIndexes]
		for element in self.elements:
			element._setRegion(self)