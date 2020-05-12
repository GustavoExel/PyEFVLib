import numpy as np

class DirichletBoundaryCondition:
	def __init__(self, boundary, value, handle):
		self.boundary = boundary
		self.value = value
		self.handle = handle

	def getValue(self, index): # Index argument for future variable boundaryCondition
		return self.value

class NeumannBoundaryCondition:
	def __init__(self, boundary, value, handle):
		self.boundary = boundary
		self.value = value
		self.handle = handle

	def getValue(self, index): # Index argument for future variable boundaryCondition
		return self.value
