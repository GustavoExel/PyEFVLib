import numpy as np

class DirichletBoundaryCondition:
	__type__="DIRICHLET"
	def __init__(self, boundary, value, handle):
		self.boundary = boundary
		self.value = value
		self.handle = handle

	def getValue(self, index): # Index argument for future variable boundaryCondition
		return self.value

class NeumannBoundaryCondition:
	__type__="NEUMANN"
	def __init__(self, boundary, value, handle):
		self.boundary = boundary
		self.value = value
		self.handle = handle

	def getValue(self, index): # Index argument for future variable boundaryCondition
		return self.value
