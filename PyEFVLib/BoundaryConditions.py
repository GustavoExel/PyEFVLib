import numpy as np
from numpy import pi, sin, cos, tan, arcsin, arccos, arctan, sinh, cosh, tanh, arcsinh, arccosh, arctanh, sqrt, e , log, exp, inf, mod, floor

def getFunction(expr):
	def function(x,y,z,t):
		return eval( expr.replace('x',str(x)).replace('y',str(y)).replace('z',str(z)).replace('t',str(t)) )
	return function

class _BoundaryConditionBase:
	def __init__(self, grid, boundary, value, handle, expression=False):
		self.grid = grid
		self.boundary = boundary
		self.value = value
		self.handle = handle
		self.expression = expression

		if self.expression:
			self.function = getFunction(value)

	def getValue(self, index, time=0.0): # Index argument for future variable boundaryCondition
		if self.expression:
			x,y,z = self.grid.vertices[index].getCoordinates()
			return self.function(x,y,z,time)
		else:
			return self.value

class DirichletBoundaryCondition(_BoundaryConditionBase):
	__type__="DIRICHLET"

class NeumannBoundaryCondition(_BoundaryConditionBase):
	__type__="NEUMANN"
