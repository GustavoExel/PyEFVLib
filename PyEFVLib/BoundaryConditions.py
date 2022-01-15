import numpy as np
from numpy import pi, sin, cos, tan, arcsin, arccos, arctan, sinh, cosh, tanh, arcsinh, arccosh, arctanh, sqrt, e , log, exp, inf, mod, floor

def getFunction(expr):
	def function(x,y,z,t):
		return eval( expr.replace('x',str(x)).replace('y',str(y)).replace('z',str(z)).replace('t',str(t)) )
	return function

class _BoundaryConditionBase:
	"""
	The purpose of this class is to assign a boundary condition from ProblemData
	with the possibility of assigning a variable boundary condition.

	When implementing the boundary conditions, getValue(vertex.handle, time=time)
	returns a constant value or evaluates an expression which can be function of
	x, y, z and t.

	DirichletBoundaryCondition prescribes a fixed value for the field at the boundary
	NeumannBoundaryCondition prescribes the flux flowing in and out of the boundary

	Attributes
	----------
	boundary	: PyEFVLib.Boundary
		The boundary that the boundary condition object refers to.
	value		: None
		stores a constant float `if not expression else` stores an expression string
	handle		: int
		Indexes all boundary conditions belonging to problemData.
	expression	: bool
		Keeps track of whether value contains a constant float or a string expression
	function	: function
		Is the function that evaluates the expression
	__type__	: str
		Can be either 'DIRICHLET' or 'NEUMANN'

	Methods
	----------
	getValue(index: int, time: float = 0.0) -> float:
		Returns the value of the boundary condition at the grid.vertices[index] location
		which can be a constant value or a function of x, y, z, t.
	"""
	def __init__(self, grid, boundary, value, handle: int, expression: bool = False):
		self.grid = grid
		self.boundary = boundary
		self.value = value
		self.handle = handle
		self.expression = expression

		if self.expression:
			self.function = getFunction(value)

	def getValue(self, index: int, time: float = 0.0) -> float:	
		"""
		Args:
			index (int)  : refers to the vertex.handle from which x, y, z will be extracted
			time (float) : refers to the time 't' at which the expression will be evaluated
		Returns:
			if self.expression:
				vertex = self.grid.vertices[index]
				return self.function( vertex.x, vertex.y, vertex.z, time )
			else:
				return self.value
		"""
		if self.expression:
			vertex = self.grid.vertices[index]
			return self.function( vertex.x, vertex.y, vertex.z, time )
		else:
			return self.value

class DirichletBoundaryCondition(_BoundaryConditionBase):
	__type__="DIRICHLET"

class NeumannBoundaryCondition(_BoundaryConditionBase):
	__type__="NEUMANN"
