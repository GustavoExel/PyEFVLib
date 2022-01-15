from PyEFVLib.Point import Point
import numpy as np

class OuterFace:
	"""
	The outer faces are used when applying the boundary conditions, and we need to
	use the area vector of the face on which we are applying the boundary condition.
	When applying Neumann boundary conditions we also need the shape functions and
	their derivatives evaluated at the outer face centroid.

	Attributes
	----------
	local	 : int
		Local index on facet.outerFaces
	vertex	 : PyEFVLib.Vertex
		Neighboring vertex
	facet	 : PyEFVLib.Facet
		Facet to which outerFace belongs
	centroid : PyEFVLib.Point
		Outer face's centroid
	area	 : PyEFVLib.Point
		Outer face's area vector
	"""
	def __init__(self, vertex, facet):
		self.vertex = vertex
		self.facet = facet

	def _build(self):
		self.__computeCentroid()
		self.__computeAreaVector()
		self.__computeGlobalDerivatives()


	def __computeCentroid(self):
		shapeFunctionValues = self.facet.element.shape.outerFaceShapeFunctionValues[self.facet.elementLocalIndex][self.local]
		elementVerticesCoords = [vertex.getCoordinates()[:self.facet.element.shape.dimension] for vertex in self.facet.element.vertices]

		self.centroid = Point(*np.dot(shapeFunctionValues, elementVerticesCoords))

	def __computeAreaVector(self):
		self.area = self.facet.area / len(self.facet.vertices)

	def __computeGlobalDerivatives(self):
		if self.facet.element.region.grid.dimension == 2:
			localDerivatives = self.facet.element.shape.vertexShapeFunctionDerivatives[ self.vertex.getLocal(self.facet.element) ]
			self.globalDerivatives = self.facet.element.getGlobalDerivatives(localDerivatives)
		else:
			pass
			# raise Exception("Need to implement 3D vertexShapeFunctionDerivatives")

	def getShapeFunctions(self):
		"""
		Returns the shape functions evaluated at the outerFace centroid. The order in which
		they apper refers to the order of the vertices in the facet.element.vertices list.
		"""
		return self.facet.element.shape.outerFaceShapeFunctionValues[self.facet.elementLocalIndex][self.local]