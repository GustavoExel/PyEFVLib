import numpy as np
from PyEFVLib.Point import Point
from PyEFVLib.InnerFace import InnerFace

class Element:
	"""
	The element entity stores its shape, vertices, inner and outer faces, and can be
	used to calculate its subvolumes, and interpolations witin itself.

	It is defined by its vertices and its shape, and then is subdivided into m parts
	where m is its number of vertices. The faces that divide the element into subelements
	are called inner faces, and they also make up the control surface.

	The Element-based Finite Volume Method uses the elements to interpolate fields inside
	the element, which can then be used to interpolate gradients and spatial derivatives
	of those fields. The way it's done is by using shape functions.

	Each element has a local coordinate system (xi, eta, zeta), on which the locations of
	the elment vertices is known. Over those coordinates we define the shape functions, 
	which are functions N_k(xi, eta, zeta) corresponding to the vertex v_k. Then, we can
	calculate the change of coordinates in the following way

	x(xi, eta, zeta) = x_1 N_1(xi, eta, zeta) + ... + x_m N_m(xi, eta, zeta)

	In the matrix notation we would call the vector with the shape functions at some 
	location N, and the matrix with the vertices coordinates in its columns C, such that

	(x,y,z) = C.T @ N

	In a very analogous way we can interpolate fields inside the element. Let U be a vector
	with the values of the field u(x,y,z) at the vertices of the element. Then we can 
	interpolate the field u(x,y,z) in the location where N was computed by evaluating

	u(xi, eta, zeta) = u_1 N_1(xi, eta, zeta) + ... + u_m N_m(xi, eta, zeta)
	u = dot( U, N )

	If we differentiate the first interpolation equation and apply the chain rule, we can 
	also interpolate the gradient of the field u inside the element. Let D be the matrix
	containing the partial derivative of the shape functions with respect to the local
	coordinates, such that the row k contains the derivatives of N_k with respect to
	(xi, eta, zeta), and each column has the derivative with respect with each local
	coordinate. That way we define the Jacobian matrix J as being

	J = C.T @ D

	And the gradient of u(x,y,z) is given by

	grad(u) = inv(J.T) @ D.T @ U
	grad(u) = G @ U

	G = inv(J.T) @ D.T

	In order to compute all these matrices using PyEFVLib, we are going to use mainly the 
	`Element` and `Shape` classes. The `Shape` classes contains the matrices N (...ShapeFunctionValues),
	and D (...ShapeFunctionDerivatives). The class `Element` computes J (getTransposedJacobian)
	and G (getGlobalDerivatives) using D as input.

	Finally, for the computation of the subelement volumes, the determinant of the jacobian
	matrix is evaluated at the centroid of each subelement, and it's used as a conversion factor
	between the already known local volume stored in the `Shape` class and the actual volume.

	Attributes
	----------
	handle				: int
	 	Indexes all elements belonging to grid.
	region				: PyEFVLib.Region
	 	The region to which element belongs.
	vertices			: list[PyEFVLib.Vertex]
	 	A list of the vertices of element.
	innerFaces			: list[PyEFVLib.InnerFace]
	 	A list of the element's innerFaces.
	outerFaces			: list[PyEFVLib.OuterFace]
	 	A list of the element's outerFaces.
	faces				: list[PyEFVLib.Face]
	 	A list of the element's innerFaces and outerFaces.
	shape				: PyEFVLib.Shape
	 	The element shape.
	volume				: float
	 	The sum of all subElementVolumes (control volume section inside element).
	subElementVolumes	: list[float]
	 	The list with each vertex subElementVolume.
	centroid			: PyEFVLib.Point
	 	The element centroid.
	numberOfVertices	: int
	 	The number of vertices of the element.

	Methods
	----------
	getTransposedJacobian(shapeFunctionDerivatives: np.ndarray[m, d]) -> np.ndarray[d, d]:
		Returns the Jacobian matrix transposed in the location given by the shapeFunctionDerivatives
	getGlobalDerivatives(derivatives: np.ndarray[m, d]) -> np.ndarray[m, d]:
		Returns the gradient operator G	given the local derivatives	

	Where m = numberOfVertices and d = dimension
	"""
	def __init__(self, grid, verticesIndexes, handle):
		self.handle = handle
		self.grid = grid
		self.vertices = [grid.vertices[vertexIndex] for vertexIndex in verticesIndexes]

		for vertex in self.vertices:
			vertex._addElement(self)
		
		self.innerFaces = []
		self.outerFaces = []
		self.faces      = []

		self.__tellShape()
		self.__buildSubelement()
		self.__buildInnerFaces()

	def __tellShape(self):
		for shape in self.grid.shapes:
			if shape._is(self):
				self.shape = shape
				return
		raise Exception("This element has not been registered in Grid yet")

	def __buildInnerFaces(self):
		for i in range(self.shape.numberOfInnerFaces):
			innerFace = InnerFace(self, self.grid.innerFaceCounter, i)
			innerFace.area = self.shape.getInnerFaceAreaVector(i, self.centroid, self.vertices)

			self.innerFaces.append(innerFace)
			self.faces.append(innerFace)

		self.grid.innerFaceCounter += self.shape.numberOfInnerFaces

	def __buildSubelement(self):
		self.subelementVolumes = []
		self.volume = 0.0
		for local in range(self.numberOfVertices):
			shapeFunctionDerivatives = self.shape.subelementShapeFunctionDerivatives[local]
			subelementVolume = self.shape.subelementTransformedVolumes[local] * np.linalg.det(self.getTransposedJacobian(shapeFunctionDerivatives))

			if subelementVolume < 0.0:
				self.grid.correctForNegativeVolume = True
				self.grid.gridData.elementsConnectivities[self.handle] = self.grid.gridData.elementsConnectivities[self.handle][::-1]
			if subelementVolume == 0.0:
				print(f"Null volume detected at the element {self.handle}, subelement {local}")

			self.volume += subelementVolume
			self.vertices[local].volume += subelementVolume
			self.subelementVolumes.append(subelementVolume)

	def _setRegion(self, region):
		self.region = region

	def _setOuterFace(self, outerFace):
		self.outerFaces.append(outerFace)
		self.faces.append(outerFace)

	@property
	def centroid(self) -> Point:
		if hasattr(self, "_centroid"):
			return self._centroid
		else:
			self._centroid = sum( self.vertices ) / self.numberOfVertices
			return self._centroid

	@property
	def numberOfVertices(self) -> int:
		return len(self.vertices)

	def getTransposedJacobian(self, shapeFunctionDerivatives: np.ndarray) -> np.ndarray:
		"""
		Returns the gradient operator G	given the local derivatives	

		>>> help(PyEFVLib.Element) # for more info
		"""
		coords = np.array([vertex.getCoordinates()[:self.shape.dimension] for vertex in self.vertices])
		return shapeFunctionDerivatives.T @ coords

	def getGlobalDerivatives(self, shapeFunctionDerivatives: np.ndarray) -> np.ndarray:
		"""
		Returns the Jacobian matrix transposed in the location given by the shapeFunctionDerivatives

		>>> help(PyEFVLib.Element) # for more info
		"""
		return np.linalg.inv(self.getTransposedJacobian(shapeFunctionDerivatives)) @ shapeFunctionDerivatives.T