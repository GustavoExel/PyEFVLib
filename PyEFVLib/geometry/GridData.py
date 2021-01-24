class GridData:
	def __init__(self, path):
		self.path = path

	def setDimension(self, gridDimension):
		self.dimension = gridDimension

	def setVertices(self, verticesCoordinates):
		self.verticesCoordinates = verticesCoordinates

	def setElementConnectivity(self, elementsConnectivities):
		self.elementsConnectivities = elementsConnectivities

	def setRegions(self, regionsNames, regionsElementsIndexes):
		self.regionsNames = regionsNames
		self.regionsElementsIndexes = regionsElementsIndexes

	def setBoundaries(self, boundariesNames, boundariesIndexes, boundariesConnectivities):
		self.boundariesConnectivities = boundariesConnectivities
		self.boundariesNames = boundariesNames
		self.boundariesIndexes = boundariesIndexes

	def setShapes(self, shapes):
		self.lines, self.triangles, self.quadrilaterals, self.tetrahedrons, self.hexahedrons, self.prisms, self.pyramids = shapes
