class GridData:
	def __init__(self, path):
		self.path = path
		
	def setVertices(self, vertices):
		self.vertices = vertices

	def setElementConnectivity(self, connectivity):
		self.elemConnectivity = connectivity

	def setRegions(self, regionNames, regionElements):
		self.regionNames = regionNames
		self.regionElements = regionElements

	def setBoundaries(self, boundaryNames, boundaryElements, boundaryElementsConnectivity):
		self.boundaryElementsConnectivity = boundaryElementsConnectivity
		self.boundaryNames = boundaryNames
		self.boundaryElements = boundaryElements

	def setShapes(self, shapes):
		self.lines, self.triangles, self.quadrilaterals, self.tetrahedrons, self.hexahedrons, self.prisms, self.pyramids = shapes