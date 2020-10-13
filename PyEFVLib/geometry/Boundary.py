from PyEFVLib.geometry.Facet import Facet
from PyEFVLib.geometry.Point import Point
from PyEFVLib.geometry.OuterFace import OuterFace
import numpy as np

class Boundary:
	def __init__(self, name, handleOfFirstOuterFace, handle):
		self.name = name
		self.handleOfFirstOuterFace = handleOfFirstOuterFace
		self.handle = handle
		self.vertices = np.array([])
		self.facets = np.array([])

	def addVertex(self, vertex):
		self.vertices = np.append(self.vertices, vertex)

	def addFacet(self, facet):
		self.facets = np.append(self.facets, facet)

class BoundaryData:
	def __init__(self, name, connectivity, handle):
		self.name = name
		self.handle = handle
		self.connectivity = connectivity
		self.vertices = list(set(sum(connectivity,[])))

class BoundaryBuilder:
	def __init__(self, grid):
		self.grid = grid
		self.buildBoundaryData()
		self.buildBoundaries()

	def buildBoundaryData(self):
		self.boundaries = []

		names = self.grid.gridData.boundariesNames
		connectivities = self.grid.gridData.boundariesConnectivities
		boundaries = self.grid.gridData.boundariesIndexes

		i = 0
		for name, boundary in zip(names, boundaries):
			boundaryConnectivity = [ connectivities[b] for b in boundary ]
			self.boundaries.append( BoundaryData(name, boundaryConnectivity, i) )
			i += 1

	def buildBoundaries(self):
		self.facetHandle = 0
		self.handleOfFirstOuterFace = 0

		for boundaryData in self.boundaries:
			boundary = Boundary(boundaryData.name, self.handleOfFirstOuterFace, boundaryData.handle)
			# The boundary recieves its facets, which recieve its outer faces
			# The facets have an element, vertices, an area, a global and a (element) local index
			# The outer faces have an area, a centroid and a vertex.
			self.addBoundaryFacets(boundary, boundaryData)
			
			# The boundary gets its vertices
			for handle in boundaryData.vertices:
				boundary.addVertex(self.grid.vertices[handle])
			
			self.grid.boundaries = np.append(self.grid.boundaries, boundary)

	def addBoundaryFacets(self, boundary, boundaryData):
		# Keeps track of which elements share faces with the boundary
		self.scanElements(boundaryData)
		for facetConnectivity in boundaryData.connectivity:
			# The facet recieves its vertices and its element.
			facet = self.buildFacet(facetConnectivity)

			facet.handleOfFirstOuterFace = self.handleOfFirstOuterFace
			facet.area = self.computeFacetAreaVector(facet.vertices)
			# The outer face gets an area, its centroid, and a vertex, and the facet gets its outer faces
			self.buildOuterFaces(facet)
			boundary.addFacet(facet)

			self.handleOfFirstOuterFace += facet.vertices.size


	def scanElements(self, boundaryData):
		# Keeps track of which elements share faces with the boundary
		self.boundariesIndexesVertices = []
		self.boundariesIndexes = []
		for element, elementVertices in zip(self.grid.elements, self.grid.gridData.elementsConnectivities):
			if len(set(elementVertices).intersection(boundaryData.vertices)) >= element.shape.dimension:
				self.boundariesIndexesVertices.append(elementVertices)
				self.boundariesIndexes.append(element)

	def buildFacet(self, facetConnectivity):
		# boundaryElement 	 : Element which contains the facet
		# localFacetVertices : facet's vertices local indices at the element
		# elemFacetIndex	 : facet's local index of the element
		for boundaryElement, boundaryElementVertices in zip(self.boundariesIndexes, self.boundariesIndexesVertices):
			# If both facet nodes belong to boundaryElementVertices, it's its element
			if set(facetConnectivity).issubset(boundaryElementVertices):
				localFacetVertices = [boundaryElementVertices.index(globalHandle) for globalHandle in facetConnectivity]
				
				for elemFacetIndex in range(boundaryElement.shape.numberOfFacets):
					# Check which facet of the element contains the boundary
					if set(localFacetVertices) == set(boundaryElement.shape.facetVerticesIndices[elemFacetIndex]):
						facet = Facet(boundaryElement, elemFacetIndex, self.facetHandle)
						
						for local in localFacetVertices:
							facet.addVertex(boundaryElement.vertices[local])
						
						self.facetHandle += 1
						return facet
				raise Exception("Unknown error")
		raise Exception("Boundary belongs to no element")

	def computeFacetAreaVector(self, vertices):
		if vertices.size == 2:
			return Point( vertices[0].y-vertices[1].y , vertices[1].x-vertices[0].x, 0.0 ) 

		if vertices.size == 3:
			v0,v1,v2=vertices[0].getCoordinates(), vertices[1].getCoordinates(), vertices[2].getCoordinates(), 
			return Point( *np.cross(v1-v0,v2-v0)/2 )

		if vertices.size == 4:
		    CM = vertices[2] - vertices[0]
		    LR = vertices[3] - vertices[1]
		    x = 0.5 * (CM.y*LR.z - CM.z*LR.y)
		    y = 0.5 * (CM.z*LR.x - CM.x*LR.z)
		    z = 0.5 * (CM.x*LR.y - CM.y*LR.x)
		    return Point(x, y, z)

	def buildOuterFaces(self, facet):
		for o in range(facet.vertices.size):
			outerFace = OuterFace(facet.vertices[o], facet, o, facet.handleOfFirstOuterFace + o)

			weights = facet.element.shape.outerFaceShapeFunctionValues[facet.elementLocalIndex][o]
			centroidCoord = np.zeros(3)
			for elemVertex, weight in zip(facet.element.vertices, weights):
				centroidCoord += elemVertex.getCoordinates() * weight
			outerFace.centroid = Point(*centroidCoord)

			outerFace.area = Point(*(facet.area.getCoordinates() / facet.vertices.size))
			facet.addOuterFace(outerFace)
