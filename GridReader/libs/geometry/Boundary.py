from libs.geometry.Facet import Facet
from libs.geometry.Point import Point
from libs.geometry.OuterFace import OuterFace
import numpy as np

class BoundaryData:
	def __init__(self, name, facetsConnectivity, handle):
		self.name = name
		self.handle = handle
		self.facetsConnectivity = facetsConnectivity
		self.vertices = list(set(sum(facetsConnectivity,[])))


class BoundaryBuilder:
	def __init__(self, grid):
		self.grid = grid
		self.buildBoundaryData()
		self.buildBoundaries()

	def buildBoundaryData(self):
		self.boundaries = []

		names = self.grid.gridData.boundaryNames
		connectivities = self.grid.gridData.boundaryElementsConnectivity
		facets = self.grid.gridData.boundaryElements

		i = 0
		for name, facet in zip(names, facets):
			facetConnectivity = [ connectivities[f] for f in facet ]
			self.boundaries.append( BoundaryData(name, facetConnectivity, i) )
			i += 1

	def buildBoundaries(self):
		self.facetHandle = 0
		self.handleOfFirstOuterFace = 0
		for boundaryData in self.boundaries:
			boundary = Boundary(boundaryData.name, self.handleOfFirstOuterFace, boundaryData.handle)
			self.addBoundaryFacets(boundary, boundaryData)
			for handle in boundaryData.vertices:
				boundary.addVertex(self.grid.vertices[handle])
			self.grid.boundaries = np.append(self.grid.boundaries, boundary)

	def addBoundaryFacets(self, boundary, boundaryData):
		self.scanElements(boundaryData)
		for facetConnectivity in boundaryData.facetsConnectivity:
			facet = self.buildFacet(facetConnectivity)

			facet.handleOfFirstOuterFace = self.handleOfFirstOuterFace
			facet.area = self.computeFacetAreaVector(facet.vertices)
			self.buildOuterFaces(facet)
			boundary.addFacet(facet)

			self.handleOfFirstOuterFace += facet.vertices.size


	def scanElements(self, boundaryData):
		# Keeps track of which elements share faces with the boundary
		self.boundaryElementsVertices = []
		self.boundaryElements = []
		for element, elementVertices in zip(self.grid.elements, self.grid.gridData.elemConnectivity):
			if len(set(elementVertices).intersection(boundaryData.vertices)) >= element.shape.dimension:
				self.boundaryElementsVertices.append(elementVertices)
				self.boundaryElements.append(element)

	def buildFacet(self, facetConnectivity):
		for boundaryElement, boundaryElementVertices in zip(self.boundaryElements, self.boundaryElementsVertices):
			if set(facetConnectivity).issubset(boundaryElementVertices):
				localFacetVertices = [boundaryElementVertices.index(globalHandle) for globalHandle in facetConnectivity]
				for elemFacetIndex in range(boundaryElement.shape.numberOfFacets):
					if set(localFacetVertices) == set(boundaryElement.shape.facetVerticesIndices[elemFacetIndex]):
						facet = Facet(boundaryElement, elemFacetIndex, self.facetHandle)
						# Can also specify the shape of the facet
						for local in boundaryElement.shape.facetVerticesIndices[elemFacetIndex]:
							facet.addVertex(boundaryElement.vertices[local])
						return facet

	def computeFacetAreaVector(self, vertices):
		if vertices.size == 2:
			return Point( vertices[0].y-vertices[1].y , vertices[1].x-vertices[0].x, 0.0 ) 

		if vertices.size == 3:
		    d10 = Point(vertices[1].x-vertices[0].x, vertices[1].y-vertices[0].y, vertices[1].z-vertices[0].z)
		    d20 = Point(vertices[2].x-vertices[0].x, vertices[2].y-vertices[0].y, vertices[2].z-vertices[0].z)
		    x = (d10.y*d20.z - d20.y*d10.z) / 2.0
		    y = (d10.z*d20.x - d20.z*d10.x) / 2.0
		    z = (d10.x*d20.y - d20.x*d10.y) / 2.0
		    return Point(x, y, z);

		if vertices.size == 4:
		    CM = Point(vertices[2].x-vertices[0].x, vertices[2].y-vertices[0].y, vertices[2].z-vertices[0].z)
		    LR = Point(vertices[3].x-vertices[1].x, vertices[3].y-vertices[1].y, vertices[3].z-vertices[1].z)
		    x = 0.5 * (CM.y()*LR.z() - CM.z()*LR.y())
		    y = 0.5 * (CM.z()*LR.x() - CM.x()*LR.z())
		    z = 0.5 * (CM.x()*LR.y() - CM.y()*LR.x())
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
			facet.outerFaces = np.append(facet.outerFaces, outerFace)

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