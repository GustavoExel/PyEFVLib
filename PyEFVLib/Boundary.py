from PyEFVLib.Facet import Facet
from PyEFVLib.Point import Point
from PyEFVLib.OuterFace import OuterFace
import numpy as np

class Boundary:
	"""
	Contains the boundary objects which stores all the vertices, facets and outer facets.
	It's important when applying the boundary conditions, for getting the vertices belonging
	to the boundary, and computing the outer faces areas.

	Attributes
	----------
	name	 : str
		Identifies the boundary in order to apply the boundary conditions.
	handle	 : int
		Indexes all boundaries belonging to grid.
	vertices : list[PyEFVLib.Vertex]
		A list with the boundary vertices
	facets	 : list[PyEFVLib.Facet]
		A list with the boundary facets
	"""
	def __init__(self, grid, facetsIdexes: list, name: str, handle: int):
		self.grid = grid
		self.name = name
		self.handle = handle

		self.facets = [grid.facets[facetIndex] for facetIndex in facetsIdexes]

		self.__setVertices()
		self.__matchElementsToFacets()
		self.__setFacets()
		self.__setOuterFaces()

	def __setVertices(self):
		vertices = []
		for facet in self.facets:
			for vertex in facet.vertices:
				if vertex not in vertices:
					vertices.append(vertex)
		self.vertices = vertices

	def __matchElementsToFacets(self):
		# Find all elements that share vertices with the boundary facets
		boundaryElements = []
		verticesIndexes  = [ vertex.handle for vertex in self.vertices ]
		for element, elementVertices in zip(self.grid.elements, self.grid.gridData.elementsConnectivities):
			if len(set(elementVertices).intersection(verticesIndexes)) >= self.grid.dimension:
				boundaryElements.append(element)

		for facet in self.facets:
			# Find the element that shares vertices with the facet
			facetElement = None
			for element in boundaryElements:
				if set(facet.vertices).issubset(element.vertices):
					facetElement = element
					break
			# Find the local index of the facet within the facetElement
			if facetElement != None:
				localFacetVertices = [vertex.getLocal(facetElement) for vertex in facet.vertices]
				for local in range(facetElement.shape.numberOfFacets):
					if set(localFacetVertices) == set(facetElement.shape.facetVerticesIndexes[local]):
						elementLocalIndex = local
						break

			facet._setElement(facetElement, elementLocalIndex)

	def __setFacets(self):
		for local, facet in enumerate(self.facets):
			facet._setBoundary(self, local)

			# Correct area
			innerVertex = [vertex for vertex in facet.element.vertices if vertex not in facet.vertices][0]
			dot = np.dot( facet.area.getCoordinates(), (facet.centroid - innerVertex).getCoordinates() )
			if dot < 0:
				facet.area *= -1

	def __setOuterFaces(self):
		for facet in self.facets:
			verticesLocalsInElement = list(facet.element.shape.facetVerticesIndexes[facet.elementLocalIndex])
			for outerFace in facet.outerFaces:
				outerFace.local = verticesLocalsInElement.index(outerFace.vertex.getLocal(facet.element))
				outerFace._build()

				facet.element._setOuterFace( outerFace )

				outerFace.handle = self.grid.outerFaceCounter
				self.grid.outerFaceCounter += 1