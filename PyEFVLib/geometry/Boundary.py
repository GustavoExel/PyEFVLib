from PyEFVLib.geometry.Facet import Facet
from PyEFVLib.geometry.Point import Point
from PyEFVLib.geometry.OuterFace import OuterFace
import numpy as np

class Boundary:
	def __init__(self, grid, facetsIdexes, name, handle):
		self.grid = grid
		self.name = name
		self.handle = handle

		self.facets = np.array([grid.facets[facetIndex] for facetIndex in facetsIdexes])

		self.setVertices()
		self.matchElementsToFacets()
		self.setFacets()
		self.setOuterFaces()

	def setVertices(self):
		vertices = []
		for facet in self.facets:
			for vertex in facet.vertices:
				if vertex not in vertices:
					vertices.append(vertex)
		self.vertices = np.array(vertices)

	def matchElementsToFacets(self):
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

			facet.setElement(facetElement, elementLocalIndex)

	def setFacets(self):
		for local, facet in enumerate(self.facets):
			facet.setBoundary(self, local)

			# Correct area
			innerVertex = [vertex for vertex in facet.element.vertices if vertex not in facet.vertices][0]
			dot = np.dot( facet.area.getCoordinates(), (facet.centroid - innerVertex).getCoordinates() )
			if dot < 0:
				facet.area *= -1

	def setOuterFaces(self):
		for facet in self.facets:
			verticesLocalsInElement = list(facet.element.shape.facetVerticesIndexes[facet.elementLocalIndex])
			for outerFace in facet.outerFaces:
				outerFace.setLocal( verticesLocalsInElement.index(outerFace.vertex.getLocal(facet.element)) )
				
				outerFace.computeCentroid()
				outerFace.computeAreaVector()

				facet.element.setOuterFace( outerFace )

				outerFace.handle = self.grid.outerFaceCounter
				self.grid.outerFaceCounter += 1