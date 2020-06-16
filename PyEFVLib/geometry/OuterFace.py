import numpy as np

class OuterFace:
	def __init__(self, vertex, facet, local, handle):
		self.vertex = vertex
		self.facet = facet
		self.local = local
		self.handle = handle
		self.vertexLocalIndex = vertex.getLocal( facet.element )

