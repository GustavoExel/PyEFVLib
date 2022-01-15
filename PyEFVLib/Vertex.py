import numpy as np
from PyEFVLib.Point import Point

class Vertex(Point):
	"""
	The vertex is a really imporant entity in the EbFVM, since it defines
	each control volume. In it, are stored the elements to which it belongs,
	as well as the control volume's volume.

	Attributes
	----------
	handle		: int
		The vertex index in the grid
	elements	: list[ PyEFVLib.Element ]
		A list of all elements to which vertex belongs
	volume		: float
		The volume of the control volume at which vertex is centered

	Methods
	----------
	getLocal(element : PyEFVLib.Elment) -> int:
		Returns the local index in element.vertices
	getInnerFaces() -> list[ PyEFVLib.InnerFace ]:
		Returns the innerFaces that make up its control surface
	getSubElementVolume(element : PyEFVLib.Elment) -> float:
		Return the volume of the subelement of element that contains vertex

	"""
	def __init__(self, coordinates, handle):
		Point.__init__(self, *coordinates)
		self.handle = handle
		self.elements = []
		self.volume = 0.0

	def _addElement(self, element):
		self.elements.append(element)

	def getLocal(self, element) -> int:
		"""
		Returns the local index in element.vertices
		"""
		return list(element.vertices).index(self)

	def getInnerFaces(self) -> list:
		"""
		Returns the innerFaces that make up its control surface
		"""
		return [ innerFace for element in self.elements for innerFace in element.innerFaces if self in innerFace.getNeighborVertices() ]

	def getSubElementVolume(self, element) -> int:
		"""
		Returns the volume of the subelement of element that contains vertex
		"""
		return element.subelementVolumes[ self.getLocal(element) ]