import numpy as np
from PyEFVLib.geometry.Point import Point

def areCoplanar(p1,p2,p3,p4):
	return np.dot( (p2-p1), np.cross((p3-p1), (p4-p1)) ) == 0

class Triangle:
	dimension						   = 2
	numberOfInnerFaces				   = 3
	numberOfFacets					   = 3
	subelementTransformedVolumes	   = np.array([1.0 / 6.0, 1.0 / 6.0, 1.0 / 6.0])
	innerFaceShapeFunctionValues 	   = np.array([[5.0/12.0, 5.0/12.0, 1.0/6.0],[1.0/6.0, 5.0/12.0, 5.0/12.0],[5.0/12.0, 1.0/6.0, 5.0/12.0]])
	innerFaceShapeFunctionDerivatives  = np.array([[[-1.0/1.0, -1.0/1.0], [1.0/1.0, 0.0/1.0], [0.0/1.0, 1.0/1.0]],[[-1.0/1.0, -1.0/1.0], [1.0/1.0, 0.0/1.0], [0.0/1.0, 1.0/1.0]],[[-1.0/1.0, -1.0/1.0], [1.0/1.0, 0.0/1.0], [0.0/1.0, 1.0/1.0]]])
	innerFaceNeighborVertices		   = np.array([[0, 1],[1, 2],[2, 0]])
	subelementShapeFunctionValues 	   = np.array([[[-1.0/1.0, -1.0/1.0], [1.0/1.0, 0.0/1.0], [0.0/1.0, 1.0/1.0]],[[-1.0/1.0, -1.0/1.0], [1.0/1.0, 0.0/1.0], [0.0/1.0, 1.0/1.0]],[[-1.0/1.0, -1.0/1.0], [1.0/1.0, 0.0/1.0], [0.0/1.0, 1.0/1.0]]])
	subelementShapeFunctionDerivatives = np.array([[[-1.0/1.0, -1.0/1.0], [1.0/1.0, 0.0/1.0], [0.0/1.0, 1.0/1.0]],[[-1.0/1.0, -1.0/1.0], [1.0/1.0, 0.0/1.0], [0.0/1.0, 1.0/1.0]],[[-1.0/1.0, -1.0/1.0], [1.0/1.0, 0.0/1.0], [0.0/1.0, 1.0/1.0]]])
	facetVerticesIndices 			   = np.array([[1, 0],[2, 1],[0, 2]])
	outerFaceShapeFunctionValues 	   = np.array([[[1.0/4.0, 3.0/4.0, 0.0/1.0],[3.0/4.0, 1.0/4.0, 0.0/1.0]],[[0.0/1.0, 1.0/4.0, 3.0/4.0],[0.0/1.0, 3.0/4.0, 1.0/4.0]],[[3.0/4.0, 0.0/1.0, 1.0/4.0],[1.0/4.0, 0.0/1.0, 3.0/4.0]]])
	vertexShapeFunctionDerivatives	   = np.array([[[-1.0,-1.0],[1.0,0.0],[0.0,1.0]],[[-1.0,-1.0],[1.0,0.0],[0.0,1.0]],[[-1.0,-1.0],[1.0,0.0],[0.0,1.0]]])

	def __init__(self, element):
		self.element = element

	@staticmethod
	def _is(elem):
		if len(elem.vertices) == 3:
			return True
		else:
			return False

	@staticmethod
	def getInnerFaceAreaVector(local, elementCentroid, elementVertices):
		vertex1 = elementVertices[Triangle.innerFaceNeighborVertices[local][0]]
		vertex2 = elementVertices[Triangle.innerFaceNeighborVertices[local][1]]
		areaVectorCoords = sum([w*vertex.getCoordinates() for w,vertex in zip([1,-0.5,-0.5],[elementCentroid, vertex1, vertex2])])
		return Point(areaVectorCoords[1], -areaVectorCoords[0], 0.0)

class Quadrilateral:
	dimension						   = 2
	numberOfInnerFaces				   = 4
	numberOfFacets 					   = 4
	subelementTransformedVolumes	   = np.array([1.0/4.0, 1.0/4.0, 1.0/4.0, 1.0/4.0])
	innerFaceShapeFunctionValues	   = np.array([[3.0/8.0, 3.0/8.0, 1.0/8.0, 1.0/8.0], [1.0/8.0, 3.0/8.0, 3.0/8.0, 1.0/8.0], [1.0/8.0, 1.0/8.0, 3.0/8.0, 3.0/8.0], [3.0/8.0, 1.0/8.0, 1.0/8.0, 3.0/8.0]])
	innerFaceShapeFunctionDerivatives  = np.array([[[-3.0/4.0, -1.0/2.0], [3.0/4.0, -1.0/2.0], [1.0/4.0, 1.0/2.0], [-1.0/4.0, 1.0/2.0]],[[-1.0/2.0, -1.0/4.0], [1.0/2.0, -3.0/4.0], [1.0/2.0, 3.0/4.0], [-1.0/2.0, 1.0/4.0]],[[-1.0/4.0, -1.0/2.0], [1.0/4.0, -1.0/2.0], [3.0/4.0, 1.0/2.0], [-3.0/4.0, 1.0/2.0]],[[-1.0/2.0, -3.0/4.0], [1.0/2.0, -1.0/4.0], [1.0/2.0, 1.0/4.0], [-1.0/2.0, 3.0/4.0]]])
	innerFaceNeighborVertices		   = np.array([[0, 1],[1, 2],[2, 3],[3, 0]])
	subelementShapeFunctionValues	   = np.array([[9.0/16.0, 3.0/16.0, 1.0/16.0, 3.0/16.0],[3.0/16.0, 9.0/16.0, 3.0/16.0, 1.0/16.0],[1.0/16.0, 3.0/16.0, 9.0/16.0, 3.0/16.0],[3.0/16.0, 1.0/16.0, 3.0/16.0, 9.0/16.0]])
	subelementShapeFunctionDerivatives = np.array([[[-3.0/4.0, -3.0/4.0], [3.0/4.0, -1.0/4.0], [1.0/4.0, 1.0/4.0], [-1.0/4.0, 3.0/4.0]],[[-3.0/4.0, -1.0/4.0], [3.0/4.0, -3.0/4.0], [1.0/4.0, 3.0/4.0], [-1.0/4.0, 1.0/4.0]],[[-1.0/4.0, -1.0/4.0], [1.0/4.0, -3.0/4.0], [3.0/4.0, 3.0/4.0], [-3.0/4.0, 1.0/4.0]],[[-1.0/4.0, -3.0/4.0], [1.0/4.0, -1.0/4.0], [3.0/4.0, 1.0/4.0], [-3.0/4.0, 3.0/4.0]]])
	facetVerticesIndices 			   = np.array([[1, 0],[2, 1],[3, 2],[0, 3]])
	outerFaceShapeFunctionValues	   = np.array([[[1.0/4.0, 3.0/4.0, 0.0/1.0, 0.0/1.0],[3.0/4.0, 1.0/4.0, 0.0/1.0, 0.0/1.0]],[[0.0/1.0, 1.0/4.0, 3.0/4.0, 0.0/1.0],[0.0/1.0, 3.0/4.0, 1.0/4.0, 0.0/1.0]],[[0.0/1.0, 0.0/1.0, 1.0/4.0, 3.0/4.0],[0.0/1.0, 0.0/1.0, 3.0/4.0, 1.0/4.0]],[[3.0/4.0, 0.0/1.0, 0.0/1.0, 1.0/4.0],[1.0/4.0, 0.0/1.0, 0.0/1.0, 3.0/4.0]]])
	vertexShapeFunctionDerivatives	   = np.array([[[-1.0,-1.0],[1.0,0.0],[0.0,0.0],[0.0,1.0]],[[-1.0,0.0],[1.0,0.0],[0.0,1.0],[0.0,0.0]],[[0.0,0.0],[0.0,-1.0],[1.0,1.0],[-1.0,0.0]],[[0.0,-1.0],[0.0,0.0],[1.0,0.0],[0.0,1.0]]])

	def __init__(self, element):
		self.element = element			

	@staticmethod
	def _is(elem):
		if len(elem.vertices) == 4 and areCoplanar(*[v.getCoordinates() for v in elem.vertices]):
			return True
		else:
			return False

	@staticmethod
	def getInnerFaceAreaVector(local, elementCentroid, elementVertices):
		vertex1 = elementVertices[Quadrilateral.innerFaceNeighborVertices[local][0]]
		vertex2 = elementVertices[Quadrilateral.innerFaceNeighborVertices[local][1]]
		areaVectorCoords = sum([w*vertex.getCoordinates() for w,vertex in zip([1,-0.5,-0.5],[elementCentroid, vertex1, vertex2])])
		return Point(*areaVectorCoords)

class Tetrahedron:
	dimension						   = 3
	numberOfInnerFaces				   = 6
	numberOfFacets 					   = 4
	subelementTransformedVolumes	   = np.array([1.0/24.0, 1.0/24.0, 1.0/24.0, 1.0/24.0])
	innerFaceShapeFunctionValues	   = np.array([[17.0/48.0, 17.0/48.0, 7.0/48.0, 7.0/48.0],[7.0/48.0, 17.0/48.0, 17.0/48.0, 7.0/48.0],[17.0/48.0, 7.0/48.0, 17.0/48.0, 7.0/48.0],[17.0/48.0, 7.0/48.0, 7.0/48.0, 17.0/48.0],[7.0/48.0, 7.0/48.0, 17.0/48.0, 17.0/48.0],[7.0/48.0, 17.0/48.0, 7.0/48.0, 17.0/48.0]])
	innerFaceShapeFunctionDerivatives  = np.array([[[-1.0/1.0, -1.0/1.0, -1.0/1.0], [1.0/1.0, 0.0/1.0, 0.0/1.0], [0.0/1.0, 1.0/1.0, 0.0/1.0], [0.0/1.0, 0.0/1.0, 1.0/1.0]],[[-1.0/1.0, -1.0/1.0, -1.0/1.0], [1.0/1.0, 0.0/1.0, 0.0/1.0], [0.0/1.0, 1.0/1.0, 0.0/1.0], [0.0/1.0, 0.0/1.0, 1.0/1.0]],[[-1.0/1.0, -1.0/1.0, -1.0/1.0], [1.0/1.0, 0.0/1.0, 0.0/1.0], [0.0/1.0, 1.0/1.0, 0.0/1.0], [0.0/1.0, 0.0/1.0, 1.0/1.0]],[[-1.0/1.0, -1.0/1.0, -1.0/1.0], [1.0/1.0, 0.0/1.0, 0.0/1.0], [0.0/1.0, 1.0/1.0, 0.0/1.0], [0.0/1.0, 0.0/1.0, 1.0/1.0]],[[-1.0/1.0, -1.0/1.0, -1.0/1.0], [1.0/1.0, 0.0/1.0, 0.0/1.0], [0.0/1.0, 1.0/1.0, 0.0/1.0], [0.0/1.0, 0.0/1.0, 1.0/1.0]],[[-1.0/1.0, -1.0/1.0, -1.0/1.0], [1.0/1.0, 0.0/1.0, 0.0/1.0], [0.0/1.0, 1.0/1.0, 0.0/1.0], [0.0/1.0, 0.0/1.0, 1.0/1.0]]])
	innerFaceNeighborVertices		   = np.array([[0, 1, 3, 2],[1, 2, 3, 0],[2, 0, 3, 1],[0, 3, 2, 1],[1, 3, 0, 2],[2, 3, 1, 0]])
	subelementShapeFunctionValues	   = np.array([[15.0/32.0, 17.0/96.0, 17.0/96.0, 17.0/96.0],[17.0/96.0, 15.0/32.0, 17.0/96.0, 17.0/96.0],[17.0/96.0, 17.0/96.0, 15.0/32.0, 17.0/96.0],[17.0/96.0, 17.0/96.0, 17.0/96.0, 15.0/32.0]])
	subelementShapeFunctionDerivatives = np.array([[[-1.0/1.0, -1.0/1.0, -1.0/1.0], [1.0/1.0, 0.0/1.0, 0.0/1.0], [0.0/1.0, 1.0/1.0, 0.0/1.0], [0.0/1.0, 0.0/1.0, 1.0/1.0]],[[-1.0/1.0, -1.0/1.0, -1.0/1.0], [1.0/1.0, 0.0/1.0, 0.0/1.0], [0.0/1.0, 1.0/1.0, 0.0/1.0], [0.0/1.0, 0.0/1.0, 1.0/1.0]],[[-1.0/1.0, -1.0/1.0, -1.0/1.0], [1.0/1.0, 0.0/1.0, 0.0/1.0], [0.0/1.0, 1.0/1.0, 0.0/1.0], [0.0/1.0, 0.0/1.0, 1.0/1.0]],[[-1.0/1.0, -1.0/1.0, -1.0/1.0], [1.0/1.0, 0.0/1.0, 0.0/1.0], [0.0/1.0, 1.0/1.0, 0.0/1.0], [0.0/1.0, 0.0/1.0, 1.0/1.0]]])
	facetVerticesIndices 			   = np.array([[0, 2, 1],[0, 3, 2],[0, 1, 3],[1, 2, 3]])
	outerFaceShapeFunctionValues	   = np.array([[[7.0/12.0, 5.0/24.0, 5.0/24.0, 0.0/1.0],[5.0/24.0, 5.0/24.0, 7.0/12.0, 0.0/1.0],[5.0/24.0, 7.0/12.0, 5.0/24.0, 0.0/1.0]],[[7.0/12.0, 0.0/1.0, 5.0/24.0, 5.0/24.0],[5.0/24.0, 0.0/1.0, 5.0/24.0, 7.0/12.0],[5.0/24.0, 0.0/1.0, 7.0/12.0, 5.0/24.0]],[[7.0/12.0, 5.0/24.0, 0.0/1.0, 5.0/24.0],[5.0/24.0, 7.0/12.0, 0.0/1.0, 5.0/24.0],[5.0/24.0, 5.0/24.0, 0.0/1.0, 7.0/12.0]],[[0.0/1.0, 7.0/12.0, 5.0/24.0, 5.0/24.0],[0.0/1.0, 5.0/24.0, 7.0/12.0, 5.0/24.0],[0.0/1.0, 5.0/24.0, 5.0/24.0, 7.0/12.0]]])
	# vertexShapeFunctionDerivatives	   = np.array()

	def __init__(self, element):
		self.element = element			

	@staticmethod
	def _is(elem):
		if len(elem.vertices) == 4 and not areCoplanar(*[v.getCoordinates() for v in elem.vertices]):
			return True
		else:
			return False

	def getInnerFaceAreaVector(self, local, elementCentroid, elementVertices):
		f,b,n1,n2 = [ elementVertices[index].getCoordinates() for index in self.innerFaceNeighborVertices[local] ]
		fc1=(f+b+n1)/3
		fc2=(f+b+n2)/3
		mp=(f+b)/2
		ec=elementCentroid.getCoordinates()
		t1=np.cross(n1-mp, n2-mp)/2
		t2=np.cross(n1-ec, n2-ec)/2
		p=Point(*(t1+t2))
		return p