from libs.geometry.GridData import GridData

# Msh format encodes elements in the following manner:
# i x x p e v1 v2 v3 - i[index], xx[shape code], p[physical index], e[elementary index], vi[vertex]
#

class MSHReader:
	def __init__(self, path):
		self.path = path
		self.read()


	def init(self):
		self.nodes 			 				= []
		self.elements 			 			= []
		self.regionNames 					= dict()
		self.regionElements 	 			= dict()
		self.boundaryElements 			 	= []
		self.boundaryNames 					= dict()
		self.boundaries						= dict()
		self.lines							= []
		self.triangles 						= []
		self.quadrilaterals 				= []
		self.tetrahedrons 					= []
		self.hexahedrons					= []
		self.prisms							= []
		self.pyramids						= []

		self.shapeCodes = {'1 2' : 'line', '2 2' : 'triangle', '3 2' : 'quadrilateral', '4 2' : 'tetrahedron'}
		self.dimensionDict = {'line' : 1, 'triangle' : 2, 'quadrilateral' : 2, 'tetrahedron' : 3}

		with open(self.path, "r") as f:
			self.fileLines = f.readlines()
		self.fields = {'physicalNames' : self.fileLines.index("$PhysicalNames\n"), 'nodes' : self.fileLines.index("$Nodes\n"), 'elements' : self.fileLines.index("$Elements\n")}

	def defineDimension(self):
		self.dimension = max([ int(line.split()[0]) for line in self.fileLines[ self.fields['physicalNames']+2 : self.fields['nodes']-1 ] ])

	def read(self):
		self.init()
		self.defineDimension()
		self.readPhysicalEntities()
		self.readNodes()
		self.readElements()

	def readPhysicalEntities(self):
		for line in self.fileLines[ self.fields['physicalNames']+2 : self.fields['nodes']-1 ]:
			index = int(line.split()[1])
			if int(line.split()[0]) == self.dimension:
				self.regionNames[index] = line.split()[-1][1:-1]
				self.regionElements[index] = []
			else:
				self.boundaryNames[index] =  line.split()[-1][1:-1]
				self.boundaries[index] = []

	def readNodes(self):
		for line in self.fileLines[ self.fields['nodes']+2 : self.fields['elements']-1 ]:
			idx,x,y,z = line.split()
			self.nodes.append([float(x), float(y), float(z)])

	def readElements(self):
		for line in self.fileLines[ self.fields['elements']+2 : -1 ]:
			code = ' '.join(line.split()[1:3])

			if self.dimensionDict[ self.shapeCodes[code] ] == self.dimension:
				# This is an element of the grid. If grid is 3D then this element is 3D and the same for 2D.
				# It will be appended to the elements list of connectivity, and its index of elements list will be stored in region elements in the right region key.
				self.elements.append([int(v)-1 for v in line.split()[5:]])

				regionId = int(line.split()[3])
				elementId = len(self.elements)-1
				self.regionElements[regionId].append(elementId)

			else:
				self.boundaryElements.append([int(v)-1 for v in line.split()[5:]])
				boundaryId = int(line.split()[3])
				elementId = len(self.boundaryElements)-1
				self.boundaries[boundaryId].append(elementId)

			if self.dimension == 2:
				if self.shapeCodes[code] == 'triangle':
					self.triangles.append( len(self.elements)-1 )

				if self.shapeCodes[code] == 'quadrilateral':
					self.quadrilaterals.append( len(self.elements)-1 )

			elif self.dimension == 3:
				if self.shapeCodes[code] == 'triangle':
					self.triangles.append( len(self.boundaryElements)-1 )

				if self.shapeCodes[code] == 'quadrilateral':
					self.quadrilaterals.append( len(self.boundaryElements)-1 )

				if self.shapeCodes[code] == 'tetrahedron':
					self.tetrahedrons.append( len(self.elements)-1 )

				if self.shapeCodes[code] == 'hexahedron':
					self.hexahedrons.append( len(self.elements)-1 )

				if self.shapeCodes[code] == 'prism':
					self.prisms.append( len(self.elements)-1 )

				if self.shapeCodes[code] == 'pyramid':
					self.pyramids.append( len(self.elements)-1 )

	def getData(self):
		regionNames = [self.regionNames[key] for key in self.regionNames]
		regionElements = [self.regionElements[key] for key in self.regionElements]

		boundaryNames = [self.boundaryNames[key] for key in self.boundaryNames]
		boundaryElements = [self.boundaries[key] for key in self.boundaries]

		gridData = GridData()
		gridData.setVertices(self.nodes)
		gridData.setElementConnectivity(self.elements)
		gridData.setRegions(regionNames, regionElements)
		gridData.setBoundaries(boundaryNames, boundaryElements, self.boundaryElements)
		gridData.setShapes([self.lines, self.triangles, self.quadrilaterals, self.tetrahedrons, self.hexahedrons, self.prisms, self.pyramids])

		return gridData