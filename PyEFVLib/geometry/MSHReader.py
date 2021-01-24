from PyEFVLib.geometry.GridData import GridData

# Msh format encodes elements in the following manner:
# i x x p e v1 v2 v3 - i[index], xx[shape code], p[physical index], e[elementary index], vi[vertex]

class MSHReader:
	def __init__(self, path):
		self.path = path

		self.open()
		self.checkFileVersion()
		self.read()

	def open(self):
		with open(self.path, 'r') as f:
			fileLines = f.read().split('$')
		self.fileData = [ [ line.split() for line in fileLines[i].split('\n')[1:-1] ] for i in range(1,8,2) ]

	def checkFileVersion(self):
		fileVersion = float(self.fileData[0][0][0])
		if fileVersion < 2.0 or fileVersion > 2.2:
			raise( Exception("MSH version must be 2.2") )

	def read(self):
		# This numberOf... are indicated on top of each msh file section
		self.numberOfSections = int( self.fileData[1][0][0] )
		self.numberOfVertices = int( self.fileData[2][0][0] )
		self.numberOfConnectivities = int( self.fileData[3][0][0] )

		# This fileDatas store raw information from file, only converting from string
		# DIMENSION, INDEX, NAME => REGIONS / BOUNDARIES
		self.sectionsFileData = [ ( int(dimension), int(index) , name[1:-1] ) for dimension, index, name in self.fileData[1][1:] ]

		# X,Y,Z COORDINATES
		self.verticesFileData = [ (float(x), float(y), float(z)) for idx,x,y,z in self.fileData[2][1:] ]

		# SHAPE CODE, SECTION INDEX, [[ VERTICES INDEXES ]]
		self.connectivitiesFileData = [ ( ''.join([c1,c2]), int(p_id), [ int(v)-1 for v in nodes ] ) for idx, c1, c2, p_id, e_id, *nodes in self.fileData[3][1:] ]

		# Sorts elements by sections
		self.sectionsElements = [ [ e[2] for e in self.connectivitiesFileData if e[1] == section[1] ] for section in self.sectionsFileData]

		shapeCodes = {"line":"12", "triangle":"22", "quadrilateral":"32", "tetrahedron":"42", "pyramid":"72", "prism":"62", "hexagon":"52"}
		self.shapes = { shape : [ idx for idx, e in enumerate(self.connectivitiesFileData) if e[0] == shapeCodes[shape] ] for shape in shapeCodes.keys() }

		self.gridDimension = max(self.sectionsFileData, key=lambda p:p[0])[0]

	def getData(self):
		elementsConnectivities, regionsNames, regionsElementsIndexes = [], [], []
		boundariesConnectivities, boundariesNames, boundariesIndexes = [], [], []

		maxDimension = max([dimension for dimension, index, name in self.sectionsFileData])
		for [dimension, index, name], sectionElements in zip( self.sectionsFileData, self.sectionsElements ):
			if dimension == maxDimension:
				indexOfFirstConnectivity = len(elementsConnectivities)

				elementsConnectivities += sectionElements
				regionsNames.append(name)
				regionsElementsIndexes.append( list(range(indexOfFirstConnectivity, indexOfFirstConnectivity+len(sectionElements))) )

			else:
				indexOfFirstConnectivity = len(boundariesConnectivities)

				boundariesConnectivities += sectionElements
				boundariesNames.append(name)
				boundariesIndexes.append( list(range(indexOfFirstConnectivity, indexOfFirstConnectivity+len(sectionElements))) )

		gridData = GridData(self.path)
		gridData.setDimension(self.gridDimension)
		gridData.setVertices(self.verticesFileData)
		gridData.setElementConnectivity(elementsConnectivities)
		gridData.setRegions(regionsNames, regionsElementsIndexes)
		gridData.setBoundaries(boundariesNames, boundariesIndexes, boundariesConnectivities)
		gridData.setShapes([ self.shapes['line'], self.shapes['triangle'], self.shapes['quadrilateral'], self.shapes['tetrahedron'], [], [], [] ])

		return gridData
