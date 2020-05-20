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
			raise( Exception("MSH version must be ") )

	def read(self):
		self.numberOfSections = int( self.fileData[1][0][0] )
		self.numberOfVertices = int( self.fileData[2][0][0] )
		self.numberOfConnectivities = int( self.fileData[3][0][0] )

		self.sectionsFileData = [ ( int(dimension), int(index) , name[1:-1] ) for dimension, index, name in self.fileData[1][1:] ]
		self.verticesFileData = [ (float(x), float(y), float(z)) for idx,x,y,z in self.fileData[2][1:] ]
		self.connectivitiesFileData = [ ( ''.join([c1,c2]), int(p_id), [ int(v)-1 for v in nodes ] ) for idx, c1, c2, p_id, e_id, *nodes in self.fileData[3][1:] ]

		self.sectionElements = [ [ e[2] for e in self.connectivitiesFileData if e[1] == section[1] ] for section in self.sectionsFileData]

		shapeCodes = {"line":"12", "triangle":"22", "quadrilateral":"32", "tetrahedron":"42"}
		self.shapes = { shape : [ (e[1],e[2]) for e in self.connectivitiesFileData if e[0] == shapeCodes[shape] ] for shape in shapeCodes.keys() }


	def getData(self):
		elementsConnectivities, regionNames, regionsElementsIndexes = [], [], []
		boundariesConnectivities, boundaryNames, boundariesIndexes = [], [], []

		maxDimension = max([dimension for dimension, index, name in self.sectionsFileData])
		for [dimension, index, name], sectionElements in zip( self.sectionsFileData, self.sectionElements ):
			if dimension == maxDimension:
				indexOfFirstConnectivity = len(elementsConnectivities)
				
				elementsConnectivities += sectionElements
				regionNames.append(name)
				regionsElementsIndexes.append( list(range(indexOfFirstConnectivity, indexOfFirstConnectivity+len(sectionElements))) )

			else:
				indexOfFirstConnectivity = len(boundariesConnectivities)
				
				boundariesConnectivities += sectionElements
				boundaryNames.append(name)
				boundariesIndexes.append( list(range(indexOfFirstConnectivity, indexOfFirstConnectivity+len(sectionElements))) )

		gridData = GridData(self.path)
		gridData.setVertices(self.verticesFileData)
		gridData.setElementConnectivity(elementsConnectivities)
		gridData.setRegions(regionNames, regionsElementsIndexes)
		gridData.setBoundaries(boundaryNames, boundariesIndexes, boundariesConnectivities)
		gridData.setShapes([ self.shapes['line'], self.shapes['triangle'], self.shapes['quadrilateral'], self.shapes['tetrahedron'], [], [], [] ])

		return gridData