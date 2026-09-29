from PyEFVLib.GridData import GridData
import meshio
import numpy
import xmltodict


class XDMFReader:
	def __init__(self, directory, filename, boundariesFilename=False, subdomainsFilename=False):
		self.directory = directory
		self.filename = filename
		self.boundariesFilename = boundariesFilename
		self.subdomainsFilename = subdomainsFilename
		self.path = "{}/{}".format(self.directory, self.filename)

	def readZoneList(self, filename):
		f = open("{}/{}".format(self.directory, filename), "r")
		self.zoneList = xmltodict.parse(f.read())

	def setFacetData(self, string):
		self.facetData = string

	def setSubdomainData(self, string):
		self.subdomainData = string

	def readFile(self, directory, filename):
		return meshio.read("{}/{}".format(directory, filename))

	def getDimension(self, mesh):
		dimensions = {"vertex": 0, "line": 1, "triangle": 2, "quad": 2,
			"tetra": 3, "hexahedron": 3, "wedge": 3, "pyramid": 3}
		return max(dimensions[block.type] for block in mesh.cells)

	def getVertices(self, mesh, gridDimension):
		verticesFileData = []
		for point in mesh.points:
			if gridDimension == 2:
				currPoint = (point[0], point[1], 0.0)
			elif gridDimension == 3:
				currPoint = (point[0], point[1], point[2])
			verticesFileData.append(currPoint)
		return verticesFileData

	def getEntityConnectivities(self, setobj):
		connectivities = []
		for cellBlock in setobj.cells:
			cellType = cellBlock.type
			cells = cellBlock.data
			for cell in cells:
				connectivities.append(cell.tolist())
		return connectivities

	def getRegionNames(self, setobj, setname):
		setdict = {}
		for i,index in enumerate(numpy.concatenate(setobj.cell_data[setname])):
			if index not in setdict.keys():
				setdict[index] = []
			setdict[index].append(i)
		regionsNames = []
		for key in setdict:
			regionsNames.append(self.zoneList['ZoneList']['Zone'+ str(key)]['@name'])
		return regionsNames

	def getRegionIndexes(self, setobj, setname):
		setdict = {}
		for i,index in enumerate(numpy.concatenate(setobj.cell_data[setname])):
			if index not in setdict.keys():
				setdict[index] = []
			setdict[index].append(i)
		regionIndexes = []
		for key in setdict:
			regionIndexes.append(setdict[key])
		return regionIndexes

	def getShapes(self):
		shapes = {}
		shapes['line'] = []
		shapes['triangle'] = []
		shapes['quadrilateral'] = []
		shapes['tetrahedron'] = []
		shapes['hexahedron'] = []
		shapes['prism'] = []
		shapes['pyramid'] = []
		return shapes

	def updateShapes(self, shapes, setobj, setname):
		for i,cellBlock in enumerate(setobj.cells):
			aliases = {"quad": "quadrilateral", "tetra": "tetrahedron", "hexahedron": "hexahedron", "wedge": "prism"}
			cellType = aliases.get(cellBlock.type, cellBlock.type)
			cells = cellBlock.data
			cellData = setobj.cell_data[setname][i]
			for j,data in enumerate(cellData):
				shapes[cellType].append((data, cells[j].tolist()))
		return shapes
	
	def read(self):
		self.mesh = self.readFile(self.directory, self.filename)
		if self.boundariesFilename:
			self.boundaries = self.readFile(self.directory, self.boundariesFilename)
		else:
			self.boundaries = self.mesh
		if self.subdomainsFilename:
			self.subdomains = self.readFile(self.directory, self.subdomainsFilename)
		else:
			self.subdomains = self.mesh
		self.gridDimension = self.getDimension(self.mesh)
		self.verticesFileData = self.getVertices(self.mesh, self.gridDimension)
		self.boundaryConnectivities = self.getEntityConnectivities(self.boundaries)
		self.boundaryNames = self.getRegionNames(self.boundaries, self.facetData)
		self.boundaryIndexes = self.getRegionIndexes(self.boundaries, self.facetData)
		self.elementsConnectivities = self.getEntityConnectivities(self.subdomains)
		self.subdomainNames = self.getRegionNames(self.subdomains, self.subdomainData)
		self.subdomainIndexes = self.getRegionIndexes(self.subdomains, self.subdomainData)
		self.shapes = self.getShapes()
		self.shapes = self.updateShapes(self.shapes, self.boundaries, self.facetData)
		self.shapes = self.updateShapes(self.shapes, self.subdomains, self.subdomainData)

	def getData(self):
		gridData = GridData(self.path)
		gridData.setDimension(self.gridDimension)
		gridData.setVertices(self.verticesFileData)
		gridData.setElementConnectivity(self.elementsConnectivities)
		gridData.setRegions(self.subdomainNames, self.subdomainIndexes)
		gridData.setBoundaries(self.boundaryNames, self.boundaryIndexes, self.boundaryConnectivities)
		gridData.setShapes([self.shapes['line'], self.shapes['triangle'], self.shapes['quadrilateral'], self.shapes['tetrahedron'], self.shapes['hexahedron'], self.shapes['prism'], self.shapes['pyramid']])
		return gridData
