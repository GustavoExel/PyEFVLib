from libs.geometry.MSHReader import MSHReader
from libs.geometry.Grid import Grid


class TestCase:
	def __init__(self, path = "/home/gustavoe/Documents/Sinmec/GTRelated/GridReader/meshes/test.msh"):
		mshReader = MSHReader(path)
		self.grid = Grid(mshReader.getData())

		self.runTests()
		self.printTests()

	def runTests(self):
		# Note that every test function must begin with 'test', therefore:
		self.results = []
		for attrName in dir(self):
			if "test" in attrName:
				getattr(self, attrName)()

	def printTests(self):
		print(f"Ran {len(self.results)} tests. {self.results.count(True)} concluded, {self.results.count(False)} failed")

	def testNodes(self):
		self.results.append( self.grid.vertices.size == 58 + 106 ) 			# Gmsh tell there are 60 nodes on lines, and 106 nodes on surfaces

	def testTriangles(self):
		self.results.append( len(self.grid.gridData.triangles) == 124 ) 	# Gmsh tell there are 124 triangles

	def testQuadrilaterals(self):
		self.results.append( len(self.grid.gridData.quadrilaterals) == 72 ) # Gmsh tell there are 72 quadrilaterals (quadrangles)

	def testRegions(self):
		self.results.append( self.grid.regions.size == 2 ) 					# Gmsh tell there are 2 regions (surfaces)

	def testBoundaries(self):
		self.results.append( self.grid.boundaries.size == 8 )				# Gmsh tell there are 8 boundaries (lines)

	def testVolumes(self):
		verticesVolumes = sum([v.volume for v in self.grid.vertices])
		elementsVolumes = sum([e.volume for e in self.grid.elements])
		self.results.append( abs(verticesVolumes - elementsVolumes) < 1e-14 )

if __name__ == "__main__":
	t = TestCase()