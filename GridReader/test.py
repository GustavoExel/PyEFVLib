from libs.geometry.MSHReader import MSHReader
from libs.geometry.Grid import Grid


class TestCase:
	def __init__(self, grid, reader):
		self.grid = grid
		self.reader = reader
		self.results = dict()
		self.runTests()
		self.printResults()


	def runTests(self):
		for attrName in dir(self):
			if "test" in attrName:
				getattr(self, attrName)()

	def printResults(self):
		positive = len([key for key in self.results.keys() if self.results[key]])
		negative = len([key for key in self.results.keys() if not self.results[key]])
		print(f"Ran {positive+negative} tests\n{negative} faliure(s) found!")

	def testNumberOfNodes(self):
		lines = self.reader.lines
		numberOfNodesIndex = lines.index("$Nodes\n") + 1
		numberOfNodesRead = int(lines[numberOfNodesIndex])
		numberOfNodesGet = self.grid.vertices.size

		if numberOfNodesGet == numberOfNodesRead:
			self.results["testNumberOfNodes"] = 1
		else:
			self.results["testNumberOfNodes"] = 0

	def testNumberOfElements(self):
		lines = self.reader.lines
		numberOfElementsIndex = lines.index("$Elements\n") + 1
		numberOfElementsRead = int(lines[numberOfElementsIndex])
		numberOfElementsGet = self.grid.elements.size

		if numberOfElementsGet == numberOfElementsRead:
			self.results["testNumberOfElements"] = 1
		else:
			self.results["testNumberOfElements"] = 0

	def testNumberOfRegions(self):
		lines = self.reader.lines
		numberOfRegionNamesIndex = lines.index("$PhysicalNames\n") + 1
		numberOfRegionNames = int(lines[numberOfRegionNamesIndex])
		numberOfRegions = self.grid.regions.size

		if numberOfRegions == numberOfRegionNames:
			self.results["testNumberOfRegions"] = 1
		else:
			self.results["testNumberOfRegions"] = 0

	def testVolumes(self):
		vertexVolumes = sum([vertex.volume for vertex in self.grid.vertices])
		elementVolumes = sum([element.volume for element in self.grid.elements])

		if vertexVolumes - elementVolumes < 1e-15:	#	There was a difference of 2e-16
			self.results["testVolumes"] = 1
		else:
			self.results["testVolumes"] = 0


# if __name__ == "__main__":
# 	r = MSHReader("/home/gustavoe/Documents/Sinmec/GTRelated/GridReader/meshes/Square.msh")#QuadPlate.msh")
# 	g = Grid(r.getData())
# 	t = TestCase(g,r)


from libs.simulation.ProblemData2D import ProblemData2D
if __name__ == "__main__":
	r = MSHReader("/home/gustavoe/Documents/Sinmec/GTRelated/GridReader/meshes/Square.msh")#QuadPlate.msh")
	g = Grid(r.getData())
	pD = ProblemData2D(g)
