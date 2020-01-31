import json

# Future alteration, Problem data reads a Script.json file and passes the NumericalSettings, BoundaryConditions, PropertyData and mesh paths
class ProblemData2D:
	def __init__(self, grid):
		self.readPaths()
		self.grid = grid
		self.propertyData = PropertyData(grid, self.paths["Property"])
		self.numericalSettings = NumericalSettings(self.paths["Numerical"])
		self.boundaryConditions = BoundaryConditions(grid, self.paths["Boundary"])

	def readPaths(self):
		self.scriptPath = "/home/gustavoe/Documents/Sinmec/GTRelated/GridReader/config/Script.json"
		with open(self.scriptPath, "r") as f:
			self.paths = json.load(f)

# Remember to associate each region with its property, by name. If that is not possible, match by order
# SUPER ENROLATION!!!! VEJA O READDATA
class PropertyData:
	def __init__(self, grid, path):
		self.path = path
		self.grid = grid
		self.readData()

	def readData(self):
		self.density, self.heatCapacity, self.conductivity, self.internalHeatGeneration = [], [], [], []

		with open(self.path, "r") as f:
			data = json.load(f)
		jsonRegionNames = data.keys()
		gridRegionNames = [r.name for r in self.grid.regions]

		if len(jsonRegionNames) < len(gridRegionNames):
			raise Exception(f"Not enougth regions in PropertyData.json. There are {len(self.grid.regions)} regions in the mesh")

		if set(gridRegionNames).issubset(set(jsonRegionNames)):
			jsonRegionNames = gridRegionNames

		for name in jsonRegionNames:
			self.density.append(data[name]["Density"]["scalar"])
			self.heatCapacity.append(data[name]["HeatCapacity"]["scalar"])
			self.conductivity.append(data[name]["Conductivity"]["scalar"])
			self.internalHeatGeneration.append(data[name]["InternalHeatGeneration"]["scalar"])


class NumericalSettings:
	def __init__(self, path):
		self.path = path
		self.readData()

	def readData(self):
		with open(self.path, "r") as f:
			data = json.load(f)
		self.timeStep = data["TimeStep"]
		self.finalTime = data["FinalTime"]
		self.tolerance = data["Tolerance"]

# It must contain neumannBoundaries and dirichletBoundaries. Those have neumannBoudaries objects, which contains boundaries which contains (vertices and) facets which contains outerFacets
class BoundaryConditions:
	def __init__(self, grid, path):
		self.grid = grid
		self.path = path
		self.readData()

	def readData(self):
		with open(self.path, "r") as f:
			data = json.load(f)

		self.boundaryConditionData = {key : data[key] for key in data.keys() if key != "InitialValue"}
		self.initialValue = data["InitialValue"]

	def build(self):
		pass