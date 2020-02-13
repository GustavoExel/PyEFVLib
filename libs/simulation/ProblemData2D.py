from libs.simulation.BoundaryConditions import DirichletBoundaryCondition, NeumannBoundaryCondition 
import json, os, numpy as np

class PropertyData:
	def __init__(self):
		self.readPropertyData()

	def readPropertyData(self):
		with open(self.paths["Property"], "r") as f:
			data = json.load(f)
		jsonRegionNames = data.keys()
		gridRegionNames = [r.name for r in self.grid.regions]

		if not set(gridRegionNames).issubset(set(jsonRegionNames)):
			raise Exception(f"Not enougth regions in PropertyData.json. Must contain {' '.join(gridRegionNames)}")

		self.propertyData = [ { _property : data[regionName][_property]["scalar"] for _property in data[regionName].keys()  } for regionName in gridRegionNames ]

class NumericalSettings:
	def __init__(self):
		self.readNumericalData()

	def readNumericalData(self):
		with open(self.paths["Numerical"], "r") as f:
			data = json.load(f)
		self.timeStep = data["TimeStep"]
		self.finalTime = data["FinalTime"]
		self.tolerance = data["Tolerance"]
		self.maxNumberOfIterations = data["MaximumNumberOfIterations"]

class BoundaryConditions:
	def __init__(self):
		self.readBoundaryConditionsData()
		self.build()

	def readBoundaryConditionsData(self):
		with open(self.paths["Boundary"], "r") as f:
			data = json.load(f)

		self.boundaryConditionData = {key : data[key] for key in data.keys() if key != "InitialValue"}
		self.initialValue = data["InitialValue"]

	def build(self):
		self.dirichletBoundaries = np.array([])
		self.neumannBoundaries = np.array([])

		for boundary in self.grid.boundaries:
			if boundary.name not in self.boundaryConditionData.keys():
				raise Exception(f"There is no boundary entry {boundary.name} in {self.path}") 

			bData = self.boundaryConditionData[boundary.name]
			if bData["condition"] == "DIRICHLET":
				self.dirichletBoundaries = np.append(self.dirichletBoundaries, DirichletBoundaryCondition(boundary, bData["value"]) )

			elif bData["condition"] == "NEUMANN":
				self.neumannBoundaries = np.append(self.neumannBoundaries, NeumannBoundaryCondition(boundary, bData["value"]) )

			else:
				raise Exception("Boundary condition not supported")


class ProblemData2D(PropertyData, NumericalSettings, BoundaryConditions):
	def __init__(self, simulatorName):
		self.simulatorName = simulatorName
		self.getPaths()

	def init(self):
		PropertyData.__init__(self)
		NumericalSettings.__init__(self)
		BoundaryConditions.__init__(self)

	def setGrid(self, grid):
		self.grid = grid
		self.init()

	def getPaths(self):
		self.libraryPath = '/'.join(os.path.realpath(__file__).split('/')[:-3])		# this is the GridReader path
		self.scriptPath  = f"{ self.libraryPath }/benchmark/{ self.simulatorName }/Script.json"

		with open(self.scriptPath, "r") as f:
			self.paths = json.load(f)
		self.paths = { key : self.libraryPath + '/' + self.paths[key] for key in self.paths}
