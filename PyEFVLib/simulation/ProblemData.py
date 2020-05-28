from PyEFVLib.simulation.BoundaryConditions import DirichletBoundaryCondition, NeumannBoundaryCondition 
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
			raise Exception("Not enougth regions in PropertyData.json. Must contain {}".format(' '.join(gridRegionNames)) )

		self.propertyData = [ { _property : data[regionName][_property] for _property in data[regionName].keys()  } for regionName in gridRegionNames ]

class NumericalSettings:
	def __init__(self):
		self.readNumericalData()

	def readNumericalData(self):
		with open(self.paths["Numerical"], "r") as f:
			data = json.load(f)
		try:
			self.timeStep = data["TimeStep"]
			self.finalTime = data["FinalTime"]
		except:
			self.steady=True
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

		i = 0
		for boundary in self.grid.boundaries:
			if boundary.name not in self.boundaryConditionData.keys():
				raise Exception("There is no boundary entry {} in {}".format(boundary.name, self.path))

			bData = self.boundaryConditionData[boundary.name]
			if bData["condition"] == "DIRICHLET":
				self.dirichletBoundaries = np.append(self.dirichletBoundaries, DirichletBoundaryCondition(boundary, bData["value"], i) )

			elif bData["condition"] == "NEUMANN":
				self.neumannBoundaries = np.append(self.neumannBoundaries, NeumannBoundaryCondition(boundary, bData["value"], i) )

			else:
				raise Exception("Boundary condition not supported")
			i+=1

class ProblemData(PropertyData, NumericalSettings, BoundaryConditions):
	def __init__(self, simulatorName):
		self.simulatorName = simulatorName
		self.getPaths()

	def read(self):
		PropertyData.__init__(self)
		NumericalSettings.__init__(self)
		BoundaryConditions.__init__(self)

	def setGrid(self, grid):
		self.grid = grid
		self.read()

	def getPaths(self):
		self.libraryPath = '/'.join(os.path.realpath(__file__).split('/')[:-3])		# this is the GridReader path
		self.scriptPath  = "{}/workspace/{}/Script.json".format( self.libraryPath, self.simulatorName )
		if not os.path.isfile(self.scriptPath):
			raise(Exception("File {} not found".format(self.scriptPath)))

		with open(self.scriptPath, "r") as f:
			self.paths = json.load(f)
		self.paths = { key : self.libraryPath + '/' + self.paths[key] for key in self.paths}
