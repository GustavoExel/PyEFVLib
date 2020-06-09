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
		# with open(self.paths["Boundary"], "r") as f:
		# 	data = json.load(f)
		self.boundaryConditionData = dict()
		self.initialValue = dict()
		for path in list( os.walk(self.paths["Boundary"]) )[0][2]:
			with open( os.path.join(self.paths["Boundary"], path) , "r") as f:
				data = json.load(f)

			variableName = path.strip(".json") 
			self.boundaryConditionData[variableName] = {key : data[key] for key in data.keys() if key != "InitialValue"}
			self.initialValue[variableName] = data["InitialValue"]

	def build(self):
		self.dirichletBoundaries = dict()
		self.neumannBoundaries = dict()
		boundaryConditions = { boundary.name:dict() for boundary in self.grid.boundaries }
		for variableName in self.boundaryConditionData.keys():
			dirichletBoundaries = np.array([])
			neumannBoundaries = np.array([])

			i = 0
			for boundary in self.grid.boundaries:
				if boundary.name not in self.boundaryConditionData[variableName].keys():
					raise Exception("There is no boundary entry {} in {}".format(boundary.name, self.paths["Boundary"]))

				bData = self.boundaryConditionData[variableName][boundary.name]
				if bData["condition"] == "DIRICHLET":
					dirichletBoundaries = np.append(dirichletBoundaries, DirichletBoundaryCondition(boundary, bData["value"], i) )
					boundaryConditions[boundary.name][variableName] = dirichletBoundaries[-1]#{"boundary":dirichletBoundaries[-1], 'type':"DIRICHLET"}

				elif bData["condition"] == "NEUMANN":
					neumannBoundaries = np.append(neumannBoundaries, NeumannBoundaryCondition(boundary, bData["value"], i) )
					boundaryConditions[boundary.name][variableName] = neumannBoundaries[-1]#{"boundary":neumannBoundaries[-1], 'type':"NEUMANN"}

				else:
					raise Exception("Boundary condition not supported")
				i+=1

			self.dirichletBoundaries[variableName] = dirichletBoundaries
			self.neumannBoundaries[variableName] = neumannBoundaries
		self.boundaryConditions = list( boundaryConditions.values() )

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
