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
		self.boundaryConditionData = dict()
		self.initialValues = dict()
		for path in list( os.walk(self.paths["Boundary"]) )[0][2]:
			with open( os.path.join(self.paths["Boundary"], path) , "r") as f:
				data = json.load(f)

			variableName = path.strip(".json") 
			self.boundaryConditionData[variableName] = {key : data[key] for key in data.keys() if key != "InitialValue"}
			self.initialValues[variableName] = [ data["InitialValue"] ] * self.grid.vertices.size

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
					dirichletBoundaries = np.append(dirichletBoundaries, DirichletBoundaryCondition(self.grid, boundary, bData["value"], i, expression=(bData["type"] == "VARIABLE") ) )
					boundaryConditions[boundary.name][variableName] = dirichletBoundaries[-1]

				elif bData["condition"] == "NEUMANN":
					neumannBoundaries = np.append(neumannBoundaries, NeumannBoundaryCondition(self.grid, boundary, bData["value"], i, expression=(bData["type"] == "VARIABLE") ) )
					boundaryConditions[boundary.name][variableName] = neumannBoundaries[-1]

				else:
					raise Exception("Boundary condition not supported")
				i+=1

			self.dirichletBoundaries[variableName] = dirichletBoundaries
			self.neumannBoundaries[variableName] = neumannBoundaries
		self.boundaryConditions = list( boundaryConditions.values() )

class ProblemData(PropertyData, NumericalSettings, BoundaryConditions):
	def __init__(self, simulatorName):
		self.simulatorName = os.path.join(*simulatorName.split("/"))
		self.getPaths()

	def read(self):
		PropertyData.__init__(self)
		NumericalSettings.__init__(self)
		BoundaryConditions.__init__(self)

	def setGrid(self, grid):
		self.grid = grid
		self.read()

	def getPaths(self):
		self.libraryPath = os.sep.join(os.path.realpath(__file__).split(os.sep)[:-3])		# this is the PyEFVLib path
		self.scriptPath  = os.path.join(self.libraryPath, self.simulatorName, "Script.json")
		if not os.path.isfile(self.scriptPath):
			raise(Exception("File {} not found".format(self.scriptPath)))

		with open(self.scriptPath, "r") as f:
			self.paths = json.load(f)

		self.replaceVariables()

	def replaceVariables(self):
		"""
		When a workspace folder is copyied into another, some informations must be changed all the time,
		because inside the Script.json the directory in which Script.json and the other files is explicit,
		giving more freedom. However, this is a pain every time a folder is copied. So to correct this problem
		this function will substitute every "keyword" or "variable", which will be indicated by ${var}.
		Also, for Windows Users, Script.json maintains the "/" separation syntax, and the class handles it;
		"""
		DIR = os.path.join(self.libraryPath, self.simulatorName)
		caseName = self.scriptPath.split("workspace")[1][1:].replace("/Script.json", "").replace("\\Script.json", "").replace("Script.json", "")
		variables = {"DIR" : DIR, "LIB" : self.libraryPath, "CASE": caseName}
		for key in self.paths.keys():
			
			sep = "/"
			if sep in self.paths[key]:
				self.paths[key] = os.path.join(*self.paths[key].split(sep))
			
			for variable in variables:
				var = "${%s}" % variable
				if var in self.paths[key]:
					self.paths[key] = self.paths[key].replace( var, variables[variable] )