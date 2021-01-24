import PyEFVLib
from PyEFVLib.simulation.BoundaryConditions import DirichletBoundaryCondition, NeumannBoundaryCondition 
import json, os, numpy as np

class NumericalSettings:
	def __init__(self, timeStep, finalTime=None, tolerance=None, maxNumberOfIterations=1000):
		self.timeStep = timeStep
		self.finalTime = finalTime
		self.tolerance = tolerance
		self.maxNumberOfIterations = maxNumberOfIterations

	def init(self, problemData):
		self.problemData = problemData

		self.problemData.timeStep = self.timeStep
		self.problemData.finalTime = self.finalTime
		self.problemData.tolerance = self.tolerance
		self.problemData.maxNumberOfIterations = self.maxNumberOfIterations

class PropertyData:
	def __init__(self, propertiesDict):
		self.propertiesDict = propertiesDict
		self.properties = list(list(self.propertiesDict.values())[0].keys())

	def init(self, problemData):
		self.problemData = problemData
		self.regions = self.problemData.grid.gridData.regionsNames

		for regionName in self.propertiesDict.keys():
			if not regionName in self.problemData.grid.gridData.regionsNames:
				print(f"Warning! Region {regionName} not in the grid.")
		for regionName in self.problemData.grid.gridData.regionsNames:
			if not regionName in self.propertiesDict.keys():
				raise Exception(f"Error, must specify {regionName} properties.")

	def get(self, handle, propertyName):
		if handle in range(len(self.regions)):
			if propertyName in self.properties:
				return self.propertiesDict[self.regions[handle] ][propertyName]
			else:
				raise Exception(f"Property {propertyName} not defined")
		else:
			raise Exception(f"Invalid region handle {handle}")

	def set(self, handle, propertyName, value):
		if handle in range(len(self.regions)):
			if propertyName in self.properties:
				self.propertiesDict[self.regions[handle] ][propertyName] = value
			else:
				raise Exception(f"Property {propertyName} not defined")
		else:
			raise Exception(f"Invalid region handle {handle}")

class BoundaryConditions:
	def __init__(self, boundaryConditionsDict):
		self.boundaryConditionsDict = boundaryConditionsDict

	def init(self, problemData):
		self.problemData = problemData

		variables = self.boundaryConditionsDict.keys()
		boundaryNames = self.problemData.grid.gridData.boundariesNames

		for boundaryName in list(self.boundaryConditionsDict.values())[0].keys():
			if boundaryName != "InitialValue" and not boundaryName in boundaryNames:
				print(f"Warning! Boundary {boundaryName} not in the grid.")
		for boundaryName in boundaryNames:
			if not boundaryName in list(self.boundaryConditionsDict.values())[0].keys():
				raise Exception(f"Error, must specify {boundaryName} condition.")

		self.dirichletBoundaries = dict()
		self.neumannBoundaries   = dict()
		self.boundaryConditions  = dict()
		self.initialValues		 = dict()

		# bcHandle = 0
		# for variable in variables:
		# 	self.neumannBoundaries[variable] = list()
		# 	self.dirichletBoundaries[variable] = list()
		# 	self.boundaryConditions[variable] = list()
		# 	self.initialValues[variable] = np.repeat( self.boundaryConditionsDict[variable]["InitialValue"], self.problemData.grid.vertices.size )

		# 	for boundary in self.problemData.grid.boundaries:
		# 		if self.boundaryConditionsDict[variable][boundary.name]["condition"] == PyEFVLib.Neumann:
		# 			bcValue = self.boundaryConditionsDict[variable][boundary.name]["value"]
		# 			bc = NeumannBoundaryCondition(self.problemData.grid, boundary, bcValue, bcHandle)
		# 			self.neumannBoundaries[variable].append(bc)
		# 			self.boundaryConditions[variable].append(bc)
					
		# 		elif self.boundaryConditionsDict[variable][boundary.name]["condition"] == PyEFVLib.Dirichlet:
		# 			bcValue = self.boundaryConditionsDict[variable][boundary.name]["value"]
		# 			bc = DirichletBoundaryCondition(self.problemData.grid, boundary, bcValue, bcHandle)
		# 			self.dirichletBoundaries[variable].append(bc)
		# 			self.boundaryConditions[variable].append(bc)
					
		# 		else:
		# 			raise Exception("Invalid boundary condition!")
		# 		bcHandle += 1

		self.neumannBoundaries   = { variable : [] for variable in variables }
		self.dirichletBoundaries = { variable : [] for variable in variables }
		self.boundaryConditions  = [ dict() for boundary in self.problemData.grid.boundaries ]
		self.initialValues		 = { variable : np.repeat( self.boundaryConditionsDict[variable]["InitialValue"], self.problemData.grid.vertices.size ) for variable in variables }
		
		bcHandle = 0
		for idx, boundary in enumerate(self.problemData.grid.boundaries):
			for variable in variables:
				if self.boundaryConditionsDict[variable][boundary.name]["condition"] == PyEFVLib.Neumann:
					bcValue = self.boundaryConditionsDict[variable][boundary.name]["value"]
					bc = NeumannBoundaryCondition(self.problemData.grid, boundary, bcValue, bcHandle)
					self.neumannBoundaries[variable].append(bc)
					self.boundaryConditions[idx][variable] = bc
					
				elif self.boundaryConditionsDict[variable][boundary.name]["condition"] == PyEFVLib.Dirichlet:
					bcValue = self.boundaryConditionsDict[variable][boundary.name]["value"]
					bc = DirichletBoundaryCondition(self.problemData.grid, boundary, bcValue, bcHandle)
					self.dirichletBoundaries[variable].append(bc)
					self.boundaryConditions[idx][variable] = bc
					
				else:
					raise Exception("Invalid boundary condition!")
				bcHandle += 1


		self.problemData.dirichletBoundaries = self.dirichletBoundaries
		self.problemData.neumannBoundaries   = self.neumannBoundaries
		self.problemData.boundaryConditions  = self.boundaryConditions
		self.problemData.initialValues  	 = self.initialValues

class ProblemData:
	def __init__(self, 
		meshFilePath,
		propertyData,
		boundaryConditions,
		numericalSettings=NumericalSettings(0.0),
		outputFilePath=None,
	):
		self.meshFilePath		= meshFilePath
		self.outputFilePath		= outputFilePath
		self.propertyData 		= propertyData

		self.parsePaths()
		self.createGrid()
		numericalSettings.init(self)
		propertyData.init(self)
		boundaryConditions.init(self)

	def parsePaths(self):
		self.libraryPath = os.path.realpath( os.path.join(os.path.dirname(__file__), *(2*[os.path.pardir])) )
		meshesPath = os.path.realpath( os.path.join(os.path.dirname(__file__), *(2*[os.path.pardir]), "meshes") )
		resultsPath = os.path.realpath( os.path.join(os.path.dirname(__file__), *(2*[os.path.pardir]), "results") )

		self.meshFilePath = os.path.join(*self.meshFilePath.replace("{MESHES}", meshesPath).split("/"))
		self.outputFilePath = os.path.join(*self.outputFilePath.replace("{RESULTS}",resultsPath).split("/"))

	def createGrid(self):
		self.grid = PyEFVLib.read(self.meshFilePath)

if __name__ == "__main__":
	problemData =ProblemData(
		meshFilePath  = "{MESHES}/msh/2D/Square.msh",		# required
		outputFilePath = "{RESULTS}/Results.xdmf",	# optional

		numericalSettings = NumericalSettings(
			timeStep = 1e-02,							# required
			finalTime = None,							# optional
			tolerance = 1e-06,							# optional
			maxNumberOfIterations = 1000, 				# optional, but with default
		),

		propertyData = PropertyData({
			"Body" : {
				"HeatCapacity"	: 1.0,
				"Conductivity"	: 1.0,
				"Density"		: 1.0,
			},
		}),

		boundaryConditions = BoundaryConditions({
			"temperature": {
				"InitialValue": 0.0,
				"West": {
					"condition" : PyEFVLib.Neumann,
					"type"		: "Constant", # The ideal use is PyEFVLib.Constant
					"value"		: 100.0
				},
			    "East": {
			        "condition"	: PyEFVLib.Dirichlet,
			        "type"		: "Constant", # The ideal use is PyEFVLib.Constant
			        "value"		: 0.0
			    },
			    "South": {
			        "condition"	: PyEFVLib.Neumann,
			        "type"		: "Constant", # The ideal use is PyEFVLib.Constant
			        "value"		: 0.0
			    },
			    "North": {
			        "condition"	: PyEFVLib.Neumann,
			        "type"		: "Constant", # The ideal use is PyEFVLib.Constant
			        "value"		: 0.0
			    },
			},
		}),
	)