from libs.simulation.BoundaryConditions import DirichletBoundaryCondition, NeumannBoundaryCondition 
import json, os, numpy as np

#pass #Property Data -> Reclamar se não tiver o nome da Region 
# Future alteration, Problem data reads a Script.json file and passes the NumericalSettings, BoundaryConditions, PropertyData and mesh paths
class ProblemData2D:
	def __init__(self, grid):
		self.grid = grid
		self.libPath = '/'.join(os.path.realpath(__file__).split('/')[:-3])
		self.readPaths()
		self.propertyData = PropertyData(grid, self.paths["Property"])
		self.numericalSettings = NumericalSettings(self.paths["Numerical"])
		self.boundaryConditions = BoundaryConditions(grid, self.paths["Boundary"])

	def readPaths(self):
		self.scriptPath = "/home/gustavoe/Documents/Sinmec/GTRelated/GridReader/config/Script.json"
		with open(self.scriptPath, "r") as f:
			self.paths = json.load(f)
		self.paths = { key : self.libPath + '/' + self.paths[key] for key in self.paths}

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
import os
class BoundaryConditions:
	def __init__(self, grid, path):
		self.grid = grid
		self.path = path
		self.readData()
		self.build()

	def readData(self):
		with open(self.path, "r") as f:
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
				dirichlet = DirichletBoundaryCondition(boundary, bData["value"])
				self.dirichletBoundaries = np.append(self.dirichletBoundaries, dirichlet)

			elif bData["condition"] == "NEUMANN":
				neumann = NeumannBoundaryCondition(boundary, bData["value"])
				self.neumannBoundaries = np.append(self.neumannBoundaries, neumann)

			else:
				raise Exception("Boundary condition not supported")
