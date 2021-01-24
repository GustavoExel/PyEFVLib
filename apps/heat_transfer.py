import sys,os
sys.path.append(os.path.join(os.path.dirname(__file__), os.path.pardir))
import PyEFVLib
from PyEFVLib import Solver
import numpy as np
from scipy import sparse
import scipy.sparse.linalg
import time

class HeatTransferSolver(Solver):
	def __init__(self, workspaceDirectory, **kwargs):
		# kwargs -> outputFileName, extension, transient, verbosity
		Solver.__init__(self, workspaceDirectory, **kwargs)

	def init(self):
		self.temperatureField = np.repeat(0.0, self.grid.vertices.size)
		self.prevTemperatureField = self.problemData.initialValues["temperature"]
		
		self.saver.save("temperature", self.prevTemperatureField, self.currentTime)

		self.coords,self.matrixVals = [], []
		self.difference=0.0

	def mainloop(self):
		self.assembleMatrix()
		while not self.converged:
			self.addToIndependentVector()
			self.solveLinearSystem()
			self.printIterationData()
			self.currentTime += self.timeStep
			self.saveIterationResults()
			self.checkConvergence()
			self.prevTemperatureField = self.temperatureField
			self.iteration += 1

	def add(self, i, j, val):
		self.coords.append((i,j))
		self.matrixVals.append(val)

	def assembleMatrix(self):
		def diffusionTerm():
			# Diffusion Term
			for region in self.grid.regions:
				conductivity = self.propertyData.get(region.handle, "Conductivity")
				for element in region.elements:
					for innerFace in element.innerFaces:
						diffusiveFlux = conductivity * np.matmul( np.transpose(innerFace.globalDerivatives) , innerFace.area.getCoordinates()[:self.dimension] )
						backwardVertexHandle = element.vertices[element.shape.innerFaceNeighborVertices[innerFace.local][0]].handle
						forwardVertexHandle = element.vertices[element.shape.innerFaceNeighborVertices[innerFace.local][1]].handle

						i=0
						for vertex in element.vertices:
							coefficient = -1.0 * diffusiveFlux[i]
							self.add(backwardVertexHandle, vertex.handle, coefficient)
							self.add(forwardVertexHandle, vertex.handle, -coefficient)
							i+=1

		def accumulationTerm():
			# Transient Term
			for region in self.grid.regions:
				density = self.propertyData.get(region.handle, "Density")
				heatCapacity = self.propertyData.get(region.handle, "HeatCapacity")
				accumulation = density * heatCapacity / self.timeStep

				for element in region.elements:
					local = 0
					for vertex in element.vertices:
						self.add(vertex.handle, vertex.handle, element.subelementVolumes[local] * accumulation)
						local += 1

		def dirichletBoundaryCondition():
			# Dirichlet Boundary Condition
			for bCondition in self.problemData.dirichletBoundaries["temperature"]:
				for vertex in bCondition.boundary.vertices:
					self.matrixVals, self.coords = zip( *[(val, coord) for coord, val in zip(self.coords, self.matrixVals) if coord[0] != vertex.handle] )
					self.matrixVals, self.coords = list(self.matrixVals), list(self.coords)
					self.add(vertex.handle, vertex.handle, 1.0)

		def invertMatrix():
			# Invert Matrix
			self.matrix = sparse.csc_matrix( (self.matrixVals, zip(*self.coords)), shape=(self.grid.vertices.size, self.grid.vertices.size) )
			self.inverseMatrix = sparse.linalg.inv( self.matrix )

		diffusionTerm()
		if self.transient:
			accumulationTerm()
		dirichletBoundaryCondition()
		invertMatrix()

	def addToIndependentVector(self):
		self.independent = np.zeros(self.grid.vertices.size)

		def generationTerm():
			# Generation Term
			for region in self.grid.regions:
				heatGeneration = self.propertyData.get(region.handle, "HeatGeneration")
				for element in region.elements:
					local = 0
					for vertex in element.vertices:
						self.independent[vertex.handle] += element.subelementVolumes[local] * heatGeneration
						local += 1

		def accumulationTerm():
			# Transient Term
			for region in self.grid.regions:
				density = self.propertyData.get(region.handle, "Density")
				heatCapacity = self.propertyData.get(region.handle, "HeatCapacity")
				accumulation = density * heatCapacity / self.timeStep

				for element in region.elements:
					local = 0
					for vertex in element.vertices:
						self.independent[vertex.handle] += element.subelementVolumes[local] * accumulation * self.prevTemperatureField[vertex.handle]
						local += 1

		def neumannBoundaryCondition():
			# Neumann Boundary Condition
			for bCondition in self.problemData.neumannBoundaries["temperature"]:
				for facet in bCondition.boundary.facets:
					for outerFace in facet.outerFaces:
						self.independent[outerFace.vertex.handle] += bCondition.getValue(outerFace.handle) * np.linalg.norm(outerFace.area.getCoordinates())

		def dirichletBoundaryCondition():
			# Dirichlet Boundary Condition
			for bCondition in self.problemData.dirichletBoundaries["temperature"]:
				for vertex in bCondition.boundary.vertices:
					self.independent[vertex.handle] = bCondition.getValue(vertex.handle)

		generationTerm()
		if self.transient:
			accumulationTerm()
		neumannBoundaryCondition()
		dirichletBoundaryCondition()

	def solveLinearSystem(self):
		self.temperatureField = np.matmul(self.inverseMatrix.toarray(), self.independent)

	def saveIterationResults(self):
		self.saver.save("temperature", self.temperatureField, self.currentTime)

	def checkConvergence(self):
		self.converged = False
		self.difference = max([abs(temp-oldTemp) for temp, oldTemp in zip(self.temperatureField, self.prevTemperatureField)])
		
		if self.problemData.finalTime != None and self.currentTime > self.problemData.finalTime:
			self.converged = True
			return
		if self.iteration > 0 and self.tolerance != None and self.difference < self.tolerance:
			self.converged = True
			return
		if self.problemData.maxNumberOfIterations and self.iteration >= self.problemData.maxNumberOfIterations:
			self.converged = True
			return

def heatTransfer(problemData, solve=True, extension="csv", saverType="default", transient=True, verbosity=True):
	solver = HeatTransferSolver(problemData, outputFileName="Results", extension=extension, saverType=saverType, transient=transient, verbosity=verbosity)
	if solve:
		solver.solve()
	return solver

if __name__ == "__main__":
	if "--help" in sys.argv: print("Usage: python apps/heat_transfer.py\n\t-p\t:Transient Solution\n\t-v\t:Verbosity\n\t--extension=EXT\t:Saves the results in EXT extension (msh, csv, xdmf, vtk, vtu, cgns (requires instalation), xmf, h5m, stl, obj, post, post.gz, dato, dato.gz, dat, fem, ugrid, wkt, ...)"); exit()
	extension = "xdmf" if not [1 for arg in sys.argv if "--extension" in arg] else [arg.split('=')[1] for arg in sys.argv if "--extension" in arg][0]
	saverType = "default" if not [1 for arg in sys.argv if "--saver" in arg] else [arg.split('=')[1] for arg in sys.argv if "--saver" in arg][0]

	problemData = PyEFVLib.ProblemData(
		meshFilePath = "{MESHES}/msh/2D/Square.msh",
		outputFilePath = "{RESULTS}/heat_transfer_2d/linear",
		numericalSettings = PyEFVLib.NumericalSettings( timeStep = 1e-02, finalTime = None, tolerance = 1e-06, maxNumberOfIterations = 300 ),
		propertyData = PyEFVLib.PropertyData({
			"Body" : {
				"HeatCapacity"	: 1.0,
				"Conductivity"	: 1.0,
				"Density"		: 1.0,
				"HeatGeneration": 0.0,
			},
		}),
		boundaryConditions = PyEFVLib.BoundaryConditions({
			"temperature": {
				"InitialValue": 0.0,
				"West":	 { "condition" : PyEFVLib.Dirichlet, "type" : PyEFVLib.Constant,"value" : 20.0 },
				"East":	 { "condition" : PyEFVLib.Dirichlet, "type" : PyEFVLib.Constant,"value" : 50.0 },
				"South": { "condition" : PyEFVLib.Neumann,   "type" : PyEFVLib.Constant,"value" : 0.0 },
				"North": { "condition" : PyEFVLib.Neumann,   "type" : PyEFVLib.Constant,"value" : 0.0 },
			},
		}),
	)

	heatTransfer(problemData, extension=extension, saverType=saverType, transient=not "-p" in sys.argv, verbosity="-v" in sys.argv)