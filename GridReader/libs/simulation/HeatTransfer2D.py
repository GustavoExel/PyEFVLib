import numpy as np
from libs.geometry.MSHReader import MSHReader
from libs.geometry.Grid import Grid
from libs.simulation.ProblemData2D import ProblemData2D
from libs.simulation.Timer import Timer
from libs.simulation.CgnsSaver import CgnsSaver
from libs.simulation.LinearSystemAdders import InternalGenerationAdder, TransientTermAdder, DiffusiveFluxAdder, NeumannBoundaryAdder, DirichletBoundaryAdder


class HeatTransfer2D:
	def __init__(self, path):
		self.path = path

		self.settings()
		self.run()

	def settings(self):
		reader = MSHReader("/home/gustavoe/Documents/Sinmec/GTRelated/GridReader/meshes/Square.msh")#QuadPlate.msh")
		self.grid = Grid(reader.getData())

		self.problemData = ProblemData2D(self.grid) # Make it contain property data
		self.timer = Timer(self.problemData.numericalSettings.timeStep)				# Contains dictionary with initial times for different labels: start("assembly"); stop("assembly")
		self.cgnsSaver = CgnsSaver()

		self.calculateInnerFacesGlobalDerivatives()
		self.numerical = np.zeros(self.grid.vertices.size)
		self.oldTemperature = np.zeros(self.grid.vertices.size)	# pass # problemData.getInitialTemperature() ->[Field, not scalar]

		self.matrix = np.zeros([self.grid.vertices.size, self.grid.vertices.size])
		self.difference = 0.0
		self.iteration = 0
		self.converged = False

	def run(self):
		while not self.converged:
			self.addToLinearSystem()
			self.solveLinearSystem()
			self.print()

			self.timer.incrementTime()
			CgnsSaver().save()
			self.converged = self.checkConvergence()

			self.iteration += 1    


	def calculateInnerFacesGlobalDerivatives(self):
		for element in self.grid.elements:
			for innerFace in element.innerFaces:
				derivatives = innerFace.element.shape.innerFaceShapeFunctionDerivatives[innerFace.local]
				innerFace.globalDerivatives = np.matmul(np.linalg.inv(element.getJacobian(derivatives)) , np.transpose(derivatives))
	# def addToLinearSystem(self):
	# 	numberOfVertices = self.grid.vertices.size
	# 	self.timer.start("assembly")
	# 	self.independent = np.zeros(numberOfVertices)
	# 	self.matrix = np.zeros([numberOfVertices,numberOfVertices])
	# 	InternalGenerationAdder().add() #Create another class later
	# 	TransientTermAdder().add()
	# 	DiffusiveFluxAdder().add()
	# 	NeumannBoundaryAdder().add()
	# 	## If error look here on matrices assemble
	# 	DirichletBoundaryAdder().add()
	# 	self.timer.stop("assemble")

	def addToLinearSystem(self):
		self.timer.start("assembly")
		self.independent = np.zeros(self.grid.vertices.size)

		heatGeneration = self.problemData.propertyData.internalHeatGeneration
		density = self.problemData.propertyData.density
		heatCapacity = self.problemData.propertyData.heatCapacity
		conductivity = self.problemData.propertyData.conductivity

		# Internal Heat Generation
		for region in self.grid.regions:
			for element in region.elements:
				for vertex in element.vertices:
					self.independent[vertex.handle] = vertex.volume * heatGeneration[region.handle]

		# TransientTermAdder and DiffusiveFluxAdder
		for region in self.grid.regions:
			accumulation = density[region.handle] * heatCapacity[region.handle] / self.	timer.timeStep
			for element in region.elements:
				localMatrixFlux = computeLocalMatrix(element, conductivity) # pass
				local = 0
				for vertex in element.vertices:
					index = vertex.handle
					self.independent[index] += element.subelementVolumes[local] * accumulation * self.oldTemperature[vertex.handle]
					if self.iteration == 0:
						self.matrix[index][index] += element.subelementVolumes[local] * accumulation
						qLocalIndex = 0
						for q in element.vertices:
							self.matrix[index][q.handle] = localMatrixFlux[local][qLocalIndex]
							qLocalIndex += 1
					local += 1

		# NeumannBoundaryAdder
		for bCondition in self.problemData.boundaryConditions.neumannBoundaries:
			for facet in bCondition.boundary.facets:
				for outerFace in facet.outerFaces:
					self.independent[outerFace.adjacentVertex.handle] += -1.0 * bCondition.getValue(outerFace.handle) * length(outerFace.area) # WHAT THE HECK IS THE LENGTH OF AN AREA???


		# DirichletBoundaryAdder
		for bCondition in self.problemData.boundaryConditions.dirichletBoundaries:
			for vertex in bCondition.boundary.vertices:
				self.independent[vertex.handle] = bCondition.getValue(vertex.handle)

		if iteration == 0:
			rows = np.zeros(self.problemData.boundaryConditions.dirichletBoundaries.numberOfVertices)
			for bCondition in self.problemData.boundaryConditions.dirichletBoundaries:
				for vertex in bCondition.boundary.vertices:
					self.matrix[vertex.handle] = np.zeros(self.grid.vertices.size)
					self.matrix[vertex.handle][vertex.handle] = 1.0

		self.timer.stop("assembly")


	def solveLinearSystem(self):
		self.timer.start("solve")
		# pass #self.numerical = np.linalg.solve(self.matrix, self.independent)
		self.timer.stop("solve")

	def checkConvergence(self):
		# Here oldTemperature becomes numerical (or right after this func is called)
		converged = False
		self.difference = max([abs(temp-oldTemp) for temp, oldTemp in zip(self.numerical, self.oldTemperature)])
		if self.timer.getCurrentTime() > self.problemData.numericalSettings.finalTime:
			converged = True
		elif self.iteration > 0:
			converged = self.difference < self.problemData.numericalSettings.tolerance

		return converged

	def print(self):
		if self.iteration == 0:
			print("{:>9}\t{:>14}\t{:>14}\t{:>14}".format("Iteration", "CurrentTime", "TimeStep", "Difference"))
		else:
			print("{:>9}\t{:>14e}\t{:>14e}\t{:>14e}".format(self.iteration, self.timer.getCurrentTime(), self.timer.timeStep, self.difference))


def computeLocalMatrix(element, permeability):
	numberOfVertices = element.vertices.size
	localMatrix = np.zeros([numberOfVertices, numberOfVertices])
	for innerFace in element.innerFaces:
		derivatives = innerFace.globalDerivatives
		if len(derivatives) == 2: derivatives = np.vstack([derivatives, np.zeros(derivatives[0].size)])

		diffusiveFlux = permeability[element.region.handle] * innerFace.area.getCoordinates() * np.diagonal(derivatives)

		backwardVertexLocalHandle = element.shape.innerFaceNeighbourVertices[innerFace.local][0]

		forwardVertexLocalHandle = element.shape.innerFaceNeighbourVertices[innerFace.local][1]

		for i in range(numberOfVertices):
			coefficient = -1.0 * diffusiveFlux[i]
			localMatrix[backwardVertexLocalHandle][i] += coefficient
			localMatrix[forwardVertexLocalHandle][i] -= coefficient

	return localMatrix

