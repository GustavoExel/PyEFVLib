import numpy as np

class InternalGenerationAdder:
	def __init__(self, simulator):
		self.simulator = simulator

	def add(self):
		for region in self.simulator.grid.regions:
			heatGeneration = self.simulator.problemData.propertyData[region.handle]["InternalHeatGeneration"]
			for element in region.elements:
				for vertex in element.vertices:
					self.simulator.independent[vertex.handle] = vertex.volume * heatGeneration

class HeatDiffusionAdder:
	def __init__(self, simulator):
		self.simulator = simulator

	def add(self):
		if self.simulator.iteration == 0:
			for region in self.simulator.grid.regions:
				conductivity = self.simulator.problemData.propertyData[region.handle]["Conductivity"]
				for element in region.elements:
					localMatrixFlux = computeLocalMatrix(element, conductivity,self.simulator.iteration==0)
					local = 0
					for vertex in element.vertices:
						qLocalIndex = 0
						for q in element.vertices:
							self.simulator.matrix[vertex.handle][q.handle] += localMatrixFlux[local][qLocalIndex]
							qLocalIndex += 1
						local += 1

class AccumulationAdder:
	def __init__(self, simulator):
		self.simulator = simulator

	def add(self):
		pass
		for region in self.simulator.grid.regions:
			density = self.simulator.problemData.propertyData[region.handle]["Density"]
			heatCapacity = self.simulator.problemData.propertyData[region.handle]["HeatCapacity"]
			accumulation = density * heatCapacity / self.simulator.timer.timeStep

			for element in region.elements:
				local = 0
				for vertex in element.vertices:
					index = vertex.handle
					self.simulator.independent[index] += element.subelementVolumes[local] * accumulation * self.simulator.oldTemperature[vertex.handle]
					if self.simulator.iteration == 0:
						self.simulator.matrix[index][index] += element.subelementVolumes[local] * accumulation
					local += 1



class NeumannBoundaryAdder:
	def __init__(self, simulator):
		self.simulator = simulator

	def add(self):
		for bCondition in self.simulator.problemData.neumannBoundaries:
			for facet in bCondition.boundary.facets:
				for outerFace in facet.outerFaces:
					self.simulator.independent[outerFace.vertex.handle] += -1.0 * bCondition.getValue(outerFace.handle) * np.linalg.norm(outerFace.area.getCoordinates())

class DirichletBoundaryAdder:
	def __init__(self, simulator):
		self.simulator = simulator

	def add(self):
		for bCondition in self.simulator.problemData.dirichletBoundaries:
			for vertex in bCondition.boundary.vertices:
				self.simulator.independent[vertex.handle] = bCondition.getValue(vertex.handle)

		if self.simulator.iteration == 0:
			for bCondition in self.simulator.problemData.dirichletBoundaries:
				for vertex in bCondition.boundary.vertices:
					self.simulator.matrix[vertex.handle] = np.zeros(self.simulator.grid.vertices.size)
					self.simulator.matrix[vertex.handle][vertex.handle] = 1.0



def computeLocalMatrix(element, permeability,b=False):
	numberOfVertices = element.vertices.size
	localMatrix = np.zeros([numberOfVertices, numberOfVertices])
	for innerFace in element.innerFaces:
		derivatives = innerFace.globalDerivatives
		if len(derivatives) == 2: derivatives = np.vstack([derivatives, np.zeros(derivatives[0].size)])
		diffusiveFlux = permeability * np.matmul( np.transpose(derivatives[:-1]) , innerFace.area.getCoordinates()[:-1] )

		backwardVertexLocalHandle = element.shape.innerFaceNeighbourVertices[innerFace.local][0]
		forwardVertexLocalHandle = element.shape.innerFaceNeighbourVertices[innerFace.local][1]

		for i in range(numberOfVertices):
			coefficient = -1.0 * diffusiveFlux[i]
			localMatrix[backwardVertexLocalHandle][i] += coefficient
			localMatrix[forwardVertexLocalHandle][i] -= coefficient
	return localMatrix
