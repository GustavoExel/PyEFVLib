import numpy as np
from PyEFVLib.geometry.MSHReader import MSHReader
from PyEFVLib.geometry.Grid import Grid
from PyEFVLib.simulation.ProblemData import ProblemData
from PyEFVLib.simulation.Timer import Timer
from PyEFVLib.simulation.CgnsSaver import CgnsSaver
from PyEFVLib.simulation.LinearSystemAdders import InternalGenerationAdder, HeatDiffusionAdder, AccumulationAdder, NeumannBoundaryAdder, DirichletBoundaryAdder
import pandas as pd

class HeatTransfer2D:
	def __init__(self):
		self.settings()
		self.run()
		self.finalize()

	def settings(self):
		self.problemData = ProblemData('heat_transfer_2d')
		
		reader = MSHReader(self.problemData.paths["Grid"])
		self.grid = Grid(reader.getData())
		self.problemData.setGrid(self.grid)

		self.timer = Timer(self.problemData.timeStep)				# Contains dictionary with initial times for different labels: start("assembly"); stop("assembly")
		self.cgnsSaver = CgnsSaver(self.grid, self.problemData.paths["Output"], self.problemData.libraryPath)

		self.internalGenerationAdder = InternalGenerationAdder(self)
		self.heatDiffusionAdder 	 = HeatDiffusionAdder(self)
		self.accumulationAdder 		 = AccumulationAdder(self)
		self.neumannBoundaryAdder 	 = NeumannBoundaryAdder(self)
		self.dirichletBoundaryAdder  = DirichletBoundaryAdder(self)

		self.numericalTemperature = np.zeros(self.grid.vertices.size)
		self.oldTemperature = np.repeat(self.problemData.initialValue, self.grid.vertices.size)

		self.matrix = np.zeros([self.grid.vertices.size, self.grid.vertices.size])
		self.difference = 0.0
		self.iteration = 0
		self.converged = False

	def run(self):
		while not self.converged and self.iteration < self.problemData.maxNumberOfIterations:
			self.addToLinearSystem()
			self.solveLinearSystem()
			self.print()

			self.timer.incrementTime()
			self.cgnsSaver.save(self.numericalTemperature, self.timer.getCurrentTime())
			self.converged = self.checkConvergence()

			self.iteration += 1   

	def addToLinearSystem(self):
		self.timer.start("assemble")
		self.independent = np.zeros(self.grid.vertices.size)

		self.internalGenerationAdder.add()
		self.heatDiffusionAdder.add()
		self.accumulationAdder.add()
		self.neumannBoundaryAdder.add()
		self.dirichletBoundaryAdder.add()

		self.timer.stop("assemble")

	def solveLinearSystem(self):
		self.timer.start("solve")
		self.numericalTemperature = np.linalg.solve(self.matrix, self.independent)
		self.timer.stop("solve")

	def checkConvergence(self):
		# Here oldTemperature becomes numerical (or right after this func is called)
		converged = False
		self.difference = max([abs(temp-oldTemp) for temp, oldTemp in zip(self.numericalTemperature, self.oldTemperature)])
		self.oldTemperature = self.numericalTemperature
		if self.timer.getCurrentTime() > self.problemData.finalTime:
			converged = True
		elif self.iteration > 0:
			converged = self.difference < self.problemData.tolerance

		return converged

	def startInfo(self):
		for key,path in zip( ["input", "output", "grids"] , [self.problemData.libraryPath+"/benchmark/heat_transfer_2d/" , self.problemData.paths["Output"], self.problemData.paths["Grid"]] ):
			print(f"\t\033[1;35m{key}\033[0m\n\t\t{path}\n")

		print(f"\t\033[1;35msolid\033[0m")
		for region in self.grid.regions:
			print(f"\t\t\033[36m{region.name}\033[0m")#//in blue
			for _property in self.problemData.propertyData[region.handle].keys():
				print(f"\t\t\t{_property}   : {self.problemData.propertyData[region.handle][_property]}")#//16 digits | scientific

			print("")
		print("")

	def print(self):
		if self.iteration == 0:
			self.startInfo()
			print("{:>9}\t{:>14}\t{:>14}\t{:>14}".format("Iteration", "CurrentTime", "TimeStep", "Difference"))
		else:
			print("{:>9}\t{:>14e}\t{:>14e}\t{:>14e}".format(self.iteration, self.timer.getCurrentTime(), self.timer.timeStep, self.difference))

	def finalize(self):
		self.cgnsSaver.finalize()
		print("")
		total = 0.0
		for timeLabel in self.timer.timeLabels.keys():
			total += self.timer.timeLabels[timeLabel]["elapsedTime"]
			print("\t{:<12}:{:>12.5f}s".format(timeLabel, self.timer.timeLabels[timeLabel]["elapsedTime"]))
		print("\t{:<12}:{:>12.5f}s".format("total", total))

		print("\n\t\033[1;35mresult:\033[0m", self.problemData.paths["Output"]+"Results.cgns", '\n')