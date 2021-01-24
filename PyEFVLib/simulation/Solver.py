from PyEFVLib import ( 
	Grid, ProblemData,
	CgnsSaver, CsvSaver, VtuSaver, VtmSaver, MeshioSaver,
	MSHReader,
)
import os
import time


class Solver:
	def __init__(self, problemData, extension="csv", saverType="default", transient=True, verbosity=False, **kwargs):
		self.problemData = problemData
		self.extension = extension
		self.transient = transient
		self.verbosity = verbosity
		self.saverType = saverType
		
		self.grid = self.problemData.grid

	def solve(self):
		self.settings()			# DEFINED HERE
		self.printHeader()		# DEFINED HERE
		self.init()				# USER DEFINED
		self.mainloop()			# USER DEFINED
		self.finalizeSaver()	# DEFINED HERE
		self.printFooter()

	def settings(self):
		self.initialTime = time.time()

		self.propertyData = self.problemData.propertyData
		self.outputPath = self.problemData.outputFilePath
		savers = { "cgns": CgnsSaver, "csv": CsvSaver, "vtu": VtuSaver, "vtm": VtmSaver }

		if self.saverType == "default" and self.extension in savers.keys():
			self.saver = savers[self.extension](self.grid, self.outputPath, self.problemData.libraryPath, fileName="Results")
		else:
			self.saver = MeshioSaver(self.grid, self.outputPath, self.problemData.libraryPath, extension=self.extension, fileName="Results")

		self.numberOfVertices = self.grid.vertices.size
		self.dimension = self.grid.dimension
		self.currentTime = 0.0
		self.timeStep = self.problemData.timeStep
		self.tolerance = self.problemData.tolerance

		self.iteration = 0
		self.converged = False

	def printHeader(self):
		if self.verbosity:
			for key,path in zip( ["output", "grids"] , [self.problemData.outputFilePath, self.problemData.meshFilePath] ):
				print("\t{}\n\t\t{}\n".format(key, path))
			print("\tsolid")
			for region in self.grid.regions:
				print("\t\t{}".format(region.name))
				colSize = len(max(self.problemData.propertyData.properties, key=lambda w:len(w)))
				for propertyName in self.problemData.propertyData.properties:
					print(f"\t\t{propertyName:>{colSize}} : { self.problemData.propertyData.get(region.handle, propertyName) }")
				print("")
			print("\n{:>9}\t{:>14}\t{:>14}\t{:>14}".format("Iteration", "CurrentTime", "TimeStep", "Difference"))

	def printFooter(self):
		if self.verbosity:
			print("Ended Simultaion, elapsed {:.2f}s".format(time.time()-self.initialTime))
			print("\n\tresult: ", end="")
			print(os.path.realpath(self.saver.outputPath), "\n")

	def printIterationData(self):
		if self.verbosity:
			print("{:>9}\t{:>14.2e}\t{:>14.2e}\t{:>14.2e}".format(self.iteration, self.currentTime, self.timeStep, self.difference))

	def finalizeSaver(self):
		self.saver.finalize()

"""
class SomeSolver(Solver):
	def __init__(self, workspaceDirectory, **kwargs):
		Solver.__init__(self, workspaceDirectory, **kwargs)
	def init(self);
	def mainloop(self);
	def assembleMatrix(self);
	def addToIndependentVector(self);
	def solveLinearSystem(self);
	def saveIterationResults(self);
	def checkConvergence(self);
"""