import numpy as np
import subprocess, os, sys

class Saver:
	def __init__(self, grid, outputPath, basePath, extension, fileName="Results"): 
		self.grid = grid
		self.fileName = fileName
		self.extension = extension
		self.outputPath = os.path.join( outputPath , f"{fileName}.{extension}" )
		self.basePath = basePath
		self.timeSteps  = np.array([])
		self.fields = dict()
		self.createFile()
		self.finalized = False

	def save(self, fieldName, fieldValues, currentTime):
		if not self.timeSteps.size or self.timeSteps[-1] != currentTime:
			self.timeSteps = np.append( self.timeSteps, currentTime )
		if not fieldName in self.fields.keys():
			self.fields[fieldName] = np.zeros((0, self.grid.numberOfVertices))

		self.fields[fieldName] = np.vstack([self.fields[fieldName], fieldValues])

	def __del__(self):
		if not self.finalize:
			self.finalize()

	def finalize(self):
		self.finalized = True

	def createFile(self):
		pass