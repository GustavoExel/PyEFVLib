import sys,os
sys.path.append(os.path.join(os.path.dirname(__file__), os.path.pardir))
from PyEFVLib import MSHReader, Grid, ProblemData, CgnsSaver, CsvSaver

class GUISettings:
	def __init__(self):
		pass

	def setFilePath(self, path):
		self.path = path
		self.readMesh()

	def readMesh(self):
		self.gridData = MSHReader(self.path).getData()

	def getBoundaryNames(self):
		return self.gridData.boundariesNames