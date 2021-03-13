import numpy as np
import subprocess, os, sys
from PyEFVLib.simulation.Saver import Saver

class CgnsSaver(Saver):
	def __init__(self, grid, outputPath, basePath, fileName="Results", **kwargs):
		self.binPath = os.path.join( basePath, "PyEFVLib", "simulation", "CGNS", "bin" )
		Saver.__init__(self, grid, outputPath, basePath, 'cgns', fileName)

	def finalize(self):
		with open(os.path.join(self.binPath , "data.txt"), "a") as f:
			t = '\n'.join(self.fields.keys()) + "\n";
			f.write(t)

		for fieldName in self.fields.keys():
			with open(os.path.join(self.binPath, fieldName) + ".txt", "w") as f:
				f.write( '\n'.join([' '.join([str(x) for x in field]) for field in self.fields[fieldName]]) )

		with open(os.path.join(self.binPath, "steps.txt"), "w") as f:
			f.write( ' '.join([ str(ts) for ts in self.timeSteps ]) )

		subprocess.run([os.path.join(self.binPath, "save")])
		os.remove(os.path.join(self.binPath, "data.txt"))
		os.remove(os.path.join(self.binPath, "steps.txt"))
		for fieldName in self.fields.keys():
			os.remove(os.path.join(self.binPath, fieldName) + ".txt")

		self.finalized = True

	def createFile(self):
		# Check existing Results.cgns, if so remove it
		if not os.path.isdir(os.path.dirname(self.outputPath)):
			os.makedirs(os.path.dirname(self.outputPath))

		if os.path.isfile(self.outputPath):
			os.remove(self.outputPath)
		self.export()

		# Create Results.cgns
		subprocess.run([os.path.join(self.binPath, "create")])
		os.remove(os.path.join(self.binPath, "coords.txt"))
		os.remove(os.path.join(self.binPath, "connectivity.txt"))
		os.remove(os.path.join(self.binPath, "sections.txt"))

	def export(self):
		with open(os.path.join(self.binPath, "data.txt"), "w") as f:
			t = self.outputPath +"\n"
			t += str(self.grid.numberOfVertices) + " " + str(self.grid.elements.size) + " 0\n"
			f.write(t)

		with open(os.path.join(self.binPath, "coords.txt"), "w") as f:
			f.write( '\n'.join( [ ' '.join([str(c) for c in coord]) for coord in zip(*self.grid.gridData.vertices) ] ) )

		# Later create a way of guaranteeing that regionsElementsIndexes is a range (in MSHReader)
		with open(os.path.join(self.binPath, "connectivity.txt"), "w") as f:
			t =  '\n'.join( [ ' '.join([str(e) for e in elem]) for elem in self.grid.gridData.elementsConnectivities ] )
			t += '\n'
			t += '\n'.join( [ ' '.join([str(e) for e in elem]) for elem in self.grid.gridData.boundariesConnectivities ] ) 
			f.write(t)

		with open(os.path.join(self.binPath, "sections.txt"), "w") as f:
			t=""
			for regionName, regionsElementsIndexes in zip( self.grid.gridData.regionsNames, self.grid.gridData.regionsElementsIndexes ):
				t += regionName + ' ' + str( len(self.grid.gridData.elementsConnectivities[ regionsElementsIndexes[0] ]) ) + ' ' + str(regionsElementsIndexes[0]) + ' ' + str(regionsElementsIndexes[-1]) + '\n'
			for boundaryName, boundariesIndexes in zip( self.grid.gridData.boundariesNames, self.grid.gridData.boundariesIndexes ):
				t += boundaryName + ' ' + str( len(self.grid.gridData.boundariesConnectivities[ boundariesIndexes[0] ]) ) + ' ' + str(len(self.grid.gridData.elementsConnectivities) + boundariesIndexes[0]) + ' ' + str(len(self.grid.gridData.elementsConnectivities) + boundariesIndexes[-1]) + "\n"
			f.write(t)