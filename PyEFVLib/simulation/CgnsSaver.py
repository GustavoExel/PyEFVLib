import numpy as np
import subprocess, os, sys

class CgnsSaver:
	def __init__(self, grid, outputPath, basePath, fieldsNames): 
		self.grid = grid
		self.outputPath = outputPath + "Results.cgns"
		self.basePath = basePath
		self.binPath = os.path.join( basePath, "PyEFVLib", "simulation" ) + "/"
		self.fieldsNames = fieldsNames

		self.timeSteps  = np.array([])
		self.fields = { fieldName : np.zeros((0, grid.vertices.size)) for fieldName in fieldsNames }
		self.createCgns()
		self.finalized = False


	def save(self, fieldName, fieldValues, currentTime):
		if not self.timeSteps.size or self.timeSteps[-1] != currentTime:
			self.timeSteps = np.append( self.timeSteps, currentTime )
		self.fields[fieldName] = np.vstack([self.fields[fieldName], fieldValues])

	def __del__(self):
		if not self.finalize:
			self.finalize()
	def finalize(self):
		for fieldName in self.fieldsNames:
			with open(self.binPath + f"{fieldName}.txt", "w") as f:
				f.write( '\n'.join([' '.join([str(x) for x in field]) for field in self.fields[fieldName]]) )

		with open(self.binPath + "steps.txt", "w") as f:
			f.write( ' '.join([ str(ts) for ts in self.timeSteps ]) )

		subprocess.run([self.binPath + "save"])
		os.remove(self.binPath + "data.txt")
		os.remove(self.binPath + "steps.txt")
		for fieldName in self.fieldsNames:
			os.remove(self.binPath + f"{fieldName}.txt")

		self.finalize = True

	def createCgns(self):
		# Check existing Results.cgns, if so remove it
		if os.path.isfile(self.outputPath):
			os.remove(self.outputPath)
		self.export()

		# Create Results.cgns
		subprocess.run([self.binPath + "create"])
		os.remove(self.binPath + "coords.txt")
		os.remove(self.binPath + "connectivity.txt")
		os.remove(self.binPath + "sections.txt")

	def export(self):
		with open(self.binPath + "data.txt", "w") as f:
			t = self.basePath + "/results/Results.cgns\n"
			t += str(self.grid.vertices.size) + " " + str(self.grid.elements.size) + " 0\n"
			t += '\n'.join(self.fieldsNames) + "\n";
			f.write(t)

		with open(self.binPath + "coords.txt", "w") as f:
			f.write( '\n'.join( [ ' '.join([str(c) for c in coord]) for coord in zip(*self.grid.gridData.vertices) ] ) )

		# Later create a way of guaranteeing that regionsElementsIndexes is a range (in MSHReader)
		with open(self.binPath + "connectivity.txt", "w") as f:
			t =  '\n'.join( [ ' '.join([str(e) for e in elem]) for elem in self.grid.gridData.elementsConnectivities ] )
			t += '\n'
			t += '\n'.join( [ ' '.join([str(e) for e in elem]) for elem in self.grid.gridData.boundariesConnectivities ] ) 
			f.write(t)

		with open(self.binPath + "sections.txt", "w") as f:
			t=""
			for regionName, regionsElementsIndexes in zip( self.grid.gridData.regionsNames, self.grid.gridData.regionsElementsIndexes ):
				t += regionName + ' ' + str( len(self.grid.gridData.elementsConnectivities[ regionsElementsIndexes[0] ]) ) + ' ' + str(regionsElementsIndexes[0]) + ' ' + str(regionsElementsIndexes[-1]) + '\n'
			for boundaryName, boundariesIndexes in zip( self.grid.gridData.boundariesNames, self.grid.gridData.boundariesIndexes ):
				t += boundaryName + ' ' + str( len(self.grid.gridData.boundariesConnectivities[ boundariesIndexes[0] ]) ) + ' ' + str(len(self.grid.gridData.elementsConnectivities) + boundariesIndexes[0]) + ' ' + str(len(self.grid.gridData.elementsConnectivities) + boundariesIndexes[-1]) + "\n"
			f.write(t)