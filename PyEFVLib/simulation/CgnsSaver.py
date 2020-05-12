import numpy as np
import subprocess, os, sys

class CgnsSaver:
	def __init__(self, grid, outputPath, basePath): 
		self.grid = grid
		self.outputPath = outputPath + "Results.cgns"
		self.basePath = basePath
		self.binPath = os.path.join( basePath, "PyEFVLib", "simulation" ) + "/"

		self.timeSteps  = np.	array([])
		self.fields = np.zeros((0, grid.vertices.size))
		self.createCgns()
		self.finalized = False


	def save(self, field, currentTime):
		self.timeSteps = np.append( self.timeSteps, currentTime )
		self.fields = np.vstack([self.fields, field])

	def __del__(self):
		if not self.finalize:
			self.finalize()
	def finalize(self):
		with open(self.binPath + "fields.txt", "w") as f:
			f.write( '\n'.join([' '.join([str(x) for x in field]) for field in self.fields]) )

		with open(self.binPath + "steps.txt", "w") as f:
			f.write( ' '.join([ str(ts) for ts in self.timeSteps ]) )

		subprocess.run([self.binPath + "save"])
		os.remove(self.binPath + "data.txt")
		os.remove(self.binPath + "fields.txt")
		os.remove(self.binPath + "steps.txt")

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
			t += str(self.grid.vertices.size) + "\n" + str(self.grid.elements.size) + "\n"
			f.write(t)

		with open(self.binPath + "coords.txt", "w") as f:
			f.write( '\n'.join( [ ' '.join([str(c) for c in coord]) for coord in zip(*self.grid.gridData.vertices) ] ) )

		# Later create a way of guaranteeing that regionElements is a range (in MSHReader)
		with open(self.binPath + "connectivity.txt", "w") as f:
			t =  '\n'.join( [ ' '.join([str(e) for e in elem]) for elem in self.grid.gridData.elemConnectivity ] )
			t += '\n'
			t += '\n'.join( [ ' '.join([str(e) for e in elem]) for elem in self.grid.gridData.boundaryElementsConnectivity ] ) 
			f.write(t)

		with open(self.binPath + "sections.txt", "w") as f:
			t=""
			for regionName, regionElements in zip( self.grid.gridData.regionNames, self.grid.gridData.regionElements ):
				t += regionName + ' ' + str( len(self.grid.gridData.elemConnectivity[ regionElements[0] ]) ) + ' ' + str(regionElements[0]) + ' ' + str(regionElements[-1]) + '\n'
			for boundaryName, boundaryElements in zip( self.grid.gridData.boundaryNames, self.grid.gridData.boundaryElements ):
				t += boundaryName + ' ' + str( len(self.grid.gridData.boundaryElementsConnectivity[ boundaryElements[0] ]) ) + ' ' + str(len(self.grid.gridData.elemConnectivity) + boundaryElements[0]) + ' ' + str(len(self.grid.gridData.elemConnectivity) + boundaryElements[-1]) + "\n"
			f.write(t)