import numpy as np
import h5py
import subprocess

class CgnsSaver:
	def __init__(self, timer, grid, outputPath="/home/gustavoe/Documents/Sinmec/GTRelated/GridReader/results/", basePath="/home/gustavoe/Documents/Sinmec/GTRelated/GridReader"):
		self.timer = timer
		self.grid = grid
		self.outputPath = outputPath + "Results.cgns"
		self.blankCgnsPath = basePath + "/meshes/cgns/blank.cgns"
		self.timeSteps  = np.	array([])
		self.timeFields = np.zeros((0, grid.vertices.size))

		self.configCgnsFile()

	def save(self, numericalField, currentTime):
		self.timeSteps	= np.append(self.timeSteps,  currentTime)
		self.timeFields = np.vstack([self.timeFields, numericalField])

	def finalize(self):
		del self.file["/BASE/TimeIterativeValues"]
		self.create_group("/BASE", "TimeIterativeValues", "I4", label="BaseIterativeData_t", data=np.array([ len(self.timeSteps)-1 ]))
		self.create_group("/BASE/TimeIterativeValues", "TimeValues", "R8", label="DataArray_t", data=self.timeSteps[:-1] )

		for count, result in enumerate(self.timeFields[1:],1):
			self.create_group("/BASE/ZONE", f"TimeStep{count}", "MT", label="FlowSolution_t")
			self.create_group(f"/BASE/ZONE/TimeStep{count}", "numerical temperature", "R8", label="DataArray_t", data=result)

	def configCgnsFile(self):
		self.createFile()
		self.setZone()

	def createFile(self):
		r=subprocess.run(f"yes | cp {self.blankCgnsPath} {self.outputPath}", shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
		print(r.stderr.decode("utf-8"))
		self.file = h5py.File(self.outputPath, "r+")

	def createZone(self):
		self.create_group("/BASE", "ZONE", "I4", data=np.array([1]))
		self.create_group("/BASE/ZONE", "ZoneType", "C1", data=np.array(list(b"Unstructured"), dtype=np.int16))

	def setZone(self):
		self.createZone()
		self.setCoordinates()
		self.setRegions()
		self.setBoundaries()


	def setCoordinates(self):
		self.create_group("/BASE/ZONE", "GridCoordinates", "MT")

		X,Y,Z = zip(*self.grid.gridData.vertices)
		self.create_group("/BASE/ZONE/GridCoordinates", "CoordinateX", "R8", data=X)
		self.create_group("/BASE/ZONE/GridCoordinates", "CoordinateY", "R8", data=Y)
		self.create_group("/BASE/ZONE/GridCoordinates", "CoordinateZ", "R8", data=Z)

		self.numberOfVertices = len(X)

	def setRegions(self):
		i=1
		for rName, rElements in zip(self.grid.gridData.regionNames, self.grid.gridData.regionElements):

			elements = [self.grid.gridData.elemConnectivity[index] for index in rElements]
			elements = [[nI+1 for nI in e] for e in elements]
			lower_r = len(sum(([[0]]+self.grid.gridData.regionElements)[:i],[]))
			upper_r = lower_r+len(elements)-1

			self.create_group("/BASE/ZONE", rName, "I4", label = "Elements_t", data=np.array([5,0], dtype=np.int32))
			
			self.create_group(f"/BASE/ZONE/{rName}", "ElementRange", "I4", label = "IndexRange_t", data=np.array([lower_r, upper_r], dtype=np.int32))
			self.create_group(f"/BASE/ZONE/{rName}", "ElementConnectivity", "I4", label = "DataArray_t", data=np.array(sum(elements,[]), dtype=np.int32))

			i+=1

	def setBoundaries(self):
		i=1
		for bName, bElements in zip(self.grid.gridData.boundaryNames, self.grid.gridData.boundaryElements):

			elements = [self.grid.gridData.boundaryElementsConnectivity[index] for index in bElements]
			lower_r = len(sum(([[0]]+self.grid.gridData.boundaryElements)[:i],[])) + len(sum(self.grid.gridData.regionElements,[]))
			upper_r = lower_r+len(elements)-1

			self.create_group("/BASE/ZONE", bName, "I4", label = "Elements_t", data=np.array([3,0], dtype=np.int32))
			self.create_group(f"/BASE/ZONE/{bName}", "ElementRange", "I4", label = "IndexRange_t", data=np.array([lower_r, upper_r], dtype=np.int32))
			self.create_group(f"/BASE/ZONE/{bName}", "ElementConnectivity", "I4", label = "DataArray_t", data=np.array(sum(elements,[]), dtype=np.int32))

			i+=1

	def create_group(self, parent_path, group_name, _type, label = None, label_lgth = 15, data=None):
		if label == None: label = group_name.capitalize()+"_t"
		if len(label) > label_lgth: label_lgth = 33
		parent = self.file[parent_path]
		g = parent.create_group(group_name)
		g.attrs["label"] = np.string_(label.ljust(label_lgth,"\x00"))
		g.attrs["type"]  = np.string_(_type.ljust(3, "\x00"))
		g.attrs["flags"] = np.array([1], dtype=np.int32)
		if data.__class__ != None.__class__:
			g.create_dataset(" data", data=data)


if __name__ == "__main__":
	import sys
	sys.path.append("/home/gustavoe/Documents/Sinmec/GTRelated/GridReader")
	from libs.simulation.HeatTransfer2D import HeatTransfer2D
	from libs.geometry.Grid import Grid
	from libs.geometry.MSHReader import MSHReader
	from libs.simulation.ProblemData2D import ProblemData2D

	grid = Grid(MSHReader(ProblemData2D('heat_transfer_2d').paths["Grid"]).getData())
	c = CgnsSaver('', grid)
	c.finalize()