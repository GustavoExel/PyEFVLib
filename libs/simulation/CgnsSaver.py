import numpy as np
import h5py
import subprocess

class CgnsSaver:
	def __init__(self, timer, grid, outputPath, basePath):
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
		self.attribute("/BASE/TimeIterativeValues", np.array([ len(self.timeSteps) ]) )
		self.attribute("/BASE/TimeIterativeValues/TimeValues", self.timeSteps)  # Make sure that self.timeSteps is a numpy.Array()

		for count, result in enumerate(self.timeFields[1:],1):
			self.copy("/BASE/ZONE", "TimeStep", "TimeStep"+str(count))
			self.attribute(f"/BASE/ZONE/TimeStep{str(count)}/numerical temperature", result)
			t = self.file[f"/BASE/ZONE/TimeStep{str(count)}"]
			t.attrs["label"] = np.string_("FlowSolution_t".ljust(33,"\x00"))
			# t.attrs["name"] = np.string_(f"TimeStep{str(count)}".ljust(33,"\x00"))
			t.attrs["type"] = np.string_("MT".ljust(3,"\x00"))
			
			# self.file["/BASE/ZONE"].create_group(f"TimeStep{count}")
			# self.file["/BASE/ZONE"][f"TimeStep{count}"].create_group("numerical temperature")
			# self.file["/BASE/ZONE"][f"TimeStep{count}/numerical temperature"].create_dataset(" data", data=result)
			# self.file["/BASE/ZONE"][f"TimeStep{count}/numerical temperature"].attrs["label"] = np.string_('DataArray_t'.ljust(33,'\x00'))
			# # self.file["/BASE/ZONE"][f"TimeStep{count}/numerical temperature"].attrs["name"] = np.string_('numerical temperature'.ljust(33,'\x00'))
			# self.file["/BASE/ZONE"][f"TimeStep{count}/numerical temperature"].attrs["type"] = np.string_('R8\x00')
			# self.file["/BASE/ZONE"][f"TimeStep{count}/numerical temperature"].attrs["flags"] = np.array([1])
			# self.file["/BASE/ZONE"][f"TimeStep{count}"].attrs["label"] = np.string_("FlowSolution_t".ljust(33, '\x00'))
			# # self.file["/BASE/ZONE"][f"TimeStep{count}"].attrs["name"]  = np.string_(f"TimeStep{count}".ljust(33, '\x00'))
			# self.file["/BASE/ZONE"][f"TimeStep{count}"].attrs["type"]  = np.string_("MT".ljust(3, '\x00'))
			# self.file["/BASE/ZONE"][f"TimeStep{count}"].attrs["flags"] = np.array([1])

		del self.file["/BASE/ZONE/TimeStep"]

		self.file.close()
		del self.timeFields
		del self.timeSteps

	def configCgnsFile(self):
		self.createFile()
		self.setZone()

	def createFile(self):
		r=subprocess.run(f"yes | cp {self.blankCgnsPath} {self.outputPath}", shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
		print(r.stderr.decode("utf-8"))
		self.file = h5py.File(self.outputPath, "r+")

	def setZone(self):
		self.setCoordinates()
		self.setRegions()
		self.setBoundaries()

		self.attribute("/BASE/ZONE", np.array([[len(self.grid.gridData.vertices)],[	len(sum(self.grid.gridData.regionElements,[]))],[0]], dtype=np.int32))
		
		del self.file["/BASE/ZONE/PHYSICAL_NAME"]

	def setCoordinates(self):
		X,Y,Z = zip(*self.grid.gridData.vertices)
		self.attribute("/BASE/ZONE/GridCoordinates/CoordinateX", np.array(X, dtype=np.float64))
		self.attribute("/BASE/ZONE/GridCoordinates/CoordinateY", np.array(Y, dtype=np.float64))
		self.attribute("/BASE/ZONE/GridCoordinates/CoordinateZ", np.array(Z, dtype=np.float64))
		self.numberOfVertices = len(X)

	def setRegions(self):
		i=1
		for rName, rElements in zip(self.grid.gridData.regionNames, self.grid.gridData.regionElements):
			self.copy("/BASE/ZONE", "PHYSICAL_NAME", rName )

			elements = [self.grid.gridData.elemConnectivity[index] for index in rElements]
			elements = [[nI+1 for nI in e] for e in elements]
			lower_r = len(sum(([[0]]+self.grid.gridData.regionElements)[:i],[]))
			upper_r = lower_r+len(elements)-1

			self.attribute(f"/BASE/ZONE/{rName}/ElementRange", np.array([lower_r, upper_r], dtype=np.int32) )
			self.attribute(f"/BASE/ZONE/{rName}/ElementConnectivity", np.array(sum(elements,[]), dtype=np.int32) )

			i+=1


	def setBoundaries(self):
		i=1
		for bName, bElements in zip(self.grid.gridData.boundaryNames, self.grid.gridData.boundaryElements):
			self.copy("/BASE/ZONE", "PHYSICAL_NAME", bName )
			self.attribute("/BASE/ZONE/"+bName, np.array([3,0], dtype=np.int32))

			elements = [self.grid.gridData.boundaryElementsConnectivity[index] for index in bElements]
			lower_r = len(sum(([[0]]+self.grid.gridData.boundaryElements)[:i],[])) + len(sum(self.grid.gridData.regionElements,[]))
			upper_r = lower_r+len(elements)-1

			self.attribute(f"/BASE/ZONE/{bName}/ElementRange", np.array([lower_r, upper_r], dtype=np.int32) )
			self.attribute(f"/BASE/ZONE/{bName}/ElementConnectivity", np.array(sum(elements,[]), dtype=np.int32) )

			i+=1

	def attribute(self, path, data):
		del self.file[path+"/ data"]
		self.file[path].create_dataset(" data", data=data)

	def copy(self, root, group, name):
		self.file[root].create_group("temp")
		self.file[root+"/temp/"+name] = self.file[root+"/"+group]
		self.file[root].copy("temp/"+name, self.file[root])
		del self.file[root+"/temp"]