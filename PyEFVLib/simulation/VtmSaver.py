import numpy as np
import subprocess, os, sys, shutil
from PyEFVLib.geometry.Shape import Triangle, Quadrilateral, Tetrahedron, Hexahedron, Prism, Pyramid
from PyEFVLib.simulation.Saver import Saver

class VtmSaver(Saver):
	def __init__(self, grid, outputPath, basePath, fileName="Results", **kwargs): 
		Saver.__init__(self, grid, outputPath, basePath, 'vtu', fileName)
		self.outputPath = '.'.join(self.outputPath.split('.')[:-1])

	def finalize(self):
		shapesDict	 = { Triangle:5,Quadrilateral:9,Tetrahedron:10,Hexahedron:12,Prism:13,Pyramid:14 }
		connectivity = [ [ vertex.handle for vertex in element.vertices ] for element in self.grid.elements ]
		shapes		 = [ shapesDict[element.shape] for element in self.grid.elements ]
		offsets		 = [ len(conn) for conn in connectivity ]
		offsets		 = [ sum(offsets[:i+1]) for i in range( len(offsets) ) ]

		if os.path.isdir(self.outputPath):
			shutil.rmtree(self.outputPath)
		os.makedirs(self.outputPath)
		os.makedirs( os.path.join(self.outputPath, "timedata") )

		for i, timeStep in enumerate(self.timeSteps):
			idx = i+1
			with open(os.path.join(self.outputPath, f"results_{idx}.vtm"), "w") as vtmFile:
				vtmFile.write( '<VTKFile type="vtkMultiBlockDataSet" version="1.0" byte_order="LittleEndian" header_type="UInt64">\n\t<vtkMultiBlockDataSet>\n' )
				vtmFile.write( f'\t\t<Block name="BASE"><DataSet name="ZONE" file="timedata/results_{idx}.vtu"/></Block>\n\t</vtkMultiBlockDataSet>\n' )
				vtmFile.write( f'\t<FieldData><DataArray type="Float64" Name="TimeValue" format="ascii"> {timeStep} </DataArray></FieldData>\n</VTKFile>\n' )

			with open(os.path.join(self.outputPath, "timedata", f"results_{idx}.vtu"), "w") as vtuFile:	
				vtuFile.write( '<VTKFile type="UnstructuredGrid" version="1.0" byte_order="LittleEndian" header_type="UInt64">\n\t<UnstructuredGrid>\n\t\t<FieldData>\n' )
				vtuFile.write( f'\t\t\t<DataArray type="Float64" Name="TimeValue" NumberOfTuples="1" format="ascii"> {timeStep} </DataArray>\n\t\t</FieldData>\n' )
				vtuFile.write( f'\t\t<Piece NumberOfPoints="{self.grid.numberOfVertices}" NumberOfCells="{self.grid.elements.size}">\n\t\t\t<PointData>\n' )
				for fieldName in self.fields.keys():
					vtuFile.write( f'\t\t\t\t<DataArray type="Float64" Name="{fieldName}" format="ascii">\n' )
					vtuFile.write("\t\t\t\t\t" + "".join([ f"{d:.15f} " for d in self.fields[fieldName][i] ]) + "\n")

					vtuFile.write( '\t\t\t\t</DataArray>\n' )
				vtuFile.write( '\t\t\t</PointData>\n\t\t\t<CellData></CellData>\n\t\t\t<Points>\n\t\t\t\t<DataArray type="Float64" Name="Points" NumberOfComponents="3" format="ascii">\n' )
				vtuFile.write( "".join( [ f"\t\t\t\t\t{c:.8f}\n" for vertex in self.grid.vertices for c in vertex.getCoordinates() ] ) )
				vtuFile.write( '\t\t\t\t\t<InformationKey name="L2_NORM_RANGE" location="vtkDataArray" length="2">\n\t\t\t\t\t\t<Value index="0"> 0 </Value>\n\t\t\t\t\t\t<Value index="1"> 1.7320508076 </Value>\n\t\t\t\t\t</InformationKey>\n\t\t\t\t</DataArray>\n\t\t\t</Points>\n\t\t\t<Cells>\n\t\t\t\t<DataArray type="Int64" Name="connectivity" format="ascii">\n' )
				vtuFile.write( "\t\t\t\t\t" + "".join( [ f"{v:.0f} " for c in connectivity for v in c ] ) + "\n")
				vtuFile.write( '\t\t\t\t</DataArray>\n\t\t\t\t<DataArray type="Int64" Name="offsets" format="ascii">\n' )
				vtuFile.write( "\t\t\t\t\t" + "".join( [ f"{o:.0f} " for o in offsets ] ) + "\n")
				vtuFile.write( '\t\t\t\t</DataArray>\n\t\t\t\t<DataArray type="UInt8" Name="types" format="ascii">\n' )
				vtuFile.write( "\t\t\t\t\t" + "".join( [ f"{s:.0f} " for s in shapes ] ) + "\n" )
				vtuFile.write( '\t\t\t\t</DataArray>\n\t\t\t</Cells>\n\t\t</Piece>\n\t</UnstructuredGrid>\n</VTKFile>\n' )

		self.finalized = True