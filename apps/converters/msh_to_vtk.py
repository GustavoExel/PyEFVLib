"""
<?xml version="1.0"?>
<VTKFile type="UnstructuredGrid" version="0.1" byte_order="LittleEndian">
	<UnstructuredGrid>
		<Piece>
			<Points>
				<DataArray type="Float64" Name="Points" NumberOfComponents="3" format="ascii">

				</DataArray>
			</Points>

			<Cells>
				<DataArray type="Int32" Name="connectivity" format="ascii"></DataArray>
				<DataArray type="Int32" Name="offsets" format="ascii"></DataArray>
				<DataArray type="Int32" Name="types" format="ascii"></DataArray>
			</Cells>

			<PointData>
				<DataArray type="Float64" Name="temperature" format="ascii"></DataArray>
			</PointData>

		</Piece>
	</UnstructuredGrid>
</VTKFile>

Obs.:
	Types:
		Line:3, Tri:5, Quad:9, Tetra:10, Hexa:12, Prism:13, Pyramid:14
		{"line":"12", "triangle":"22", "quadrilateral":"32", "tetrahedron":"42", "pyramid":"72", "prism":"62", "hexagon":"52"}
		shapesDict = {(1,2):3,(2,2):5,(3,2):9,(4,2):10,(7,2):14,(6,2):13,(5,2):12}

"""
import meshio
import sys
import pandas as pd
import numpy as np

if "--help" in sys.argv:
	print("Usage:\npython msh_to_vtk.py [MSH file path] [CSV data file path] [Output file path]")

try:
	mshFilePath = sys.argv[1]
	csvDataFilePath = sys.argv[2]
	outputFilePath = sys.argv[3]
except:
	print("Usage:\npython msh_to_vtk.py [MSH file path] [Output file path]")

mesh = meshio.read(mshFilePath)
dataFrame = pd.read_csv(csvDataFilePath)
data = np.array( dataFrame[ dataFrame.columns[-1] ] )

with open(outputFilePath, "w") as file:
	shapesDict = { (1,2):3, (2,2):5, (3,2):9, (4,2):10, (7,2):14, (6,2):13, (5,2):12 }
	dimensionDict = {3:1,5:2,9:2,10:3,14:3,13:3,12:3}

	with open(mshFilePath, "r") as mshFile:
		mshRawData = mshFile.readlines()

	elementsRawData = [[int(d) for i,d in enumerate(line.split()) if i in [1,2] or i>=5] for line in mshRawData[ mshRawData.index("$Elements\n")+2 : mshRawData.index("$EndElements\n") ]]
	shapes, connectivity = zip(*[ ( shapesDict[(s1,s2)], np.array(conn)-1 ) for s1,s2,*conn in elementsRawData ])
	dim = max([dimensionDict[s] for s in shapes])
	shapes, connectivity = zip(*[ (s,c) for s,c in zip(shapes, connectivity) if dimensionDict[s]==dim ])
	offsets = [ len(conn) for conn in connectivity ]
	offsets = [ sum(offsets[:i+1]) for i in range( len(offsets) ) ]

	file.write("<?xml version=\"1.0\"?>\n<VTKFile type=\"UnstructuredGrid\" version=\"0.1\" byte_order=\"LittleEndian\">\n\t<UnstructuredGrid>\n")
	file.write(f"\t\t<Piece NumberOfPoints=\"{ len(mesh.points) }\" NumberOfCells=\"{ len(shapes) }\">\n\t\t\t<Points>\n")
	file.write("\t\t\t\t<DataArray type=\"Float64\" Name=\"Points\" NumberOfComponents=\"3\" format=\"ascii\">\n")
	
	# WRITE THE COORDINATES HERE
	file.write( "".join( [ f"\t\t\t\t\t{c:.16e}\n" for c in mesh.points.reshape( mesh.points.size ) ] ) )
	
	file.write("\t\t\t\t</DataArray>\n\t\t\t</Points>\n\t\t\t<Cells>\n")
	file.write("\t\t\t\t<DataArray type=\"Int32\" Name=\"connectivity\" format=\"ascii\">\n")
	# WRITE connectivity HERE
	file.write( "\t\t\t\t\t" + "".join( [ f"{v:.0f} " for c in connectivity for v in c ] ) + "\n")

	file.write("\t\t\t\t</DataArray>\n")
	file.write("\t\t\t\t<DataArray type=\"Int32\" Name=\"offsets\" format=\"ascii\">\n")
	# WRITE offsets HERE
	file.write( "\t\t\t\t\t" + "".join( [ f"{o:.0f} " for o in offsets ] ) + "\n")

	file.write("\t\t\t\t</DataArray>\n")
	file.write("\t\t\t\t<DataArray type=\"Int32\" Name=\"types\" format=\"ascii\">\n")
	# WRITE types HERE
	file.write( "\t\t\t\t\t" + "".join( [ f"{s:.0f} " for s in shapes ] ) + "\n" )

	file.write("\t\t\t\t</DataArray>\n")
	file.write("\t\t\t</Cells>\n")
	file.write("\t\t\t<PointData>\n")
	file.write("\t\t\t\t<DataArray type=\"Float64\" Name=\"TEMPERATURE\" format=\"ascii\">\n")
	# WRITE FIELD DATA HERE
	file.write("\t\t\t\t\t" + "".join([f"{d:.0f} " for d in data]) + "\n")

	file.write("\t\t\t\t</DataArray>\n")
	file.write("\t\t\t</PointData>\n\t\t</Piece>\n\t</UnstructuredGrid>\n</VTKFile>\n")