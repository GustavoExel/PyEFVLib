import sys, os
import meshio
import pandas as pd
import numpy as np
pyefvlibPath = os.path.join(os.path.dirname(__file__), os.path.pardir, os.path.pardir)
sys.path.append(pyefvlibPath)
import PyEFVLib

# py csv_msh_to_xdmf.py ../../meshes/msh/2D/column6.msh ../../results/geomechanics/Results.csv
if len(sys.argv)<3 or (len(sys.argv)>1 and sys.argv[1] == "--help"):
	print("Usage:\npy csv_msh_to_xdmf.py csv_filePath msh_filePath")
	exit()

if len(sys.argv) >= 3:
	csvFilePath = sys.argv[1]
	mshFilePath = sys.argv[2]
	outputPath = (sys.argv[3] if len(sys.argv)==4 else f"{pyefvlibPath}/results") + "/output.xdmf"

	df = pd.read_csv(csvFilePath)

	fieldNames = [ colName.split(" - ")[-1] for colName in df.columns[3:] ]
	fieldNames = [ field for idx,field in enumerate(fieldNames) if field not in fieldNames[:idx] ]
	
	lastTimeStep = int(df.columns[-1].split(" - ")[0].replace("TimeStep",""))
	timeSteps = range(1, lastTimeStep+1)

	fields = dict()

	for fieldName in fieldNames:
		fields[fieldName] = []
		for timeStep in timeSteps:
			fields[fieldName].append( np.array(df[f"TimeStep{timeStep} - {fieldName}"]) )


	grid = PyEFVLib.read(mshFilePath)

	points = np.array( [v.getCoordinates() for v in grid.vertices] )

	meshioShapes   = ["triangle", "quad", "tetra", "pyramid", "wedge", "hexahedron"]
	pyEFVLibShapes = [PyEFVLib.Triangle, PyEFVLib.Quadrilateral, PyEFVLib.Tetrahedron, PyEFVLib.Pyramid, PyEFVLib.Prism, PyEFVLib.Hexahedron]

	cells = [ ( shape , np.array([[vertex.handle for vertex in element.vertices] for element in grid.elements if element.shape == shapeClass], dtype=np.uint64) ) for shape, shapeClass in zip(meshioShapes, pyEFVLibShapes) ]
	cells = [ cell for cell in cells if cell[1].size ]

	with meshio.xdmf.TimeSeriesWriter(outputPath) as writer:
		writer.write_points_cells(points, cells)
		for idx, timeStep, in enumerate(timeSteps):
			data  = { fieldName : fields[fieldName][idx] for fieldName in fields }
			writer.write_data(timeStep, point_data=data)