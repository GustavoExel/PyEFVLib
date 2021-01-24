import sys,os
sys.path.append(os.path.join(os.path.dirname(__file__), os.path.pardir))
from PyEFVLib import MSHReader, CgnsSaver, Grid
import pandas as pd

if "--help" in sys.argv:
	print("Usage: python apps/csv_to_cgns.py [mesh path] [csv path]")
if len(sys.argv) < 3:
	raise Exception("Must insert msh and csv paths")
for arg in sys.argv:
	if ".msh" in arg:
		mshPath = arg
	elif ".csv" in arg:
		csvPath = arg
try:
	mshPath = os.path.realpath( mshPath )
	csvPath = os.path.realpath( csvPath )
	outputPath = os.path.sep.join( csvPath.split(os.path.sep)[:-1] )
except:
	raise Exception("Invalid input path(s)")

data 	  = pd.read_csv(csvPath) 
reader 	  = MSHReader(mshPath)
grid   	  = Grid(reader.getData())
saver = CgnsSaver(grid, outputPath, os.path.join(os.path.dirname(__file__),os.path.pardir))

# First we'll assume that the data is ordered acording to the msh file
for columnLabel in data.columns[3:]:
	columnData = data[columnLabel]
	fieldName = " - ".join( columnLabel.split(" - ")[1:] )
	timeStep = int( columnLabel.split(" - ")[0].replace("TimeStep", "") )
	saver.save(fieldName, data[columnLabel], timeStep)
saver.finalize()

print(f"Output Path: {saver.outputPath}")