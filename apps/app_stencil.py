import sys,os
sys.path.append(os.path.join(os.path.dirname(__file__), os.path.pardir))
from PyEFVLib import MSHReader, Grid, ProblemData, CgnsSaver, CsvSaver


if __name__ == "__main__":

	model = "workspace/heat_transfer_2d/linear"
	problemData = ProblemData(model)

	print(problemData.paths["Grid"])
	reader = MSHReader(problemData.paths["Grid"])
	grid = Grid(reader.getData())

	grid.buildStencil()
	print(grid.stencil)
