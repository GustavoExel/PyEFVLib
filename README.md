# PyEFVLib

This package intends to support the solution of PDEs using the Element-based Finite Volume Method (EbFVM). The input mesh may be \*.msh or \*.xdmf files, and the output may be \*.csv or \*.cgns.

## Dependencies & Installation

- [Python 3](https://www.python.org/downloads/) (3.8.2);
- [matplotlib](https://matplotlib.org/) (3.3.2);
- [meshio](https://pypi.org/project/meshio/) (4.0.15);
- [numpy](https://numpy.org/) (1.17.4);
- [pandas](https://pandas.pydata.org/)(1.1.3);
- [petsc4py](https://pypi.org/project/petsc4py/) (3.12.0);
- [scipy](https://www.scipy.org/) (1.5.3);
- [xmltodict](https://pypi.org/project/xmltodict/) (0.12.0).

## Usage

```python
import PyEFVLib

grid = PyEFVLib.read("meshes/msh/2D/Square.msh")

totalVolume = 0.0
for element in grid.elements:
	for vertex in element:
		totalVolume += vertex.volume

print(totalVolume)
```

## Latest changes
- ProblemData doesn't read from workspace anymore, instead you can inform the settings like so
```python
import PyEFVLib

problemData = PyEFVLib.ProblemData(
	meshFilePath = "{MESHES}/msh/2D/Square.msh",
	outputFilePath = "{RESULTS}/heat_transfer_2d/linear",
	numericalSettings = PyEFVLib.NumericalSettings( timeStep = 1e-02, tolerance = 1e-06, maxNumberOfIterations = 300 ),
	propertyData = PyEFVLib.PropertyData({
		"Body" : {
			"HeatCapacity"	: 1.0,
			"Conductivity"	: 1.0,
			"Density"		: 1.0,
			"HeatGeneration": 0.0,
		},
	}),
	boundaryConditions = PyEFVLib.BoundaryConditions({
		"temperature": {
			"InitialValue": 0.0,
			"West":	 { "condition" : PyEFVLib.Dirichlet, "type" : PyEFVLib.Constant,"value" : 20.0 },
			"East":	 { "condition" : PyEFVLib.Dirichlet, "type" : PyEFVLib.Constant,"value" : 50.0 },
			"South": { "condition" : PyEFVLib.Neumann,   "type" : PyEFVLib.Constant,"value" : 0.0 },
			"North": { "condition" : PyEFVLib.Neumann,   "type" : PyEFVLib.Constant,"value" : 0.0 },
		},
	}),
)
```

- Getting the region property now must be done using the get method
```python
# Before
problemData.propertyData[region.handle]["PropertyName"]
# Now
problemData.propertyData.get(region.handle, "PropertyName")
```

- Creating a grid object now is made easier with problemData
```python
# Before
problemData = ProblemData(model)
reader = MSHReader(problemData.meshFilePath)
grid = Grid(reader.getData())
problemData.setGrid(grid)
problemData.read()

# Now
problemData = PyEFVLib.ProblemData(meshFilePath, outputFilePath, numericalSettings, propertyData, boundaryConditions)
grid = problemData.grid
```
