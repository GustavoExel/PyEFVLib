# GridReader
Library for reading (currently) .msh mesh files, and create a Grid object with different objects and methods that facititates aplying the Elements-based Finite Volumes Method.

## Prerequisites

This project is writen in Python 3
You may install numpy

```bash
pip install numpy
```

## Usage

```python
from libs.geometry.MSHReader import MSHReader
from libs.geometry.Grid import Grid


reader = MSHReader("/home/gustavoe/Documents/Sinmec/GTRelated/GridReader/meshes/Square.msh")
grid = Grid(r.getData())

totalVolume = 0.0
for element in grid.elements:
	for vertex in element:
		totalVolume += vertex.volume

print(totalVolume)
```
