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

Also, for CGNS writing, C++ is used, and two libraries are required: [CGNS](https://cgns.github.io) and [Boost](https://www.boost.org).
<!-- After installing them, configure their install directories in PyEFVLib > simulation > CGNS > CMakeLists.txt, and compile using
```bash
./install.sh
``` -->

## Usage

```python
from PyEFVLib import MSHReader, Grid, Point
import os, numpy as np

path = os.path.join(*[os.path.dirname(__file__), os.path.pardir, "meshes", "Square.msh"])

reader = MSHReader(path)
grid   = Grid(reader.getData())

totalVolume = 0.0
for element in grid.elements:
	for vertex in element:
		totalVolume += vertex.volume

print(totalVolume)
```