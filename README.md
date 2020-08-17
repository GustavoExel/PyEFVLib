# PyEFVLib
Library for implementing the Element-based Finite Volumes Method. The input mesh is a .msh file, and the output is a .cgns file. 

## Prerequisites

This project is writen in Python 3
You may install numpy, scipy, pandas and matplotlib

```bash
pip install numpy
pip install scipy
pip install pandas
pip install matplotlib
```

Also, for CGNS writing C++ is used, and two libraries are required: [CGNS](https://cgns.github.io) and [Boost](https://www.boost.org).
After installing them, configure their install directories in PyEFVLib > simulation > CGNS > CMakeLists.txt, and compile using
```bash
./install.sh
```

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