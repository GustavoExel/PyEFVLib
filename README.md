# EbFVM
PyEFVLib is a Python library that provides tools to support the solution of systems of partial differential equations (PDEs) using the Element-based Finite Volumes Method (EbFVM)

## Summary
1. **Instalation**
2. **Numerical Method**
3. **PyEFVLib Classes - Geometrical Entities**
4. **PyEFVLib Classes - Non Geometrical Entities**
5. **Tutorial**
6. **3rd-party Softwares**
7. **2D vs 3D considerations**
8. **bellbird**

## 1. Instalation
Install PyEFVLib using pip, by typing in the command line:

```bash
pip install PyEFVLib
```

## 2. Numerical Method

### 2.1 - Partial Differential Equations
- A differential equation is an equation that relates one or more functions and their derivatives.
- A partial differential equation (PDE) is an equation that relates one or more functions dependent on two or more variables, to the partial derivatives of these functions in relation to the independent variables. 
- In physics, differential equations relate physical properties (e.g., pressure, temperature, displacement, turbulent kinetic energy, speed) to their spatial and temporal derivatives, and model natural and everyday phenomena. 

---
### 2.2 - Example: Integrating a PDE in space

![control_volume](https://user-images.githubusercontent.com/57679731/118341693-6b4b8d00-b4f6-11eb-915e-357af75ecda8.png)

Take the [heat equation](https://en.wikipedia.org/wiki/Heat_equation#Non-uniform_isotropic_medium) as an example

![eq1](https://latex.codecogs.com/png.latex?%5Cdpi%7B120%7D%20%5Cbg_white%20%5Clarge%20%5Cnabla%20%5Ccdotp%20%28%20k%5Cnabla%20T%29%20&plus;q%27%27%27%3D%5Crho%20c_%7Bp%7D%5Ctfrac%7B%5Cpartial%20T%7D%7B%5Cpartial%20t%7D)

Integrating in the control volume (CV)

![eq2](https://latex.codecogs.com/png.latex?%5Cdpi%7B120%7D%20%5Cbg_white%20%5Clarge%20%5Cint%20_%7B%7B%7B%5COmega%7D_%7Bi%7D%7D%7D%5B%20%5Cnabla%20%5Ccdotp%20%28%20k%5Cnabla%20T%29%20&plus;q%27%27%27%5D%20%5C%20%7Bd%7B%5COmega%7D_%7Bi%7D%7D%20%3D%5Cint%20_%7B%7B%7B%5COmega%7D_%7Bi%7D%7D%7D%20%5Crho%20c_%7Bp%7D%5Ctfrac%7B%5Cpartial%20T%7D%7B%5Cpartial%20t%7D%20d%7B%5COmega%7D_%7Bi%7D)

And applying the [divergence theorem](https://en.wikipedia.org/wiki/Divergence_theorem)

![eq3](https://latex.codecogs.com/png.latex?%5Cdpi%7B120%7D%20%5Cbg_white%20%5Clarge%20%5Cint%20_%7B%5CGamma%20%7B_%7Bi%7D%7D%7D%28%20k%5Cnabla%20T%29%20%5Ccdotp%20%5Chat%7Bn%7D%20%5C%20%7Bd%5CGamma%20_%7Bi%7D%7D%20&plus;%5Cint%20_%7B%7B%7B%5COmega%7D_%7Bi%7D%7D%7D%20q%27%27%27%7Bd%7B%5COmega%7D_%7Bi%7D%7D%20%3D%5Cint%20_%7B%7B%7B%5COmega%7D_%7Bi%7D%7D%7D%20%5Crho%20c_%7Bp%7D%5Ctfrac%7B%5Cpartial%20T%7D%7B%5Cpartial%20t%7D%20d%7B%5COmega%7D_%7Bi%7D)

![eq4](https://latex.codecogs.com/png.latex?%5Cdpi%7B120%7D%20%5Cbg_white%20%5Clarge%20%5Csum%20_%7Bf%5Cin%20%5CGamma%20_%7Bi%7D%7D%28%20k%5Cnabla%20T%29%20%5Ccdotp%20%5Coverrightarrow%7B%5CDelta%20s%7D_%7Bf%7D%20&plus;q%27%27%27%5C%20%5CDelta%20%5COmega%20_%7Bi%7D%20%3D%5Crho%20c_%7Bp%7D%5Ctfrac%7B%5Cpartial%20T_%7Bi%7D%7D%7B%5Cpartial%20t%7D%20%5CDelta%20%5COmega%20_%7Bi%7D)

![eq5](https://latex.codecogs.com/png.latex?%5Cdpi%7B120%7D%20%5Cbg_white%20%5Clarge%20%5Csum%20_%7Bf%5Cin%20%5CGamma%20_%7Bi%7D%7D%5Csum%20_%7Bpi%5Cin%20f%7D%28%20k%5Cnabla%20T%29%20%5Ccdotp%20%5Coverrightarrow%7B%5CDelta%20s%7D_%7Bpi%7D%20&plus;q%27%27%27%5CDelta%20%5COmega%20_%7Bi%7D%20%3D%5Crho%20c_%7Bp%7D%5Ctfrac%7B%5CDelta%20T_%7Bi%7D%7D%7B%5CDelta%20t%7D%20%5CDelta%20%5COmega%20_%7Bi%7D)

---
### 2.3 - Domain Discretization

![mesh discretization](https://user-images.githubusercontent.com/57679731/118343663-9f2bb000-b500-11eb-8483-1d3d12c9f37a.png)


Previously it was shown that a differential equation in the differential form can be integrated in any control volume, but to achieve a detailed solution along the domain of interest, the finite volume method discretizes the domain in several control volumes forming a mesh. The more refined the mesh, the better the solution to the problem will be.
> Note: The black bordered triangles shown in the figure are actually elements, the control volumes are outlined by a blue lines.

---
### 2.4 - Field Approximations
As it is the nature of PDEs, the equations involve spatial derivatives, but since we are storing a finite number of variables, we need a way to aproximate the field continuously to evaluate the spatial derivatives, and EbFVM does this using the shape functions.

We'll start by evaluating the field within an element of our mesh. It'll be expressed by a linear combination of the property values at the vertices of the element, and the weight of the property at each vertex is given by the shape functions at each point inside the element, as in the figure below: 

![shape functions](https://user-images.githubusercontent.com/57679731/118343196-dba9dc80-b4fd-11eb-9e02-69758af70f80.png)

And as we know the entire field inside the element, we can evaluate its partial derivatives, as shown below: 

![eq4](https://latex.codecogs.com/png.latex?%5Cdpi%7B120%7D%20%5Cbg_white%20%5Clarge%20p%28%20r%29%20%3D%5Csum%20_%7Bk%7D%5Cmathcal%7BN%7D_%7Bk%7D%28%20r%29%20%5Ccdotp%20p_%7Bk%7D%5C%5C%20%5Cfrac%7B%5Cpartial%20p%28%20r%29%7D%7B%5Cpartial%20x_%7Bi%7D%7D%20%3D%5Csum%20_%7Bk%7D%5Cfrac%7B%5Cpartial%20%5Cmathcal%7BN%7D_%7Bk%7D%28%20r%29%7D%7B%5Cpartial%20x_%7Bi%7D%7D%20%5Ccdotp%20p_%7Bk%7D)

---
### 2.5 - Problem Solution
Previously, the differential equations expressed the values of the unknowns as a function of their partial derivatives, but as we can express the partial derivatives as a function of the values of the field at the neighboring vertices, so the solution of the problem becomes a system of algebraic equations. 

![eq6](https://latex.codecogs.com/png.latex?%5Cdpi%7B120%7D%20%5Cbg_white%20%5Clarge%20%5Cbegin%7Bcases%7D%20f_%7B1%7D%28%20p_%7B1%7D%20%2Cp_%7B2%7D%20%2C...%2Cp_%7BN%7D%29%20%3D0%5C%5C%20f_%7B2%7D%28%20p_%7B1%7D%20%2Cp_%7B2%7D%20%2C...%2Cp_%7BN%7D%29%20%3D0%5C%5C%20%5C%20%5C%20%5C%20%5Cvdots%20%5C%5C%20f_%7BN%7D%28%20p_%7B1%7D%20%2Cp_%7B2%7D%20%2C...%2Cp_%7BN%7D%29%20%3D0%20%5Cend%7Bcases%7D)

And in the case of a linear system of equations 

![eq6](https://latex.codecogs.com/png.latex?%5Cdpi%7B120%7D%20%5Cbg_white%20%5Clarge%20%5Cbegin%7Bbmatrix%7D%20%5Csquare%20%26%20%5Csquare%20%26%20%5Csquare%20%26%20%5Cvdots%20%26%20%5Csquare%20%5C%5C%20%5Csquare%20%26%20%5Csquare%20%26%20%5Csquare%20%26%20%5Cvdots%20%26%20%5Csquare%20%5C%5C%20%5Csquare%20%26%20%5Csquare%20%26%20%5Csquare%20%26%20%5Cvdots%20%26%20%5Csquare%20%5C%5C%20%5Cdotsc%20%26%20%5Cdotsc%20%26%20%5Cdotsc%20%26%20%5Cddots%20%26%20%5Csquare%20%5C%5C%20%5Csquare%20%26%20%5Csquare%20%26%20%5Csquare%20%26%20%5Csquare%20%26%20%5Csquare%20%5Cend%7Bbmatrix%7D%5Cbegin%7Bbmatrix%7D%20p_%7B1%7D%5C%5C%20p_%7B2%7D%5C%5C%20p_%7B3%7D%5C%5C%20%5Cvdots%20%5C%5C%20p_%7BN%7D%20%5Cend%7Bbmatrix%7D%20%3D%5Cbegin%7Bbmatrix%7D%20i_%7B1%7D%5C%5C%20i_%7B2%7D%5C%5C%20i_%7B3%7D%5C%5C%20%5Cvdots%20%5C%5C%20i_%7BN%7D%20%5Cend%7Bbmatrix%7D)

---
## 3. PyEFVLib Classes - Geometrical Entities

### 3.1 - Grid

| Attribute 	    | Description			                                      |
| :-------------- | :---------------------------------------------------- |
| vertices  	    | a list of the *grid*'s *vertices*	                      |
| elements  	    | a list of the *grid*'s *elements*                       |
| regions  	      | a list of the *grid*'s *regions*                        |
| boundaries  	  | a list of the *grid*'s *boundaries*                     |
| dimension  	    | the dimension of the mesh (2 or 3)                  |
| gridData  	    | a data structure that holds the information provided by mesh file |


---
### 3.2 - Vertex

| Attribute 	  | Description			                                               |
| :------------ | :------------------------------------------------------------- |
| x     		    | *vertex* x coordinate                                          |
| y		          | *vertex* y coordinate		                                       |
| z  	    	    | *vertex* z coordinate	                                         |
| handle  	    | the *vertex* index in the *grid*                               |
| elements      | a list of all *elements* to which *vertex* belongs             |
| volume  	    | the volume of the control volume at which *vertex* is centered |

| Method                  | Parameters            | Return type        | Description                                                                       |
| :---------------------- | :-------------------- | :----------------- | :-------------------------------------------------------------------------------- |
| getLocal                | element               | int                | returns the local index of the *vertex* within the *element*                      |
| getInnerFaces           |                       | list\[InnerFace\]  | returns the innerFaces next to *vertex*                                           |
| getSubElementVolume     | element               | float              | returns the volume of the section of *element* that is inside *vertex* control volume |

> Vertex is a child class of the class Point

---
### 3.3 - Region

Region is a class used to perform simulations where there are multiple regions with different properties, so the properties are specific to each region.
If the domain is homogeneous (has constant properties) one single region will be created containing all the elements. The region names are set in the mesh file.

![region](https://user-images.githubusercontent.com/57679731/118344923-765ae900-b507-11eb-87a6-611a64501889.png)

| Attribute 	      | Description			                             |
| :---------------- | :-------------------------------------       |
| name     		      | the name of the *region*		                 |
| handle  	        | the *region* index in the *grid*             |
| elements          | a list of all elements belonging to *region* |


---
### 3.4 - Element

![element](https://user-images.githubusercontent.com/57679731/118381990-e1b7c000-b5c6-11eb-9661-0a8912f84cc8.png)

| Attribute 	           | Description			                                                                            |
| :--------------------- | :------------------------------------------------------------------------------------------- |
| handle                 | the *element* index in the *grid*                                                            |
| region                 | the *region* to which *element* belongs                                                      |
| vertices               | a list of the *vertices* of *element*                                                      |
| innerFaces             | a list of the *element*'s *innerFaces*                                                     |
| outerFaces             | a list of the *element*'s *outerFaces*                                                     |
| faces                  | a list of the *element*'s *innerFaces* and *outerFaces*                                    |
| shape                  | the *element* *Shape*                                                                        |
| volume                 | the sum of all subElementVolumes (control volume section inside *element*)                   |
| subElementVolumes      | the list with each *vertex* subElementVolume                                                 |

| Method                  | Parameters               | Return type        | Description                                                                       |
| :---------------------- | :----------------------- | :----------------- | :-------------------------------------------------------------------------------- |
| getTransposedJacobian   | shapeFunctionDerivatives | numpy.array\[d,m\] | returns the transposed of the jacobian of shapeFunctionDerivatives                |
| getGlobalDerivatives    | derivatives              | numpy.array\[d,d\] | returns the derivatives of the *vertices* shape functions with respect to x,y,z   |

> Where d is the grid dimension and m is the number of *vertices* in the *element*


---
### 3.5 - InnerFace

The innerFaces are inside the element, but they surround the control volumes. The reason we need the shape functions and its derivatives evaluated at the inner faces centroids is because whenever there's a surface integral in our differential equation (over the control surface), it'll be summed over "integration points" over the control surface, which happen to be at the centroids of the inner faces.

| Attribute 	           | Description			                                                                            |
| :--------------------- | :------------------------------------------------------------------------------------------- |
| handle                 | the *innerFace* index in the *grid*                                                          |
| local                  | the *innerFace* local index in the *element*                                                 |
| element                | the *innerFace* element                                                                      |
| area                   | the *innerFace* area (*Point* object)                                                    |
| centroid               | the *innerFace* centroid (*Point* object)                                                    |
| globalDerivatives      | a numpy.array\[d,m\] containing the gradient operator. Multiply by a m-dimensional array containing the property at the vertices of its element to get the gradient of the property at the *innerFace*                                                          |

| Method                     | Parameters  | Return type        | Description                                                                       |
| :------------------------- | :---------: | :----------------- | :-------------------------------------------------------------------------------- |
| getShapeFunctions          |      -      | numpy.array\[d,m\] | The shape functions of the *element*'s *vertices*' shapeFunctions evaluated at the *innerFace* |
| getNeighborVertices        |      -      | (Vertex,Vertex)    | The surrounding *vertices* of *innerFace*                                         |
| getNeighborVerticesLocals  |      -      | (int, int)         | The surrounding *vertices*' local indices in the *element* of *innerFace*         |
| getNeighborVerticesHandles |      -      | (int, int)         | The surrounding *vertices*' handles of *innerFace*                                |


---
### 3.6 - Boundary

| Attribute 	           | Description			                                                                            |
| :--------------------- | :------------------------------------------------------------------------------------------- |
| name                   | it used to apply the *boundary conditions* on it                                             |
| handle                 | the *boundary* index on the grid                                                             |
| vertices               | a list with the *vertices* that make the *boundary*                                          |
| facets                 | a list with the *facets* that make the *boundary*                                            |


---
### 3.7 - Facet

![facet](https://user-images.githubusercontent.com/57679731/118346725-d5265f80-b513-11eb-8ae0-d99ba40b0b33.png)

| Attribute 	           | Description			                                                                            |
| :--------------------- | :------------------------------------------------------------------------------------------- |
| handle                 | the *facet* index on the grid                                                                |
| boundary               | the *boundary* to which *facet* belongs                                                      |
| element                | the *element* to which *facet* belongs                                                       |
| vertices               | the *vertices* that lie on facet                                                             |
| outerFaces             | the *outerFaces* that make *facet*                                                           |
| area                   | *facet*'s area (*Point* object)                                                              |
| centroid               | *facet*'s centroid (*Point* object)                                                          |
| elementLocalIndex      | its local index on *element*                                                                 |
| boundaryLocalIndex     | its local index on *boundary*                                                                |


---
### 3.8 - OuterFace

| Attribute 	           | Description			                                                                          |
| :--------------------- | :------------------------------------------------------------------------------------------|
| local                  | local index on *facet*                                                                     |
| vertex                 | neighboring *vertex*                                                                       |
| facet                  | *facet* to which *outerFace* belongs                                                       |
| centroid               | its centroid (*Point* object)                                                              |
| area                   | its area (*Point* object)                                                                  |

| Method            | Parameters  | Return type        | Description                                                                                    |
| :-----------------| :---------: | :----------------- | :--------------------------------------------------------------------------------------------- |
| getShapeFunctions |      -      | numpy.array\[d,m\] | The shape functions of the *element*'s *vertices*' shapeFunctions evaluated at the *outerFace* |


---
### 3.9 - Shape

The Shape subclasses are:
- Triangle
- Quadrilateral
- Tetrahedron
- Hexahedron
- Prism
- Pyramid

All of them share the same attributes and methods:

| Attribute 	                         | Description			                                                                            |
| :----------------------------------- | :------------------------------------------------------------------------------------------- |
| dimension                            | the *shape* dimension                                                                        |
| numberOfInnerFaces                   | the number of *inner faces*                                                                  |
| numberOfFacets                       | the number of *facets*                                                                       |
| subelementTransformedVolumes         | the volumes of the subelements in the transformed space                                      |
| innerFaceShapeFunctionValues         | the shape function values evaluated at the *innerFace* centroid                              |
| innerFaceShapeFunctionDerivatives    | the shape function derivatives evaluated at the *innerFace* centroid                         |
| innerFaceNeighborVertices            | the *element* local indices of the *vertices* that surround each *innerFace*                 |
| subelementShapeFunctionValues        | the shape function values evaluated at the subelement centroid                               |
| subelementShapeFunctionDerivatives   | the shape function derivatives evaluated at the subelement centroid                          |
| facetVerticesIndexes                 | the *element* local indices of the group of *vertices* that make each *facet*                |
| outerFaceShapeFunctionValues         | the shape function values evaluated at the *outerFace* centroid                              |
| vertexShapeFunctionDerivatives       | the shape function derivatives evaluated at the *vertices*                                   |


---
### 3.10 - BoundaryConditions
The two BoundaryConditions classes are
- DirichletBoundaryCondition
- NeumannBoundaryCondition

| Attribute 	           | Description			                                                                            |
| :--------------------- | :------------------------------------------------------------------------------------------- |
| handle                 | the *boundaryCondition* index on the *grid*                                                  |
| value                  | the prescribed condition value                                                               |
| expression             | True if *value* is a string and the expression must be parsed, otherwise False               |
| boundary               | the *boundary* to which the condition is being applied                                       |

| Method   | Parameters            | Return type | Description                                                                       |
| :------- | :-------------------- | :---------- | :-------------------------------------------------------------------------------- |
| getValue | index, time=0.0       | float       | if *expression*==False, returns *value*. Otherwise, will parse the *value* string |

If `expression==True`, the `value` string is an expression that will be parsed using the buit-in method `eval`.
It may contain the variables `x`, `y` and `z`, which are equal respectively to `grid.vertices[index].x`, `grid.vertices[index].y` and `grid.vertices[index].z`, as well as `t`, which is equal to the argument *time*.
The following functions/constants can also be used: `pi, sin, cos, tan, arcsin, arccos, arctan, sinh, cosh, tanh, arcsinh, arccosh, arctanh, sqrt, e , log, exp, inf, mod, floor` (which have been imported from the module [numpy](https://numpy.org/)).

The Dirichlet boundary condition means that the unknown variable is being prescribed at the boundary, and the Neumann boundary condition means that the unknown variable directional gradient (in the direction of the normal vector to the boundary) is being prescribed.


---
## 4. PyEFVLib Classes - Non Geometrical Entities

Besides grid entities generation that facilitates the solution of PDEs, PyEFVLib also offers some I/O classes, some of which are integrated with other I/O libraries such as [meshio](https://github.com/nschloe/meshio).

### 4.1 - ProblemData

The ProblemData class purpose is to provide the solver every information needed to solve the problem, in such a way that one can write a code in the following way:
```python
def solver(problemData):
  ...
  
if __name__ == "__main__":
  problemData = PyEFVLib.ProblemData(...)
  solution = solver(problemData)
```

This way in order to keep everything simple ProblemData only needs 5 arguments to initialize

| Arguments    	      | Arg Type   		     | Description                                                                    |
| :------------------ | :----------------- | :----------------------------------------------------------------------------- |
| meshFilePath        | string             | provides the path to the input mesh                                            |
| outputFilePath      | string			       | provides the directory where the results will be written                       |
| numericalSettings   | NumericalSettings  | provides all the numerical method details                                      |
| propertyData        | PropertyData       | provides all the property info of the domain                                   |
| boundaryConditions  | BoundaryConditions | provides all the boundary conditions of the problem                            |

- **4.1.1 - NumericalSettings**

| Arguments    	        | Arg Type | Default  | Description                                                          |
| :------------------   | :------- | :------: | -------------------------------------------------------------------- |
| timeStep              | float    |    -     | step in time to be increased                                          |
| finalTime             | float  	 |   inf    | when the simulation reaches this time it'll stop                      |
| tolerance             | float    |   1e-4   | can be a steady state criteria or an iterative tolerance (it is up to who writes the solver and can be used in both ways)   |
| maxNumberOfIterations | int      |   1000   | when the simulation reach this number of iterations it'll stop        |

- **4.1.2 - PropertyData**

Usage:
```python
numericalSettings = PyEFVLib.NumericalSettings({
  "RegionName1": {
    "property1": value1,
    "property2": value2,
    "property3": value3,
  },
  "RegionName2": {
    "property1": value4,
    "property2": value5,
    "property3": value6,
  },
})
```

- **4.1.3 - BoundaryConditions**

Usage:
```python
boundaryConditions = PyEFVLib.BoundaryConditions({
  "variable1": {
    "InitialValue": initialValue1,
    "boundaryName1": { "condition": PyEFVLib.Dirichlet, "type": PyEFVLib.Constant, "value": 0.0 },
    "boundaryName2": { "condition": PyEFVLib.Dirichlet, "type": PyEFVLib.Variable, "value": "100*sin(x*t)" },
    "boundaryName3": { "condition": PyEFVLib.Neumann, "type": PyEFVLib.Constant, "value": 100.0 },
  },
  "variable2": {
    "InitialValue": initialValue2,
    "boundaryName1": { "condition": PyEFVLib.Neumann, "type": PyEFVLib.Variable, "value": "0.0 if y<0.5 else 50.0" },
    "boundaryName2": { "condition": PyEFVLib.Dirichlet, "type": PyEFVLib.Constant, "value": -50.0 },
    "boundaryName3": { "condition": PyEFVLib.Neumann, "type": PyEFVLib.Constant, "value": 0.0 },
  },
})
```

---
### 4.2 - Readers
The current readers are:
- MSHReader
- XDMFReader

(Soon - MeshioReader)

The only argument that the readers take is the mesh path.
To initialize a grid object from a reader do as such:
```python
grid = PyEFVLib.Grid( MSHReader(filePath).getData() )
```
or
```python
grid = PyEFVLib.read(filePath)
```

> To create your own reader, all it is needed is to end your reader by creating a GridData object (which only takes the mesh file path as argument), and populate it by using the methods: setDimension, setVertices, setElementConnectivity, setRegions, setBoundaries, setShapes

---
### 4.3 - Savers
The current savers are:
- CsvSaver (saves the results in a .csv spreadsheet-like file)
- VtuSaver (saves the results in only one moment in a .vtu file)
- VtmSaver (saves the results in multiple .vtm files from all time steps)
- MeshioSaver (mostly saves the results in one moment, except by xdmf) \[recommended\]

> MeshioSaver supports the following extensions: msh, mdpa, ply, stl, vtk, vtu, xdmf, xmf, h5m, med, inp, mesh, meshb, bdf, fem, nas, obj, off, post, post.gz, dato, dato.gz, su2, svg, dat, tec, ugrid, wkt. Works best with xdmf

| Arguments    	      | Arg Type   		     | Description                                                                    |
| :------------------ | :----------------- | :----------------------------------------------------------------------------- |
| grid                | Grid               | the grid object                                                                |
| outputFilePath      | str    		         | provided by problemData.outputFilePath                                         |
| libraryPath         | str                | provided by problemData.libraryPath                                            |
| extension           | str                | exclusive of MeshioSaver, insert one of the listed above                       |

| Method                  | Parameters                       | Description                                                                       |
| :---------------------- | :------------------------------- | :-------------------------------------------------------------------------------- |
| save                    | variableName, field, currentTime | saves the field                                                                   |
| finalize                | -                                | MUST be called at the end of the simulation in order for the results to be saved  |

---
## 5. Tutorial
In order to demonstrate the library classes we'll sove the heat transfer equation. Let's remember the equations:
![eq7](https://latex.codecogs.com/gif.latex?%5Cdpi%7B120%7D%20%5Cbg_white%20%5Clarge%20%5Cnabla%20%5Ccdotp%20%28%20k%5Cnabla%20T%29%20&plus;q%27%27%27%3D%5Crho%20c_%7Bp%7D%5Ctfrac%7B%5Cpartial%20T%7D%7B%5Cpartial%20t%7D%5C%5C%20%5Cint%20_%7B%5COmega%20_%7Bi%7D%7D%5B%20%5Cnabla%20%5Ccdotp%20%28%20k%5Cnabla%20T%29%20&plus;q%27%27%27%5D%20%5C%20d%5COmega%20_%7Bi%7D%20%3D%5Cint%20_%7B%5COmega%20_%7Bi%7D%7D%20%5Crho%20c_%7Bp%7D%5Ctfrac%7B%5Cpartial%20T%7D%7B%5Cpartial%20t%7D%20d%5COmega%20_%7Bi%7D%5C%5C%20%5Cint%20_%7B%5CGamma%20_%7Bi%7D%7D%28%20k%5Cnabla%20T%29%20%5Ccdotp%20%5Chat%7Bn%7D%20%5C%20d%5CGamma%20_%7Bi%7D%20&plus;%5Cint%20_%7B%5COmega%20_%7Bi%7D%7D%20q%27%27%27%5C%20d%5COmega%20_%7Bi%7D%20%3D%5Cint%20_%7B%5COmega%20_%7Bi%7D%7D%20%5Crho%20c_%7Bp%7D%5Ctfrac%7B%5Cpartial%20T%7D%7B%5Cpartial%20t%7D%20d%5COmega%20_%7Bi%7D%5C%5C%20%5Csum%20_%7Bf%5C%20%5Cin%20%5C%20%5CGamma%20_%7Bi%7D%7D%5Csum%20_%7Bip%5C%20%5Cin%20%5C%20f%7D%28%20k%5Cnabla%20T%29_%7Bip%7D%20%5Ccdotp%20%5Coverrightarrow%7B%5CDelta%20s_%7Bip%7D%7D%20&plus;q%27%27%27%5CDelta%20%5COmega%20_%7Bi%7D%20%3D%5Crho%20c_%7Bp%7D%5Ctfrac%7B%5Cpartial%20T%7D%7B%5Cpartial%20t%7D%20%5CDelta%20%5COmega%20_%7Bi%7D)

![eq8](https://latex.codecogs.com/gif.latex?%5Cdpi%7B120%7D%20%5Cbg_white%20%5Clarge%20%5Cint%5Climits%20_%7Bt%7D%5E%7Bt&plus;%5CDelta%20t%7D%5Cleft%5B%5Csum%20_%7Bf%5C%20%5Cin%20%5C%20%5CGamma%20_%7Bi%7D%7D%5Csum%20_%7Bip%5C%20%5Cin%20%5C%20f%7D%28%20k%5Cnabla%20T%29_%7Bip%7D%20%5Ccdotp%20%5Coverrightarrow%7B%5CDelta%20s_%7Bip%7D%7D%5Cright%5D%20dt&plus;%5Cint%5Climits%20_%7Bt%7D%5E%7Bt&plus;%5CDelta%20t%7D%20q%27%27%27%5CDelta%20%5COmega%20_%7Bi%7D%20%5C%20dt%3D%5Cint%5Climits%20_%7Bt%7D%5E%7Bt&plus;%5CDelta%20t%7D%20%5Crho%20c_%7Bp%7D%5Ctfrac%7B%5Cpartial%20T%7D%7B%5Cpartial%20t%7D%20%5CDelta%20%5COmega%20_%7Bi%7D%20%5C%20dt%5C%5C%20%5Csum%20_%7Bf%5C%20%5Cin%20%5C%20%5CGamma%20_%7Bi%7D%7D%5Csum%20_%7Bip%5C%20%5Cin%20%5C%20f%7D%28%20k%5Cnabla%20T%29_%7Bip%7D%20%5Ccdotp%20%5Coverrightarrow%7B%5CDelta%20s_%7Bip%7D%7D%20%5C%20%5CDelta%20t&plus;q%27%27%27%5CDelta%20%5COmega%20_%7Bi%7D%20%5CDelta%20t%3D%5Crho%20c_%7Bp%7D%5Cleft%28%20T-T%5E%7Bold%7D%5Cright%29%20%5CDelta%20%5COmega%20_%7Bi%7D%5C%5C%20%5Crho%20c_%7Bp%7D%5Ctfrac%7B1%7D%7B%5CDelta%20t%7D%20T%5CDelta%20%5COmega%20_%7Bi%7D%20-%5Csum%20_%7Bf%5C%20%5Cin%20%5C%20%5CGamma%20_%7Bi%7D%7D%5Csum%20_%7Bip%5C%20%5Cin%20%5C%20f%7D%28%20k%5Cnabla%20T%29_%7Bip%7D%20%5Ccdotp%20%5Coverrightarrow%7B%5CDelta%20s_%7Bip%7D%7D%20%3Dq%27%27%27%5CDelta%20%5COmega%20_%7Bi%7D%20&plus;%5Crho%20c_%7Bp%7D%5Ctfrac%7B1%7D%7B%5CDelta%20t%7D%20T%5E%7Bold%7D%20%5CDelta%20%5COmega%20_%7Bi%7D%5C%5C%20%5Crho%20c_%7Bp%7D%5Ctfrac%7B1%7D%7B%5CDelta%20t%7D%20T%5CDelta%20%5COmega%20_%7Bi%7D%20-%5Csum%20_%7Bf%5C%20%5Cin%20%5C%20%5CGamma%20_%7Bi%7D%7D%5Csum%20_%7Bip%5C%20%5Cin%20%5C%20f%7D%20k%5Coverrightarrow%7B%5CDelta%20s_%7Bip%7D%5E%7BT%7D%7D%5Cmathbf%7BB%7D_%7Bip%7D%5Coverrightarrow%7BT%5E%7Belem%7D%7D%20%3Dq%27%27%27%5CDelta%20%5COmega%20_%7Bi%7D%20&plus;%5Crho%20c_%7Bp%7D%5Ctfrac%7B1%7D%7B%5CDelta%20t%7D%20T%5E%7Bold%7D%20%5CDelta%20%5COmega%20_%7Bi%7D)

We'll start by importing the library and initializing some variables

```python
import PyEFVLib
import numpy as np

def heatTransfer(problemData):
    propertyData     = problemData.propertyData
    timeStep         = problemData.timeStep
    grid             = problemData.grid
    numberOfVertices = problemData.grid.numberOfVertices
    dimension        = problemData.grid.dimension
  
    saver = PyEFVLib.MeshioSaver(grid, problemData.outputFilePath, problemData.libraryPath, extension="xdmf")
  
    temperatureField = np.repeat(0.0, numberOfVertices)
    oldTemperatureField = problemData.initialValues["temperature"].copy()
```

The first 5 lines where we call some problemData variables are optional, but since they are used very often, it is handy to declare them like that. The saver must be initialized before starting to save the results, so here is a good place to do it. The temperatureField is our unknown to be found, and the oldTemperatureField will be updated after each increment in time, and will be used to evaluate the temporal derivative.

Notice that in our discretized equation the left hand side contains only terms involving the temperature (which are going to be added to the matrix) and the right hand side contains only terms independent of the temperature(and are going to be added to the independent vector). The linear system matrix will have N rows and N columns, where the i-th row corresponds to the equation applied to the i-th control volume (or vertex), and the columns correspond to the temperature in the i-th vertex. So, to implement the matrix of our linear system we'll do the following:

```python
    ...
    oldTemperatureField = problemData.initialValues["temperature"].copy()
    
    def assembleMatrix():
        matrix = np.zeros((numberOfVertices, numberOfVertices))
    
        for region in grid.regions:
            k     = propertyData.get(region.handle, "conductivity")
            q     = propertyData.get(region.handle, "heatGeneration")
            rho = propertyData.get(region.handle, "density")
            cp    = propertyData.get(region.handle, "specificHeat")
            
            for element in region.elements:
                # Accumulation term (rho * cp * T/Δt)
                for vertex in element.vertices:
                    matrix[vertex.handle][vertex.handle] += vertex.getSubElementVolume(element) * rho * cp / timeStep
                    
                # Diffusion term ( -k*grad(T) Δs <-> -k * Δs^T * B_ip * T_elem )
                for innerFace in element.innerFaces:
                    area = innerFace.area.getCoordinates()[:dimension]
                    backwardsVertexHandle, forwardVertexHandle = innerFace.getNeighborVerticesHandles()
                    
                    coefficients = -k * np.matmul(area.T, innerFace.globalDerivatives)
                    for local, vertex in enumerate(element.vertices):
                        matrix[backwardsVertexHandle][vertex.handle] += coefficients[local]
                        matrix[forwardVertexHandle][vertex.handle] -= coefficients[local]

```

Remember that the `region` object contains the elements with the same properties, that's why we start the `assembleMatrix` function looipng over the regions and calling the properties at that region. Then we loop over the elements, since the innerFaces belong inside the elements, and for the vertices loop we need the subElementVolume.

The ![term1](https://latex.codecogs.com/gif.latex?%5Cbg_white%20%5Crho%20c_%7Bp%7D%5Ctfrac%7B1%7D%7B%5CDelta%20t%7D%20T%5CDelta%20%5COmega%20_%7Bi%7D) term comes from the i-th equation relating the i-th vertex temperature, that's why we add the term to `matrix[vertex.handle][vertex.handle]`

Also note that from our surface integral over the controle surface, we got the summation over the centroids of the innerFaces (or integration points, denoted in the equations as *ip*). So the loop over the innerFaces of element will add this term to our linear system of equations. The ![Δs_ip](https://latex.codecogs.com/gif.latex?%5Cbg_white%20%5Coverrightarrow%7B%5CDelta%20s_%7Bip%7D%7D) (or `innerFace.area`) is the area vector of the innerFaces, which has modulus equal to the area of the innerFace and points in the outwards normal direction to the innerFace. So in our summation we add `coefficient` to backwardsVertex because `innerFace.area` is pointing outwards from its control volume, and we subtract it from forwardVertex because `-innerFace.area` is poining outwards from its control volume.

Since `innerFace.area` is a Point object (and not a numpy.array) we convert it into a 3x1 numpy.array using the `getCoordinates` method, and because `innerFace.globalDerivatives` has `dimension` rows and `element.vertices.size` columns, we slice it into a `dimension`x1 numpy.array (removing the z coordinate if the simulation is 2D).

The ![B_ip](https://latex.codecogs.com/gif.latex?%5Cbg_white%20%5Cmathbf%7BB%7D_%7Bip%7D) term is called gradient operator and it is stored in `innerFace.globalDerivatives`. It's defined as
![B_ip definition](https://latex.codecogs.com/gif.latex?%5Cbg_white%20%5Cmathbf%7BB%7D_%7Bip%7D%20%3D%5Cbegin%7Bbmatrix%7D%20%5Ctfrac%7B%5Cpartial%20%7D%7B%5Cpartial%20x%7D%5Cmathcal%7BN%7D_%7B1%7D%28%20x_%7Bip%7D%20%2Cy_%7Bip%7D%29%20%26%20%5Ctfrac%7B%5Cpartial%20%7D%7B%5Cpartial%20x%7D%5Cmathcal%7BN%7D_%7B2%7D%28%20x_%7Bip%7D%20%2Cy_%7Bip%7D%29%20%26%20%5Ctfrac%7B%5Cpartial%20%7D%7B%5Cpartial%20x%7D%5Cmathcal%7BN%7D_%7B3%7D%28%20x_%7Bip%7D%20%2Cy_%7Bip%7D%29%5C%5C%20%5Ctfrac%7B%5Cpartial%20%7D%7B%5Cpartial%20y%7D%5Cmathcal%7BN%7D_%7B1%7D%28%20x_%7Bip%7D%20%2Cy_%7Bip%7D%29%20%26%20%5Ctfrac%7B%5Cpartial%20%7D%7B%5Cpartial%20y%7D%5Cmathcal%7BN%7D_%7B2%7D%28%20x_%7Bip%7D%20%2Cy_%7Bip%7D%29%20%26%20%5Ctfrac%7B%5Cpartial%20%7D%7B%5Cpartial%20y%7D%5Cmathcal%7BN%7D_%7B3%7D%28%20x_%7Bip%7D%20%2Cy_%7Bip%7D%29%20%5Cend%7Bbmatrix%7D)

and it has the following property:
![B_ip property](https://latex.codecogs.com/gif.latex?%5Cbg_white%20%5Cnabla%20p%28%20x_%7Br%7D%20%2Cy_%7Br%7D%29%20%3D%5Cmathbf%7BB%7D_%7Br%7D%5Cbegin%7Bbmatrix%7D%20p_%7B1%7D%20%26%20p_%7B2%7D%20%26%20p_%7B3%7D%20%5Cend%7Bbmatrix%7D%5E%7BT%7D).

<br/><br/>
Now we need to turn our attention to the boundary conditions. Previously we simplified our equation by discretizing it only in the interior of the domain, without making considerations about its boundaries. Note that we can express Control Surface = Control Surface (inside the domain) + Control Surface (where the Neumann boundary condition is applied) + Control Surface (where the Dirichlet boundary condition is applied).
![eq9](https://latex.codecogs.com/gif.latex?%5Cdpi%7B120%7D%20%5Cbg_white%20%5Clarge%20%5CGamma%20_%7Bi%7D%20%3D%5CGamma%20_%7Bi%7D%20%5Ccup%20%5CGamma%20_%7Bi%2C_%7BNeumann%7D%7D%20%5Ccup%20%5CGamma%20_%7Bi%2C_%7BDirichlet%7D%7D)

So in the ![eq10](https://latex.codecogs.com/gif.latex?%5Cbg_white%20%5Cint%20_%7B%5CGamma%20_%7Bi%7D%7D%28%20k%5Cnabla%20T%29%20%5Ccdotp%20%5Chat%7Bn%7D%20%5C%20d%5CGamma%20_%7Bi%7D%20&plus;%5Cint%20_%7B%5COmega%20_%7Bi%7D%7D%20q%27%27%27%5C%20d%5COmega%20_%7Bi%7D%20%3D%5Cint%20_%7B%5COmega%20_%7Bi%7D%7D%20%5Crho%20c_%7Bp%7D%5Ctfrac%7B%5Cpartial%20T%7D%7B%5Cpartial%20t%7D%20d%5COmega%20_%7Bi%7D) step we could have made the following substitution:

![eq11](https://latex.codecogs.com/gif.latex?%5Cbg_white%20%5Cint%20_%7B%5CGamma%20_%7Bi%2C%7D%7D%28%20k%5Cnabla%20T%29%20%5Ccdotp%20%5Chat%7Bn%7D%20%5C%20d%5CGamma%20_%7Bi%7D%20%3D%5Cint%20_%7B%5CGamma%20_%7Bi%2C%5Ctext%7Binner%7D%7D%7D%28%20k%5Cnabla%20T%29%20%5Ccdotp%20%5Chat%7Bn%7D%20%5C%20d%5CGamma%20_%7Bi%7D%20&plus;%5Cint%20_%7B%5CGamma%20_%7Bi%2C%5Ctext%7BNeumann%7D%7D%7D%28%20k%5Cnabla%20T%29%20%5Ccdotp%20%5Chat%7Bn%7D%20%5C%20d%5CGamma%20_%7Bi%7D%20&plus;%5Cint%20_%7B%5CGamma%20_%7Bi%2C%5Ctext%7BDirichlet%7D%7D%7D%28%20k%5Cnabla%20T%29%20%5Ccdotp%20%5Chat%7Bn%7D%20%5C%20d%5CGamma%20_%7Bi%7D)

Where

![eq12](https://latex.codecogs.com/gif.latex?%5Cbg_white%20%5Cint%20_%7B%5CGamma%20_%7Bi%2C%5Ctext%7Binner%7D%7D%7D%28%20k%5Cnabla%20T%29%20%5Ccdotp%20%5Chat%7Bn%7D%20%5C%20d%5CGamma%20_%7Bi%7D%20%3D%5Csum%20_%7Bf%5C%20%5Cin%20%5C%20%5CGamma%20_%7Bi%7D%7D%5Csum%20_%7Bip%5C%20%5Cin%20%5C%20f%7D%28%20k%5Cnabla%20T%29_%7Bip%7D%20%5Ccdotp%20%5Coverrightarrow%7B%5CDelta%20s_%7Bip%7D%7D)

And

![eq13](https://latex.codecogs.com/gif.latex?%5Cbg_white%20%5Cint%20_%7B%5CGamma%20_%7Bi%2C%5Ctext%7BNeumann%7D%7D%7D%28%20k%5Cnabla%20T%29%20%5Ccdotp%20%5Chat%7Bn%7D%20%5C%20d%5CGamma%20_%7Bi%7D%20%3D%5Coverrightarrow%7Bq%27%27_%7B_%7B%5Ctext%7BNeumann%7D%7D%7D%7D%20%5Ccdotp%20%5Coverrightarrow%7B%5CDelta%20s_%7Bip%7D%7D)

The dirichlet boundary condition is very simple, if it is applied to the i-th vertex, all the elements on the i-th row of the matrix must be 0, except for the i-th element of the i-th row, which is equal to 1. Meanwhile the i-th element of the independent vector must be the prescribed temperature.

![eq14](https://latex.codecogs.com/gif.latex?%5Cbg_white%20T_%7B1%7D%20%5Ccdotp%200&plus;T_%7B2%7D%20%5Ccdotp%200&plus;...&plus;T_%7Bi%7D%20%5Ccdotp%201&plus;...&plus;T_%7BN%7D%20%5Ccdotp%200%3DT_%7Bi%2C%5Ctext%7Bprescribed%7D%7D)

That being said, the way we implement the boundary conditions in our matrix is the following

(still in the `assembleMatrix` function)

```python
...
def heatTransfer(problemData):
    ...
    def assembleMatrix():
        ...
                        matrix[forwardVertexHandle][vertex.handle] += k * np.matmul(area.T, innerFace.globalDerivatives)
                        
        # Dirichlet Boundary Condition
        for boundaryCondition in problemData.dirichletBoundaries["temperature"]:
            for vertex in boundaryCondition.boundary.vertices:
                matrix[vertex.handle] = np.zeros(numberOfVertices)
                matrix[vertex.handle][vertex.handle] = 1.0
        
        return matrix
```

Now we need an `assembleIndependentVector` function:

```python
...
def heatTransfer(problemData):
    ...
    def assembleIndependentVector():
        independent = np.zeros(numberOfVertices)
    
        for region in grid.regions:
            k   = propertyData.get(region.handle, "conductivity")
            q   = propertyData.get(region.handle, "heatGeneration")
            rho = propertyData.get(region.handle, "density")
            cp  = propertyData.get(region.handle, "specificHeat")
            
            for element in region.elements:
                for vertex in element.vertices:
                    # Accumulation term (rho * cp * T_old/Δt)
                    independent[vertex.handle] += vertex.getSubElementVolume(element) * rho * cp * oldTemperatureField[vertex.handle] / timeStep
                    
                    # Generation term (q''')
                    independent[vertex.handle] += vertex.getSubElementVolume(element) * q
                    
            # Neumann Boundary Condition
            for boundaryCondition in problemData.neumannBoundaries["temperature"]:
                for facet in boundaryCondition.boundary.facets:
                    for outerFace in facet.outerFaces:
                        independent[outerFace.vertex.handle] += boundaryCondition.getValue(outerFace.vertex.handle, currentTime) * np.linalg.norm(outerFace.area.getCoordinates())
                        
            # Dirichlet Boundary Condition
            for boundaryCondition in problemData.dirichletBoundaries["temperature"]:
                for vertex in boundaryCondition.boundary.vertices:
                    independent[vertex.handle] = boundaryCondition.getValue(vertex.handle, currentTime)
                    
        return independent
```

Because the `assembleMatrix` function and the boundary condition details have already been explained, this function is more straightforward.

The accumulation term is added to the independent vector, and it calls the `oldTemperatureField` vector, which must be updated after each step in time.

Both the Neumann boundary condition and the Dirichlet boundary condition are applied to the facets of the boundaries, which are in a way different from the vertices. However, since the Dirichlet boundary condition prescribes directly the temperature (or whichever property) at the vertices of the boundary, the temperature will be set at those points, that's why it is applied after the Neumann boundary condition, because otherwise the temperature would be mistakenly set.

On the `boundaryCondition.getValue(outerFace.vertex.handle)`, remember that if the boundaryCondition is of the type `PyEFVLib.Constant`, it will return the prescribed boundaryCondition anyways, however if it is of the type `PyEFVLib.Variable`, it can depend of the x,y and z coordinates of the vertex and also of the `currentTime`. So we must pass the handle of the `outerFace.vertex` so it can check the x,y,z coordinates of it.

<br/><br/>
Next we start the loop to solve our problem:

```python
    steadyStateCriterion = problemData.tolerance
    difference = 2*steadyStateCriterion
    stepInTime = 0
    currentTime = 0.0
    convergedToSteadyState = False
    
    matrix = assembleMatrix()
    while not convergedToSteadyState:
        independent = assembleIndependentVector()
        temperatureField = np.linalg.solve(matrix, independent)
        difference = max( abs(temperatureField - oldTemperatureField) )
        
        oldTemperatureField = temperatureField.copy()
        currentTime += timeStep
        stepInTime += 1
        
        saver.save("temperatureField", temperatureField, currentTime)
        print("{:>9}        {:>14.2e}     {:>14.2e}     {:>14.2e}".format(stepInTime, currentTime, timeStep, difference))
        
        convergedToSteadyState = ( difference <= steadyStateCriterion ) or ( currentTime >= problemData.finalTime ) or ( stepInTime >= problemData.maxNumberOfIterations )
            
    saver.finalize()
```

Here we declare some handy variables at the start of the loop, and notice that we assemble the matrix before starting the loop. That's because none of the terms added in matrix change over time. If we wanted for example to have a variable timeStep, then we'd have to assemble the matrix inside the transient loop every instant of time. In this example, other detail could have been changed for the simulation to perform better. Instead of solving the linear system of equations every single step in time, we could have inverted the matrix before the transient loop (right after assembling the matrix) and simply multiplied the inverted matrix by the independent vector. Here it has not been done for didactic purposes.

The `difference` variable is meant to track how the `temperatureField` stop changing after a while, and to stop the simulation after the changes are negligible. The `.copy()` method is needed to update `oldTemperatureField` because of how the mutable list class behaves. After solving the problem it is needed to call the `saver.finalize()` method, to write the output file with the results.

<br/><br/>
Finally, we need to call our `heatTransfer` function, like so:

```python
if __name__ == "__main__":
    problemData = PyEFVLib.ProblemData(
        meshFilePath = "mesh.msh",
        outputFilePath = "results/",
        numericalSettings = PyEFVLib.NumericalSettings(timeStep = 1000.0, tolerance = 1e-4, maxNumberOfIterations = 500),
        propertyData = PyEFVLib.PropertyData({
            "Body": {
                "conductivity"   : 22.0,   # [W/m.K]
                "density"        : 8960.0, # [kg/m³]
                "specificHeat"   : 377.0,  # [J/kg.K]
                "heatGeneration" : 5000.0, # [W/m³]
            },
        }),
        boundaryConditions = PyEFVLib.BoundaryConditions({
            "temperature": {
                "InitialValue": 300.0,
                "West": { "condition" : PyEFVLib.Dirichlet, "type" : PyEFVLib.Constant, "value" : 300.0 },
                "East": { "condition" : PyEFVLib.Dirichlet, "type" : PyEFVLib.Constant, "value" : 400.0 },
                "South": { "condition" : PyEFVLib.Neumann, "type" : PyEFVLib.Constant, "value" : 0.0 },
                "North": { "condition" : PyEFVLib.Neumann, "type" : PyEFVLib.Constant, "value" : 0.0 },                
            }
        }),
    )
    
    heatTransfer(problemData)
```

---
## 6. 3rd-party Softwares
### 6.1 - Gmsh
Notice that in the workflow presented, the mesh is simply imported by passing its path, but it need to be generated. In order to do so, the recommended software is [Gmsh](https://gmsh.info/), which is an open-source mesh generator with built-in CAD tools. The only (and minor) downside of the mesh importing at the moment is that the mesh must be from the 2.2 version of .msh. But that's very simple. Just download the latest version of Gmsh, follow the step-by-step:

1. Open the Help menu
2. Choose Current Options and Workspace
3. Search 'version' in the search bar, or scroll down until you find the parameter 'Mesh.MshFileVersion'
4. Double click it, type 2.2, and click Ok
5. Optionally open the File menu, and choose the option 'Save Options as Default', so you don't have to go through this step-by-step every time

![gmsh tut](https://user-images.githubusercontent.com/57679731/118505884-352f2880-b703-11eb-9676-270d726f7fc5.png)


### 6.2 - Paraview
After solving the simulation, you'll want to view and analyize the results you generated. One option is to save it as a .csv file (using the CsvSaver) and open it using [pandas](https://pandas.pydata.org/), and you could visualize it using [matplotlib](https://matplotlib.org/). But if you want to view your results, a much better way of doing this is with [Paraview](https://www.paraview.org/).

Paraview is an open-source data analysis and visualization application, and by saving the results in formats such as .xdmf (recommended), .vtu, .vtm, you can visualize it using paraview.

---
## 7. 2D vs 3D considerations
One great advantage of using the EbFVM, is that the difference between 2D and 3D simulations is very little. The major changes are in the field and gradient aproximations using different elements and therefore different shape functions. However PyEFVLib handles it without the user having to pay attention. For example in our heat transfer application show as an example in section 5, the only change one need to make to simulate a 3D domain is simply to change the mesh file path (and generate a 3D mesh), and add the boundary conditions according to the new mesh.

---
## 8. bellbird
Although PyEFVLib makes the use of EbFVM quite simple for the user, it is still required for the user to know the numerical method, discretize the equations, and when solving bigger systems of differential equations the coupling of the equation in the linear system can become a little tricky, and can demand more attention to little details, while also requiring a lot more lines of code. For that reason, the [bellbird](https://github.com/GustavoExel/bellbird) library makes it for the user way more simple to solve systems of PDEs. The user informs the equations as string, and the library integrates them, discretizes them, and writes a script very similar to the shown in the tutorial using PyEFVLib.

Then, if the user wants, they can change the code, but with a huge headstart.

At the moment, the library is still under development, but it can already solve some system of equations and offer a basis for anyone who wants to solve their own problem.
