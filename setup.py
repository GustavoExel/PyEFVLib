from setuptools import setup

with open("README.md", "r") as f:
	long_description = f.read()

setup(
	name="PyEFVLib",
	version="1.0.11",
	description="A library that facilitates the implementation of the Element-based Finite Volumes Method (EbFVM) to solve partial differential equations (PDEs)",
	py_modules=["PyEFVLib/__init__", "PyEFVLib/Vertex", "PyEFVLib/OuterFace", "PyEFVLib/GridData", "PyEFVLib/InnerFace", "PyEFVLib/Point", "PyEFVLib/Shape", "PyEFVLib/Region", "PyEFVLib/Facet", "PyEFVLib/MSHReader", "PyEFVLib/XDMFReader", "PyEFVLib/Boundary", "PyEFVLib/Element", "PyEFVLib/Grid", "PyEFVLib/BoundaryConditions", "PyEFVLib/Saver", "PyEFVLib/CsvSaver", "PyEFVLib/VtuSaver", "PyEFVLib/VtmSaver", "PyEFVLib/MeshioSaver", "PyEFVLib/ProblemData", "PyEFVLib/Solver", "PyEFVLib/LinearSystem"],
	# package_dir={"": "PyEFVLib"},
	classifiers=[
		"Programming Language :: Python :: 3",
		"Programming Language :: Python :: 3.6",
		"Programming Language :: Python :: 3.7",
		"Programming Language :: Python :: 3.8",
		"License :: OSI Approved :: MIT License",
		"Operating System :: OS Independent"
	],
	long_description=long_description,
	long_description_content_type="text/markdown",
	install_requires=[
		"numpy ~= 1.19.2",
		"pandas ~= 1.1.2",
		"scipy ~= 1.5.2",
		"meshio ~= 4.3.4",
		"h5py ~= 3.2.1",
		"xmltodict ~= 0.12.0",
	],
	url="https://github.com/GustavoExel/PyEFVLib",
	author="Gustavo Exel",
	author_email="gustavoexelgpe@gmail.com",

)