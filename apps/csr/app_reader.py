import sys,os
sys.path.append(os.path.join(os.path.dirname(__file__), os.path.pardir))

from PyEFVLib import MSHReader, Grid, Point
import os, numpy as np

# Tests mesh generated with Gmsh
path = [os.path.dirname(__file__), os.path.pardir, "meshes", "msh", "2D", "Square.msh"]
reader = MSHReader(os.path.join(*path))
grid = Grid(reader.getData())
elementVolumes = sum([element.volume for element in grid.elements])
controlVolumes = sum([vertex.volume for vertex in grid.vertices])
print("*.msh file:")
print("Sum of elements volume:", elementVolumes)
print("Sum of CV volume:", controlVolumes)

from PyEFVLib import XDMFReader

# Tests mesh generated with Gmsh and converted to XDMF with Meshio
directory = [os.path.dirname(__file__), os.path.pardir, "meshes/xdmf"]
directory = os.path.join(*directory)
filename = "Square.xdmf"
boundariesFilename = "Square_facets.xdmf"
subdomainsFilename = "Square_physical_region.xdmf"
zoneList = "Square_zone_list.xml"
reader = XDMFReader(directory=directory, filename=filename, boundariesFilename=boundariesFilename, subdomainsFilename=subdomainsFilename)
reader.readZoneList(zoneList)
reader.setFacetData('gmsh:physical')
reader.setSubdomainData('gmsh:physical')
reader.read()
grid = Grid(reader.getData())
elementVolumes = sum([element.volume for element in grid.elements])
controlVolumes = sum([vertex.volume for vertex in grid.vertices])
print("*.xdmf file (from Gmsh):")
print("Sum of elements volume:", elementVolumes)
print("Sum of CV volume:", controlVolumes)

# Tests mesh generated with Ansys-ICEM and converted to XDMF with MSHtoXDMF
filename = "Subdomains.xdmf"
boundariesFilename = "Subdomains_facets.xdmf"
subdomainsFilename = "Subdomains_physical_region.xdmf"
zoneList = "Subdomains_zone_list.xml"
reader = XDMFReader(directory=directory, filename=filename, boundariesFilename=boundariesFilename, subdomainsFilename=subdomainsFilename)
reader.readZoneList(zoneList)
reader.setFacetData('boundaries')
reader.setSubdomainData('subdomains')
reader.read()
grid = Grid(reader.getData())
elementVolumes = sum([element.volume for element in grid.elements])
controlVolumes = sum([vertex.volume for vertex in grid.vertices])
print("*.xdmf file (from Ansys):")
print("Sum of elements volume:", elementVolumes)
print("Sum of CV volume:", controlVolumes)