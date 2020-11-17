import sys,os
sys.path.append(os.path.join(os.path.dirname(__file__), os.path.pardir))
from PyEFVLib import MSHReader, XDMFReader, Grid, ProblemData, CgnsSaver, CsvSaver


if __name__ == "__main__":

	# Tests stencil from mesh generated with Gmsh
	model = "workspace/heat_transfer_2d/linear"
	problemData = ProblemData(model)
	print(problemData.paths["Grid"])
	reader = MSHReader(problemData.paths["Grid"])
	grid = Grid(reader.getData())
	grid.buildStencil()
	print(grid.stencil)

	# Tests stencil from mesh generated with Gmsh and converted to XDMF with Meshio
	directory = [os.path.dirname(__file__), os.path.pardir, "meshes/xdmf"]
	directory = os.path.join(*directory)
	filename = "3x3 tri.xdmf"
	boundariesFilename = "3x3 tri_facets.xdmf"
	subdomainsFilename = "3x3 tri_physical_region.xdmf"
	zoneList = "3x3 tri_zone_list.xml"
	print('{}/{}'.format(directory, filename))
	reader = XDMFReader(directory=directory, filename=filename, boundariesFilename=boundariesFilename, subdomainsFilename=subdomainsFilename)
	reader.readZoneList(zoneList)
	reader.setFacetData('gmsh:physical')
	reader.setSubdomainData('gmsh:physical')
	reader.read()
	grid = Grid(reader.getData())
	grid.buildStencil()
	print(grid.stencil)

	# Tests stencil from mesh generated with Ansys-ICEM and converted to XDMF with MSHtoXDMF
	directory = [os.path.dirname(__file__), os.path.pardir, "meshes/xdmf"]
	directory = os.path.join(*directory)
	filename = "Subdomains.xdmf"
	boundariesFilename = "Subdomains_facets.xdmf"
	subdomainsFilename = "Subdomains_physical_region.xdmf"
	zoneList = "Subdomains_zone_list.xml"
	print('{}/{}'.format(directory, filename))
	reader = XDMFReader(directory=directory, filename=filename, boundariesFilename=boundariesFilename, subdomainsFilename=subdomainsFilename)
	reader.readZoneList(zoneList)
	reader.setFacetData('boundaries')
	reader.setSubdomainData('subdomains')
	reader.read()
	grid = Grid(reader.getData())
	grid.buildStencil()
	print(grid.stencil)
