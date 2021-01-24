import xmltodict
import meshio
import sys


filename = sys.argv[1]
msh = meshio.read(filename + ".msh")
for cellblock in msh.cells:
	if cellblock.type == "triangle":
		triangle_cells = cellblock.data
	elif  cellblock.type == "line":
		line_cells = cellblock.data
for key in msh.cell_data_dict["gmsh:physical"].keys():
	if key == "line":
		line_data = msh.cell_data_dict["gmsh:physical"][key]
	elif key == "triangle":
		triangle_data = msh.cell_data_dict["gmsh:physical"][key]
main_mesh = meshio.Mesh(points=msh.points, cells={"triangle": triangle_cells})
boundary_mesh =meshio.Mesh(points=msh.points, cells={"line": line_cells}, cell_data={"gmsh:physical": [line_data]})
subdomains_mesh = meshio.Mesh(points=msh.points, cells={"triangle": triangle_cells}, cell_data={"gmsh:physical": [triangle_data]})
meshio.write("../../xdmf/" + filename + ".xdmf", main_mesh)
meshio.write("../../xdmf/" + filename + "_facets.xdmf", boundary_mesh)
meshio.write("../../xdmf/" + filename + "_physical_region.xdmf", subdomains_mesh)
zone_list = {'ZoneList': {}}
for key in msh.field_data.keys():
	tag = msh.field_data[key][0]
	dim = msh.field_data[key][1]
	if dim == 2:
		zone_type = 'fluid'
	else:
		zone_type = 'wall'
	zone_list['ZoneList']['Zone' + str(tag)] = {'@type': zone_type, '@name': key}
file = "../../xdmf/{}_zone_list.xml".format(filename)
f = open(file, "w")
f.write(xmltodict.unparse(zone_list, pretty=True))
f.close()
