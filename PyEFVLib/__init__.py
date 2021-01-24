from PyEFVLib.geometry.Vertex import *
from PyEFVLib.geometry.OuterFace import *
from PyEFVLib.geometry.GridData import *
from PyEFVLib.geometry.InnerFace import *
from PyEFVLib.geometry.Point import *
from PyEFVLib.geometry.Shape import *
from PyEFVLib.geometry.Region import *
from PyEFVLib.geometry.Facet import *
from PyEFVLib.geometry.MSHReader import *
from PyEFVLib.geometry.XDMFReader import *
from PyEFVLib.geometry.Boundary import *
from PyEFVLib.geometry.Element import *
from PyEFVLib.geometry.Grid import *
from PyEFVLib.simulation.BoundaryConditions import *
from PyEFVLib.simulation.CgnsSaver import *
from PyEFVLib.simulation.CsvSaver import *
from PyEFVLib.simulation.VtuSaver import *
from PyEFVLib.simulation.VtmSaver import *
from PyEFVLib.simulation.MeshioSaver import *
from PyEFVLib.simulation.ProblemData import *
from PyEFVLib.simulation.Solver import *
from PyEFVLib.simulation.LinearSystem import *

READERS_DICTIONARY = {
	"msh": MSHReader,
	# "xdmf": XDMFReader,	# As classes são bem diferentes, tem que conversar pra fazer bater
}

Neumann   = "NEUMANN_BOUNDARY_CONDITION"
Dirichlet = "DIRICHLET_BOUNDARY_CONDITION"
Constant  = "BOUNDARY_CONDITION_CONSTANT_VALUE" 
Variable  = "BOUNDARY_CONDITION_VARIABLE_VALUE" 

def read(filePath):
	extension = filePath.split('.')[-1]

	if extension not in READERS_DICTIONARY.keys():
		raise Exception("File extension not supported yet! Input your extension sugestion at https://github.com/GustavoExel/PyEFVLib")
	else:
		return Grid( READERS_DICTIONARY[ extension ]( filePath ).getData() )