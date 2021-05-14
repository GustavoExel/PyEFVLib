from PyEFVLib.Vertex import *
from PyEFVLib.OuterFace import *
from PyEFVLib.GridData import *
from PyEFVLib.InnerFace import *
from PyEFVLib.Point import *
from PyEFVLib.Shape import *
from PyEFVLib.Region import *
from PyEFVLib.Facet import *
from PyEFVLib.MSHReader import *
from PyEFVLib.XDMFReader import *
from PyEFVLib.Boundary import *
from PyEFVLib.Element import *
from PyEFVLib.Grid import *
from PyEFVLib.BoundaryConditions import *
from PyEFVLib.CsvSaver import *
from PyEFVLib.VtuSaver import *
from PyEFVLib.VtmSaver import *
from PyEFVLib.MeshioSaver import *
from PyEFVLib.ProblemData import *
from PyEFVLib.Solver import *
from PyEFVLib.LinearSystem import *

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