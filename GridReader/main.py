from libs.geometry.MSHReader import MSHReader
from libs.geometry.Grid import Grid
import time
from libs.simulation.HeatTransfer2D import HeatTransfer2D

r = MSHReader("/home/gustavoe/Documents/Sinmec/GTRelated/GridReader/meshes/Square.msh")#QuadPlate.msh")
g = Grid(r.getData())

# s = HeatTransfer2D("/home/gustavoe/Documents/Sinmec/GTRelated/GridReader/meshes/Square.msh")

