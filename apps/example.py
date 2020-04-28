import sys,os
sys.path.append(os.path.join(os.path.dirname(__file__), os.path.pardir))

from libs.geometry.MSHReader import MSHReader
from libs.geometry.Grid import Grid
from libs.geometry.Vertex import Point
import os, numpy as np

path = [os.path.dirname(__file__), os.path.pardir, "meshes", "Square.msh"]

reader = MSHReader(os.path.join(*path))
grid   = Grid(reader.getData())

elementVolumes = sum([element.volume for element in grid.elements])
controlVolumes = sum([vertex.volume for vertex in grid.vertices])

print("A soma dos volumes dos elementos é", elementVolumes)
print("Já a soma dos volumes de controle é", controlVolumes)