import numpy as np

class Point:
	def __init__(self, x, y, z=0.0):
		self.x = x
		self.y = y
		self.z = z
		
	def getCoordinates(self):
		return np.array([self.x, self.y, self.z])