import numpy as np

class Point:
	def __init__(self, x, y, z=0.0):
		self.x = x
		self.y = y
		self.z = z

	def __add__(self, p2):
		return Point(self.x+p2.x, self.y+p2.y, self.z+p2.z)

	def __sub__(self, p2):
		return Point(self.x-p2.x, self.y-p2.y, self.z-p2.z)

	def __mul__(self, p2):
		if p2.__class__ == self.__class__:
			return Point(self.x*p2.x, self.y*p2.y, self.z*p2.z)
		else:
			return Point(self.x*p2, self.y*p2, self.z*p2)

	def __truediv__(self, n):
		return Point(self.x/n, self.y/n, self.z/n)

	def __repr__(self):
		return "Point({}, {}, {})".format(self.x, self.y, self.z)

	def getCoordinates(self):
		return np.array([self.x, self.y, self.z])
