import numpy as np
from scipy.sparse import csr_matrix

class LinearSystem(object):
	def __init__(self, stencil, nDOF):
		self.stencil = stencil
		self.nDOF = nDOF
		self.nVertices = len(stencil)
		self.size = self.nDOF*self.nVertices

	def initialize(self):
		self.rhs = np.zeros(self.size)
		self.matrix = np.zeros((self.size, self.size))

	def addValueToMatrix(self, row, col, value):
		self.matrix[row][col] += value

	def setValueToMatrix(self, row, col, value):
		self.matrix[row][col] = value

	def addValueToRHS(self, row, value):
		self.rhs[row] += value

	def setValueToRHS(self, row, value):
		self.rhs[row] = value

	def restartRHS(self):
		self.rhs = np.zeros(self.size)




class LinearSystemCSR(LinearSystem):
	def initialize(self):
		self.rhs = np.zeros(self.nDOF*self.nVertices)
		indPtr = [0]
		indices = []
		for i in range(self.nDOF):
			for vertex_stencil in self.stencil:
				indPtr.append(indPtr[-1] + self.nDOF*len(vertex_stencil))
				for j in range(self.nDOF):
					for vertex in vertex_stencil:
						indices.append(vertex + j*self.nVertices)
		indPtr = np.array(indPtr)
		indices = np.array(indices)
		data = np.zeros(indices.size)
		self.matrix = csr_matrix((data, indices, indPtr), shape=(self.size, self.size))

	def __getIndex(self, row, col):
		pos1 = self.matrix.indptr[row]
		pos2 = self.matrix.indptr[row+1]
		aux_index = np.where(self.matrix.indices[pos1: pos2] == col)
		if aux_index[0].size != 0:
			return pos1 + aux_index[0][0]
		else:
			print("New entries not allowed.")

	def addValueToMatrix(self, row, col, value):
		index = self.__getIndex(row, col)
		self.matrix.data[index] += value

	def setValueToMatrix(self, row, col, value):
		index = self.__getIndex(row, col)
		self.matrix.data[index] = value

	def matZeroRow(self, row, diagonal_value):
		for index in range(self.matrix.indptr[row], self.matrix.indptr[row+1]):
			self.matrix.data[index] = 0.0
		self.setValueToMatrix(row, row, diagonal_value)






if __name__ == '__main__':
	stencil = [	[0, 1, 5, 4], [0, 1, 5, 2, 6], [1, 2, 6, 3, 7], [2, 3, 7], [4, 5, 9, 0, 8], [0, 1, 5, 4, 9, 6, 10], [1, 2, 6, 5, 10, 7, 11],
				[2, 3, 7, 6, 11], [8, 9, 13, 4, 12], [4, 5, 9, 8, 13, 10, 14], [5, 6, 10, 9, 14, 11, 15], [6, 7, 11, 10, 15], [13, 12, 8],
				[8, 9, 13, 12, 14], [9, 10, 14, 13, 15], [10, 11, 15, 14]]

	stencil = [[0, 1], [0, 2, 3], [], [0, 4], [], [2, 3, 5]]
	nDOF = 1
	ls_csr = LinearSystemCSR(stencil, nDOF)
	ls_csr.initialize()


	i = 0
	for v_stencil in stencil:
		for j in v_stencil:
			ls_csr.addValueToMatrix(i, j, -3)
		i += 1

	ls_csr.addValueToMatrix(3, 0, -12)
	ls_csr.addValueToMatrix(3, 0, -12)
	ls_csr.setValueToMatrix(3, 4, -5)
	print(ls_csr.matrix.todense())
	ls_csr.matZeroRow(5, -20)
	print(ls_csr.matrix.todense())
