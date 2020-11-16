from petsc4py import PETSc
import numpy as np
from scipy.sparse import csr_matrix


def createPETScVector(size):
	vector = PETSc.Vec()
	vector.create()
	vector.setSizes(size)
	vector.setUp()
	vector.setValues(range(size), np.zeros(size))
	return vector

def createPETScMatrix(size, mat_type='dense'):
	matrix = PETSc.Mat()
	matrix.create()
	matrix.setSizes([size, size])
	matrix.setType(mat_type)
	matrix.setUp()	
	return matrix

class LinearSystem(object):
	def __init__(self, stencil, nDOF, PETSc_backend=False):
		self.stencil = stencil
		self.nDOF = nDOF
		self.nVertices = len(stencil)
		self.size = self.nDOF*self.nVertices
		if PETSc_backend:
			self.backend = 'petsc'
		else:
			self.backend = 'numpy/scipy'

	def initialize(self):
		if self.backend == 'numpy/scipy':
			self.rhs = np.zeros(self.size)
			self.matrix = np.zeros((self.size, self.size))
		elif self.backend == 'petsc':
			self.rhs = createPETScVector(self.size)
			self.matrix = createPETScMatrix(self.size)

	def addValueToMatrix(self, row, col, value):
		if self.backend == 'numpy/scipy':
			self.matrix[row][col] += value
		elif self.backend == 'petsc':
			self.matrix.setValue(row=row, col=col, value=value, addv=True)

	def setValueToMatrix(self, row, col, value):
		if self.backend == 'numpy/scipy':
			self.matrix[row][col] = value
		elif self.backend == 'petsc':
			self.matrix.setValue(row=row, col=col, value=value, addv=False)

	def getValueFromMatrix(self, row, col):
		if self.backend == 'numpy/scipy':
			return self.matrix[row][col]
		elif self.backend == 'petsc':
			return self.matrix.getValue(row=row, col=col)

	def addValueToRHS(self, row, value):
		if self.backend == 'numpy/scipy':
			self.rhs[row] += value
		elif self.backend == 'petsc':
			self.rhs.setValue(index=row, value=value, addv=True)

	def setValueToRHS(self, row, value):
		if self.backend == 'numpy/scipy':
			self.rhs[row] = value
		elif self.backend == 'petsc':
			self.rhs.setValue(index=row, value=value, addv=False)

	def getValueFromRHS(self, row):
		if self.backend == 'numpy/scipy':
			return self.rhs[row]
		elif self.backend == 'petsc':
			return self.rhs.getValue(index=row)

	def restartRHS(self):
		if self.backend == 'numpy/scipy':
			self.rhs = np.zeros(self.size)
		elif self.backend == 'petsc':
			self.rhs = createPETScVector(self.size)

	def matZeroRow(self, row, diagonal_value):
		for col in range(self.size):
			self.setValueToMatrix(row, col, 0)
		self.setValueToMatrix(row, row, diagonal_value)

	def getDense(self):
		if self.backend == 'numpy/scipy':
			return self.matrix
		elif self.backend == 'petsc':
			return self.matrix.convert('dense').getDenseArray()

	def getRHS(self):
		if self.backend == 'numpy/scipy':
			return self.rhs
		elif self.backend == 'petsc':
			return self.rhs.getArray()

	def assembly(self):
		if self.backend == 'petsc':
			self.matrix.assemblyBegin()
			self.matrix.assemblyEnd()
			self.rhs.assemblyBegin()
			self.rhs.assemblyEnd()
		else:
			return

class LinearSystemCSR(LinearSystem):
	def initialize(self):
		if self.backend == 'numpy/scipy':
			self.rhs = np.zeros(self.size)
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
		elif self.backend == 'petsc':
			self.rhs = createPETScVector(self.size)
			self.matrix = createPETScMatrix(self.size, mat_type='aij')

	def __getIndex(self, row, col):
		pos1 = self.matrix.indptr[row]
		pos2 = self.matrix.indptr[row+1]
		aux_index = np.where(self.matrix.indices[pos1: pos2] == col)
		if aux_index[0].size != 0:
			return pos1 + aux_index[0][0]
		else:
			print("New entries not allowed.")

	def addValueToMatrix(self, row, col, value):
		if self.backend == 'numpy/scipy':
			index = self.__getIndex(row, col)
			self.matrix.data[index] += value
		elif self.backend == 'petsc':
			self.matrix.setValue(row=row, col=col, value=value, addv=True)

	def setValueToMatrix(self, row, col, value):
		if self.backend == 'numpy/scipy':
			index = self.__getIndex(row, col)
			self.matrix.data[index] = value
		elif self.backend == 'petsc':
			self.matrix.setValue(row=row, col=col, value=value, addv=False)

	def getValueFromMatrix(self, row, col):
		if self.backend == 'numpy/scipy':
			index = self.__getIndex(row, col)
			return self.matrix.data[index]
		elif self.backend == 'petsc':
			return self.matrix.getValue(row=row, col=col)

	def matZeroRow(self, row, diagonal_value):
		if self.backend == 'numpy/scipy':
			for index in range(self.matrix.indptr[row], self.matrix.indptr[row+1]):
				self.matrix.data[index] = 0.0
			self.setValueToMatrix(row, row, diagonal_value)
		elif self.backend == 'petsc':
			for col in range(self.size):
				self.setValueToMatrix(row, col, 0)
			self.setValueToMatrix(row, row, diagonal_value)

	def getDense(self):
		if self.backend == 'numpy/scipy':
			return self.matrix.todense()
		elif self.backend == 'petsc':
			return self.matrix.convert('dense').getDenseArray()

	def getSparse(self):
		if self.backend == 'numpy/scipy':
			return (self.matrix.indptr, self.matrix.indices, self.matrix.data)
		elif self.backend == 'petsc':
			return self.matrix.getValuesCSR()


if __name__ == '__main__':
	stencil = [	[0, 1, 5, 4], [0, 1, 5, 2, 6], [1, 2, 6, 3, 7], [2, 3, 7], [4, 5, 9, 0, 8], [0, 1, 5, 4, 9, 6, 10], [1, 2, 6, 5, 10, 7, 11],
				[2, 3, 7, 6, 11], [8, 9, 13, 4, 12], [4, 5, 9, 8, 13, 10, 14], [5, 6, 10, 9, 14, 11, 15], [6, 7, 11, 10, 15], [13, 12, 8],
				[8, 9, 13, 12, 14], [9, 10, 14, 13, 15], [10, 11, 15, 14]]
	stencil = [[0, 1], [0, 2, 3], [], [0, 4], [], [2, 3, 5]]
	nDOF = 1
	ls = LinearSystemCSR(stencil, nDOF)
	ls.initialize()
	i = 0
	for v_stencil in stencil:
		for j in v_stencil:
			ls.addValueToMatrix(i, j, -3)
		i += 1
	ls.addValueToMatrix(3, 0, -12)
	ls.addValueToMatrix(3, 0, -12)
	ls.setValueToMatrix(3, 4, -5)
	ls.assembly()
	print(ls.getDense())
	ls.matZeroRow(5, -20)
	print(ls.getDense())