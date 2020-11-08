#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Nov  8 19:00:15 2020

@author: herminio
"""

import sys,os
sys.path.append(os.path.join(os.path.dirname(__file__), os.path.pardir))
from PyEFVLib.simulation import LinearSystem as ls
import matplotlib.pyplot as plt

stencil = [	[0, 1, 5, 4], [0, 1, 5, 2, 6], [1, 2, 6, 3, 7], [2, 3, 7], [4, 5, 9, 0, 8], [0, 1, 5, 4, 9, 6, 10], [1, 2, 6, 5, 10, 7, 11],
				[2, 3, 7, 6, 11], [8, 9, 13, 4, 12], [4, 5, 9, 8, 13, 10, 14], [5, 6, 10, 9, 14, 11, 15], [6, 7, 11, 10, 15], [13, 12, 8],
				[8, 9, 13, 12, 14], [9, 10, 14, 13, 15], [10, 11, 15, 14]]

nDOF = 3
nVertices = len(stencil)
size = (nDOF-1)* nVertices
ls_csr = ls.LinearSystemCSR(stencil, nDOF)
ls_csr.initialize()

for i, v_stencil in enumerate(stencil):
	for j in v_stencil:
		for dof_i in range(nDOF):
			for dof_j in range(nDOF):
				ls_csr.addValueToMatrix(i + dof_i*nVertices, j + dof_j*nVertices, -3)

ls_csr.matZeroRow(nDOF*nVertices-1, 1)

M = ls_csr.matrix.todense()

f, ax = plt.subplots(nrows=1, ncols=1, figsize=(10,8))
ax.spy(M, markersize=5)

plt.tight_layout()
plt.show()
