import numpy as np
from scipy import sparse
from matplotlib import pyplot as plt

def belongsToPolygon(x, y, X, Y):
	"""
	Lema: Um ponto está dentro de um polígono se a quantidade de interseções de uma linha horizontal sobre sua coordenada Y for ímpar em ambos os lados
	Algoritmo:
	Se um dos Y está abaixo do ponto e outro acima ou em cima:
		Se pelo menos um dos pontos estiver à esquerda:
			Se a interseção estier à esquerda:
				Conta um vértice
	"""
	numberOfCorners = len(X)
	oddNodes=False

	for i in range(numberOfCorners):
		# horizontal
		if abs(Y[i]-Y[i-1]) < 1e-10 and abs(Y[i-1]-y) < 1e-10 and min(X[i],X[i-1]) <= x <= max(X[i],X[i-1]):
			return 2
		# vertical
		if abs(X[i]-X[i-1]) < 1e-10 and abs(X[i-1]-x) < 1e-10 and min(Y[i],Y[i-1]) <= y <= max(Y[i],Y[i-1]):
			return 2
		# diagonal
		if abs(Y[i-1]-Y[i]) > 1e-10 and abs(X[i-1]-X[i]) > 1e-10 and abs( X[i]+(y-Y[i])*(X[i-1]-X[i])/(Y[i-1]-Y[i]) - x ) < 1e-10:
			return 2
		
		if (Y[i]< y and Y[i-1]>=y or Y[i-1]< y and Y[i]>=y) and (X[i]<=x or X[i-1]<=x):
			oddNodes = ( oddNodes != bool( X[i]+(y-Y[i])/(Y[i-1]-Y[i])*(X[i-1]-X[i]) < x ) )
	return oddNodes

def generateMesh(X,Y,dx,dy,boundaryNames,outputPath):
	# Creata the Grid
	nx = int( (max(X)-min(X)) / dx ) + 1
	ny = int( (max(Y)-min(Y)) / dy ) + 1
	gridNodes = [(x,y) for y in np.linspace(0,(max(Y)-min(Y)),ny) for x in np.linspace(0,(max(X)-min(X)),nx)]
	gridSquares = [ (vIdx, vIdx+1, vIdx+1+nx, vIdx+nx) for i, vIdx in enumerate([vIdx for vIdx in range(nx*ny-nx-1) if (vIdx+1)%nx!=0]) ]

	# Check which gridNodes are inside the polygon
	gridNodesInMesh = [ int(belongsToPolygon(*node, X, Y)) for node in gridNodes ]
	meshNodes = [node for node, isIn in zip(gridNodes,gridNodesInMesh) if isIn]
	gridToMeshDict = { idx:sum( [bool(b) for b in gridNodesInMesh[:idx]] ) for idx, (node, isIn) in enumerate( zip(gridNodes,gridNodesInMesh) ) if isIn }

	# Store gridSquares that belong to the polygon
	squaresInMesh = []
	squaresCut = []
	squaresCutIntersection = []
	for idx, square in enumerate(gridSquares):
		squareNodesIn = [ node for node in square if gridNodesInMesh[node] ]
		squareNodesInside = [ node for node in square if gridNodesInMesh[node] == 1 ]
		if len(squareNodesIn) == 4:
			squaresInMesh.append(idx)
		elif len(squareNodesIn) > 1:
			squaresCut.append( idx )
			squaresCutIntersection.append(squareNodesIn)

	# Square connectivities
	mesh2DElements = [[ gridToMeshDict[gridNodeIdx] for gridNodeIdx in gridSquares[squareIdx] ] for squareIdx in squaresInMesh]

	# Create cut gridSquares
	for squareIdx, squareNodesInMesh in zip( squaresCut, squaresCutIntersection ):
		squareNodes = gridSquares[squareIdx]

		squareIntersections = []
		for i in range(4):
			n2 = gridNodes[ squareNodes[(i+1)%4] ]
			n1 = gridNodes[ squareNodes[i] ]
			x1, y1 = n1[0], n1[1]
			x2, y2 = n2[0], n2[1]
			for j in range(len(X)):
				nB = ( X[(j+1)%len(X)], Y[(j+1)%len(X)] )
				nA = ( X[j], Y[j] )
				xA, yA = nA[0], nA[1]
				xB, yB = nB[0], nB[1]

				denominator = (x2-x1)*(yB-yA)-(xB-xA)*(y2-y1)
				if denominator != 0:
					xi = ((x2*y1-x1*y2)*(xB-xA)-(xB*yA-xA*yB)*(x2-x1))/denominator
					yi = ((x2*y1-x1*y2)*(yB-yA)-(xB*yA-xA*yB)*(y2-y1))/denominator

					if min(x1,x2) < xi < max(x1,x2) or min(y1,y2) < yi < max(y1,y2):
						if not [1 for sn in gridSquares[squareIdx] if abs(gridNodes[sn][0]-xi)<1e-10 and abs(gridNodes[sn][1]-yi)<1e-10]:
							squareIntersections.append((xi,yi))

		# Create cutted polygons
		if 3 <= len(squareNodesInMesh)+len(squareIntersections) <= 4:
			poly = []
			totalCount = squareNodeCount = intersectionCount = 0
			while True:
				if squareNodeCount < len(squareNodesInMesh) and gridSquares[squareIdx][totalCount] == squareNodesInMesh[squareNodeCount]:
					poly.append( gridToMeshDict[ squareNodesInMesh[squareNodeCount] ] )
					squareNodeCount += 1
				elif intersectionCount < len(squareIntersections):
					xi, yi = squareIntersections[intersectionCount]
					nIdx = [ i for i,(xn,yn) in enumerate(meshNodes) if abs(xn-xi)<1e-10 and abs(yn-yi)<1e-10 ]
					if nIdx:
						poly.append( nIdx[0] )						
					else:
						meshNodes.append( (xi,yi) )
						poly.append( len(meshNodes)-1 )
					intersectionCount += 1
				totalCount += 1

				if squareNodeCount >= len(squareNodesInMesh) and intersectionCount >= len(squareIntersections):
					break
			
			valid = True
			for i in range(len(poly)):
				if not belongsToPolygon( 0.5*meshNodes[ poly[i] ][0] + 0.5*meshNodes[ poly[i-1] ][0], 0.5*meshNodes[ poly[i] ][1] + 0.5*meshNodes[ poly[i-1] ][1], X, Y ):
					valid = False
			if valid:
				mesh2DElements.append(poly)

		elif False and len(squareNodesInMesh)+len(squareIntersections) == 5:
			midpoint = ( (squareIntersections[0][0]+squareIntersections[1][0])/2, (squareIntersections[0][1]+squareIntersections[1][1])/2 )
			meshNodes.append(midpoint)
			midpointIdx = len(meshNodes)-1
			
			totalCount = squareNodeCount = intersectionCount = 0
			for i in range(2):
				poly = []
				for j in range(3):
					if gridSquares[squareIdx][totalCount] == squareNodesInMesh[squareNodeCount]:
						poly.append(squareNodesInMesh[squareNodeCount])
						squareNodeCount += 1
					else:
						# TEM QUE CONSERTAR DENOVO SE NÃO FICA VÉRTICE DUPLICADO
						meshNodes.append( squareIntersections[intersectionCount] )
						poly.append( len(meshNodes)-1 )
						intersectionCount += 1
					totalCount += 1

					poly.append(midpointIdx)
					mesh2DElements.append(poly)
					squareNodeCount -= 1
					totalCount -= 1

	# Mesh Boundaries
	boundaries = []
	for idx, name in enumerate(boundaryNames):
		xA, yA = X[idx], Y[idx]
		xB, yB = X[(idx+1)%len(X)], Y[(idx+1)%len(X)]

		bNodes = sorted([ idx for idx, node in enumerate(meshNodes) if belongsToPolygon(*node, (xA,xB), (yA,yB)) ])
		boundaries.append(bNodes)

	elements = []
	for bIdx, boundary in enumerate(boundaries):
		elements += [[len(elements)+i+1, *(1, 2), *(bIdx+2, bIdx+2), *( boundary[i]+1, boundary[i+1]+1 )] for i in range(len(boundary)-1) ]
	elements += [[len(elements)+i+1, *(len(element)-1, 2), *(1,1), *[nIdx+1 for nIdx in element]] for i, element in enumerate(mesh2DElements) ]

	physicalNames = {"Body":2}
	physicalNames.update({ bName:1 for bName in boundaryNames })

	# Write MSH file
	text  = "$MeshFormat\n2.2 0 8\n$EndMeshFormat\n$PhysicalNames\n"
	text += str( len(physicalNames) ) + "\n"
	text += "".join([ "{} {} \"{}\"\n".format(physicalNames[name], idx+1, name) for idx, name in enumerate(physicalNames.keys()) ])
	text += "$EndPhysicalNames\n$Nodes\n"
	text += str( len(meshNodes) ) + "\n"
	text += "".join([ "{} {} {} 0.0\n".format(idx+1, x, y) for idx,(x,y) in enumerate(meshNodes) ])
	text += "$EndNodes\n$Elements\n"
	text += str( len(elements) ) + "\n"
	text += "\n".join([" ".join([str(ev) for ev in e]) for e in elements])
	text += "\n$EndElements\n"

	with open(outputPath, "w") as f:
		f.write(text)

if __name__ == "__main__":
	# X = [0,2,3]
	# Y = [0,1,0]
	X = [0,1.9,3]
	Y = [0,0.9,0]
	# X = [0,2,2,0]
	# Y = [0,0,2,2]
	# X = [0,1,1,0.8,0.8,0]
	# Y = [0,0,1,1,0.2,0.2]
	namez = ["A","B","C"]
	# namez = ["A","B","C","D"]	
	# namez = ["A", "B","C","D","E","F"]
	dx = 0.05
	dy = 0.05

	generateMesh(X,Y,dx,dy,namez,"a.msh")