import numpy as np
import matplotlib.pyplot as plt

X, Y = (0.0,1.0,1.0,0.7,0.7,0.0), (0.0,0.0,1.0,1.0,0.3,0.3)
dx = 0.2
dy = 0.1
outputPath = "teemp.msh"

def belongsToPolygon(x, y, X, Y):
	numberOfCorners = len(X)
	oddNodes=False
	for i in range(numberOfCorners):
		verticalLine = abs(Y[i]-Y[i-1]) < 1e-10
		horizontalLine = abs(X[i]-X[i-1]) < 1e-10
		inBetweenX = min(X[i],X[i-1]) <= x <= max(X[i],X[i-1])
		inBetweenY = min(Y[i],Y[i-1]) <= y <= max(Y[i],Y[i-1])
		if verticalLine and abs(Y[i-1]-y) < 1e-10 and inBetweenX:
			return 2
		if horizontalLine and abs(X[i-1]-x) < 1e-10 and inBetweenY:
			return 2
		if not verticalLine and not horizontalLine and inBetweenX and inBetweenY and abs( X[i]+(y-Y[i])*(X[i-1]-X[i])/(Y[i-1]-Y[i]) - x ) < 1e-10:
			return 2
		
		if (Y[i]< y and Y[i-1]>=y or Y[i-1]< y and Y[i]>=y) and (X[i]<=x or X[i-1]<=x):
			oddNodes = ( oddNodes != bool( X[i]+(y-Y[i])/(Y[i-1]-Y[i])*(X[i-1]-X[i]) < x ) )
	return oddNodes
def generateMesh(X,Y,dx,dy,boundaryNames,outputPath):
	preview = True		
	# Create the Grid
	nx = int( (max(X)-min(X)) / dx ) + 1
	ny = int( (max(Y)-min(Y)) / dy ) + 1
	gridNodes = [(x,y) for y in np.linspace(0,(max(Y)-min(Y)),ny) for x in np.linspace(0,(max(X)-min(X)),nx)]
	gridSquares = [ (vIdx, vIdx+1, vIdx+1+nx, vIdx+nx) for i, vIdx in enumerate([vIdx for vIdx in range(nx*ny-nx-1) if (vIdx+1)%nx!=0]) ]

	# Check which grid nodes are inside the polygon
	gridNodesInMesh = [ int(belongsToPolygon(*node, X, Y)) for node in gridNodes ]
	meshNodes = [node for node, isIn in zip(gridNodes,gridNodesInMesh) if isIn]
	gridToMeshDict = { idx:sum( [bool(b) for b in gridNodesInMesh[:idx]] ) for idx, (node, isIn) in enumerate( zip(gridNodes,gridNodesInMesh) ) if isIn }

	# Check which grid squares belong to the polygon, and the ones that pass through its boundaries
	squaresInMesh = []
	squaresCut = []
	squaresCutIntersection = []
	for idx, square in enumerate(gridSquares):
		squareNodesIn = [ node for node in square if gridNodesInMesh[node] ]
		squareNodesInside = [ node for node in square if gridNodesInMesh[node] == 1 ]
		if len(squareNodesIn) == 4:
			squaresInMesh.append(idx)
		elif len(squareNodesIn) >= 1:
			squaresCut.append( idx )
			squaresCutIntersection.append(squareNodesIn)

	if preview:
		for square in gridSquares:
			sX,sY = zip(*[gridNodes[vtxIdx] for vtxIdx in square])
			sX,sY = sX+(sX[0],), sY+(sY[0],)
			plt.plot(sX,sY,color="b")
		plt.plot(X+(X[0],),Y+(Y[0],),color="k")

		for square in squaresInMesh:
			sX,sY = zip(*[gridNodes[v] for v in gridSquares[square]])
			plt.plot(sX+(sX[0],),sY+(sY[0],),color="orange", linewidth=2)

		innerNodes = meshNodes.copy()


	# Mesh element connectivities
	mesh2DElements = [[ gridToMeshDict[gridNodeIdx] for gridNodeIdx in gridSquares[squareIdx] ] for squareIdx in squaresInMesh]
	originalOne = mesh2DElements.copy()

	# Create 2D elements near the boundary
	for _, (squareIdx, squareNodesInMesh) in enumerate(zip( squaresCut, squaresCutIntersection )):
		if _ != 3:
			continue
		squareNodes = gridSquares[squareIdx]

		squareIntersections = []
		squareIntersectionsLocalIndex = []

		# Compute intersections between the square and the polygon
		for i in range(4):
			# Coordinates of square nodes
			n1, n2 = gridNodes[ squareNodes[i] ], gridNodes[ squareNodes[(i+1)%4] ]
			x1, y1, x2, y2 = n1[0], n1[1], n2[0], n2[1]
			for j in range(len(X)):
				# Coordinates of polygon edge
				nA, nB = ( X[j], Y[j] ), ( X[(j+1)%len(X)], Y[(j+1)%len(X)] )
				xA, yA, xB, yB = nA[0], nA[1], nB[0], nB[1]

				denominator = (x2-x1)*(yB-yA)-(xB-xA)*(y2-y1)
				if denominator != 0.0:
					xi = ((x2*y1-x1*y2)*(xB-xA)-(xB*yA-xA*yB)*(x2-x1))/denominator
					yi = ((x2*y1-x1*y2)*(yB-yA)-(xB*yA-xA*yB)*(y2-y1))/denominator

					if ( min(x1,x2) < xi < max(x1,x2) or min(y1,y2) < yi < max(y1,y2) ) and ( min(xA,xB) <= xi <= max(xA,xB) or min(yA,yB) <= yi <= max(yA,yB) ) :
						if not [1 for sn in gridSquares[squareIdx] if abs(gridNodes[sn][0]-xi)<1e-10 and abs(gridNodes[sn][1]-yi)<1e-10]:
							# Check if intersection has already been counted
							if not [n for n in squareIntersections if abs(n[0]-xi)<1e-10 and abs(n[1]-yi)<1e-10]:
								squareIntersections.append((xi,yi))
								squareIntersectionsLocalIndex.append(i)

		if 3 <= len(squareNodesInMesh)+len(squareIntersections) <= 4:
			poly = []
			totalCount = squareNodeCount = intersectionCount = 0
			while True:
				if squareNodeCount < len(squareNodesInMesh) and gridSquares[squareIdx][totalCount] == squareNodesInMesh[squareNodeCount]:
					poly.append( gridToMeshDict[ squareNodesInMesh[squareNodeCount] ] )
					squareNodeCount += 1

				if intersectionCount < len(squareIntersections) and totalCount == squareIntersectionsLocalIndex[intersectionCount]:
					xi, yi = squareIntersections[intersectionCount]
					nIdx = [ i for i,(xn,yn) in enumerate(meshNodes) if abs(xn-xi)<1e-10 and abs(yn-yi)<1e-10 ]
					if nIdx:
						poly.append( nIdx[0] )						
					else:
						meshNodes.append( (xi,yi) )
						poly.append( len(meshNodes)-1 )
					intersectionCount += 1

				if squareNodeCount >= len(squareNodesInMesh) and intersectionCount >= len(squareIntersections):
					break
				totalCount += 1

			
			valid = True
			for i in range(len(poly)):
				if not belongsToPolygon( 0.5*meshNodes[ poly[i] ][0] + 0.5*meshNodes[ poly[i-1] ][0], 0.5*meshNodes[ poly[i] ][1] + 0.5*meshNodes[ poly[i-1] ][1], X, Y ):
					valid = False
			if valid:
				mesh2DElements.append(poly)


		elif len(squareNodesInMesh)+len(squareIntersections) == 5:
			# Check if there is a boundary vertex inside element.
			squareX, squareY = zip(*[gridNodes[n] for n in squareNodes])
			midpoints = [(x,y) for x,y in zip(X,Y) if min(squareX) < x < max(squareX) and min(squareY) < y < max(squareY)]
			if midpoints:
				midpoint = midpoints[0]
			else:
				midpoint = ( (squareIntersections[0][0]+squareIntersections[1][0])/2, (squareIntersections[0][1]+squareIntersections[1][1])/2 )
			meshNodes.append(midpoint)
			midpointIdx = len(meshNodes)-1

			intersectionIdxs = []
			for i in range(2):
				xi, yi = squareIntersections[i]
				nIdx = [ i for i,(xn,yn) in enumerate(meshNodes) if abs(xn-xi)<1e-10 and abs(yn-yi)<1e-10 ]
				if nIdx:
					intersectionIdxs.append( nIdx[0] )						
				else:
					meshNodes.append( (xi,yi) )
					intersectionIdxs.append( len(meshNodes)-1 )

			if gridSquares[squareIdx][0] != squareNodesInMesh[0]:
				hexa = [midpointIdx, intersectionIdxs[0], *[gridToMeshDict[gridIdx] for gridIdx in squareNodesInMesh], intersectionIdxs[1]]
				mpIdx = 0
			else:
				hexa = []
				squareNodeCount = intersectionCount = 0
				for totalCount in range(4):
					if squareNodeCount < len(squareNodesInMesh) and gridSquares[squareIdx][totalCount] == squareNodesInMesh[squareNodeCount]:
						hexa.append( gridToMeshDict[ squareNodesInMesh[squareNodeCount] ] )
						squareNodeCount += 1
					elif intersectionCount < len(squareIntersections):
						hexa += [ intersectionIdxs[0], midpointIdx, intersectionIdxs[1] ]
						mpIdx = totalCount+1
						intersectionCount += 1

			poly1 = hexa[mpIdx-3:mpIdx+1] if hexa[mpIdx-3:mpIdx+1] else hexa[mpIdx-3:] + hexa[:mpIdx+1]
			poly2 = hexa[mpIdx:mpIdx-2]   if hexa[mpIdx:mpIdx-2]   else hexa[mpIdx:]   + hexa[:mpIdx-2]

			mesh2DElements += [poly1, poly2]


	# Create boundaries elements
	boundaries = []
	for idx, name in enumerate(boundaryNames):
		xA, yA, xB, yB = X[idx], Y[idx], X[(idx+1)%len(X)], Y[(idx+1)%len(X)]
		bNodes = [ idx for idx, node in enumerate(meshNodes) if belongsToPolygon(*node, (xA,xB), (yA,yB)) ]
		bNodes = sorted(bNodes, key=lambda i:meshNodes[i][0])
		boundaries.append(bNodes)


	if preview:
		polys = [p for p in mesh2DElements if not p in originalOne ]
		for poly in polys:
			pX,pY = zip(*[meshNodes[v] for v in poly])
			plt.plot(pX+(pX[0],), pY+(pY[0],),color="red",linewidth=2)
		nX,nY = zip(*meshNodes)
		plt.scatter(nX, nY)
		nX,nY = zip(*innerNodes)
		plt.scatter(nX, nY)
		plt.scatter(X,Y,color="k")
		plt.show()



	# Set the elements list in the MSH format
	elements = []
	for bIdx, boundary in enumerate(boundaries):
		elements += [[len(elements)+i+1, *(1, 2), *(bIdx+2, bIdx+2), *( boundary[i]+1, boundary[i+1]+1 )] for i in range(len(boundary)-1) ]
	elements += [[len(elements)+i+1, *(len(element)-1, 2), *(1,1), *[nIdx+1 for nIdx in element]] for i, element in enumerate(mesh2DElements) ]

	# Set the physical names ( which are region and boundaries names )
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



generateMesh(X,Y,dx,dy,["a", "b", "c", "d", "e"],outputPath)








##################################
