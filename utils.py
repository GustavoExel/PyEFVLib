def c():
	os.system("clear")

def MM():
	bi=[vertex.handle for vertex in boundary.vertices]
	mm=[[c for j,c in enumerate(line) if j in bi] for i,line in enumerate(matrix) if i in bi]
	return np.array(mm)

def print2DArray(A):
	t=""
	l=[]
	for line in A:
		l.append(', '.join([f'{ln:{">10.3e" if ln!=int(ln) else ">10"}}' for ln in line]))
		# l.append(', '.join([str(round(ln,3)) for ln in line]))
		l[-1] = "[%s]" % l[-1]
	t=',\n'.join(l)
	t="[%s]" % t
	print(t)

def f():
	II=0
	boundary=bc["u"].boundary
	uBoundaryType = bc["u"].__type__
	vBoundaryType = bc["v"].__type__
	if uBoundaryType != vBoundaryType:
		dirichletLabel = "u" if uBoundaryType=="DIRICHLET" else "v"
		neumannLabel = "u" if uBoundaryType=="NEUMANN" else "v"
		dirichletIndex = 0 if uBoundaryType=="DIRICHLET" else 1
		neumannIndex = 0 if uBoundaryType=="NEUMANN" else 1
		D=lambda vertex: int( vertex.handle + numberOfVertices * dirichletIndex )
		N=lambda vertex: int( vertex.handle + numberOfVertices * neumannIndex )
		for vertex in boundary.vertices:
			matrix[D(vertex)] = np.zeros(2*numberOfVertices)
			matrix[D(vertex)][D(vertex)] = 1.0
			independent[D(vertex)] = bc[dirichletLabel].getValue(vertex.handle)
			matrix[N(vertex)] = np.zeros(2*numberOfVertices)
			independent[N(vertex)] = 0.0
		print(II, "Must be empty")
		print2DArray(MM())
		print("\n\n\n")
		for facet in boundary.facets:
			shearModulus = problemData.propertyData[region.handle]["ShearModulus"]
			constitutiveMatrix = getConstitutiveMatrix(facet.element.region)
			for outerFace in facet.outerFaces:
				globalDerivatives = getOuterFaceGlobalDerivatives(outerFace)
				Nx, Ny = globalDerivatives
				Sx,Sy,_ = outerFace.area.getCoordinates()
				normalArea = [Sx,Sy][neumannIndex]
				shearArea = [Sx,Sy][dirichletIndex]
				transposedVoigtArea = getTransposedVoigtArea(outerFace)
				voigtGradientOperator = getVoigtGradientOperator(globalDerivatives)
				coefficient = np.einsum("ij,jk,kmn->imn", transposedVoigtArea, constitutiveMatrix, voigtGradientOperator)
				local=0
				for vertex in facet.element.vertices:
					print(f"Adding {coefficient[neumannIndex][0][local]} and subtracting {shearArea * shearModulus * Ny[local]}")
					matrix[N(outerFace.vertex)][U(vertex.handle)] += coefficient[neumannIndex][0][local]
					matrix[N(outerFace.vertex)][V(vertex.handle)] += coefficient[neumannIndex][1][local]
					matrix[N(outerFace.vertex)][U(vertex.handle)] -= shearArea * shearModulus * Ny[local]
					matrix[N(outerFace.vertex)][V(vertex.handle)] -= shearArea * shearModulus * Nx[local]
					independent[N(outerFace.vertex)] += normalArea * bc[neumannLabel].getValue( outerFace.vertex.handle )
					local+=1
					print(f"II={II}, local={local}")
					print2DArray(MM())
					print("\n\n---\n")