import numpy as np
import meshio

load_width = 1
width  = 16
height = 8
sdx = 2
sdy = 2

# This values are the values of the pixels from an image
# Since the image was of a FVM, the elements are being defined between the volumes
volXs = [132,137,141,145,150,154,158,162,167,171,175,180,184,188,192.5,197,201,205,210,214,219,224,230,240,252,267,284,304,327.5,352,379,412,459,526,620,727,834,937,987]
volXs = np.array([x for x0,x1 in zip(volXs[:-1],volXs[1:]) for x in np.linspace(x0,x1,sdx)])
volXs -= volXs[0]; volXs *= width/volXs[-1]
mdPtsX = 0.5*(volXs[:-1]+volXs[1:])[1:-1]
x = np.array( [volXs[0]] + list(mdPtsX) + [volXs[-1]] )

volYs = [11,18,24.5,31.5,38,45,52,58,65,73,81,90,98,107,115,123.5,132,140,149,157,166,174,182,191,199.5,208,216,225,234,245,256,269,284,301,320,343,373,411,453,495,537,579,622,664,685]
volYs = np.array([y for y0,y1 in zip(volYs[:-1],volYs[1:]) for y in np.linspace(y0,y1,sdy)])
volYs -= volYs[-1]; volYs = volYs[::-1]; volYs *= height/volYs[-1]
mdPtsY = 0.5*(volYs[:-1]+volYs[1:])[1:-1]
y = np.array( [volYs[0]] + list(mdPtsY) + [volYs[-1]] )

nx = len(x) - 1
ny = len(y) - 1
nl = min(enumerate(abs((volXs-load_width))),key=lambda t:t[1])[0]

xv,yv = np.meshgrid(x,y)
xv.resize(xv.size); yv.resize(yv.size)

points = np.array(list(zip(xv,yv)))
points = np.concatenate((points.T,[np.zeros(len(points))])).T

idx = lambda i,j: i+(nx+1)*j
quad_cells = np.array([[idx(i,j),idx(i+1,j),idx(i+1,j+1),idx(i,j+1)] for j in range(ny) for i in range(nx)])
line_cells = np.array([ [idx(*(p,q)[ ::(1,-1)[k] ]), idx(*(p+1,q)[ ::(1,-1)[k] ])] for k in range(2) for q in (0,(ny,nx)[k]) for p in range((nx,ny)[k]) ])

quad_data = np.ones(len(quad_cells), dtype=np.int32)
line_data = np.array(nx*[2]+nl*[3]+(nx-nl)*[4]+ny*[5]+ny*[6], dtype=np.int32)

cells = [("line",line_cells), ("quad",quad_cells),]
cell_data = {"gmsh:physical":[line_data, quad_data],"gmsh:geometrical":[line_data, quad_data]}

field_data = {"Body":[1,2], "South":[2,1], "Load":[3,1], "North":[4,1], "West":[5,1], "East":[6,1]}

mesh = meshio.Mesh(points, cells, cell_data=cell_data, field_data=field_data)

mesh.write("output.msh", file_format="gmsh22", binary=False)