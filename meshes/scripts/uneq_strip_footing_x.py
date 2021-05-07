import numpy as np
import meshio

load_width = 1
width  = 16
height = 8

# This values are the values of the pixels from an image
# Since the image was of a FVM, the elements are being defined between the volumes
volXs = np.array([36,44,59,74,88,103,117,132,147,163,183,205,227,251,278,308,339,374,410,456,520,590,623], dtype=np.float64)
volXs -= volXs[0]; volXs *= width/volXs[-1]
mdPts = 0.5*(volXs[:-1]+volXs[1:])[1:-1]

x = np.array( [volXs[0]] + list(mdPts) + [volXs[-1]] )

# import matplotlib.pyplot as plt
# plt.scatter(volXs, np.ones(volXs.size), color='k', marker='.')
# plt.scatter(mdPts, np.ones(mdPts.size), color='r', marker='.')
# plt.show()

nx = len(x) - 1
ny = 80
nl = min(enumerate(abs((volXs-load_width))),key=lambda t:t[1])[0]

y = np.linspace(0,height,ny+1)

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