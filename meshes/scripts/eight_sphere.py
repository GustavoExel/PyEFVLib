import numpy as np
import meshio

R = 1
nr = 10
na = 10

rs = np.linspace(0,R,nr+1)[1:]
ts = np.linspace(0,np.pi/2,na)
ps = np.linspace(0,np.pi/2,na)

xv = np.array([r * np.cos(p) * np.sin(t) for r in rs for p in ps for t in ts])
yv = np.array([r * np.sin(p) * np.sin(t) for r in rs for p in ps for t in ts])
zv = np.array([r * np.cos(t) for r in rs for p in ps for t in ts])

points = np.array(list(zip(xv,yv,zv)))

idx = lambda ri,pi,ti: ti + na*pi + na*na*ri
hexa_cells = np.array([[idx(ri,pi,ti),idx(ri,pi,ti+1),idx(ri,pi+1,ti+1),idx(ri,pi+1,ti),idx(ri+1,pi,ti),idx(ri+1,pi,ti+1),idx(ri+1,pi+1,ti+1),idx(ri+1,pi+1,ti)] for ri in range(nr-1) for pi in range(na-1) for ti in range(na-1)])
hexa_data = np.ones(len(hexa_cells), dtype=np.int32)

cells = [("hexahedron",hexa_cells),]
cell_data = {"gmsh:physical":[hexa_data],"gmsh:geometrical":[hexa_data]}

mesh = meshio.Mesh(points, cells, cell_data=cell_data)
mesh.write("output.msh", file_format="gmsh22", binary=False)


# idx = lambda i,j: i+(nx+1)*j
# quad_cells = np.array([[idx(i,j),idx(i+1,j),idx(i+1,j+1),idx(i,j+1)] for j in range(ny) for i in range(nx)])
# line_cells = np.array([ [idx(*(p,q)[ ::(1,-1)[k] ]), idx(*(p+1,q)[ ::(1,-1)[k] ])] for k in range(2) for q in (0,(ny,nx)[k]) for p in range((nx,ny)[k]) ])

# quad_data = np.ones(len(quad_cells), dtype=np.int32)
# line_data = np.array(nx*[2]+nx*[3]+ny*[4]+ny*[5], dtype=np.int32)

# cells = [("line",line_cells), ("quad",quad_cells),]
# cell_data = {"gmsh:physical":[line_data, quad_data],"gmsh:geometrical":[line_data, quad_data]}

# field_data = {"Body":[1,2], "South":[2,1], "North":[3,1], "West":[4,1], "East":[5,1]}

# mesh = meshio.Mesh(points, cells, cell_data=cell_data, field_data=field_data)

# mesh.write("output.msh", file_format="gmsh22", binary=False)