for boundary in boundaries:
    if boundary["u"].type == "Dirichlet" && boundary["v"].type == "Dirichlet":
        for vertex in boundary.vertices:
            matrix["u"][vertex] = [0] * (2*numberOfVertices)
            matrix["v"][vertex] = [0] * (2*numberOfVertices)
            matrix["u"][vertex]["u"][vertex] = 1.0
            matrix["v"][vertex]["v"][vertex] = 1.0
            independent["u"][vertex] = u_boundaryCondition
            independent["v"][vertex] = v_boundaryCondition

    elif boundary["u"].type == "Neumann" && boundary["v"].type == "Neumann":
        for vertex in boundary.vertices:
            independent["u"][vertex]=0.0
            independent["v"][vertex]=0.0
            matrix["u"][vertex]=[0] * (2*numberOfVertices)
            matrix["v"][vertex]=[0] * (2*numberOfVertices)

        for outerFace in boundary.outerFaces:
            xNormalStress=sx*(Tx*sx+Ty*sy)/(sx**2+sy**2)
            yNormalStress=sy*(Tx*sx+Ty*sy)/(sx**2+sy**2)
            shearStress=sqrt((Tx-xNormalStress)**2+(Ty-yNormalStress)**2)

            independent["u"] += sx*xNormalStress+sy*shearStress
            independent["v"] += sx*shearStress+sy*yNormalStress



    elif boundary["u"].type == "Dirichlet" && boundary["v"].type == "Neumann":
    elif boundary["u"].type == "Neumann" && boundary["v"].type == "Dirichlet":


for boundary in boundaries:
    if boundary["u"].type == "Dirichlet" && boundary["v"].type == "Dirichlet":
        for vertex in boundary.vertices:
            matrix["u"][vertex] = [0] * (2*numberOfVertices)
            matrix["v"][vertex] = [0] * (2*numberOfVertices)
            independent["u"][vertex] = u_boundaryCondition
            independent["v"][vertex] = v_boundaryCondition

    Sx, Sy = outerFace.area
    elif boundary["u"].type == "Neumann" && boundary["v"].type == "Neumann":
        for vertex in boundary.vertices:
            independent["u"][vertex]=0.0
            independent["v"][vertex]=0.0
        for outerFace in boundary.outerFaces:
            xNormalStress=Sx*(Tx*Sx+Ty*Sy)/(Sx**2+Sy**2)
            yNormalStress=Sy*(Tx*Sx+Ty*Sy)/(Sx**2+Sy**2)
            shearStress=sqrt((Tx-xNormalStress)**2+(Ty-yNormalStress)**2)

            independent["u"] += Sx*xNormalStress+Sy*shearStress
            independent["v"] += Sx*shearStress+Sy*yNormalStress

    elif boundary["u"].type == "Dirichlet" && boundary["v"].type == "Neumann":
        for vertex in boundary.vertices:
            matrix["u"][vertex] = [0] * (2*numberOfVertices)
            independent["u"][vertex] = u_boundaryCondition
            independent["v"][vertex] = 0.0

        for outerFace in boundary.outerFaces:
            yNormalStress=Sy*(Tx*Sx+Ty*Sy)/(Sx**2+Sy**2)
            independent["v"] += Sy*yNormalStress

        for vertex in outerFace.element.vertices:
            Nx,Ny=outerFace.globalShapeFunctionDerivatives[vertex]
            matrix["u"][outerFace][vertex] -= Sy*shearModulus*Ny
            matrix["v"][outerFace][vertex] -= Sy*shearModulus*Nx

    elif boundary["u"].type == "Neumann" && boundary["v"].type == "Dirichlet":
        for vertex in boundary.vertices:
            matrix["v"][vertex] = [0] * (2*numberOfVertices)
            independent["v"][vertex] = v_boundaryCondition
            independent["u"][vertex] = 0.0

        for outerFace in boundary.outerFaces:
            xNormalStress=Sx*(Tx*Sx+Ty*Sy)/(Sx**2+Sy**2)
            independent["u"] += Sx*xNormalStress
        for vertex in outerFace.element.vertices:
            Nx,Ny=outerFace.globalShapeFunctionDerivatives[vertex]
            matrix["u"][outerFace][vertex] -= Sy*shearModulus*Ny
            matrix["v"][outerFace][vertex] -= Sy*shearModulus*Nx