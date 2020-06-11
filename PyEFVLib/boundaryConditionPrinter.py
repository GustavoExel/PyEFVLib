def stressEquilibriumBoundaryConditionsPrinter(bcs):
	boundaryConditions={ bc["u"].boundary.name:bc for bc in bcs }
	values=[
		boundaryConditions["NORTH"]["u"].value,
		boundaryConditions["SOUTH"]["u"].value,
		boundaryConditions["WEST"]["u"].value,
		boundaryConditions["EAST"]["u"].value,
		boundaryConditions["NORTH"]["v"].value,
		boundaryConditions["SOUTH"]["v"].value,
		boundaryConditions["WEST"]["v"].value,
		boundaryConditions["EAST"]["v"].value
	]
	units=[
		"mm" if boundaryConditions["NORTH"]["u"].__type__ == "DIRICHLET" else "MPa",
		"mm" if boundaryConditions["SOUTH"]["u"].__type__ == "DIRICHLET" else "MPa",
		"mm" if boundaryConditions["WEST"]["u"].__type__ == "DIRICHLET" else "MPa",
		"mm" if boundaryConditions["EAST"]["u"].__type__ == "DIRICHLET" else "MPa",
		"mm" if boundaryConditions["NORTH"]["v"].__type__ == "DIRICHLET" else "MPa",
		"mm" if boundaryConditions["SOUTH"]["v"].__type__ == "DIRICHLET" else "MPa",
		"mm" if boundaryConditions["WEST"]["v"].__type__ == "DIRICHLET" else "MPa",
		"mm" if boundaryConditions["EAST"]["v"].__type__ == "DIRICHLET" else "MPa"
	]
	factors=[
		1e+3 if boundaryConditions["NORTH"]["u"].__type__ == "DIRICHLET" else 1e-6,
		1e+3 if boundaryConditions["SOUTH"]["u"].__type__ == "DIRICHLET" else 1e-6,
		1e+3 if boundaryConditions["WEST"]["u"].__type__ == "DIRICHLET" else 1e-6,
		1e+3 if boundaryConditions["EAST"]["u"].__type__ == "DIRICHLET" else 1e-6,
		1e+3 if boundaryConditions["NORTH"]["v"].__type__ == "DIRICHLET" else 1e-6,
		1e+3 if boundaryConditions["SOUTH"]["v"].__type__ == "DIRICHLET" else 1e-6,
		1e+3 if boundaryConditions["WEST"]["v"].__type__ == "DIRICHLET" else 1e-6,
		1e+3 if boundaryConditions["EAST"]["v"].__type__ == "DIRICHLET" else 1e-6
	]
	txt="             {4:>13}  |                            \n                            |     {0:<13}          \n                            v<-----                      \n                +---------------------+                  \n                |                     |                  \n   {6:>10}  ||                     || {7:<9}       \n               ||                     ||                 \n               v|                     |v                 \n{2:>9}  ---->|                     |<----   {3:<10}\n                |                     |                  \n                |                     |                  \n                |                     |                  \n                |                     |                  \n                |                     |                  \n                +---------------------+                  \n                      ----> ^                            \n          {1:>13}     |                            \n                            | {5:<12}               \n"
	print(txt.format(*[ str(round(fac*val,4))+unit for val,fac,unit in zip(values,factors,units) ]))