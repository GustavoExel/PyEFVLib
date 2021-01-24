def stressEquilibriumBoundaryConditionsPrinter(bcs):
	boundaryConditions={ bc["u"].boundary.name:bc for bc in bcs }
	values=[
		-boundaryConditions["North"]["u"].value,
		+boundaryConditions["South"]["u"].value,
		+boundaryConditions["West"]["u"].value,
		-boundaryConditions["East"]["u"].value,
		-boundaryConditions["North"]["v"].value,
		+boundaryConditions["South"]["v"].value,
		-boundaryConditions["West"]["v"].value,
		-boundaryConditions["East"]["v"].value
	]
	units=[
		"mm" if boundaryConditions["North"]["u"].__type__ == "DIRICHLET" else "MPa",
		"mm" if boundaryConditions["South"]["u"].__type__ == "DIRICHLET" else "MPa",
		"mm" if boundaryConditions["West"]["u"].__type__ == "DIRICHLET" else "MPa",
		"mm" if boundaryConditions["East"]["u"].__type__ == "DIRICHLET" else "MPa",
		"mm" if boundaryConditions["North"]["v"].__type__ == "DIRICHLET" else "MPa",
		"mm" if boundaryConditions["South"]["v"].__type__ == "DIRICHLET" else "MPa",
		"mm" if boundaryConditions["West"]["v"].__type__ == "DIRICHLET" else "MPa",
		"mm" if boundaryConditions["East"]["v"].__type__ == "DIRICHLET" else "MPa"
	]
	factors=[
		1e+3 if boundaryConditions["North"]["u"].__type__ == "DIRICHLET" else 1e-6,
		1e+3 if boundaryConditions["South"]["u"].__type__ == "DIRICHLET" else 1e-6,
		1e+3 if boundaryConditions["West"]["u"].__type__ == "DIRICHLET" else 1e-6,
		1e+3 if boundaryConditions["East"]["u"].__type__ == "DIRICHLET" else 1e-6,
		1e+3 if boundaryConditions["North"]["v"].__type__ == "DIRICHLET" else 1e-6,
		1e+3 if boundaryConditions["South"]["v"].__type__ == "DIRICHLET" else 1e-6,
		1e+3 if boundaryConditions["West"]["v"].__type__ == "DIRICHLET" else 1e-6,
		1e+3 if boundaryConditions["East"]["v"].__type__ == "DIRICHLET" else 1e-6
	]
	txt="             {4:>13}  |                            \n                            |     {0:<13}          \n                            v<-----                      \n                +---------------------+                  \n                |                     |                  \n   {6:>10}  ||                     || {7:<9}       \n               ||                     ||                 \n               v|                     |v                 \n{2:>9}  ---->|                     |<----   {3:<10}\n                |                     |                  \n                |                     |                  \n                |                     |                  \n                |                     |                  \n                |                     |                  \n                +---------------------+                  \n                      ----> ^                            \n          {1:>13}     |                            \n                            | {5:<12}               \n"
	print(txt.format(*[ str(round(fac*val,4))+unit for val,fac,unit in zip(values,factors,units) ]))