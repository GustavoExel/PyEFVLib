import PyEFVLib
from PyEFVLib import ProblemData, NumericalSettings, PropertyData, BoundaryConditions

problemData =ProblemData(
	meshFilePath  = "{MESHES}/msh/2D/Square.msh",		# required
	outputFilePath = "{RESULTS}/Results.xdmf",	# optional

	numericalSettings = NumericalSettings(
		timeStep = 1e-02,							# required
		finalTime = None,							# optional
		tolerance = 1e-06,							# optional
		maxNumberOfIterations = 1000, 				# optional, but with default
	),

	propertyData = PropertyData({
		"Body" : {
			"HeatCapacity"	: 1.0,
			"Conductivity"	: 1.0,
			"Density"		: 1.0,
		},
	}),

	boundaryConditions = BoundaryConditions({
		"temperature": {
			"InitialValue": 0.0,
			"West": {
				"condition" : PyEFVLib.Neumann,
				"type"		: PyEFVLib.Constant,
				"value"		: 100.0
			},
		    "East": {
		        "condition"	: PyEFVLib.Dirichlet,
		        "type"		: PyEFVLib.Constant,
		        "value"		: 0.0
		    },
		    "South": {
		        "condition"	: PyEFVLib.Neumann,
		        "type"		: PyEFVLib.Constant,
		        "value"		: 0.0
		    },
		    "North": {
		        "condition"	: PyEFVLib.Neumann,
		        "type"		: PyEFVLib.Constant,
		        "value"		: 0.0
		    },
		},
	}),
)