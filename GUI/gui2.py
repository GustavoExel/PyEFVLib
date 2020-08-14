import sys,os
sys.path.append(os.path.join(os.path.dirname(__file__), os.path.pardir))
from PyEFVLib import MSHReader, Grid, ProblemData, CgnsSaver, CsvSaver, NeumannBoundaryCondition, DirichletBoundaryCondition
from apps.heat_transfer import heatTransfer

import tkinter as tk
from tkinter import filedialog, ttk

import numpy as np

class Application:
	def __init__(self):
		self.root = tk.Tk()
		self.root.title("PyEFVLib GUI")
		self.root.bind("<Key>", lambda key: self.root.destroy() if key.char=="\x17" else 0) # Close window if Ctrl+W is pressed

		self.HEIGHT, self.WIDTH = (500, 600)

		self.mainPage 		= MainPage(self, self.root)
		self.bcPage 		= BCPage(self, self.root)
		self.propertiesPage = PropertiesPage(self, self.root)
		self.postPage 		= PostPage(self, self.root)

		# Later it'll be mainPage
		self.bcPage.show()

		self.root.mainloop()

	def getData(self):
		self.unitConversions = {
			"K" : lambda T: T,
			"°C": lambda T: T+273,
			"°F": lambda T: (T+459.67)*5/9,
			"K/m": lambda D: D,
			"kg/m³": lambda d: d,
			"g/cm³": lambda d: 1000.0*d,
			"J/kg.K": lambda Cp: Cp,
			"W/m.K": lambda k: k,
			"K/m³": lambda q: q
		}

		# IT IS STILL NECESSARY TO CONVERT ALL UNITS

		self.bcData = dict()
		# Boundary Conditions
		for bName, entry, unit, bType in zip(self.bcPage.boundariesNames, self.bcPage.boundaryValueEntries, self.bcPage.boundaryUnitVars, self.bcPage.boundaryTypeVars):
			self.bcData[bName] = {
				"condition": ["NEUMANN", "DIRICHLET"][bType.get()-1],
				"type": "CONSTANT",
				"value": float( entry.get() )
			}

		# Property Data
		self.propertyValues = dict()
		for region in self.propertiesPage.propertyEntries.keys():
			self.propertyValues[region] = dict()
			for _property in self.propertiesPage.propertyEntries[region].keys():
				self.propertyValues[region][_property] = float( self.propertiesPage.propertyEntries[region][_property].get() )

		# print(self.propertiesPage.propertyUnitVars)

	def dumpData(self):
		pass

	def runSimulation(self):
		grid = Grid( MSHReader(self.bcPage.fileLabel["text"]).getData() )

		heatTransfer(
			libraryPath = os.path.join(os.path.dirname(__file__), os.path.pardir),
			outputPath = os.path.join(os.path.dirname(__file__), os.path.pardir, "results", "gui"),
			extension = "csv",
			
			grid 	  = grid,
			propertyData = [self.propertyValues["Body"]],#################
			
			initialValues = {"temperature": np.zeros(grid.vertices.size)},
			neumannBoundaries = {"temperature":[NeumannBoundaryCondition(grid, boundary, self.bcData[boundary.name]["value"], handle) for handle, boundary in enumerate(grid.boundaries) if self.bcData[boundary.name]["condition"] == "NEUMANN"]},
			dirichletBoundaries = {"temperature":[DirichletBoundaryCondition(grid, boundary, self.bcData[boundary.name]["value"], handle) for handle, boundary in enumerate(grid.boundaries) if self.bcData[boundary.name]["condition"] == "DIRICHLET"]},

			timeStep  = 10000.0,		#####################
			finalTime = 1e+06,		#####################
			maxNumberOfIterations = 1e+04,		#####################
			tolerance = 1e-06,		#####################
			
			# transient = not "-p" in sys.argv,
			# verbosity = not "-s" in sys.argv
		)

class Page:
	def __init__(self, app, root):
		self.app = app
		self.root = root

		self.HEIGHT, self.WIDTH = self.app.HEIGHT, self.app.WIDTH

		self.populated = False

	def populate(self):
		self.populated = True

	def show(self):
		if not self.populated:
			self.populate()
		else:
			self.canvas.pack(side="top", fill="both", expand="yes")

class MainPage(Page):
	def __init__(self, app, root):
		Page.__init__(self, app, root)

	def populate(self):
		pass

	def show(self):
		pass

class BCPage(Page):
	def __init__(self, app, root):
		Page.__init__(self, app, root)

	def populate(self):
		self.canvas = tk.Canvas(self.root, height=self.HEIGHT, width=self.WIDTH)
		self.canvas.pack(side="top", fill="both", expand="yes")

		self.openFrame = tk.LabelFrame(self.canvas, text="Open Mesh File")
		self.openFrame.place(relx=0.02, rely=0.02, relheight=0.18, relwidth=0.96, anchor="nw")

		self.openButton = tk.Button(self.openFrame, text="Open Mesh File", command=self.openFile)
		self.openButton.place(relx=0.02, rely=0.5, relheight=0.50, relwidth=0.25, anchor="w")

		self.fileLabel = tk.Label(self.openFrame, text="", bg="white", anchor="e")
		self.fileLabel.place(relx=0.30, rely=0.5, relheight=0.50, relwidth=0.65, anchor="w")

		self.BCFrame = tk.LabelFrame(self.canvas, text="Boundary Conditions Settings")
		self.BCFrame.place(relx=0.02, rely=0.22, relheight=0.63, relwidth=0.96, anchor="nw")

		# Footer
		self.bottomFrame = tk.Frame(self.canvas)
		self.bottomFrame.place(relx=0.02, rely=0.87, relheight=0.1, relwidth=0.96, anchor="nw" )

		self.nextButton = tk.Button(self.bottomFrame, text="Next", command=self.next)
		self.nextButton.place(relx=1.0, rely=0.50, relheight=0.75, relwidth=0.20, anchor="e")

		self.prevButton = tk.Button(self.bottomFrame, text="Prev", command=self.prev)
		self.prevButton.place(relx=0.78, rely=0.50, relheight=0.75, relwidth=0.20, anchor="e")

		self.populated = True

	def placeBCForms(self):
		self.boundariesNames = MSHReader(self.fileLabel["text"]).getData().boundariesNames
				
		BCCanvas = tk.Canvas(self.BCFrame)
		BCCanvas.place(relx=0.0, rely=0.0, relwidth=1.0, relheight=1.0, anchor="nw")

		scrollbar = ttk.Scrollbar(self.BCFrame, orient="vertical", command=BCCanvas.yview)
		scrollbar.place(relx=1.0, rely=0.0, relheight=1.0, anchor="ne")

		# with Windows OS
		self.root.bind("<MouseWheel>", lambda e: BCCanvas.yview_scroll(-1*int(e.delta/120), "units"))
		# with Linux OS
		self.root.bind("<Button-4>", lambda e: BCCanvas.yview_scroll(-1, "units"))
		self.root.bind("<Button-5>", lambda e: BCCanvas.yview_scroll(+1, "units"))

		BCWindow = ttk.Frame(BCCanvas)
		BCWindow.bind("<Configure>", lambda event: BCCanvas.configure(scrollregion=BCCanvas.bbox("all")))

		BCCanvas.create_window((0, 0), window=BCWindow, anchor="nw")
		BCCanvas.configure(yscrollcommand=scrollbar.set)

		BCWindow.columnconfigure(index=0, pad=5)
		BCWindow.columnconfigure(index=1, pad=5)
		BCWindow.columnconfigure(index=2, pad=5)
		BCWindow.columnconfigure(index=3, pad=5)
		BCWindow.rowconfigure(index=0, pad=5)
		BCWindow.rowconfigure(index=1, pad=5)
		BCWindow.rowconfigure(index=2, pad=5)
		BCWindow.rowconfigure(index=3, pad=5)

		ttk.Label(BCWindow, text="Boundary Name").grid(row=0, column=0)
		ttk.Label(BCWindow, text="Value").grid(row=0, column=1)
		ttk.Label(BCWindow, text="Unit").grid(row=0, column=2)
		ttk.Label(BCWindow, text="  Neumann BC  ").grid(row=0, column=3)
		ttk.Label(BCWindow, text="  Dirichlet BC  ").grid(row=0, column=4)

		# Obs.: If 3D then neumann units are [Temp]/[Distance^2]
		self.neumannUnits = ["K/m"]
		self.dirichletUnits = ["K", "°C", "°F"]

		self.boundaryUnitVars		  = []
		self.boundaryTypeVars		  = []
		self.boundaryValueEntries	  = []
		self.boundaryUnitMenus		  = []

		i=1
		for boundaryName in self.boundariesNames:
			options = []

			unitVar = tk.StringVar(BCWindow)
			unitVar.set("")
			bTypeVar = tk.IntVar(BCWindow)

			nameLabel = ttk.Label(BCWindow, text=boundaryName)
			nameLabel.grid(row=i, column=0, padx=5, sticky="W")

			valEntry = tk.Entry(BCWindow)
			valEntry.grid(row=i,column=1)

			unitMenu = ttk.OptionMenu(BCWindow, unitVar, *options)
			unitMenu.grid(row=i, column=2)

			neumannButton = tk.Radiobutton(BCWindow, variable=bTypeVar, value=1, command=self.setNeumann(i-1))
			neumannButton.grid(row=i, column=3)

			dirichletButton = tk.Radiobutton(BCWindow, variable=bTypeVar, value=2, command=self.setDirichlet(i-1))
			dirichletButton.grid(row=i, column=4)

			self.boundaryUnitVars.append(unitVar)
			self.boundaryTypeVars.append(bTypeVar)
			self.boundaryValueEntries.append(valEntry)
			self.boundaryUnitMenus.append(unitMenu)
			i+=1

	def setNeumann(self, i):
		def setNeumannFunc():
			unitVar = self.boundaryUnitVars[i]
			optionMenu = self.boundaryUnitMenus[i]

			unitVar.set("")
			optionMenu["menu"].delete(0,"end") # remove full list 
			for opt in self.neumannUnits: 
				optionMenu['menu'].add_command(label=opt, command=tk._setit(unitVar, opt))
			unitVar.set(self.neumannUnits[0]) # default value set 

		return setNeumannFunc

	def setDirichlet(self, i):
		def setDirichletFunc():
			unitVar = self.boundaryUnitVars[i]
			optionMenu = self.boundaryUnitMenus[i]

			unitVar.set("")
			optionMenu["menu"].delete(0,"end") # remove full list 
			for opt in self.dirichletUnits: 
				optionMenu['menu'].add_command(label=opt, command=tk._setit(unitVar, opt))
			unitVar.set(self.dirichletUnits[0]) # default value set 
		return setDirichletFunc

	def openFile(self):
		initialdir = os.path.join( os.path.dirname(__file__), os.path.pardir, "meshes" )
		self.fileLabel["text"] = tk.filedialog.askopenfilename(initialdir=initialdir, title="Select a mesh file", filetypes=(("MSH files", "*.msh"),("All files", "*")))
		self.placeBCForms()

	def next(self):
		self.canvas.pack_forget()

		self.app.propertiesPage.show()

	def prev(self):
		pass

class PropertiesPage(Page):
	def __init__(self, app, root):
		Page.__init__(self, app, root)

	def populate(self):
		filePath = self.app.bcPage.fileLabel["text"]
		regionNames = MSHReader(filePath).getData().regionsNames
		print(regionNames)

		self.canvas = tk.Canvas(self.root, height=self.HEIGHT, width=self.WIDTH)
		self.canvas.pack(side="top", fill="both", expand="yes")

		self.bottomFrame = tk.Frame(self.canvas)
		self.bottomFrame.place(relx=0.02, rely=0.87, relheight=0.1, relwidth=0.96, anchor="nw" )

		self.nextButton = tk.Button(self.bottomFrame, text="RUN", command=self.next)
		self.nextButton.place(relx=1.0, rely=0.50, relheight=0.75, relwidth=0.20, anchor="e")

		self.prevButton = tk.Button(self.bottomFrame, text="Prev", command=self.prev)
		self.prevButton.place(relx=0.78, rely=0.50, relheight=0.75, relwidth=0.20, anchor="e")

		self.properties = ["Density","HeatCapacity","Conductivity","HeatGeneration"]
		self.propertyUnits = {
			"Density":       ["kg/m³", "g/cm³"],
			"HeatCapacity":   ["J/kg.K"],
			"Conductivity":   ["W/m.K"],
			"HeatGeneration": ["K/m³"]
		}

		self.propertyEntries = dict()
		self.propertyUnitVars = dict()
		self.propertiesFrames = []
		self.currentRegion = 0

		def nextRegion():
			if self.currentRegion < len(regionNames)-1:
				self.propertiesFrames[self.currentRegion].place_forget()
				self.propertiesFrames[self.currentRegion+1].place(relx=0.02, rely=0.02, relheight=0.81, relwidth=0.96, anchor="nw")
				self.currentRegion += 1

		def prevRegion():
			if self.currentRegion > 0:
				self.propertiesFrames[self.currentRegion].place_forget()
				self.propertiesFrames[self.currentRegion-1].place(relx=0.02, rely=0.02, relheight=0.81, relwidth=0.96, anchor="nw")
				self.currentRegion -= 1

		for regionCount, region in enumerate(regionNames):
			propertiesFrame = tk.LabelFrame(self.canvas, text="Material Properties")
			self.propertiesFrames.append(propertiesFrame)
			centerFrame = tk.Frame(propertiesFrame)

			self.propertyEntries[region] = dict()
			self.propertyUnitVars[region] = dict()

			self.regionCountLabel = tk.Label(centerFrame, text=f"{region}\t[{regionCount+1}/{len(regionNames)}]")
			self.regionCountLabel.grid(row=0, column=0, pady=5, sticky="W")

			i=1
			for propertyName in self.properties:
				options = self.propertyUnits[propertyName]

				unitVar = tk.StringVar(centerFrame)
				unitVar.set(self.propertyUnits[propertyName][0])
				bTypeVar = tk.IntVar(centerFrame)

				nameLabel = ttk.Label(centerFrame, text=propertyName)
				nameLabel.grid(row=i, column=0, sticky="W", padx=5)

				valEntry = ttk.Entry(centerFrame)
				valEntry.grid(row=i,column=1)

				unitMenu = tk.OptionMenu(centerFrame, unitVar, *options)
				unitMenu.grid(row=i, column=2, sticky="E")


				self.propertyEntries[region][propertyName] = valEntry
				self.propertyUnitVars[region][propertyName] = unitVar
				i+=1

			prevButton = tk.Button(centerFrame, text="<", command=prevRegion)
			prevButton.grid(row=i, column=2, sticky="E")

			nextButton = tk.Button(centerFrame, text=">", command=nextRegion)
			nextButton.grid(row=i, column=3)

			centerFrame.place(relx=0.5, rely=0.0, anchor="n")

		self.propertiesFrames[0].place(relx=0.02, rely=0.02, relheight=0.81, relwidth=0.96, anchor="nw")
		
		self.populated = True

	def next(self):
		print("Running Simulation...")
		self.app.getData()
		self.app.runSimulation()

	def prev(self):
		self.canvas.pack_forget()
		self.app.bcPage.show()

class PostPage(Page):
	def __init__(self, app, root):
		Page.__init__(self, app, root)

	def populate(self):
		pass

	def show(self):
		pass


if __name__ == "__main__":
	app = Application()
	# app.root.mainloop()