import sys,os
sys.path.append(os.path.join(os.path.dirname(__file__), os.path.pardir))
from PyEFVLib import MSHReader, Grid, ProblemData, CgnsSaver, CsvSaver, NeumannBoundaryCondition, DirichletBoundaryCondition
from apps.heat_transfer import heatTransfer

import tkinter as tk
from tkinter import filedialog, ttk, messagebox

import numpy as np

class Application:
	def __init__(self):
		self.root = tk.Tk()
		self.root.title("PyEFVLib GUI")
		self.root.bind("<Key>", lambda key: self.root.destroy() if key.char=="\x17" else 0) # Close window if Ctrl+W is pressed

		self.HEIGHT, self.WIDTH = (500, 600)

		self.heatTransferApplication = HeatTransferApplication(self.root)
		self.solidMechanicsApplication = SolidMechanicsApplication(self.root)

		self.root.mainloop()


class SolidMechanicsApplication:
	def __init__(self, root):
		self.root = root

class HeatTransferApplication:
	def __init__(self, root):
		self.root = root

		self.settings()

		self.populatePage1()
		self.showPage1()

	def settings(self):
		self.neumannUnits = {"temperature": ["K/m"]}
		self.dirichletUnits = {"temperature": ["K", "°C", "°F"]}
		
		self.properties = ["Density","HeatCapacity","Conductivity","HeatGeneration"]
		self.propertyUnits = {
			"Density":       ["kg/m³", "g/cm³"],
			"HeatCapacity":   ["J/kg.K"],
			"Conductivity":   ["W/m.K"],
			"HeatGeneration": ["K/m³"]
		}


		self.meshFileName = ""

	def populatePage1(self):
		self.page1 = tk.Canvas(self.root, height=500, width=600)

		self.openFrame = tk.LabelFrame(self.page1, text="Open Mesh File")
		self.openFrame.place(relx=0.02, rely=0.02, relheight=0.18, relwidth=0.96, anchor="nw")

		self.openButton = tk.Button(self.openFrame, text="Open Mesh File", command=self.openFile)
		self.openButton.place(relx=0.02, rely=0.5, relheight=0.50, relwidth=0.25, anchor="w")

		self.fileLabel = tk.Label(self.openFrame, text="", bg="white", anchor="e")
		self.fileLabel.place(relx=0.30, rely=0.5, relheight=0.50, relwidth=0.65, anchor="w")

		self.BCFrame = tk.LabelFrame(self.page1, text="Boundary Conditions Settings")
		self.BCFrame.place(relx=0.02, rely=0.22, relheight=0.63, relwidth=0.96, anchor="nw")

		# Footer
		self.bottomFrame = tk.Frame(self.page1)
		self.bottomFrame.place(relx=0.02, rely=0.87, relheight=0.1, relwidth=0.96, anchor="nw" )

		self.nextButton = tk.Button(self.bottomFrame, text="Next", command=self.page1Next)
		self.nextButton.place(relx=1.0, rely=0.50, relheight=0.75, relwidth=0.20, anchor="e")

		self.prevButton = tk.Button(self.bottomFrame, text="Prev", command=self.page1Prev)
		self.prevButton.place(relx=0.78, rely=0.50, relheight=0.75, relwidth=0.20, anchor="e")

	def showPage1(self):
		self.page1.pack(side="top", fill="both", expand="yes")

	def populateBCFrame(self):
		self.boundariesNames = self.meshData.boundariesNames
		def chooseNeumannRadioButton(i): #
			def chooseNeumannRadioButtonFunc():
				unitVar = self.boundaryUnitVars[i]
				optionMenu = self.boundaryUnitMenus[i]

				unitVar.set("")
				optionMenu["menu"].delete(0,"end") # remove full list 
				for opt in self.neumannUnits["temperature"]: 
					optionMenu['menu'].add_command(label=opt, command=tk._setit(unitVar, opt))
				unitVar.set(self.neumannUnits["temperature"][0]) # default value set 

			return chooseNeumannRadioButtonFunc
		def chooseDirichletRadioButton(i): #
			def chooseDirichletRadioButtonFunc():
				unitVar = self.boundaryUnitVars[i]
				optionMenu = self.boundaryUnitMenus[i]

				unitVar.set("")
				optionMenu["menu"].delete(0,"end") # remove full list 
				for opt in self.dirichletUnits["temperature"]: 
					optionMenu['menu'].add_command(label=opt, command=tk._setit(unitVar, opt))
				unitVar.set(self.dirichletUnits["temperature"][0]) # default value set 
			return chooseDirichletRadioButtonFunc
		def placeStatic():
			BCCanvas = tk.Canvas(self.BCFrame)
			BCCanvas.place(relx=0.0, rely=0.0, relwidth=1.0, relheight=1.0, anchor="nw")

			scrollbar = ttk.Scrollbar(self.BCFrame, orient="vertical", command=BCCanvas.yview)
			scrollbar.place(relx=1.0, rely=0.0, relheight=1.0, anchor="ne")

			# with Windows OS
			self.root.bind("<MouseWheel>", lambda e: BCCanvas.yview_scroll(-1*int(e.delta/120), "units"))
			# with Linux OS
			self.root.bind("<Button-4>", lambda e: BCCanvas.yview_scroll(-1, "units"))
			self.root.bind("<Button-5>", lambda e: BCCanvas.yview_scroll(+1, "units"))

			self.BCWindow = ttk.Frame(BCCanvas)
			self.BCWindow.bind("<Configure>", lambda event: BCCanvas.configure(scrollregion=BCCanvas.bbox("all")))

			BCCanvas.create_window((0, 0), window=self.BCWindow, anchor="nw")
			BCCanvas.configure(yscrollcommand=scrollbar.set)

			self.BCWindow.columnconfigure(index=0, pad=5)
			self.BCWindow.columnconfigure(index=1, pad=5)
			self.BCWindow.columnconfigure(index=2, pad=5)
			self.BCWindow.columnconfigure(index=3, pad=5)
			self.BCWindow.rowconfigure(index=0, pad=5)
			self.BCWindow.rowconfigure(index=1, pad=5)
			self.BCWindow.rowconfigure(index=2, pad=5)
			self.BCWindow.rowconfigure(index=3, pad=5)

			ttk.Label(self.BCWindow, text="Boundary Name").grid(row=0, column=0)
			ttk.Label(self.BCWindow, text="Value").grid(row=0, column=1)
			ttk.Label(self.BCWindow, text="Unit").grid(row=0, column=2)
			ttk.Label(self.BCWindow, text="  Neumann BC  ").grid(row=0, column=3)
			ttk.Label(self.BCWindow, text="  Dirichlet BC  ").grid(row=0, column=4)
		def placeInputs():
			# Obs.: If 3D then neumann units are [Temp]/[Distance^2]
			self.boundaryUnitVars		  = []
			self.boundaryTypeVars		  = []
			self.boundaryValueEntries	  = []
			self.boundaryUnitMenus		  = []

			i=1
			for boundaryName in self.boundariesNames:
				options = []

				unitVar = tk.StringVar(self.BCWindow)
				unitVar.set("")
				bTypeVar = tk.IntVar(self.BCWindow)

				nameLabel = ttk.Label(self.BCWindow, text=boundaryName)
				nameLabel.grid(row=i, column=0, padx=5, sticky="W")

				valEntry = tk.Entry(self.BCWindow)
				valEntry.grid(row=i,column=1)

				unitMenu = ttk.OptionMenu(self.BCWindow, unitVar, *options)
				unitMenu.grid(row=i, column=2)

				neumannButton = tk.Radiobutton(self.BCWindow, variable=bTypeVar, value=1, command=chooseNeumannRadioButton(i-1))
				neumannButton.grid(row=i, column=3)

				dirichletButton = tk.Radiobutton(self.BCWindow, variable=bTypeVar, value=2, command=chooseDirichletRadioButton(i-1))
				dirichletButton.grid(row=i, column=4)

				self.boundaryUnitVars.append(unitVar)
				self.boundaryTypeVars.append(bTypeVar)
				self.boundaryValueEntries.append(valEntry)
				self.boundaryUnitMenus.append(unitMenu)

				i+=1
		
		placeStatic()
		placeInputs()

	def populatePage2(self): #
		self.page2 = tk.Canvas(self.root, height=500, width=600)

		self.bottomFrame = tk.Frame(self.page2)
		self.bottomFrame.place(relx=0.02, rely=0.87, relheight=0.1, relwidth=0.96, anchor="nw" )

		self.nextButton = tk.Button(self.bottomFrame, text="RUN", command=self.page2Next)
		self.nextButton.place(relx=1.0, rely=0.50, relheight=0.75, relwidth=0.20, anchor="e")

		self.prevButton = tk.Button(self.bottomFrame, text="Prev", command=self.page2Prev)
		self.prevButton.place(relx=0.78, rely=0.50, relheight=0.75, relwidth=0.20, anchor="e")

		self.propertyEntries = dict()
		self.propertyUnitVars = dict()
		self.propertiesFrames = []
		self.currentRegion = 0

		def nextRegion():
			if self.currentRegion < len(self.meshData.regionsNames)-1:
				self.propertiesFrames[self.currentRegion].place_forget()
				self.propertiesFrames[self.currentRegion+1].place(relx=0.02, rely=0.02, relheight=0.81, relwidth=0.96, anchor="nw")
				self.currentRegion += 1

		def prevRegion():
			if self.currentRegion > 0:
				self.propertiesFrames[self.currentRegion].place_forget()
				self.propertiesFrames[self.currentRegion-1].place(relx=0.02, rely=0.02, relheight=0.81, relwidth=0.96, anchor="nw")
				self.currentRegion -= 1

		for regionCount, region in enumerate(self.meshData.regionsNames):
			propertiesFrame = tk.LabelFrame(self.page2, text="Material Properties")
			self.propertiesFrames.append(propertiesFrame)
			centerFrame = tk.Frame(propertiesFrame)

			self.propertyEntries[region] = dict()
			self.propertyUnitVars[region] = dict()

			self.regionCountLabel = tk.Label(centerFrame, text=f"{region}\t[{regionCount+1}/{len(self.meshData.regionsNames)}]")
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

	def showPage2(self): #
		self.page2.pack(side="top", fill="both", expand="yes")

	def openFile(self):
		initialdir = os.path.join( os.path.dirname(__file__), os.path.pardir, "meshes" )
		self.meshFileName = tk.filedialog.askopenfilename(initialdir=initialdir, title="Select a mesh file", filetypes=(("MSH files", "*.msh"),("All files", "*")))

		if self.meshFileName: # Prevents against not choosing a file
			try:
				self.meshData = MSHReader(self.meshFileName).getData()
			except:
				messagebox.showwarning("Warning","Invalid mesh file")
				raise Exception("Invalid mesh file")
			if not ".msh" in self.meshFileName:
				messagebox.showwarning("Warning","Must be a .msh file")
				raise Exception("Must be a .msh file")
			

			self.fileLabel["text"] = self.meshFileName
			self.populateBCFrame()
			self.populatePage2()

	def page1Prev(self):
		pass
	
	def page1Next(self):
		self.checkPage1Data()
		self.page1.pack_forget()
		self.showPage2()
	
	def page2Prev(self):
		self.page2.pack_forget()
		self.showPage1()
	
	def page2Next(self):
		self.checkPage2Data()
		self.runSimulation()
	
	def checkPage1Data(self):
		if not self.meshFileName:
			messagebox.showwarning("Warning","Must Select a mesh File")
			# raise Exception("Must select a mesh file")
		for boundaryName, entry, bType in zip(self.boundariesNames, self.boundaryValueEntries, self.boundaryTypeVars):
			try:
				float( entry.get() )
			except:
				messagebox.showwarning("Warning","Invalid value in \"{}\" field".format(boundaryName))
				raise Exception("Invalid value in \"{}\" field".format(boundaryName))
			if bType.get() == 0:
				messagebox.showwarning("Warning","Must select either Neumann or Dirichlet Boundary Condition")
				raise Exception("Must select either Neumann or Dirichlet Boundary Condition")
	
	def checkPage2Data(self):
		for region in self.meshData.regionsNames:
			for propertyName in self.properties:
					propertyValue = self.propertyEntries[region][propertyName].get()
					try:
						float( propertyValue )
					except:
						messagebox.showwarning("Warning", "Invalid value in \"{}\" {} field".format(region, propertyName))
						raise Exception("Invalid value in \"{}\" {} field".format(region, propertyName))
	
	def runSimulation(self):
		grid = Grid( MSHReader(self.meshFileName).getData() )

		# Boundary Conditions
		self.boundaryConditionsData = dict()
		for bName, entry, unit, bType in zip(self.boundariesNames, self.boundaryValueEntries, self.boundaryUnitVars, self.boundaryTypeVars):
			self.boundaryConditionsData[bName] = {
				"condition": ["NEUMANN", "DIRICHLET"][bType.get()-1],
				"type": "CONSTANT",
				"value": float( entry.get() )
			}

		# Property Data
		self.propertyData = dict()
		for region in self.propertyEntries.keys():
			self.propertyData[region] = dict()
			for _property in self.propertyEntries[region].keys():
				self.propertyData[region][_property] = float( self.propertyEntries[region][_property].get() )


		heatTransfer(
			libraryPath = os.path.join(os.path.dirname(__file__), os.path.pardir),
			outputPath = os.path.join(os.path.dirname(__file__), os.path.pardir, "results", "gui"),
			extension = "csv",
			
			grid 	  = grid,
			propertyData = [self.propertyData["Body"]],#################
			
			initialValues = {"temperature": np.zeros(grid.vertices.size)},
			neumannBoundaries = {"temperature":[NeumannBoundaryCondition(grid, boundary, self.boundaryConditionsData[boundary.name]["value"], handle) for handle, boundary in enumerate(grid.boundaries) if self.boundaryConditionsData[boundary.name]["condition"] == "NEUMANN"]},
			dirichletBoundaries = {"temperature":[DirichletBoundaryCondition(grid, boundary, self.boundaryConditionsData[boundary.name]["value"], handle) for handle, boundary in enumerate(grid.boundaries) if self.boundaryConditionsData[boundary.name]["condition"] == "DIRICHLET"]},
 
			timeStep  = 10000.0,		#####################
			finalTime = 1e+06,		#####################
			maxNumberOfIterations = 1e+04,		#####################
			tolerance = 1e-06,		#####################
			
			transient = True,	#####################
			verbosity = True
		)
		print(grid)

if __name__ == "__main__":
	app = Application()
