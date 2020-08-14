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

		self.numericalSettingsOp = {	# "option": [entry, unit, checkbox]
			"Time Step"					:{"option":[1,1,0,1],"units":["s", "h", "days", "weeks"]},
			"Final time"				:{"option":[1,1,1,0],"units":["s", "h", "days", "weeks"]},
			"Max number of iterations"	:{"option":[1,0,1,0],"units":[""]},
			"Tolerance"					:{"option":[1,1,1,1],"units":["K"]},
			"Transient"					:{"option":[0,0,1,1],"units":[""]}
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
		bottomFrame = tk.Frame(self.page1)
		bottomFrame.place(relx=0.02, rely=0.87, relheight=0.1, relwidth=0.96, anchor="nw" )

		nextButton = tk.Button(bottomFrame, text="Next", command=self.page1Next)
		nextButton.place(relx=1.0, rely=0.50, relheight=0.75, relwidth=0.20, anchor="e")

		prevButton = tk.Button(bottomFrame, text="Prev", command=self.page1Prev, state="disabled")
		prevButton.place(relx=0.78, rely=0.50, relheight=0.75, relwidth=0.20, anchor="e")

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

		# Bottom Frame
		bottomFrame = tk.Frame(self.page2)
		bottomFrame.place(relx=0.02, rely=0.87, relheight=0.1, relwidth=0.96, anchor="nw" )

		nextButton = tk.Button(bottomFrame, text="RUN", command=self.page2Next)
		nextButton.place(relx=1.0, rely=0.50, relheight=0.75, relwidth=0.20, anchor="e")

		prevButton = tk.Button(bottomFrame, text="Prev", command=self.page2Prev)
		prevButton.place(relx=0.78, rely=0.50, relheight=0.75, relwidth=0.20, anchor="e")

		# Property Frame
		self.propertyEntries = dict()
		self.propertyUnitVars = dict()
		self.propertiesFrames = []
		self.currentRegion = 0

		numberOfRegions = len(self.meshData.regionsNames)

		def nextRegion():
			if self.currentRegion < numberOfRegions-1:
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

			self.regionCountLabel = tk.Label(centerFrame, text=f"{region}\t[{regionCount+1}/{numberOfRegions}]")
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

			prevButton = tk.Button(centerFrame, text="<", command=prevRegion, state="disabled" if numberOfRegions==1 or regionCount==0 else "normal")
			prevButton.grid(row=i, column=2, sticky="E")

			nextButton = tk.Button(centerFrame, text=">", command=nextRegion, state="disabled" if numberOfRegions==1 or regionCount==numberOfRegions-1 else "normal")
			nextButton.grid(row=i, column=3)

			centerFrame.place(relx=0.5, rely=0.0, anchor="n")

		self.propertiesFrames[0].place(relx=0.02, rely=0.02, relheight=0.45, relwidth=0.96, anchor="nw")
		
		# Numerical Frame
		numericalFrame = tk.LabelFrame(self.page2, text="Numerical Settings")

		self.numericalSettingsEntries = []
		self.numericalSettingsUnits = []
		self.numericalSettingsUnitMenus = []
		self.numericalSettingsBools = []
		self.numericalSettingsCheckboxes = []

		def toggleCheckbox(i):
			def toggleCheckboxFunc():
				setting = list(self.numericalSettingsOp.keys())[i]
				state = self.numericalSettingsBools[i].get()

				if setting != "Transient":
					# Entry
					if self.numericalSettingsOp[setting]["option"][0]:
						if state:
							self.numericalSettingsEntries[i].grid(row=i, column=1, padx=5, pady=5, sticky="W")
						else:
							self.numericalSettingsEntries[i].grid_forget()
					# Unit
					if self.numericalSettingsOp[setting]["option"][1]:
						if state:
							self.numericalSettingsUnitMenus[i].grid(row=i, column=2, sticky="E")
						else:
							self.numericalSettingsUnitMenus[i].grid_forget()
				else:
					if state:
						j=0
						for loopSetting, entry, unitMenu, boolVar, checkbox in zip( self.numericalSettingsOp.keys(), self.numericalSettingsEntries, self.numericalSettingsUnitMenus, self.numericalSettingsBools, self.numericalSettingsCheckboxes ):
							if self.numericalSettingsOp[loopSetting]["option"][0] and self.numericalSettingsOp[loopSetting]["option"][3]:
								entry.grid(row=j, column=1, padx=5, pady=5, sticky="W")
							
							if self.numericalSettingsOp[loopSetting]["option"][1] and self.numericalSettingsOp[loopSetting]["option"][3]:
								unitMenu.grid(row=j, column=2, sticky="E")
							boolVar.set(self.numericalSettingsOp[loopSetting]["option"][3])
							if loopSetting != "Time Step":
								checkbox.configure(state="normal")
							j+=1
					else:
						for loopSetting, entry, unitMenu, boolVar, checkbox in zip( self.numericalSettingsOp.keys(), self.numericalSettingsEntries, self.numericalSettingsUnitMenus, self.numericalSettingsBools, self.numericalSettingsCheckboxes ):
							if loopSetting != "Transient":
								entry.grid_forget()
								unitMenu.grid_forget()
								boolVar.set(False)
								checkbox.configure(state="disabled")

			return toggleCheckboxFunc

		i=0
		for setting in self.numericalSettingsOp.keys():
			# Label
			label = tk.Label(numericalFrame, text=setting)
			label.grid(row=i, column=0, padx=5, pady=5, sticky="W")

			# Entry
			entry = tk.Entry(numericalFrame)
			if self.numericalSettingsOp[setting]["option"][0] and self.numericalSettingsOp[setting]["option"][3]:
				entry.grid(row=i, column=1, padx=5, pady=5, sticky="W")

			self.numericalSettingsEntries.append(entry)

			# Unit
			unitVar = tk.StringVar(numericalFrame)
			unitVar.set(self.numericalSettingsOp[setting]["units"][0])

			unitMenu = tk.OptionMenu(numericalFrame, unitVar, *self.numericalSettingsOp[setting]["units"])
			if self.numericalSettingsOp[setting]["option"][1] and self.numericalSettingsOp[setting]["option"][3]:
				unitMenu.grid(row=i, column=2, sticky="E")

			self.numericalSettingsUnits.append(unitVar)
			self.numericalSettingsUnitMenus.append(unitMenu)

			# Checkbox
			boolVar = tk.BooleanVar()
			boolVar.set(self.numericalSettingsOp[setting]["option"][3])

			checkbox = tk.Checkbutton(numericalFrame, variable=boolVar, command=toggleCheckbox(i), state= "normal" if self.numericalSettingsOp[setting]["option"][2] else "disabled")
			checkbox.grid(row=i, column=3)

			self.numericalSettingsBools.append(boolVar)
			self.numericalSettingsCheckboxes.append(checkbox)

			i+=1

		numericalFrame.place(relx=0.02, rely=0.47, relheight=0.38, relwidth=0.96, anchor="nw")


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
	
		transientIndex = list(self.numericalSettingsOp.keys()).index("Transient")
		finalTimeIndex = list(self.numericalSettingsOp.keys()).index("Final time")
		maxNOfItIndex = list(self.numericalSettingsOp.keys()).index("Max number of iterations")
		toleranceIndex = list(self.numericalSettingsOp.keys()).index("Tolerance")
		
		if self.numericalSettingsBools[transientIndex].get():
			if not self.numericalSettingsBools[finalTimeIndex].get() and not self.numericalSettingsBools[maxNOfItIndex].get() and not self.numericalSettingsBools[toleranceIndex].get():
				messagebox.showwarning("Warning", "There is no way of reaching convergence. Set at least one of final time, max number of iterations or tolerance fields")
				raise Exception("There is no way of reaching convergence. Set at least one of final time, max number of iterations or tolerance fields")

		i=0
		for numericalSetting in self.numericalSettingsOp.keys():
			if self.numericalSettingsOp[numericalSetting]["option"][0] and self.numericalSettingsBools[i].get():
				try:
					float( self.numericalSettingsEntries[i].get() )
				except:
					messagebox.showwarning("Warning", "Invalid input in {} field".format(numericalSetting))
					raise Exception("Invalid input in {} field".format(numericalSetting))
			i+=1

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
		# TEM QUE TOMAR CUIDADO, PQ PROPERTY DATA DEVE ESTAR NA ORDEM DAS REGIÕES NA GRID. MAS ATUALMENTE ESTÁ NA ORDEM DAS REGIÕES
		# DO ARQUIVO MSH. TORCER PRA SER SEMPRE O MESMO, MAS MUDA SAPORRAÍ
		self.propertyData = []
		for region in self.propertyEntries.keys():
			self.propertyData.append( dict() )
			for _property in self.propertyEntries[region].keys():
				self.propertyData[-1][_property] = float( self.propertyEntries[region][_property].get() )

		# Numerical Settings
			if self.numericalSettingsBools[list(self.numericalSettingsOp.keys()).index("Time Step")].get():
				timeStep  = float(self.numericalSettingsEntries[list(self.numericalSettingsOp.keys()).index("Time Step")].get())
			else:
				timeStep = 0.0

			if self.numericalSettingsBools[list(self.numericalSettingsOp.keys()).index("Final time")].get():
				finalTime = float(self.numericalSettingsEntries[list(self.numericalSettingsOp.keys()).index("Final time")].get())
			else:
				finalTime = None

			if self.numericalSettingsBools[list(self.numericalSettingsOp.keys()).index("Max number of iterations")].get():
				maxNumberOfIterations = float(self.numericalSettingsEntries[list(self.numericalSettingsOp.keys()).index("Max number of iterations")].get())
			else:
				maxNumberOfIterations = None

			if self.numericalSettingsBools[list(self.numericalSettingsOp.keys()).index("Tolerance")].get():
				tolerance = float(self.numericalSettingsEntries[list(self.numericalSettingsOp.keys()).index("Tolerance")].get())
			else:
				tolerance = None

			
			transient = self.numericalSettingsBools[list(self.numericalSettingsOp.keys()).index("Transient")].get()



		heatTransfer(
			libraryPath = os.path.join(os.path.dirname(__file__), os.path.pardir),
			outputPath = os.path.join(os.path.dirname(__file__), os.path.pardir, "results", "gui"),
			extension = "csv",
			
			grid 	  = grid,
			propertyData = self.propertyData,
			
			initialValues = {"temperature": np.zeros(grid.vertices.size)},
			neumannBoundaries = {"temperature":[NeumannBoundaryCondition(grid, boundary, self.boundaryConditionsData[boundary.name]["value"], handle) for handle, boundary in enumerate(grid.boundaries) if self.boundaryConditionsData[boundary.name]["condition"] == "NEUMANN"]},
			dirichletBoundaries = {"temperature":[DirichletBoundaryCondition(grid, boundary, self.boundaryConditionsData[boundary.name]["value"], handle) for handle, boundary in enumerate(grid.boundaries) if self.boundaryConditionsData[boundary.name]["condition"] == "DIRICHLET"]},
 
			timeStep  = timeStep ,
			finalTime = finalTime,
			maxNumberOfIterations = maxNumberOfIterations,
			tolerance = tolerance,
			
			transient = transient,
			verbosity = True
		)
		self.root.destroy()

class HeatTransferApplication2(HeatTransferApplication):
	def __init__(self, root):
		self.root = root

		self.meshFileName = "/home/gustavoe/Documents/Sinmec/HTRelated/PyEFVLib/meshes/Square.msh"
		# self.meshFileName = "/home/gustavoe/Documents/Sinmec/HTRelated/PyEFVLib/meshes/dualRegion.msh"
		self.meshData = MSHReader(self.meshFileName).getData()

		self.settings()

		self.populatePage2()
		self.showPage2()

if __name__ == "__main__":
	app = Application()
