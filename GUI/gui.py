import sys,os
sys.path.append(os.path.join(os.path.dirname(__file__), os.path.pardir))
from PyEFVLib import MSHReader, Grid, ProblemData, CgnsSaver, CsvSaver, NeumannBoundaryCondition, DirichletBoundaryCondition
from apps.heat_transfer import heatTransfer
from apps.stress_equilibrium import stressEquilibrium

import tkinter as tk
from tkinter import filedialog, ttk, messagebox

import numpy as np

class PyEFVLibGUI:
	def __init__(self):
		self.root = tk.Tk()
		self.root.title("PyEFVLib GUI")
		self.root.bind("<Key>", lambda key: self.root.destroy() if key.char=="\x17" else 0) # Close window if Ctrl+W is pressed

		self.HEIGHT, self.WIDTH = (500, 600)

		self.mainMenu = MainMenu(self.root, self)
		self.heatTransferApplication = HeatTransferApplication(self.root, self)
		self.solidMechanicsApplication = SolidMechanicsApplication(self.root, self)

		self.mainMenu.init()

		self.root.mainloop()

class MainMenu:
	def __init__(self, root, application):
		self.root = root
		self.app = application

	def init(self):
		self.populate()
		self.show()

	def populate(self):
		self.page = tk.Frame(self.root, width=600, height=500)

		self.heatTransferButton = tk.Button(self.page, text="Heat Transfer Application", font="Arial 24", command=self.openHeatTransfer)
		self.heatTransferButton.place(relx=0.5, rely=0.5, relwidth=0.75, relheight=0.25, anchor="s")

		self.solidMechanicsButton = tk.Button(self.page, text="Solid Mechanics Application", font="Arial 24", command=self.openSolidMechanics)
		self.solidMechanicsButton.place(relx=0.5, rely=0.5, relwidth=0.75, relheight=0.25, anchor="n")

	def show(self):
		self.page.pack(side="top", fill="both", expand="yes")

	def openHeatTransfer(self):
		self.page.destroy()
		self.app.heatTransferApplication.init()

	def openSolidMechanics(self):
		self.page.destroy()
		self.app.solidMechanicsApplication.init()


class Application:
	def __init__(self, root, application):
		self.root = root
		self.app = application
		self.settings()

	def init(self):
		self.populatePage1()
		self.showPage1()

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

		prevButton = tk.Button(bottomFrame, text="Prev", command=self.page1Prev)
		prevButton.place(relx=0.78, rely=0.50, relheight=0.75, relwidth=0.20, anchor="e")

	def showPage1(self):
		self.page1.pack(side="top", fill="both", expand="yes")

	def populateBCFrame(self):
		self.boundariesNames = self.meshData.boundariesNames
		self.currentField = 0
		fieldsBCCanvas = []
		fieldsIValueFrames = []
		fieldsBCWindows = dict()

		self.boundaryUnitVars		= dict()
		self.boundaryTypeVars		= dict()
		self.boundaryValueEntries	= dict()
		self.boundaryUnitMenus		= dict()

		self.initialValueEntries	= dict()
		self.initialValueUnitVars	= dict()

		def chooseNeumannRadioButton(i, field): #
			def chooseNeumannRadioButtonFunc():
				unitVar = self.boundaryUnitVars[field][i]
				optionMenu = self.boundaryUnitMenus[field][i]

				unitVar.set("")
				optionMenu["menu"].delete(0,"end") # remove full list 
				for opt in self.neumannUnits[field]: 
					optionMenu['menu'].add_command(label=opt, command=tk._setit(unitVar, opt))
				unitVar.set(self.neumannUnits[field][0]) # default value set 

			return chooseNeumannRadioButtonFunc
		def chooseDirichletRadioButton(i, field): #
			def chooseDirichletRadioButtonFunc():
				unitVar = self.boundaryUnitVars[field][i]
				optionMenu = self.boundaryUnitMenus[field][i]

				unitVar.set("")
				optionMenu["menu"].delete(0,"end") # remove full list 
				for opt in self.dirichletUnits[field]: 
					optionMenu['menu'].add_command(label=opt, command=tk._setit(unitVar, opt))
				unitVar.set(self.dirichletUnits[field][0]) # default value set 
			return chooseDirichletRadioButtonFunc
		def prevField():
			if self.currentField > 0:
				fieldsBCCanvas[self.currentField].place_forget()
				fieldsIValueFrames[self.currentField].place_forget()
				self.currentField -= 1
				fieldsBCCanvas[self.currentField].place(relx=0.0, rely=0.0, relwidth=1.0, relheight=0.85, anchor="nw")
				fieldsIValueFrames[self.currentField].place(relx=0.0, rely=0.85, relwidth=1.0, relheight=0.15, anchor="nw")
				self.BCFrame.configure(text="Boundary Conditions Settings - {} [{}/{}]".format( self.fields[self.currentField] , self.currentField+1, len(self.fields)))
		def nextField():
			if self.currentField < len(self.fields)-1:
				fieldsBCCanvas[self.currentField].place_forget()
				fieldsIValueFrames[self.currentField].place_forget()
				self.currentField += 1
				fieldsBCCanvas[self.currentField].place(relx=0.0, rely=0.0, relwidth=1.0, relheight=0.85, anchor="nw")
				fieldsIValueFrames[self.currentField].place(relx=0.0, rely=0.85, relwidth=1.0, relheight=0.15, anchor="nw")
				self.BCFrame.configure(text="Boundary Conditions Settings - {} [{}/{}]".format( self.fields[self.currentField] , self.currentField+1, len(self.fields)))
		def placeStatic(field):
			BCCanvas = tk.Canvas(self.BCFrame)
			fieldsBCCanvas.append(BCCanvas)
			# BCCanvas.place(relx=0.0, rely=0.0, relwidth=1.0, relheight=0.85, anchor="nw")

			scrollbar = ttk.Scrollbar(self.BCFrame, orient="vertical", command=BCCanvas.yview)
			scrollbar.place(relx=1.0, rely=0.0, relheight=0.85, anchor="ne")

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
			
			fieldsBCWindows[field] = BCWindow
		def placeInputs(field):
			# Obs.: If 3D then neumann units are [Temp]/[Distance^2]
			BCWindow = fieldsBCWindows[field]
			self.boundaryUnitVars[field]		  = []
			self.boundaryTypeVars[field]		  = []
			self.boundaryValueEntries[field]	  = []
			self.boundaryUnitMenus[field]		  = []

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
				unitMenu.grid(row=i, column=2, sticky="E")

				neumannButton = tk.Radiobutton(BCWindow, variable=bTypeVar, value=1, command=chooseNeumannRadioButton(i-1, field))
				neumannButton.grid(row=i, column=3)

				dirichletButton = tk.Radiobutton(BCWindow, variable=bTypeVar, value=2, command=chooseDirichletRadioButton(i-1, field))
				dirichletButton.grid(row=i, column=4)

				self.boundaryUnitVars[field].append(unitVar)
				self.boundaryTypeVars[field].append(bTypeVar)
				self.boundaryValueEntries[field].append(valEntry)
				self.boundaryUnitMenus[field].append(unitMenu)

				i+=1
		def placeInitialValue(field):
			initialValueFrame = tk.LabelFrame(self.BCFrame)
			fieldsIValueFrames.append(initialValueFrame)
			# initialValueFrame.place(relx=0.0, rely=0.85, relwidth=1.0, relheight=0.15, anchor="nw")

			initialValueLabel = tk.Label(initialValueFrame, text="Frame's\ninitial value")
			initialValueLabel.place(relx=0.0, rely=0.5, anchor="w")

			self.initialValueEntries[field] = tk.Entry(initialValueFrame)
			self.initialValueEntries[field].place(x=100, rely=0.5, anchor="w")

			self.initialValueUnitVars[field] = tk.StringVar(initialValueFrame)
			self.initialValueUnitVars[field].set("")

			unitMenu = ttk.OptionMenu(initialValueFrame, self.initialValueUnitVars[field], *self.dirichletUnits[field])
			unitMenu.place(x=270, rely=0.5, anchor="w")

			prevButton = tk.Button(initialValueFrame, text="<", command=prevField, state="disabled" if field == self.fields[0] else "normal")
			prevButton.place(relx=0.9, rely=0.5, anchor="e")

			nextButton = tk.Button(initialValueFrame, text=">", command=nextField, state="disabled" if field == self.fields[-1] else "normal")
			nextButton.place(relx=0.9, rely=0.5, anchor="w")

		for field in self.fields:
			placeStatic(field)
			placeInputs(field)
			placeInitialValue(field)

		fieldsBCCanvas[0].place(relx=0.0, rely=0.0, relwidth=1.0, relheight=0.85, anchor="nw")
		fieldsIValueFrames[0].place(relx=0.0, rely=0.85, relwidth=1.0, relheight=0.15, anchor="nw")
		self.BCFrame.configure(text="Boundary Conditions Settings - {} [{}/{}]".format( self.fields[self.currentField] , self.currentField+1, len(self.fields)))

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
		self.page1.destroy()
		self.app.mainMenu.init()
	
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
		for field in self.fields:
			for boundaryName, entry, bType in zip(self.boundariesNames, self.boundaryValueEntries[field], self.boundaryTypeVars[field]):
				try:
					float( entry.get() )
				except:
					messagebox.showwarning("Warning","Invalid value in \"{} - {}\" field".format(field, boundaryName))
					raise Exception("Invalid value in \"{} - {}\" field".format(field, boundaryName))
				if bType.get() == 0:
					messagebox.showwarning("Warning","Must select either Neumann or Dirichlet Boundary Condition")
					raise Exception("Must select either Neumann or Dirichlet Boundary Condition")
			try:
				float( self.initialValueEntries[field].get() )
			except:
				messagebox.showwarning("Warning", "Invalid Value in Initial Value field")
				raise Exception("Invalid Value in Initial Value field")
	
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


class SolidMechanicsApplication(Application):
	def settings(self):
		self.fields = ["u", "v"]

		self.neumannUnits = {"u": ["Pa", "kPa", "MPa", "GPa", "kgf/m²", "psi"], "v": ["mm/mm"]}
		self.dirichletUnits = {"u": ["m", "mm", "cm", "μm", "inch"], "v": ["m", "mm", "cm", "μm", "inch"]}

		self.properties = ["Density", "PoissonsRatio", "ShearModulus", "Gravity"]
		self.propertyUnits = {
			"Density"		: ["kg/m³", "g/cm³"], 		
			"PoissonsRatio" : [""], 
			"ShearModulus"	: ["MPa", "GPa", "kPa", "Pa", "psi"], 
			"Gravity"		: ["m/s²", "mm/s²", "km/h²"]

		}

		self.numericalSettingsOp = {	# "option": [entry, unit, checkbox]
			"Time Step"					:{"option":[1,1,0,0],"units":["s", "min", "h", "days", "weeks"]},
			"Final time"				:{"option":[1,1,0,0],"units":["s", "min", "h", "days", "weeks"]},
			"Max number of iterations"	:{"option":[1,0,0,0],"units":[""]},
			"Tolerance"					:{"option":[1,1,0,0],"units":["K"]},
			"Transient"					:{"option":[0,0,0,0],"units":[""]}
		}

		self.unitsConvert = {
			"m"		: lambda x: x,
			"mm"		: lambda x: 0.001*x,
			"cm"		: lambda x: 0.01*x,
			"μm"		: lambda x: 1e-06*x,
			"inch"		: lambda x: 0.0254*x,
			"kg/m³"		: lambda rho: rho,
			"g/cm³"		: lambda rho: 1000.0*rho,
			""			: lambda _: _,
			"MPa"		: lambda P: 1e+06*P,
			"kgf/m²"	: lambda P: 9.807*P,
			"GPa"		: lambda P: 1e+09*P,
			"kPa"		: lambda P: 1e+03*P,
			"Pa"		: lambda P: P,
			"psi"		: lambda P: 6894.757*P,
			"m/s²"		: lambda a: a,
			"mm/s²"		: lambda a: 0.001*a,
			"km/h²"		: lambda a: 7.716e-05*a,
			"s"			: lambda t: t,
			"min"		: lambda t: 60.0*t,
			"h"			: lambda t: 3600*t,
			"days"		: lambda t: 86400.0*t,
			"weeks"		: lambda t: 604800.0*t,
		}

		self.meshFileName = ""

	def runSimulation(self):
		grid = Grid( self.meshData )

		# Boundary Conditions
		initialValues = { field : self.unitsConvert[ self.initialValueUnitVars[field].get() ](float( self.initialValueEntries[field].get() )) for field in self.fields}
		
		neumannBoundaries = dict()
		dirichletBoundaries = dict()
		boundaryConditionsDict = {bName : dict() for bName in self.boundariesNames}

		for field in self.fields:
			neumannBoundaries[field] = []
			dirichletBoundaries[field] = []
			handle = 0

			for bName, boundary, entry, unit, bType in zip(self.boundariesNames, grid.boundaries, self.boundaryValueEntries[field], self.boundaryUnitVars[field], self.boundaryTypeVars[field]):
				value = self.unitsConvert[unit.get()]( float( entry.get() ) )
				if bType.get() == 1:
					bc = NeumannBoundaryCondition(grid, boundary, value, handle)
					neumannBoundaries[field].append(bc)
				if bType.get() == 2:
					bc = DirichletBoundaryCondition(grid, boundary, value, handle)
					dirichletBoundaries[field].append(bc)
				boundaryConditionsDict[bName][field] = bc
				handle += 1

		# Property Data
		self.propertyData = []
		for region in grid.regions:
			self.propertyData.append( dict() )
			for propertyName in self.propertyEntries[region.name].keys():
				self.propertyData[-1][propertyName] = float( self.propertyEntries[region.name][propertyName].get() )
				self.propertyData[-1][propertyName] = self.unitsConvert[ self.propertyUnitVars[region.name][propertyName].get() ]( self.propertyData[-1][propertyName] )

		stressEquilibrium(
			libraryPath = os.path.join(os.path.dirname(__file__), os.path.pardir),
			outputPath = os.path.join(os.path.dirname(__file__), os.path.pardir, "results", "gui"),
			extension = "csv",
			
			grid 	  = grid,
			propertyData = self.propertyData,

			initialValues = initialValues,
			neumannBoundaries = neumannBoundaries,
			dirichletBoundaries = dirichletBoundaries,
			boundaryConditions = list(boundaryConditionsDict.values()),

			verbosity=True 
		)
		self.root.destroy()

class HeatTransferApplication(Application):
	def settings(self):
		self.fields = ["temperature"]

		self.neumannUnits = {"temperature": ["K/m"]}
		self.dirichletUnits = {"temperature": ["K", "°C", "°F"]}
		
		self.properties = ["Density","HeatCapacity","Conductivity","HeatGeneration"]
		self.propertyUnits = {
			"Density":        ["kg/m³", "g/cm³"],
			"HeatCapacity":   ["J/kg.K"],
			"Conductivity":   ["W/m.K"],
			"HeatGeneration": ["K/m³"]
		}

		self.numericalSettingsOp = {	# "option": [entry, unit, checkbox]
			"Time Step"					:{"option":[1,1,0,1],"units":["s", "min", "h", "days", "weeks"]},
			"Final time"				:{"option":[1,1,1,0],"units":["s", "min", "h", "days", "weeks"]},
			"Max number of iterations"	:{"option":[1,0,1,0],"units":[""]},
			"Tolerance"					:{"option":[1,1,1,1],"units":["K"]},
			"Transient"					:{"option":[0,0,1,1],"units":[""]}
		}

		self.unitsConvert = {
			"K"		: lambda T: T,
			"°C"	: lambda T: T+273.15,
			"°F"	: lambda T: 5.0*(T-32.0)/9.0+273.15,
			"K/m"	: lambda d: d,
			"kg/m³"	: lambda rho: rho,
			"g/cm³"	: lambda rho: 1000.0*rho,
			"J/kg.K": lambda Cp: Cp,
			"W/m.K"	: lambda k: k,
			"K/m³"	: lambda q: q,
			"s"		: lambda t: t,
			"min"	: lambda t: 60.0*t,
			"h"		: lambda t: 3600*t,
			"days"	: lambda t: 86400.0*t,
			"weeks"	: lambda t: 604800.0*t,
		}

		self.meshFileName = ""

	def runSimulation(self):
		grid = Grid( self.meshData )

		# Boundary Conditions
		boundaryConditionsData = dict()
		for bName, entry, unit, bType in zip(self.boundariesNames, self.boundaryValueEntries["temperature"], self.boundaryUnitVars["temperature"], self.boundaryTypeVars["temperature"]):
			boundaryConditionsData[bName] = {
				"condition": ["NEUMANN", "DIRICHLET"][bType.get()-1],
				"type": "CONSTANT",
				"value": float( entry.get() )
			}

		initialValues = {"temperature": [float(self.initialValueEntries["temperature"].get())] * grid.vertices.size},
		neumannBoundaries = {"temperature":[NeumannBoundaryCondition(grid, boundary, boundaryConditionsData[boundary.name]["value"], handle) for handle, boundary in enumerate(grid.boundaries) if boundaryConditionsData[boundary.name]["condition"] == "NEUMANN"]},
		dirichletBoundaries = {"temperature":[DirichletBoundaryCondition(grid, boundary, boundaryConditionsData[boundary.name]["value"], handle) for handle, boundary in enumerate(grid.boundaries) if boundaryConditionsData[boundary.name]["condition"] == "DIRICHLET"]}

		# Property Data
		self.propertyData = []
		for region in grid.regions:
			self.propertyData.append( dict() )
			for propertyName in self.propertyEntries[region.name].keys():
				self.propertyData[-1][propertyName] = float( self.propertyEntries[region.name][propertyName].get() )
				self.propertyData[-1][propertyName] = self.unitsConvert[ self.propertyUnitVars[region.name][propertyName].get() ]( self.propertyData[-1][propertyName] )

		# Numerical Settings
		index = list(self.numericalSettingsOp.keys()).index("Time Step")
		if self.numericalSettingsBools[index].get():
			timeStep  = float(self.numericalSettingsEntries[index].get())
			timeStep  = self.unitsConvert[ self.numericalSettingsUnits[index].get() ](timeStep) 
		else:
			timeStep = 0.0

		index = list(self.numericalSettingsOp.keys()).index("Final time")
		if self.numericalSettingsBools[index].get():
			finalTime = float(self.numericalSettingsEntries[index].get())
			finalTime  = self.unitsConvert[ self.numericalSettingsUnits[index].get() ](finalTime) 
		else:
			finalTime = None

		index = list(self.numericalSettingsOp.keys()).index("Max number of iterations")
		if self.numericalSettingsBools[index].get():
			maxNumberOfIterations = float(self.numericalSettingsEntries[index].get())
		else:
			maxNumberOfIterations = None

		index = list(self.numericalSettingsOp.keys()).index("Tolerance")
		if self.numericalSettingsBools[index].get():
			tolerance = float(self.numericalSettingsEntries[index].get())
		else:
			tolerance = None

		transient = self.numericalSettingsBools[list(self.numericalSettingsOp.keys()).index("Transient")].get()

		heatTransfer(
			libraryPath = os.path.join(os.path.dirname(__file__), os.path.pardir),
			outputPath = os.path.join(os.path.dirname(__file__), os.path.pardir, "results", "gui"),
			extension = "csv",
			
			grid 	  = grid,
			propertyData = self.propertyData,
			
			initialValues = {"temperature": [ self.unitsConvert[self.initialValueUnitVars["temperature"].get()]( float(self.initialValueEntries["temperature"].get()) ) ] * grid.vertices.size},
			neumannBoundaries = {"temperature":[NeumannBoundaryCondition(grid, boundary, boundaryConditionsData[boundary.name]["value"], handle) for handle, boundary in enumerate(grid.boundaries) if boundaryConditionsData[boundary.name]["condition"] == "NEUMANN"]},
			dirichletBoundaries = {"temperature":[DirichletBoundaryCondition(grid, boundary, boundaryConditionsData[boundary.name]["value"], handle) for handle, boundary in enumerate(grid.boundaries) if boundaryConditionsData[boundary.name]["condition"] == "DIRICHLET"]},
 
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
	app = PyEFVLibGUI()
