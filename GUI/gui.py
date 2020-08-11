import sys,os
sys.path.append(os.path.join(os.path.dirname(__file__), os.path.pardir))
from PyEFVLib import MSHReader, Grid, ProblemData, CgnsSaver, CsvSaver

import tkinter as tk
from tkinter import filedialog, ttk

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
		self.bcData = dict()
		# Boundary Conditions
		for boundaryName, entry, unit, conditionType in zip( self.bcPage.boundariesNames, self.bcPage.boundaryValueEntries, self.bcPage.boundaryUnitVars,self.bcPage.boundaryTypeVars ):
			self.bcData[boundaryName] = {
				"condition": ["NEUMANN", "DIRICHLET"][conditionType.get()-1],
				"type": "CONSTANT",
				"value": self.unitConversions[unit.get()]( float( entry.get() ) )
			}

		# Properties
		self.propertiesData = dict()
		for propertyName, entry, unit in zip( self.propertiesPage.properties, self.propertiesPage.propertyEntries, self.propertiesPage.propertyUnitVars ):
			self.propertiesData[propertyName] = self.unitConversions[unit.get()]( float( entry.get() ) )

	def dumpData(self):
		pass

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
		BCCanvas.pack(side="left", fill="both", expand=True)

		scrollbar = ttk.Scrollbar(self.BCFrame, orient="vertical", command=BCCanvas.yview)
		scrollbar.pack(side="right", fill="y")

		# with Windows OS
		self.root.bind("<MouseWheel>", lambda e: BCCanvas.yview_scroll(-1*(event.delta/120), "units"))
		# with Linux OS
		self.root.bind("<Button-4>", lambda e: BCCanvas.yview_scroll(-1, "units"))
		self.root.bind("<Button-5>", lambda e: BCCanvas.yview_scroll(+1, "units"))

		BCWindow = ttk.Frame(BCCanvas)
		BCWindow.bind("<Configure>", lambda event: BCCanvas.configure(scrollregion=BCCanvas.bbox("all")))

		BCCanvas.create_window((0, 0), window=BCWindow, anchor="nw")
		BCCanvas.configure(yscrollcommand=scrollbar.set)

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
		self.boundaryNeumannButtons	  = []
		self.boundaryDirichletButtons = []

		i=1
		for boundaryName in self.boundariesNames:
			options = []

			unitVar = tk.StringVar(BCWindow)
			unitVar.set("")
			bTypeVar = tk.IntVar(BCWindow)

			nameLabel = ttk.Label(BCWindow, text=boundaryName)
			nameLabel.grid(row=i, column=0)

			valEntry = tk.Entry(BCWindow)
			valEntry.grid(row=i,column=1)
			print(valEntry.get())

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
			self.boundaryNeumannButtons.append(neumannButton)
			self.boundaryDirichletButtons.append(dirichletButton)
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
		self.fileLabel["text"] = tk.filedialog.askopenfilename(initialdir="../meshes", title="Select a mesh file", filetypes=(("MSH files", "*.msh"),("All files", "*")))
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

		# For now only one region
		self.canvas = tk.Canvas(self.root, height=self.HEIGHT, width=self.WIDTH)
		self.canvas.pack(side="top", fill="both", expand="yes")

		self.bottomFrame = tk.Frame(self.canvas)
		self.bottomFrame.place(relx=0.02, rely=0.87, relheight=0.1, relwidth=0.96, anchor="nw" )

		self.nextButton = tk.Button(self.bottomFrame, text="RUN", command=self.next)
		self.nextButton.place(relx=1.0, rely=0.50, relheight=0.75, relwidth=0.20, anchor="e")

		self.prevButton = tk.Button(self.bottomFrame, text="Prev", command=self.prev)
		self.prevButton.place(relx=0.78, rely=0.50, relheight=0.75, relwidth=0.20, anchor="e")

		self.propertiesFrame = tk.LabelFrame(self.canvas, text="Material Properties")
		self.propertiesFrame.place(relx=0.02, rely=0.02, relheight=0.81, relwidth=0.96, anchor="nw")

		self.properties = ["Density","HeatCapacity","Conductivity","HeatGeneration"]
		self.propertyUnits = {
			"Density":       ["kg/m³", "g/cm³"],
			"HeatCapacity":   ["J/kg.K"],
			"Conductivity":   ["W/m.K"],
			"HeatGeneration": ["K/m³"]
		}
		self.propertyEntries = []
		self.propertyUnitVars = []

		i=0
		for propertyName in self.properties:
			options = self.propertyUnits[propertyName]

			unitVar = tk.StringVar(self.propertiesFrame)
			unitVar.set(self.propertyUnits[propertyName][0])
			bTypeVar = tk.IntVar(self.propertiesFrame)

			nameLabel = ttk.Label(self.propertiesFrame, text=propertyName)
			nameLabel.grid(row=i, column=0)

			valEntry = ttk.Entry(self.propertiesFrame)
			valEntry.grid(row=i,column=1)

			unitMenu = tk.OptionMenu(self.propertiesFrame, unitVar, *options)
			unitMenu.grid(row=i, column=2)

			self.propertyEntries.append(valEntry)
			self.propertyUnitVars.append(unitVar)
			i+=1

		self.populated = True

	def next(self):
		print("Running Simulation...")
		self.app.getData()

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
	app.root.mainloop()