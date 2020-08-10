import sys,os
sys.path.append(os.path.join(os.path.dirname(__file__), os.path.pardir))
from PyEFVLib import MSHReader, Grid, ProblemData, CgnsSaver, CsvSaver


import tkinter as tk
from tkinter import filedialog, ttk

class Application():
	def __init__(self):
		WIDTH, HEIGHT = 600,500

		self.root = tk.Tk()
		self.root.title("PyEFVLib GUI")
		self.root.bind("<Key>", lambda key: self.root.destroy() if key.char=="\x17" else 0) # Close window if Ctrl+W is pressed

		self.canvas = tk.Canvas(self.root, height=HEIGHT, width=WIDTH)
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


	def placeBCForms(self):
		self.boundariesNames = gridData = MSHReader(self.fileLabel["text"]).getData().boundariesNames
				
		BCCanvas = tk.Canvas(self.BCFrame)
		BCCanvas.pack(side="left", fill="both", expand=True)

		scrollbar = ttk.Scrollbar(self.BCFrame, orient="vertical", command=BCCanvas.yview)
		scrollbar.pack(side="right", fill="y")

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

			valEntry = ttk.Entry(BCWindow)
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
		self.openFrame.destroy()
		self.BCFrame.destroy()

		# For now only one region

		self.propertiesFrame = tk.LabelFrame(self.canvas, text="Material Properties")
		self.propertiesFrame.place(relx=0.02, rely=0.02, relheight=0.81, relwidth=0.96, anchor="nw")

		properties = ["Density","HeatCapacity","Conductivity","HeatGeneration"]
		propertyUnits = {"Density":       ["kg/m³", "g/cm³"],
						"HeatCapacity":   ["J/kg.K"],
						"Conductivity":   ["W/m.K"],
						"HeatGeneration": ["K/m³"]
		}
		propertyEntries = []
		propertyUnitMenus = []

		i=0
		for propertyName in properties:
			options = propertyUnits[propertyName]

			unitVar = tk.StringVar(self.propertiesFrame)
			unitVar.set("")
			bTypeVar = tk.IntVar(self.propertiesFrame)

			nameLabel = ttk.Label(self.propertiesFrame, text=propertyName)
			nameLabel.grid(row=i, column=0)

			valEntry = ttk.Entry(self.propertiesFrame)
			valEntry.grid(row=i,column=1)

			unitMenu = ttk.OptionMenu(self.propertiesFrame, unitVar, *options)
			unitMenu.grid(row=i, column=2)

			propertyEntries.append(valEntry)
			propertyUnitMenus.append(unitMenu)
			i+=1

	def prev(self):
		pass

app = Application()
app.root.mainloop()