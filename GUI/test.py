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

		self.BCScrollCanvas = tk.Canvas(self.BCFrame, bg="blue")
		self.BCScrollCanvas.pack(side="right", fill="both", expand="yes")

		self.BCScrollBar = tk.Scrollbar(self.BCFrame, orient="vertical", command=self.BCScrollCanvas.yview)
		self.BCScrollBar.place(relx=1.0,rely=0.50,relheight=1.0,width=15,anchor="e")

		self.BCScrollCanvas.configure(yscrollcommand=self.BCScrollBar.set)
		self.BCScrollCanvas.bind("<Configure>", lambda event: self.BCScrollCanvas.configure(scrollregion=self.BCScrollCanvas.bbox("all")))

		self.BCWindow = tk.Frame(self.BCScrollCanvas, bg="yellow")
		self.BCScrollCanvas.create_window((0,0), window=self.BCWindow, anchor="n")

		self.BCArea = tk.Canvas(self.BCWindow, width=WIDTH, height=HEIGHT)
		self.BCArea.pack(side="left", fill="both", expand="yes")

		self.bottomFrame = tk.Frame(self.canvas)
		self.bottomFrame.place(relx=0.02, rely=0.87, relheight=0.1, relwidth=0.96, anchor="nw" )

		self.nextButton = tk.Button(self.bottomFrame, text="Next")
		self.nextButton.place(relx=1.0, rely=0.50, relheight=0.75, relwidth=0.20, anchor="e")

		self.prevButton = tk.Button(self.bottomFrame, text="Prev")
		self.prevButton.place(relx=0.78, rely=0.50, relheight=0.75, relwidth=0.20, anchor="e")

		self.root.mainloop()

	def placeBCForms(self):
		boundariesNames = gridData = MSHReader(self.fileLabel["text"]).getData().boundariesNames
		
		i=0
		for boundaryName in boundariesNames:
			label = tk.Label(self.BCArea, text=boundaryName)
			label.grid(row=i, column=0)
			i+=1

	def openFile(self):
		self.fileLabel["text"] = tk.filedialog.askopenfilename(initialdir="../meshes", title="Select a mesh file", filetypes=(("MSH files", "*.msh"),("All files", "*")))
		self.placeBCForms()

app = Application()