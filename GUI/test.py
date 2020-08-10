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

		self.bottomFrame = tk.Frame(self.canvas)
		self.bottomFrame.place(relx=0.02, rely=0.87, relheight=0.1, relwidth=0.96, anchor="nw" )

		self.nextButton = tk.Button(self.bottomFrame, text="Next")
		self.nextButton.place(relx=1.0, rely=0.50, relheight=0.75, relwidth=0.20, anchor="e")

		self.prevButton = tk.Button(self.bottomFrame, text="Prev")
		self.prevButton.place(relx=0.78, rely=0.50, relheight=0.75, relwidth=0.20, anchor="e")

		self.root.mainloop()

	def placeBCForms(self):
		self.boundariesNames = gridData = MSHReader(self.fileLabel["text"]).getData().boundariesNames
		
		countLabel = tk.Label(self.BCFrame, text=f"1/{len(self.boundariesNames)}")
		countLabel.place(relx=0.5,rely=0.0,anchor="n")

		boundaryNameLabel = tk.Label(self.BCFrame, text=self.boundariesNames[0])
		boundaryNameLabel.place(relx=0.5,rely=0.1,anchor="nw")
		# i=0
		# for boundaryName in self.boundariesNames:
		# 	label = tk.Label(self.BCFrame, text=boundaryName)
		# 	label.grid(row=i, column=0)
		# 	i+=1

	def openFile(self):
		self.fileLabel["text"] = tk.filedialog.askopenfilename(initialdir="../meshes", title="Select a mesh file", filetypes=(("MSH files", "*.msh"),("All files", "*")))
		self.placeBCForms()

app = Application()