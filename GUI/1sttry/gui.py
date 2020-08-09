import tkinter as tk
from tkinter import filedialog, ttk
import kernel

class HeatTransferSettingsWindow:
	def __init__(self, root, guiSettings):
		self.root = root
		self.guiSettings = guiSettings

		self.HEIGHT = 500
		self.WIDTH  = 600

		self.canvas = tk.Canvas(self.root, height=self.HEIGHT, width=self.WIDTH).pack()
		self.populate()

		self.BCFrame.place(relx=0.02,rely=0.22,relwidth=0.96,relheight=0.58, anchor="nw")
		self.BCCanvas.pack(side="left", fill="both", expand="yes")
		self.scrollbar.pack(side="right", fill="y")


	def populate(self):		
		# Open File Frame
		self.openFileFrame = tk.LabelFrame(self.canvas, text="Select Mesh File") 
		self.openFileFrame.place(relx=0.02,rely=0.02,relwidth=0.96,relheight=0.18, anchor="nw")
		
		self.openFileButton = tk.Button(self.openFileFrame, text="Open mesh file", command=self.askopenfilename)
		self.openFileButton.place(relx=0.10,rely=0.5,relwidth=0.25,relheight=0.4,anchor="w")
		
		self.fileLocationLabel = tk.Label(self.openFileFrame, text="", bg="white", anchor="e")
		self.fileLocationLabel.place(relx=0.40,rely=0.5,relwidth=0.50,relheight=0.4,anchor="w")

		# Boundary Conditions Selection Frame
		self.BCFrame = tk.LabelFrame(self.canvas)
		# self.BCFrame.place(relx=0.02,rely=0.22,relwidth=0.96,relheight=0.58, anchor="nw")

		self.BCCanvas = tk.Canvas(self.BCFrame)
		# self.BCCanvas.pack(side="left", fill="both", expand="yes")

		self.scrollbar = ttk.Scrollbar(self.BCFrame, orient="vertical", command=self.BCCanvas.yview)
		# self.scrollbar.pack(side="right", fill="y")

		self.BCCanvas.configure(yscrollcommand=self.scrollbar.set)
		self.BCCanvas.bind("<Configure>", lambda event:self.BCCanvas.configure(scrollregion=self.BCCanvas.bbox("all")))

		self.BCWindow = tk.Frame(self.BCCanvas, bg="yellow")
		self.BCCanvas.create_window((100,1000), window=self.BCWindow, anchor="nw")

		for i in range(50):
			tk.Button(self.BCWindow, text=f"My Button - {i}").pack()

		# Bottom Frame
		self.bottomFrame = tk.Frame(self.canvas, bd=5) 
		self.bottomFrame.place(relx=0.0,rely=0.85,relwidth=1.0,relheight=0.15, anchor="nw")

		self.prevButton = tk.Button(self.bottomFrame, text="Prev", command=self.prev)
		self.prevButton.place(relx=0.60, rely=0.50, relwidth=0.15, relheight=0.40, anchor="w")
		self.nextButton = tk.Button(self.bottomFrame, text="Next", command=self.next)
		self.nextButton.place(relx=0.80, rely=0.50, relwidth=0.15, relheight=0.40, anchor="w")

	def askopenfilename(self):
		self.meshFileName = tk.filedialog.askopenfilename(initialdir="../meshes", title="Select a mesh file", filetypes=(("MSH files", "*.msh"),("All files", "*")))
		self.fileLocationLabel["text"] = self.meshFileName

		self.guiSettings.setFilePath(self.meshFileName)

		# Check if mesh is valid, otherwise throw a warning message

		self.populateBCFrame()

	def populateBCFrame(self):
		boundaries = []
		index=0
		y0 = 0
		for boundary in self.guiSettings.getBoundaryNames():
			rely = 0.12*index + 0.02 + y0
			boundaryLabel = tk.Label(self.BCWindow, text=boundary, bg="white")
			# boundaryLabel.pack()
			boundaryLabel.place(relx=0.10,rely=rely,relwidth=0.25,relheight=0.1,anchor="nw")
			index+=1
		
	def prev(self):
		self.root.destroy()
	def next(self):
		pass




if __name__ == "__main__":
	root = tk.Tk()
	root.title("PyEFVLib GUI")
	root.bind("<Key>", lambda key: root.destroy() if key.char=="\x17" else 0) # Close window if Ctrl+W is pressed
	
	guiSettings = kernel.GUISettings()

	htsw = HeatTransferSettingsWindow(root, guiSettings)

	root.mainloop()