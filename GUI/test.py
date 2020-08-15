# This imports all functions in tkinter module 
from tkinter import * 
from tkinter.ttk import *

# creating master window 
master = Tk() 

# This method is used to get the position 
# of the desired widget available in any 
# othet widget 
def click(event): 
	
	# Here retrieving the size of the parent 
	# widget relative to master widget 
	x = event.x_root - f.winfo_rootx() 
	y = event.y_root - f.winfo_rooty() 

	# Here grid_location() method is used to 
	# retrieve the relative position on the 
	# parent widget 
	z = f.grid_location(x, y) 

	# printing position 
	print(z) 

# Frame widget, wil work as 
# parent for buttons widget 
f = Frame(master) 
f.pack() 

# Button widgets 
b = Button(f, text = "Button") 
b.grid(row = 2, column = 3) 

c = Button(f, text = "Button2") 
c.grid(row = 1, column = 0) 

# Here binding click method with mouse 
master.bind("<Button-1>", click) 

# infinite loop 
mainloop() 
