import tkinter as tk
from tkinter import ttk

root = tk.Tk()
root.geometry("600x600")
container = ttk.Frame(root)
canvas = tk.Canvas(container)
scrollbar = ttk.Scrollbar(container, orient="vertical", command=canvas.yview)
scrollable_frame = ttk.Frame(canvas)

scrollable_frame.bind(
    "<Configure>",
    lambda e: canvas.configure(
        scrollregion=canvas.bbox("all")
    )
)

canvas.create_window((0, 0), window=scrollable_frame, anchor="nw")

canvas.configure(yscrollcommand=scrollbar.set)

for i in range(50):
	options = ["=)", ">=(", ";D", ":P", "=)"]
	unitVar = tk.StringVar(scrollable_frame)
	unitVar.set("=)")
	bVar = tk.IntVar(scrollable_frame)
	ttk.Label(scrollable_frame, text="Sample scrolling label").grid(row=i, column=0)
	ttk.Entry(scrollable_frame).grid(row=i,column=1)
	ttk.OptionMenu(scrollable_frame, unitVar, *options).grid(row=i, column=2)
	tk.Radiobutton(scrollable_frame, variable=bVar, value=1).grid(row=i, column=3)
	tk.Radiobutton(scrollable_frame, variable=bVar, value=2).grid(row=i, column=4)


container.pack(side="top", fill="both", expand=True)
canvas.pack(side="left", fill="both", expand=True)
scrollbar.pack(side="right", fill="y")

root.mainloop()
