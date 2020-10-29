import sys, os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

def show_fluxes(fileName):
	data = pd.read_csv(fileName)

	X1 = np.array(data["X1"])
	Y1 = np.array(data["Y1"])
	X2 = np.array(data["X2"])
	Y2 = np.array(data["Y2"])

	XC = [ sum(e)/4 for e in X1.reshape(X1.shape[0]//4, 4) for i in range(4)]
	YC = [ sum(e)/4 for e in Y1.reshape(Y1.shape[0]//4, 4) for i in range(4)]

	X = ((X1 + X2)/2 + XC)/2
	Y = ((Y1 + Y2)/2 + YC)/2

	U = np.array(data[ data.columns[-3] ])
	V = np.array(data[ data.columns[-2] ])

	C = np.array([np.linalg.norm((x,y)) for (x,y) in zip(X,Y)])

	fig, ax = plt.subplots()
	fig.canvas.set_window_title(fileName)
	ax.quiver(X,Y,U,V,C)
	# ax.scatter(X1,X2,marker=".",color="k")
	plt.show()
	fig.savefig(fileName.replace("csv", "png").replace("results", "images"))


if __name__ == "__main__":
	fileName = "fluxes.csv" if len(sys.argv)<2 else sys.argv[1]
	# show_fluxes(fileName)

	baseDir = "fluxos - vec\\results"
	for fileName in os.listdir(baseDir):
		show_fluxes(os.path.join(baseDir, fileName))
